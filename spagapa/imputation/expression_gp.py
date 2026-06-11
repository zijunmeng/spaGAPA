"""
Expression-informed Gaussian Process imputation for spatial APA data.

This module adds geometry-aware and expression-informed GP variants on top of
the coordinate-only GP imputer.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Literal, Optional

import numpy as np
from scipy.linalg import cho_factor, cho_solve

from .feature_builders import GeometryFeatureBuilder

logger = logging.getLogger(__name__)


Variant = Literal["spatial_radial", "additive", "product", "adaptive_additive", "layer_local"]
KernelType = Literal["rbf", "matern", "auto"]
LayerGateMode = Literal["radius", "pseudolayer"]


@dataclass
class _Standardizer:
    mean_: np.ndarray
    scale_: np.ndarray

    @classmethod
    def fit(cls, X: np.ndarray) -> "_Standardizer":
        X = np.asarray(X, dtype=float)
        mean = X.mean(axis=0)
        scale = X.std(axis=0)
        scale = np.where(scale > 0, scale, 1.0)
        return cls(mean_=mean, scale_=scale)

    def transform(self, X: np.ndarray) -> np.ndarray:
        X = np.asarray(X, dtype=float)
        return (X - self.mean_) / self.scale_


class ExpressionGPImputer:
    """
    Geometry-aware and expression-informed GP imputer.

    Parameters
    ----------
    variant : {"spatial_radial", "additive", "product", "adaptive_additive", "layer_local"}, default="additive"
        GP variant to use.
    spatial_kernel_type : {"rbf", "matern", "auto"}, default="matern"
        Kernel for the spatial component.
    expression_kernel_type : {"rbf", "matern", "auto"}, default="rbf"
        Kernel for the expression component in additive models.
    length_scale_space : float, optional
        Spatial kernel length scale. If None, estimated from training data.
    length_scale_expr : float, optional
        Expression kernel length scale. If None, estimated from training data.
    lambda_expr : float, default=1.0
        Weight of the expression kernel in expression-informed models.
    product_offset : float, default=1.0
        Baseline multiplier for product kernels:
        `K_total = K_space * (product_offset + lambda_expr * K_expr)`.
    gate_mode : {"distance"}, default="distance"
        Gating strategy for adaptive expression kernels.
    gate_tau : float, optional
        Spatial gating length scale for adaptive kernels. If None, estimated
        from the spatial training features.
    local_k : int, default=20
        Neighborhood size used to estimate the soft local gate length scale
        for layer-local kernels.
    layer_gate_mode : {"radius", "pseudolayer"}, default="radius"
        Local gate feature for layer-local kernels. Both modes are derived
        from coordinates and do not use true layer labels.
    nu : float, default=1.5
        Matérn smoothness parameter.
    alpha : float, default=1e-3
        Jitter / regularization term added to the kernel diagonal.
    use_theta : bool, default=False
        Whether to include polar angle in the spatial-radial feature set.
    normalize_y : bool, default=True
        Whether to z-score target values before fitting.
    """

    def __init__(
        self,
        variant: Variant = "additive",
        spatial_kernel_type: KernelType = "matern",
        expression_kernel_type: KernelType = "rbf",
        length_scale_space: Optional[float] = None,
        length_scale_expr: Optional[float] = None,
        lambda_expr: float = 1.0,
        product_offset: float = 1.0,
        gate_mode: str = "distance",
        gate_tau: Optional[float] = None,
        local_k: int = 20,
        layer_gate_mode: LayerGateMode = "radius",
        nu: float = 1.5,
        alpha: float = 1e-3,
        use_theta: bool = False,
        normalize_y: bool = True,
    ):
        if variant not in {"spatial_radial", "additive", "product", "adaptive_additive", "layer_local"}:
            raise ValueError(f"Unknown variant: {variant}")
        if product_offset < 0:
            raise ValueError("product_offset must be non-negative")
        if gate_mode != "distance":
            raise ValueError(f"Unsupported gate_mode: {gate_mode}")
        if gate_tau is not None and gate_tau <= 0:
            raise ValueError("gate_tau must be positive")
        if local_k <= 0:
            raise ValueError("local_k must be positive")
        if layer_gate_mode not in {"radius", "pseudolayer"}:
            raise ValueError(f"Unsupported layer_gate_mode: {layer_gate_mode}")
        self.variant = variant
        self.spatial_kernel_type = spatial_kernel_type
        self.expression_kernel_type = expression_kernel_type
        self.length_scale_space = length_scale_space
        self.length_scale_expr = length_scale_expr
        self.lambda_expr = lambda_expr
        self.product_offset = product_offset
        self.gate_mode = gate_mode
        self.gate_tau = gate_tau
        self.local_k = int(local_k)
        self.layer_gate_mode = layer_gate_mode
        self.nu = nu
        self.alpha = alpha
        self.use_theta = use_theta
        self.normalize_y = normalize_y

        self._space_standardizer: Optional[_Standardizer] = None
        self._expr_standardizer: Optional[_Standardizer] = None
        self._local_standardizer: Optional[_Standardizer] = None
        self._train_space: Optional[np.ndarray] = None
        self._train_expr: Optional[np.ndarray] = None
        self._train_local: Optional[np.ndarray] = None
        self._train_values: Optional[np.ndarray] = None
        self._cho = None
        self._alpha_vec = None
        self._y_mean = 0.0
        self._y_scale = 1.0
        self._signal_variance = 1.0
        self._gate_tau_fit = None
        self._local_gate_tau_fit = None
        self._local_center_fit: Optional[np.ndarray] = None
        self._pseudolayer_edges_fit: Optional[np.ndarray] = None

    def _kernel(self, X1: np.ndarray, X2: np.ndarray, kernel_type: KernelType, length_scale: float) -> np.ndarray:
        if kernel_type == "auto":
            kernel_type = "matern"

        dists = self._pairwise_distance(X1, X2)
        if kernel_type == "rbf":
            return np.exp(-0.5 * (dists / length_scale) ** 2)
        if kernel_type == "matern":
            if self.nu == 0.5:
                return np.exp(-dists / length_scale)
            if self.nu == 1.5:
                scaled = np.sqrt(3.0) * dists / length_scale
                return (1.0 + scaled) * np.exp(-scaled)
            if self.nu == 2.5:
                scaled = np.sqrt(5.0) * dists / length_scale
                return (1.0 + scaled + (scaled ** 2) / 3.0) * np.exp(-scaled)
            raise ValueError("Unsupported nu for matern kernel; use one of {0.5, 1.5, 2.5}")
        raise ValueError(f"Unknown kernel type: {kernel_type}")

    @staticmethod
    def _pairwise_distance(X1: np.ndarray, X2: np.ndarray) -> np.ndarray:
        diff = X1[:, None, :] - X2[None, :, :]
        return np.sqrt(np.sum(diff * diff, axis=2))

    @staticmethod
    def _estimate_length_scale(X: np.ndarray) -> float:
        if len(X) <= 1:
            return 1.0
        dists = ExpressionGPImputer._pairwise_distance(X, X)
        dists = dists[np.triu_indices_from(dists, k=1)]
        finite = dists[np.isfinite(dists) & (dists > 0)]
        if finite.size == 0:
            return 1.0
        return float(np.median(finite))

    def _prepare_spatial_features(
        self,
        coordinates: np.ndarray,
        geometry_features: dict[str, np.ndarray] | None,
    ) -> np.ndarray:
        coords = np.asarray(coordinates, dtype=float)
        if coords.ndim != 2 or coords.shape[1] != 2:
            raise ValueError("coordinates must have shape (n_spots, 2)")

        if self.variant != "spatial_radial":
            return coords

        if geometry_features is None:
            geometry_features = GeometryFeatureBuilder(use_theta=self.use_theta).build(coords)

        parts = [coords]
        radius = geometry_features.get("radius")
        if radius is None:
            raise ValueError("geometry_features must contain 'radius' for spatial_radial variant")
        parts.append(np.asarray(radius, dtype=float))

        if self.use_theta:
            theta = geometry_features.get("theta")
            if theta is None:
                raise ValueError("geometry_features must contain 'theta' when use_theta=True")
            parts.append(np.asarray(theta, dtype=float))

        return np.concatenate(parts, axis=1)

    def _uses_expression(self) -> bool:
        return self.variant in {"additive", "product", "adaptive_additive", "layer_local"}

    def _uses_adaptive_gate(self) -> bool:
        return self.variant == "adaptive_additive"

    def _uses_local_gate(self) -> bool:
        return self.variant == "layer_local"

    def _prepare_expression_features(self, expression_embedding: np.ndarray | None) -> np.ndarray:
        if not self._uses_expression():
            return np.empty((0, 0))
        if expression_embedding is None:
            raise ValueError(f"expression_embedding is required for {self.variant} variant")
        expr = np.asarray(expression_embedding, dtype=float)
        if expr.ndim != 2:
            raise ValueError("expression_embedding must be 2-dimensional")
        return expr

    def _adaptive_gate_kernel(
        self,
        X1: np.ndarray,
        X2: np.ndarray,
        length_scale: float,
    ) -> np.ndarray:
        dists = self._pairwise_distance(X1, X2)
        return np.exp(-0.5 * (dists / length_scale) ** 2)

    def _prepare_local_features(
        self,
        coordinates: np.ndarray,
        geometry_features: dict[str, np.ndarray] | None,
        fit: bool = False,
    ) -> np.ndarray:
        coords = np.asarray(coordinates, dtype=float)
        if coords.ndim != 2 or coords.shape[1] != 2:
            raise ValueError("coordinates must have shape (n_spots, 2)")

        radius = None
        if geometry_features is not None and "radius" in geometry_features:
            radius = np.asarray(geometry_features["radius"], dtype=float).reshape(-1)
            if radius.shape[0] != coords.shape[0]:
                raise ValueError("geometry_features['radius'] must match coordinates")

        if radius is None:
            if fit or self._local_center_fit is None:
                center = coords.mean(axis=0)
                if fit:
                    self._local_center_fit = center
            else:
                center = self._local_center_fit
            radius = np.linalg.norm(coords - center, axis=1)

        if self.layer_gate_mode == "radius":
            return radius[:, None]

        if fit:
            edges = np.quantile(radius, np.linspace(0.0, 1.0, 6)[1:-1])
            self._pseudolayer_edges_fit = np.unique(edges)
        if self._pseudolayer_edges_fit is None:
            raise ValueError("Pseudolayer gate is not fitted")
        pseudolayer = np.digitize(radius, self._pseudolayer_edges_fit).astype(float)
        return pseudolayer[:, None]

    @staticmethod
    def _estimate_knn_length_scale(X: np.ndarray, k: int) -> float:
        if len(X) <= 1:
            return 1.0
        dists = ExpressionGPImputer._pairwise_distance(X, X)
        np.fill_diagonal(dists, np.inf)
        kth = min(max(1, int(k)), len(X) - 1)
        kth_dists = np.partition(dists, kth - 1, axis=1)[:, kth - 1]
        finite = kth_dists[np.isfinite(kth_dists) & (kth_dists > 0)]
        if finite.size == 0:
            return ExpressionGPImputer._estimate_length_scale(X)
        return float(np.median(finite))

    def _local_gate_kernel(
        self,
        X1: np.ndarray,
        X2: np.ndarray,
        length_scale: float,
    ) -> np.ndarray:
        dists = self._pairwise_distance(X1, X2)
        return np.exp(-0.5 * (dists / length_scale) ** 2)

    def _combine_kernels(
        self,
        K_space: np.ndarray,
        K_expr: np.ndarray | None,
        K_gate: np.ndarray | None = None,
    ) -> np.ndarray:
        if self.variant == "additive":
            if K_expr is None:
                raise ValueError("K_expr is required for additive variant")
            return K_space + self.lambda_expr * K_expr
        if self.variant == "product":
            if K_expr is None:
                raise ValueError("K_expr is required for product variant")
            return K_space * (self.product_offset + self.lambda_expr * K_expr)
        if self.variant == "adaptive_additive":
            if K_expr is None:
                raise ValueError("K_expr is required for adaptive_additive variant")
            if K_gate is None:
                raise ValueError("K_gate is required for adaptive_additive variant")
            return K_space + self.lambda_expr * (K_gate * K_expr)
        if self.variant == "layer_local":
            if K_expr is None:
                raise ValueError("K_expr is required for layer_local variant")
            if K_gate is None:
                raise ValueError("K_gate is required for layer_local variant")
            return K_space * (self.product_offset + self.lambda_expr * (K_gate * K_expr))
        return K_space

    def fit(
        self,
        coordinates: np.ndarray,
        values: np.ndarray,
        mask: Optional[np.ndarray] = None,
        expression_embedding: Optional[np.ndarray] = None,
        geometry_features: Optional[dict[str, np.ndarray]] = None,
    ) -> "ExpressionGPImputer":
        values = np.asarray(values, dtype=float)
        if values.ndim != 1:
            raise ValueError("values must be 1-dimensional")

        if mask is None:
            mask = values > 0
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != values.shape:
            raise ValueError("mask must have the same shape as values")

        if mask.sum() == 0:
            raise ValueError("No training data available (all values are zero or masked)")

        space = self._prepare_spatial_features(coordinates, geometry_features)
        expr = self._prepare_expression_features(expression_embedding)
        local = (
            self._prepare_local_features(coordinates, geometry_features, fit=True)
            if self._uses_local_gate()
            else None
        )

        train_space = space[mask]
        train_values = values[mask]
        train_expr = expr[mask] if self._uses_expression() else None
        train_local = local[mask] if self._uses_local_gate() else None

        self._space_standardizer = _Standardizer.fit(train_space)
        train_space = self._space_standardizer.transform(train_space)

        if self._uses_expression():
            self._expr_standardizer = _Standardizer.fit(train_expr)
            train_expr = self._expr_standardizer.transform(train_expr)
        if self._uses_local_gate():
            self._local_standardizer = _Standardizer.fit(train_local)
            train_local = self._local_standardizer.transform(train_local)

        if self.length_scale_space is None:
            ls_space = self._estimate_length_scale(train_space)
        else:
            ls_space = self.length_scale_space
        gate_tau = self.gate_tau if self.gate_tau is not None else ls_space
        local_gate_tau = (
            self.gate_tau
            if self.gate_tau is not None
            else self._estimate_knn_length_scale(train_local, self.local_k)
            if self._uses_local_gate()
            else None
        )

        if self.normalize_y:
            self._y_mean = float(train_values.mean())
            self._y_scale = float(train_values.std())
            if self._y_scale == 0:
                self._y_scale = 1.0
            y_train = (train_values - self._y_mean) / self._y_scale
            self._signal_variance = 1.0
        else:
            self._y_mean = 0.0
            self._y_scale = 1.0
            y_train = train_values
            var = float(np.var(train_values))
            self._signal_variance = var if var > 0 else 1.0

        K_space = self._kernel(train_space, train_space, self.spatial_kernel_type, ls_space)
        K_total = K_space

        self._length_scale_space_fit = ls_space
        self._length_scale_expr_fit = None
        self._gate_tau_fit = gate_tau if self._uses_adaptive_gate() else None
        self._local_gate_tau_fit = local_gate_tau if self._uses_local_gate() else None

        if self._uses_expression():
            if self.length_scale_expr is None:
                ls_expr = self._estimate_length_scale(train_expr)
            else:
                ls_expr = self.length_scale_expr
            K_expr = self._kernel(train_expr, train_expr, self.expression_kernel_type, ls_expr)
            K_gate = None
            if self._uses_adaptive_gate():
                K_gate = self._adaptive_gate_kernel(train_space, train_space, gate_tau)
            if self._uses_local_gate():
                K_gate = self._local_gate_kernel(train_local, train_local, local_gate_tau)
            K_total = self._combine_kernels(K_space, K_expr, K_gate=K_gate)
            self._length_scale_expr_fit = ls_expr

        K_total = self._signal_variance * K_total
        K_total = K_total + self.alpha * np.eye(K_total.shape[0])

        self._cho = cho_factor(K_total, lower=True, check_finite=False)
        self._alpha_vec = cho_solve(self._cho, y_train, check_finite=False)
        self._train_space = train_space
        self._train_expr = train_expr
        self._train_local = train_local
        self._train_values = train_values
        self._mask = mask.copy()
        return self

    def predict(
        self,
        coordinates: np.ndarray,
        expression_embedding: Optional[np.ndarray] = None,
        geometry_features: Optional[dict[str, np.ndarray]] = None,
        return_std: bool = True,
    ) -> tuple[np.ndarray, Optional[np.ndarray]]:
        if self._cho is None or self._alpha_vec is None or self._train_space is None:
            raise ValueError("Model not fitted. Call fit() first.")

        test_space = self._prepare_spatial_features(coordinates, geometry_features)
        test_space = self._space_standardizer.transform(test_space)
        K_cross_space = self._kernel(
            self._train_space,
            test_space,
            self.spatial_kernel_type,
            self._length_scale_space_fit,
        )
        K_cross = K_cross_space

        if self._uses_expression():
            test_expr = self._prepare_expression_features(expression_embedding)
            test_expr = self._expr_standardizer.transform(test_expr)
            K_cross_expr = self._kernel(
                self._train_expr,
                test_expr,
                self.expression_kernel_type,
                self._length_scale_expr_fit,
            )
            K_gate_cross = None
            if self._uses_adaptive_gate():
                K_gate_cross = self._adaptive_gate_kernel(
                    self._train_space,
                    test_space,
                    self._gate_tau_fit,
                )
            if self._uses_local_gate():
                test_local = self._prepare_local_features(coordinates, geometry_features, fit=False)
                test_local = self._local_standardizer.transform(test_local)
                K_gate_cross = self._local_gate_kernel(
                    self._train_local,
                    test_local,
                    self._local_gate_tau_fit,
                )
            K_cross = self._combine_kernels(K_cross_space, K_cross_expr, K_gate=K_gate_cross)
        else:
            test_expr = None
            test_local = None

        K_cross = self._signal_variance * K_cross
        pred_norm = K_cross.T @ self._alpha_vec
        predictions = pred_norm * self._y_scale + self._y_mean

        if not return_std:
            return predictions, None

        K_self_space = self._kernel(
            test_space,
            test_space,
            self.spatial_kernel_type,
            self._length_scale_space_fit,
        )
        K_self = K_self_space
        if self._uses_expression():
            K_self_expr = self._kernel(
                test_expr,
                test_expr,
                self.expression_kernel_type,
                self._length_scale_expr_fit,
            )
            K_gate_self = None
            if self._uses_adaptive_gate():
                K_gate_self = self._adaptive_gate_kernel(
                    test_space,
                    test_space,
                    self._gate_tau_fit,
                )
            if self._uses_local_gate():
                K_gate_self = self._local_gate_kernel(
                    test_local,
                    test_local,
                    self._local_gate_tau_fit,
                )
            K_self = self._combine_kernels(K_self_space, K_self_expr, K_gate=K_gate_self)

        K_self = self._signal_variance * K_self
        v = cho_solve(self._cho, K_cross, check_finite=False)
        var = np.clip(np.diag(K_self) - np.sum(K_cross * v, axis=0), a_min=0.0, a_max=None)
        std = np.sqrt(var) * self._y_scale
        return predictions, std

    def impute(
        self,
        coordinates: np.ndarray,
        values: np.ndarray,
        mask: Optional[np.ndarray] = None,
        expression_embedding: Optional[np.ndarray] = None,
        geometry_features: Optional[dict[str, np.ndarray]] = None,
        return_uncertainty: bool = True,
    ) -> tuple[np.ndarray, Optional[np.ndarray]]:
        self.fit(
            coordinates,
            values,
            mask=mask,
            expression_embedding=expression_embedding,
            geometry_features=geometry_features,
        )
        predictions, uncertainty = self.predict(
            coordinates,
            expression_embedding=expression_embedding,
            geometry_features=geometry_features,
            return_std=return_uncertainty,
        )

        values = np.asarray(values, dtype=float)
        if mask is None:
            mask = values > 0
        mask = np.asarray(mask, dtype=bool)
        predictions[mask] = values[mask]
        if uncertainty is not None:
            uncertainty[mask] = 0.0
        return predictions, uncertainty

    def fit_batch(
        self,
        coordinates: np.ndarray,
        values: np.ndarray,
        mask: Optional[np.ndarray] = None,
        expression_embedding: Optional[np.ndarray] = None,
        geometry_features: Optional[dict[str, np.ndarray]] = None,
        n_jobs: int = 1,
        verbose: bool = True,
    ) -> "ExpressionGPImputerBatch":
        return ExpressionGPImputerBatch(
            base_imputer=self,
            coordinates=coordinates,
            values=values,
            mask=mask,
            expression_embedding=expression_embedding,
            geometry_features=geometry_features,
            n_jobs=n_jobs,
            verbose=verbose,
        )


class ExpressionGPImputerBatch:
    """
    Batch wrapper for ExpressionGPImputer.

    MVP implementation fits genes sequentially to keep feature sharing and
    memory behavior simple and predictable.
    """

    def __init__(
        self,
        base_imputer: ExpressionGPImputer,
        coordinates: np.ndarray,
        values: np.ndarray,
        mask: Optional[np.ndarray] = None,
        expression_embedding: Optional[np.ndarray] = None,
        geometry_features: Optional[dict[str, np.ndarray]] = None,
        n_jobs: int = 1,
        verbose: bool = True,
    ):
        self.base_imputer = base_imputer
        self.coordinates = np.asarray(coordinates, dtype=float)
        self.values = np.asarray(values, dtype=float)
        self.mask = mask
        self.expression_embedding = expression_embedding
        self.geometry_features = geometry_features
        self.n_jobs = n_jobs
        self.verbose = verbose
        self.imputers_: list[Optional[ExpressionGPImputer]] = []
        self._fit_batch()

    def _new_imputer(self) -> ExpressionGPImputer:
        return ExpressionGPImputer(
            variant=self.base_imputer.variant,
            spatial_kernel_type=self.base_imputer.spatial_kernel_type,
            expression_kernel_type=self.base_imputer.expression_kernel_type,
            length_scale_space=self.base_imputer.length_scale_space,
            length_scale_expr=self.base_imputer.length_scale_expr,
            lambda_expr=self.base_imputer.lambda_expr,
            product_offset=self.base_imputer.product_offset,
            gate_mode=self.base_imputer.gate_mode,
            gate_tau=self.base_imputer.gate_tau,
            local_k=self.base_imputer.local_k,
            layer_gate_mode=self.base_imputer.layer_gate_mode,
            nu=self.base_imputer.nu,
            alpha=self.base_imputer.alpha,
            use_theta=self.base_imputer.use_theta,
            normalize_y=self.base_imputer.normalize_y,
        )

    def _fit_batch(self) -> None:
        n_genes = self.values.shape[0]
        self.imputers_ = [None] * n_genes
        if self.n_jobs != 1:
            logger.info("ExpressionGPImputerBatch currently fits sequentially; ignoring n_jobs=%s", self.n_jobs)

        iterator = range(n_genes)
        for gene_idx in iterator:
            gene_values = self.values[gene_idx, :]
            gene_mask = self.mask[gene_idx, :] if self.mask is not None else None
            imputer = self._new_imputer()
            try:
                imputer.fit(
                    self.coordinates,
                    gene_values,
                    mask=gene_mask,
                    expression_embedding=self.expression_embedding,
                    geometry_features=self.geometry_features,
                )
                self.imputers_[gene_idx] = imputer
            except Exception as exc:
                logger.warning("Failed to fit expression GP for gene %s: %s", gene_idx, exc)
                self.imputers_[gene_idx] = None

    def predict(
        self,
        coordinates: Optional[np.ndarray] = None,
        expression_embedding: Optional[np.ndarray] = None,
        geometry_features: Optional[dict[str, np.ndarray]] = None,
        return_std: bool = True,
    ) -> tuple[np.ndarray, Optional[np.ndarray]]:
        if coordinates is None:
            coordinates = self.coordinates
        if expression_embedding is None:
            expression_embedding = self.expression_embedding
        if geometry_features is None:
            geometry_features = self.geometry_features

        n_genes = len(self.imputers_)
        n_spots = len(coordinates)
        predictions = np.zeros((n_genes, n_spots), dtype=float)
        uncertainty = np.zeros((n_genes, n_spots), dtype=float) if return_std else None

        for i, imputer in enumerate(self.imputers_):
            if imputer is None:
                continue
            pred, unc = imputer.predict(
                coordinates,
                expression_embedding=expression_embedding,
                geometry_features=geometry_features,
                return_std=return_std,
            )
            predictions[i, :] = pred
            if uncertainty is not None and unc is not None:
                uncertainty[i, :] = unc

        return predictions, uncertainty

    def impute(
        self,
        return_uncertainty: bool = True,
    ) -> tuple[np.ndarray, Optional[np.ndarray]]:
        predictions, uncertainty = self.predict(
            self.coordinates,
            expression_embedding=self.expression_embedding,
            geometry_features=self.geometry_features,
            return_std=return_uncertainty,
        )

        if self.mask is not None:
            predictions[self.mask] = self.values[self.mask]
            if uncertainty is not None:
                uncertainty[self.mask] = 0.0
        else:
            observed = self.values > 0
            predictions[observed] = self.values[observed]
            if uncertainty is not None:
                uncertainty[observed] = 0.0

        return predictions, uncertainty
