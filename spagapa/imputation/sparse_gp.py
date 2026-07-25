"""
Sparse Gaussian Process approximation for large-scale spatial data.

This module implements sparse GP methods (inducing points) to handle
large spatial transcriptomics datasets efficiently.
"""

import numpy as np
from typing import Optional, Tuple, Literal
from scipy.spatial import distance_matrix
from scipy.linalg import cho_solve, cho_factor, solve_triangular
from sklearn.cluster import KMeans
import logging

logger = logging.getLogger(__name__)


class SparseGPImputer:
    """
    Sparse Gaussian Process imputer using inducing points.
    
    Uses a subset of inducing points to approximate the full GP,
    reducing computational complexity from O(n³) to O(nm²) where
    m << n is the number of inducing points.
    
    Parameters
    ----------
    n_inducing : int, default=100
        Number of inducing points
    inducing_method : str, default='kmeans'
        Method for selecting inducing points: 'kmeans', 'random', or 'grid'
    length_scale : float, default=1.0
        Length scale for RBF kernel.  When left at the default AND no
        ``median_nn_dist`` is supplied to ``fit``, the caller is expected to
        scale this externally (e.g. the head-to-head benchmark passes
        ``length_scale = nn_dist * length_scale_multiplier``).
    noise_level : float, default=0.1
        Observation noise level.  Used as-is when ``local_noise`` is False;
        overridden by a per-gene estimate when ``local_noise`` is True and
        ``fit`` can compute a local-variance estimate.
    length_scale_multiplier : float or None, default=None
        When set to a float, the effective length scale is derived as
        ``median_nn_dist * length_scale_multiplier`` (the caller supplies
        ``median_nn_dist`` to ``fit`` / ``fit_batch``; the batch wrapper
        computes it automatically).  When None (the default), the explicit
        ``length_scale`` is used directly -- this preserves the historical
        behaviour for callers that pass an explicit ``length_scale``.  The
        head-to-head benchmark opts into the multiplier path; smaller values
        (e.g. 1.0-2.0) reduce over-smoothing versus the historical value of
        ~5.0.
    local_noise : bool, default=False
        When True, ``fit`` estimates a per-gene noise floor from the local
        variance of observed values instead of using the flat ``noise_level``
        constant.  This makes the predictive uncertainty track the actual
        error (a flat 0.1 noise_level yields ~constant std across tissue).
    noise_scale : float, default=1000.0
        Multiplier applied to the local-noise variance estimate when
        ``local_noise`` is True.  The RBF kernel has a fixed unit amplitude,
        so this scaling brings the per-gene noise estimate (typically ~1e-3
        for a [0,1]-bounded APA index) into the O(1) kernel-amplitude regime
        where the GP is appropriately skeptical of noisy, weakly-spatial
        data.  Empirically validated to turn Moran's-I recovery from
        negative (~-0.32) to near-zero/positive and to roughly double the
        uncertainty-error correlation.

    Attributes
    ----------
    inducing_points_ : np.ndarray
        Selected inducing point locations
    alpha_ : np.ndarray
        Fitted coefficients

    Examples
    --------
    >>> imputer = SparseGPImputer(n_inducing=100)
    >>> imputer.fit(spatial_coords, apa_counts)
    >>> imputed = imputer.predict(spatial_coords)
    """

    def __init__(
        self,
        n_inducing: int = 100,
        inducing_method: Literal['kmeans', 'random', 'grid'] = 'kmeans',
        length_scale: float = 1.0,
        noise_level: float = 0.1,
        kernel_type: Literal['rbf'] = 'rbf',
        length_scale_multiplier: Optional[float] = None,
        local_noise: bool = False,
        noise_scale: float = 1000.0,
    ):
        if kernel_type != 'rbf':
            raise ValueError("SparseGPImputer currently supports only kernel_type='rbf'")
        self.n_inducing = n_inducing
        self.inducing_method = inducing_method
        self.length_scale = length_scale
        self.noise_level = noise_level
        self.kernel_type = kernel_type
        self.length_scale_multiplier = length_scale_multiplier
        self.local_noise = local_noise
        self.noise_scale = noise_scale

        self.inducing_points_ = None
        self.alpha_ = None
        self._K_mm_inv = None
        self._train_coords = None
        self._train_values = None
        # Per-fit effective length scale (after multiplier applied) and the
        # per-gene noise estimate used by this fit.  Exposed for inspection
        # / unit tests.
        self._effective_length_scale = None
        self._effective_noise = None
    
    def _rbf_kernel(
        self,
        X1: np.ndarray,
        X2: np.ndarray,
        length_scale: Optional[float] = None
    ) -> np.ndarray:
        """
        Compute RBF (Gaussian) kernel matrix.
        
        Parameters
        ----------
        X1 : np.ndarray, shape (n1, d)
            First set of points
        X2 : np.ndarray, shape (n2, d)
            Second set of points
        length_scale : float, optional
            Length scale parameter
        
        Returns
        -------
        np.ndarray, shape (n1, n2)
            Kernel matrix
        """
        if length_scale is None:
            length_scale = self.length_scale
        
        dists = distance_matrix(X1, X2)
        return np.exp(-0.5 * (dists / length_scale) ** 2)
    
    def _select_inducing_points(
        self,
        coordinates: np.ndarray
    ) -> np.ndarray:
        """
        Select inducing points from training data.
        
        Parameters
        ----------
        coordinates : np.ndarray, shape (n, 2)
            Training coordinates
        
        Returns
        -------
        np.ndarray, shape (m, 2)
            Inducing point locations
        """
        n_points = len(coordinates)
        n_inducing = min(self.n_inducing, n_points)
        
        if self.inducing_method == 'kmeans':
            # Use k-means clustering to find representative points
            kmeans = KMeans(n_clusters=n_inducing, random_state=42, n_init=10)
            kmeans.fit(coordinates)
            inducing_points = kmeans.cluster_centers_
            logger.info(f"Selected {n_inducing} inducing points using k-means")
        
        elif self.inducing_method == 'random':
            # Random subset of training points
            indices = np.random.choice(n_points, size=n_inducing, replace=False)
            inducing_points = coordinates[indices]
            logger.info(f"Selected {n_inducing} inducing points randomly")
        
        elif self.inducing_method == 'grid':
            # Regular grid covering the data range
            x_min, x_max = coordinates[:, 0].min(), coordinates[:, 0].max()
            y_min, y_max = coordinates[:, 1].min(), coordinates[:, 1].max()
            
            # Create grid
            n_per_dim = int(np.sqrt(n_inducing))
            x_grid = np.linspace(x_min, x_max, n_per_dim)
            y_grid = np.linspace(y_min, y_max, n_per_dim)
            xx, yy = np.meshgrid(x_grid, y_grid)
            inducing_points = np.column_stack([xx.ravel(), yy.ravel()])
            logger.info(f"Selected {len(inducing_points)} inducing points on grid")
        
        else:
            raise ValueError(f"Unknown inducing method: {self.inducing_method}")
        
        return inducing_points
    
    def _estimate_local_noise(
        self,
        train_coords: np.ndarray,
        train_values: np.ndarray,
    ) -> float:
        """Estimate a per-gene noise floor from the local variance of
        observed values.

        Uses a kNN graph over the training coordinates and measures the
        median absolute deviation of each observation from the median of
        its neighbours' values (a robust local-scale estimate of the noise,
        i.e. the variability left AFTER removing the local spatial mean).
        Falls back to a fraction of the observed variance when too few
        neighbours are usable.

        The estimate is scaled by ``noise_scale`` before being returned.
        This scaling is necessary because the RBF kernel used here has a
        FIXED unit amplitude (signal variance = 1) regardless of the data
        scale, so the GP's noise-to-signal ratio is ``noise / 1.0``.  For
        a noisy, weakly-spatial field (e.g. a gene-level APA index whose
        observed variance is ~0.005) the raw local-noise variance is
        ~1e-3 -- orders of magnitude below the kernel amplitude -- and the
        GP would treat the data as an almost-noiseless smooth field,
        imposing strong spatial structure that does not exist (negative
        Moran's-I recovery).  Scaling the local estimate by ``noise_scale``
        (default 1000) brings the noise into the same O(1) regime as the
        kernel so the GP is appropriately skeptical: it predicts close to
        the local mean without over-smoothing, which turns Moran's-I
        recovery positive and makes predictive std track error.

        Parameters
        ----------
        train_coords : ndarray, shape (n_train, 2)
        train_values : ndarray, shape (n_train,)

        Returns
        -------
        float
            Estimated noise variance in kernel-amplitude units (always > 0).
        """
        n = len(train_values)
        scale = float(self.noise_scale)
        if n < 4:
            var = float(np.var(train_values)) if np.var(train_values) > 0 else self.noise_level
            return max(var * scale, 1e-3)
        try:
            from scipy.spatial import cKDTree
            k = min(8, n - 1)
            tree = cKDTree(train_coords)
            _, idx = tree.query(train_coords, k=k + 1)
            nbr_idx = idx[:, 1:]  # exclude self
            nbr_vals = train_values[nbr_idx]  # (n, k)
            # MAD of each point from its neighbours' median = robust local noise
            nbr_med = np.median(nbr_vals, axis=1)
            dev = np.abs(train_values - nbr_med)
            mad = float(np.median(dev))
            # Convert MAD to a variance-like scale (1.4826^2 * MAD^2 ~ variance
            # for Gaussian noise); scale into kernel-amplitude units; guard 0.
            noise_var = max((1.4826 * mad) ** 2 * scale, 1e-3)
            return noise_var
        except Exception:
            var = float(np.var(train_values))
            return max(var * scale, 1e-3)

    def fit(
        self,
        coordinates: np.ndarray,
        values: np.ndarray,
        mask: Optional[np.ndarray] = None,
        inducing_points: Optional[np.ndarray] = None,
        K_mm_inv: Optional[np.ndarray] = None,
        k_nm_full: Optional[np.ndarray] = None,
        median_nn_dist: Optional[float] = None,
    ):
        """
        Fit sparse GP model.

        Parameters
        ----------
        coordinates : np.ndarray, shape (n_spots, 2)
            Spatial coordinates
        values : np.ndarray, shape (n_spots,)
            Observed values
        mask : np.ndarray, optional
            Boolean mask for training data
        inducing_points : np.ndarray, optional, shape (m, 2)
            Precomputed inducing points. When supplied, the (expensive)
            KMeans selection is skipped. Inducing points depend only on
            coordinates, so a single set computed across all spots can be
            reused for every gene in a batch.
        K_mm_inv : np.ndarray, optional, shape (m, m)
            Precomputed inverse of the inducing-point kernel ``K_mm``.
            Must be consistent with ``inducing_points``. Skipping the
            per-gene ``np.linalg.inv`` is the second half of the reuse
            optimization.
        k_nm_full : np.ndarray, optional, shape (n_spots, m)
            Precomputed kernel between ALL coordinates (not just the
            mask-filtered training subset) and ``inducing_points``. When
            provided, this gene's ``K_nm`` is obtained by indexing
            ``k_nm_full`` with the training mask, skipping the per-gene
            ``distance_matrix`` evaluation. Must be consistent with
            ``inducing_points`` and the supplied ``coordinates``.
        median_nn_dist : float, optional
            Median nearest-neighbour distance of the full coordinate set.
            When supplied, the effective length scale is set to
            ``median_nn_dist * length_scale_multiplier`` (unless an explicit
            ``length_scale`` override is in effect).  This lets callers pass
            the spatial scale once and control smoothing via the multiplier.

        Returns
        -------
        self
        """
        # Filter training data
        if mask is None:
            mask = values > 0
        mask = np.asarray(mask, dtype=bool)

        train_coords = coordinates[mask]
        train_values = values[mask]

        if len(train_values) == 0:
            raise ValueError("No training data available")

        # ---- effective length scale ----
        # The multiplier path (eff_length = median_nn_dist * multiplier) is
        # opt-in: it activates only when length_scale_multiplier is set to a
        # float AND median_nn_dist is supplied.  Otherwise the explicit
        # length_scale is used (backward-compatible default).
        if (self.length_scale_multiplier is not None
                and median_nn_dist is not None and median_nn_dist > 0):
            eff_length = float(median_nn_dist) * float(self.length_scale_multiplier)
        else:
            eff_length = float(self.length_scale)
        self._effective_length_scale = eff_length

        # ---- effective noise ----
        # local_noise replaces the flat 0.1 with a per-gene estimate so the
        # predictive std tracks error instead of being ~constant across tissue.
        if self.local_noise:
            eff_noise = self._estimate_local_noise(train_coords, train_values)
        else:
            eff_noise = float(self.noise_level)
        self._effective_noise = eff_noise

        # Select inducing points (or reuse a shared precomputed set)
        if inducing_points is not None:
            self.inducing_points_ = inducing_points
        else:
            self.inducing_points_ = self._select_inducing_points(train_coords)
        m = len(self.inducing_points_)
        n = len(train_coords)

        # K_nm: kernel between training points and inducing points (n x m).
        # When the caller supplies a precomputed kernel over ALL coordinates,
        # index it by this gene's training mask instead of recomputing the
        # distance matrix -- the kernel depends only on coordinates, which
        # are shared across the whole batch.  NOTE: when median_nn_dist was
        # supplied the kernel MUST be evaluated at the effective length scale,
        # so a caller-provided k_nm_full is only valid if it was built at the
        # same eff_length (the batch wrapper handles this by precomputing once
        # at the shared eff_length).
        if k_nm_full is not None:
            train_idx = np.flatnonzero(mask)
            K_nm = np.asarray(k_nm_full, dtype=float)[train_idx]
        else:
            K_nm = self._rbf_kernel(train_coords, self.inducing_points_,
                                    length_scale=eff_length)

        # Compute K_mm and K_mm^{-1} (or reuse a shared precomputed inverse).
        # K_mm is still needed below for Sigma = K_mm + K_mn @ Lambda^{-1} @ K_nm.
        # Both paths evaluate K_mm at the effective length scale so it stays
        # consistent with K_nm (when the batch wrapper supplies a precomputed
        # K_mm_inv, it was built at the SAME shared eff_length).
        if K_mm_inv is not None:
            self._K_mm_inv = K_mm_inv
            K_mm = self._rbf_kernel(self.inducing_points_, self.inducing_points_,
                                    length_scale=eff_length)
            K_mm += 1e-6 * np.eye(m)  # Jitter, consistent with the precompute path
        else:
            K_mm = self._rbf_kernel(self.inducing_points_, self.inducing_points_,
                                    length_scale=eff_length)
            K_mm += 1e-6 * np.eye(m)  # Jitter for numerical stability
            self._K_mm_inv = np.linalg.inv(K_mm)

        # Local alias: the math below uses K_mm_inv in two places.
        K_mm_inv = self._K_mm_inv

        # Compute Lambda (diagonal correction term)
        # Lambda = diag(K_nn - Q_nn) + noise  (eff_noise = per-gene local noise
        # when local_noise=True, else the flat noise_level).
        K_nn_diag = np.ones(n)  # RBF kernel diagonal is 1
        Q_nn_diag = np.sum(K_nm @ K_mm_inv * K_nm, axis=1)
        Lambda = K_nn_diag - Q_nn_diag + eff_noise
        
        # Compute Sigma = K_mm + K_mn @ Lambda^{-1} @ K_nm
        K_mn = K_nm.T
        Sigma = K_mm + (K_mn * (1.0 / Lambda)[None, :]) @ K_nm
        
        # Cholesky decomposition for efficient solving
        L = cho_factor(Sigma, lower=True)
        
        # Compute alpha = Sigma^{-1} @ K_mn @ Lambda^{-1} @ y
        temp = K_mn @ (train_values / Lambda)
        self.alpha_ = cho_solve(L, temp)
        
        # Store training data
        self._train_coords = train_coords
        self._train_values = train_values
        
        logger.info(
            f"Fitted sparse GP with {n} training points, "
            f"{m} inducing points (compression ratio: {n/m:.1f}x)"
        )
        
        return self
    
    def predict(
        self,
        coordinates: np.ndarray,
        return_std: bool = False
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        """
        Predict at new coordinates.
        
        Parameters
        ----------
        coordinates : np.ndarray, shape (n_test, 2)
            Test coordinates
        return_std : bool, default=False
            Whether to return standard deviation
        
        Returns
        -------
        predictions : np.ndarray, shape (n_test,)
            Predicted values
        std : np.ndarray, optional
            Standard deviations
        """
        if self.inducing_points_ is None or self.alpha_ is None:
            raise ValueError("Model not fitted. Call fit() first.")

        # Use the length scale that was actually in effect at fit time
        # (may differ from self.length_scale when median_nn_dist was supplied).
        eff_length = (self._effective_length_scale
                      if self._effective_length_scale is not None
                      else self.length_scale)
        eff_noise = (self._effective_noise
                     if self._effective_noise is not None
                     else self.noise_level)

        # K_*m: kernel between test points and inducing points
        K_star_m = self._rbf_kernel(coordinates, self.inducing_points_,
                                    length_scale=eff_length)

        # Predictions: K_*m @ alpha
        predictions = K_star_m @ self.alpha_

        if return_std:
            # Compute predictive variance
            # Var = K_** - K_*m @ K_mm^{-1} @ K_m*
            K_star_star_diag = np.ones(len(coordinates))  # RBF diagonal
            Q_star_star_diag = np.sum(K_star_m @ self._K_mm_inv * K_star_m, axis=1)
            variance = K_star_star_diag - Q_star_star_diag + eff_noise
            std = np.sqrt(np.maximum(variance, 0))  # Ensure non-negative

            return predictions, std

        return predictions, None
    
    def impute(
        self,
        coordinates: np.ndarray,
        values: np.ndarray,
        mask: Optional[np.ndarray] = None,
        return_uncertainty: bool = False
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        """
        Impute missing values.
        
        Parameters
        ----------
        coordinates : np.ndarray
            Spatial coordinates
        values : np.ndarray
            Observed values
        mask : np.ndarray, optional
            Training mask
        return_uncertainty : bool, default=False
            Whether to return uncertainty
        
        Returns
        -------
        imputed_values : np.ndarray
            Imputed values
        uncertainty : np.ndarray, optional
            Uncertainty estimates
        """
        # Fit on observed data
        self.fit(coordinates, values, mask)
        
        # Predict at all coordinates
        imputed, uncertainty = self.predict(
            coordinates,
            return_std=return_uncertainty
        )
        
        # Keep observed values
        if mask is None:
            mask = values > 0
        imputed[mask] = values[mask]
        
        if uncertainty is not None:
            uncertainty[mask] = 0.0
        
        return imputed, uncertainty

    def fit_batch(
        self,
        coordinates: np.ndarray,
        values: np.ndarray,
        mask: Optional[np.ndarray] = None,
        n_jobs: int = 1,
        verbose: bool = True,
    ) -> 'SparseGPImputerBatch':
        """
        Fit sparse GP models for multiple genes.

        Parameters
        ----------
        coordinates : np.ndarray, shape (n_spots, 2)
            Spatial coordinates.
        values : np.ndarray, shape (n_genes, n_spots)
            APA values for multiple genes.
        mask : np.ndarray, optional, shape (n_genes, n_spots)
            Boolean training mask.
        n_jobs : int, default=1
            Reserved for API compatibility. Sparse batch currently runs
            sequentially to keep memory use predictable.
        verbose : bool, default=True
            Print progress every 100 genes.

        Returns
        -------
        SparseGPImputerBatch
            Fitted batch imputer.
        """
        return SparseGPImputerBatch(
            base_imputer=self,
            coordinates=coordinates,
            values=values,
            mask=mask,
            verbose=verbose,
        )


class SparseGPImputerBatch:
    """Batch wrapper with the same public shape contract as GPImputerBatch."""

    def __init__(
        self,
        base_imputer: SparseGPImputer,
        coordinates: np.ndarray,
        values: np.ndarray,
        mask: Optional[np.ndarray] = None,
        verbose: bool = True,
    ):
        if values.ndim != 2:
            raise ValueError("values must have shape (n_genes, n_spots)")
        if values.shape[1] != len(coordinates):
            raise ValueError(
                "values must have shape (n_genes, n_spots) with n_spots "
                "matching coordinates"
            )
        if mask is not None and mask.shape != values.shape:
            raise ValueError("mask must have the same shape as values")

        self.base_imputer = base_imputer
        self.coordinates = coordinates
        self.values = values
        self.mask = mask
        self.verbose = verbose
        self.imputers_ = []
        self._fit_batch()

    def _new_imputer(self) -> SparseGPImputer:
        return SparseGPImputer(
            n_inducing=self.base_imputer.n_inducing,
            inducing_method=self.base_imputer.inducing_method,
            length_scale=self.base_imputer.length_scale,
            noise_level=self.base_imputer.noise_level,
            kernel_type=self.base_imputer.kernel_type,
            length_scale_multiplier=self.base_imputer.length_scale_multiplier,
            local_noise=self.base_imputer.local_noise,
            noise_scale=self.base_imputer.noise_scale,
        )

    def _fit_batch(self):
        n_genes = self.values.shape[0]

        # Resolve the SHARED effective length scale once.  The multiplier
        # path is opt-in: it activates only when length_scale_multiplier is
        # set to a float.  In that case the median NN distance of the full
        # coordinate set defines the spatial scale and the kernel uses
        # median_nn_dist * length_scale_multiplier for every gene; this keeps
        # the precomputed K_mm and k_nm_full valid for all genes (the kernel
        # depends only on coordinates, which are shared).  When the multiplier
        # is None, the explicit length_scale is used (backward-compatible).
        base = self.base_imputer
        use_multiplier = base.length_scale_multiplier is not None
        if use_multiplier:
            try:
                from scipy.spatial import cKDTree
                nn = cKDTree(self.coordinates).query(self.coordinates, k=2)[0][:, 1]
                median_nn = float(np.median(nn))
            except Exception:
                median_nn = None
            if median_nn is not None and median_nn > 0:
                shared_eff_length = median_nn * float(base.length_scale_multiplier)
            else:
                shared_eff_length = float(base.length_scale)
        else:
            median_nn = None
            shared_eff_length = float(base.length_scale)
        self._shared_eff_length = shared_eff_length
        self._shared_median_nn = median_nn

        # Precompute inducing points ONCE on all coordinates (not the
        # gene-specific training mask). Inducing points depend only on the
        # spatial layout, which is identical for every gene, so running
        # KMeans(42k points, ~424 clusters, n_init=10) per gene is pure waste.
        # Reusing the same set + precomputed K_mm^{-1} turns an O(n_genes)
        # KMeans cost into O(1) and is the bulk of the speedup.
        shared_inducing = base._select_inducing_points(self.coordinates)
        K_mm = base._rbf_kernel(shared_inducing, shared_inducing,
                                length_scale=shared_eff_length)
        K_mm += 1e-6 * np.eye(len(shared_inducing))
        shared_K_mm_inv = np.linalg.inv(K_mm)
        # Precompute the kernel from ALL coordinates to the shared inducing
        # points once. K_nm depends only on coordinates (not gene values),
        # so each gene's training kernel is just ``k_nm_full[mask]``. This
        # eliminates one O(n*m) distance_matrix evaluation per gene.
        k_nm_full = base._rbf_kernel(self.coordinates, shared_inducing,
                                     length_scale=shared_eff_length)
        logger.info(
            f"Precomputed {len(shared_inducing)} shared inducing points "
            f"and K_nm ({k_nm_full.shape}) for {n_genes} genes "
            f"(eff_length_scale={shared_eff_length:.3f})"
        )

        for gene_idx in range(n_genes):
            if self.verbose and gene_idx > 0 and gene_idx % 100 == 0:
                logger.info(f"Fitted {gene_idx}/{n_genes} sparse GP models")

            imputer = self._new_imputer()
            gene_values = self.values[gene_idx, :]
            gene_mask = self.mask[gene_idx, :] if self.mask is not None else None
            try:
                # Pass median_nn_dist so each imputer records the same
                # _effective_length_scale (== shared_eff_length) used to build
                # k_nm_full / K_mm_inv; this keeps predict() consistent.
                imputer.fit(
                    self.coordinates, gene_values, gene_mask,
                    inducing_points=shared_inducing,
                    K_mm_inv=shared_K_mm_inv,
                    k_nm_full=k_nm_full,
                    median_nn_dist=median_nn,
                )
                self.imputers_.append(imputer)
            except Exception as exc:
                logger.warning(f"Failed to fit sparse GP for gene {gene_idx}: {exc}")
                self.imputers_.append(None)

        n_success = sum(1 for imp in self.imputers_ if imp is not None)
        logger.info(f"Fitted {n_success}/{n_genes} sparse GP models successfully")

    def predict(
        self,
        coordinates: Optional[np.ndarray] = None,
        return_std: bool = True,
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        if coordinates is None:
            coordinates = self.coordinates

        n_genes = len(self.imputers_)
        n_spots = len(coordinates)
        predictions = np.zeros((n_genes, n_spots))
        uncertainty = np.zeros((n_genes, n_spots)) if return_std else None

        for gene_idx, imputer in enumerate(self.imputers_):
            if imputer is None:
                continue
            pred, unc = imputer.predict(coordinates, return_std=return_std)
            predictions[gene_idx, :] = pred
            if uncertainty is not None:
                uncertainty[gene_idx, :] = unc

        return predictions, uncertainty

    def impute(
        self,
        return_uncertainty: bool = True,
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        predictions, uncertainty = self.predict(
            self.coordinates,
            return_std=return_uncertainty,
        )

        if self.mask is not None:
            observed = self.mask
        else:
            observed = self.values > 0

        predictions[observed] = self.values[observed]
        if uncertainty is not None:
            uncertainty[observed] = 0.0

        return predictions, uncertainty


class BlockGPImputer:
    """
    Block-wise GP imputation for very large datasets.
    
    Divides spatial domain into blocks and fits separate GP models
    for each block, with overlap for smooth transitions.
    
    Parameters
    ----------
    block_size : int, default=500
        Maximum number of points per block
    overlap : float, default=0.1
        Overlap fraction between blocks (0-1)
    base_imputer : str, default='sparse'
        Base imputer type: 'sparse' or 'full'
    **imputer_kwargs
        Additional arguments for base imputer
    
    Examples
    --------
    >>> imputer = BlockGPImputer(block_size=500, overlap=0.1)
    >>> imputed = imputer.impute(spatial_coords, apa_counts)
    """
    
    def __init__(
        self,
        block_size: int = 500,
        overlap: float = 0.1,
        base_imputer: str = 'sparse',
        **imputer_kwargs
    ):
        self.block_size = block_size
        self.overlap = overlap
        self.base_imputer = base_imputer
        self.imputer_kwargs = imputer_kwargs
        
        self.blocks_ = []
        self.imputers_ = []
    
    def _create_blocks(
        self,
        coordinates: np.ndarray
    ) -> list:
        """
        Divide coordinates into overlapping blocks.
        
        Parameters
        ----------
        coordinates : np.ndarray
            Spatial coordinates
        
        Returns
        -------
        list
            List of block indices
        """
        n_points = len(coordinates)
        
        if n_points <= self.block_size:
            # No need for blocking
            return [np.arange(n_points)]
        
        # Compute spatial extent
        x_min, x_max = coordinates[:, 0].min(), coordinates[:, 0].max()
        y_min, y_max = coordinates[:, 1].min(), coordinates[:, 1].max()
        
        # Estimate number of blocks needed
        area = (x_max - x_min) * (y_max - y_min)
        points_per_area = n_points / area
        block_area = self.block_size / points_per_area
        block_width = np.sqrt(block_area)
        
        # Create grid of blocks with overlap
        overlap_width = block_width * self.overlap
        
        x_edges = np.arange(x_min, x_max + block_width, block_width - overlap_width)
        y_edges = np.arange(y_min, y_max + block_width, block_width - overlap_width)
        
        blocks = []
        for i in range(len(x_edges) - 1):
            for j in range(len(y_edges) - 1):
                # Find points in this block
                mask = (
                    (coordinates[:, 0] >= x_edges[i]) &
                    (coordinates[:, 0] < x_edges[i+1]) &
                    (coordinates[:, 1] >= y_edges[j]) &
                    (coordinates[:, 1] < y_edges[j+1])
                )
                indices = np.where(mask)[0]
                
                if len(indices) > 0:
                    blocks.append(indices)
        
        logger.info(f"Created {len(blocks)} blocks for {n_points} points")
        
        return blocks
    
    def impute(
        self,
        coordinates: np.ndarray,
        values: np.ndarray,
        mask: Optional[np.ndarray] = None,
        return_uncertainty: bool = False
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        """
        Impute using block-wise processing.
        
        Parameters
        ----------
        coordinates : np.ndarray
            Spatial coordinates
        values : np.ndarray
            Observed values
        mask : np.ndarray, optional
            Training mask
        return_uncertainty : bool, default=False
            Whether to return uncertainty
        
        Returns
        -------
        imputed_values : np.ndarray
            Imputed values
        uncertainty : np.ndarray, optional
            Uncertainty estimates
        """
        # Create blocks
        self.blocks_ = self._create_blocks(coordinates)
        
        # Initialize output
        imputed = np.zeros_like(values)
        uncertainty = np.zeros_like(values) if return_uncertainty else None
        weights = np.zeros_like(values)
        
        # Process each block
        for block_idx, indices in enumerate(self.blocks_):
            block_coords = coordinates[indices]
            block_values = values[indices]
            block_mask = mask[indices] if mask is not None else None
            
            # Create imputer for this block
            if self.base_imputer == 'sparse':
                imputer = SparseGPImputer(**self.imputer_kwargs)
            else:
                from .gp_imputer import GPImputer
                imputer = GPImputer(**self.imputer_kwargs)
            
            # Impute block
            try:
                block_imputed, block_unc = imputer.impute(
                    block_coords,
                    block_values,
                    block_mask,
                    return_uncertainty=return_uncertainty
                )
                
                # Accumulate results (weighted average in overlap regions)
                imputed[indices] += block_imputed
                weights[indices] += 1
                
                if uncertainty is not None and block_unc is not None:
                    uncertainty[indices] += block_unc
            
            except Exception as e:
                logger.warning(f"Failed to impute block {block_idx}: {e}")
                # Use original values for failed blocks
                imputed[indices] += block_values
                weights[indices] += 1
        
        # Average overlapping regions
        weights[weights == 0] = 1  # Avoid division by zero
        imputed /= weights
        
        if uncertainty is not None:
            uncertainty /= weights
        
        logger.info(f"Completed block-wise imputation with {len(self.blocks_)} blocks")
        
        return imputed, uncertainty
