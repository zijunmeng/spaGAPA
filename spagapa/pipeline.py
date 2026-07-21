"""
Complete spaGAPA Pipeline with uncertainty-weighted analysis.

Integrates all modules: spatial validation, GP imputation (with uncertainty),
APA quantification, domain identification, differential APA, and SVAPA detection.
All downstream analyses support uncertainty weighting when available.
"""

import numpy as np
import pandas as pd
import json
from typing import Optional, Dict, List, Union, Tuple
from pathlib import Path
import warnings

from spagapa.core import APADataset
from spagapa.io import load_spatial_dataset
from spagapa.spatial import SpatialNeighbors
from spagapa.calling import SpatialValidator, QualityFilter
from spagapa.imputation import GPImputer, SparseGPImputer, ExpressionFeatureBuilder
from spagapa.bioml import (
    BioMLDomainDetector,
    GraphRegularizedAPAFactorizer,
    HighResBioMLConfig,
    MultiViewGraphBuilder,
    highres_bioml_recover,
)
from spagapa.presets import (
    VALID_ANALYSIS_PRESETS,
    profile_spatial_apa_matrix,
    resolve_analysis_preset,
)
from spagapa.quantification import APAIndexCalculator, QCReportGenerator
from spagapa.analysis import (
    DomainIdentifier,
    DifferentialAPAAnalyzer,
    SpatialPatternAnalyzer,
    GPTrendDetector,
)
from sklearn.neighbors import NearestNeighbors


class SpaGAPA:
    """
    Main spaGAPA pipeline with uncertainty-weighted analysis.

    Workflow (7 steps):
    1. Load data → APADataset
    2. Spatial validation + quality filtering
    3. GP imputation → imputed values + uncertainty
    4. APA quantification (RUD, PDUI, WUL)
    5. Spatial domain identification (uncertainty-weighted)
    6. Differential APA analysis (uncertainty-weighted)
    7. SVAPA gene detection (GP likelihood ratio + weighted Moran's I)

    All downstream analyses automatically use uncertainty weights when
    imputation uncertainty is available.

    Parameters
    ----------
    n_neighbors : int, default=6
        Spatial neighbors for graph construction.
    kernel_type : str, default='matern'
        GP kernel: 'rbf', 'matern', or 'auto'.
        gp_alpha : float, default=1e-10
        Exact GP noise/regularization parameter.
    gp_n_restarts_optimizer : int, default=1
        Number of GP hyperparameter optimizer restarts.
    gp_n_jobs : int, default=1
        Number of CPU workers for batch GP fitting.
    use_sparse_gp : bool, default=False
        Use sparse GP approximation for large datasets.
    n_inducing : int, default=100
        Inducing points for sparse GP.
    min_spatial_support : float, default=0.3
        Minimum spatial support threshold.
    min_read_count : int, default=10
        Minimum read count per site.
    min_spots : int, default=5
        Minimum spots for quality filtering.
    analysis_preset : {'auto', 'standard', 'highres_accuracy', 'highres_fast'}, default='auto'
        User-facing analysis mode. ``auto`` resolves from data shape.
    use_bioml : bool, optional
        Override preset-driven BioML multi-view graph domain recovery.
    bioml_rank : int, default=8
        Low-rank dimension for BioML APA factorization.
    bioml_domain_method : {'spectral', 'kmeans'}, default='spectral'
        Domain detector for BioML outputs.
    bioml_spatial_weight, bioml_expression_weight, bioml_apa_weight : float
        Multi-view graph fusion weights.
    verbose : bool, default=True
        Print progress.
    """

    def __init__(
        self,
        n_neighbors: int = 6,
        kernel_type: str = 'matern',
        gp_alpha: float = 1e-10,
        gp_n_restarts_optimizer: int = 1,
        gp_n_jobs: int = 1,
        use_sparse_gp: bool = False,
        n_inducing: int = 100,
        sparse_gp_inducing_method: str = 'kmeans',
        sparse_gp_length_scale: Union[str, float] = 1.0,
        sparse_gp_length_scale_multiplier: float = 1.0,
        sparse_gp_noise_level: float = 0.1,
        input_type: str = 'apa_index',
        min_spatial_support: float = 0.3,
        min_read_count: int = 10,
        min_spots: int = 5,
        analysis_preset: str = 'auto',
        use_bioml: Optional[bool] = None,
        bioml_rank: int = 8,
        bioml_lambda_graph: float = 0.5,
        bioml_lambda_l2: float = 1e-2,
        bioml_max_iter: int = 20,
        bioml_n_neighbors: int = 15,
        bioml_blend: float = 0.1,
        bioml_domains_only: bool = False,
        bioml_domain_method: str = 'spectral',
        bioml_spatial_weight: float = 0.4,
        bioml_expression_weight: float = 0.4,
        bioml_apa_weight: float = 0.2,
        highres_bioml_gp_blend: float = 0.1,
        highres_bioml_spatial_weight: Optional[float] = None,
        highres_bioml_expression_weight: Optional[float] = None,
        highres_bioml_apa_weight: Optional[float] = None,
        highres_bioml_apa_source: str = 'expression_knn',
        highres_bioml_expression_knn_k: int = 15,
        highres_bioml_neighbor_mode: str = 'adaptive',
        highres_bioml_adaptive_neighbor_scale: float = 10.0,
        highres_bioml_parent_weight: float = 0.0,
        highres_bioml_parent_neighbors: int = 8,
        expression_n_components: int = 10,
        verbose: bool = True,
    ):
        self.n_neighbors = n_neighbors
        self.kernel_type = kernel_type
        self.gp_alpha = float(gp_alpha)
        self.gp_n_restarts_optimizer = int(gp_n_restarts_optimizer)
        self.gp_n_jobs = int(gp_n_jobs)
        self.use_sparse_gp = use_sparse_gp
        self.n_inducing = n_inducing
        self.sparse_gp_inducing_method = sparse_gp_inducing_method
        self.sparse_gp_length_scale = sparse_gp_length_scale
        self.sparse_gp_length_scale_multiplier = float(sparse_gp_length_scale_multiplier)
        self.sparse_gp_noise_level = float(sparse_gp_noise_level)
        self.input_type = input_type
        self.min_spatial_support = min_spatial_support
        self.min_read_count = min_read_count
        self.min_spots = min_spots
        if analysis_preset not in VALID_ANALYSIS_PRESETS:
            valid = ", ".join(VALID_ANALYSIS_PRESETS)
            raise ValueError(f"analysis_preset must be one of: {valid}")
        self.analysis_preset = analysis_preset
        self.use_bioml = None if use_bioml is None else bool(use_bioml)
        self.bioml_rank = int(bioml_rank)
        self.bioml_lambda_graph = float(bioml_lambda_graph)
        self.bioml_lambda_l2 = float(bioml_lambda_l2)
        self.bioml_max_iter = int(bioml_max_iter)
        self.bioml_n_neighbors = int(bioml_n_neighbors)
        self.bioml_blend = float(bioml_blend)
        self.bioml_domains_only = bool(bioml_domains_only)
        self.bioml_domain_method = bioml_domain_method
        self.bioml_spatial_weight = float(bioml_spatial_weight)
        self.bioml_expression_weight = float(bioml_expression_weight)
        self.bioml_apa_weight = float(bioml_apa_weight)
        self.highres_bioml_gp_blend = float(highres_bioml_gp_blend)
        self.highres_bioml_spatial_weight = (
            None if highres_bioml_spatial_weight is None else float(highres_bioml_spatial_weight)
        )
        self.highres_bioml_expression_weight = (
            None if highres_bioml_expression_weight is None else float(highres_bioml_expression_weight)
        )
        self.highres_bioml_apa_weight = (
            None if highres_bioml_apa_weight is None else float(highres_bioml_apa_weight)
        )
        self.highres_bioml_apa_source = highres_bioml_apa_source
        self.highres_bioml_expression_knn_k = int(highres_bioml_expression_knn_k)
        self.highres_bioml_neighbor_mode = highres_bioml_neighbor_mode
        self.highres_bioml_adaptive_neighbor_scale = float(highres_bioml_adaptive_neighbor_scale)
        self.highres_bioml_parent_weight = float(highres_bioml_parent_weight)
        self.highres_bioml_parent_neighbors = int(highres_bioml_parent_neighbors)
        self.expression_n_components = int(expression_n_components)
        self.verbose = verbose

        self.dataset_: Optional[APADataset] = None
        self.imputer_ = None
        self.results_: Dict = {}
        self.resolved_analysis_preset_: Optional[str] = None
        self.dataset_profile_: Optional[Dict] = None
        self._active_bioml_params: Dict = {}
        self._active_sparse_gp_params: Dict = {}

    def _log(self, msg: str):
        if self.verbose:
            print(msg)

    def _validate_dataset_shapes(self):
        """Validate the project-wide genes x spots algorithm contract."""
        if self.dataset_ is None:
            raise ValueError("No dataset loaded")

        counts = self.dataset_.raw_counts
        coords = self.dataset_.coords

        if counts.ndim != 2:
            raise ValueError("raw_counts must have shape (n_genes, n_spots)")
        if coords.ndim != 2 or coords.shape[1] < 2:
            raise ValueError("coords must have shape (n_spots, 2)")
        expected = (self.dataset_.n_genes, self.dataset_.n_spots)
        if counts.shape != expected:
            raise ValueError(f"raw_counts shape mismatch: expected {expected}, got {counts.shape}")
        if coords.shape[0] != self.dataset_.n_spots:
            raise ValueError(
                "Coordinate count does not match spots: "
                f"{coords.shape[0]} coords vs {self.dataset_.n_spots} spots"
            )

    def _active_bioml_defaults(self, resolved_preset: str) -> Dict[str, float]:
        """Resolve run-local BioML defaults without mutating the estimator."""
        params = {
            'spatial_weight': self.bioml_spatial_weight,
            'expression_weight': self.bioml_expression_weight,
            'apa_weight': self.bioml_apa_weight,
            'blend': self.bioml_blend,
            'highres_gp_blend': self.highres_bioml_gp_blend,
        }
        if resolved_preset not in {'highres_accuracy', 'highres_fast'}:
            return params

        default_weights = (
            abs(self.bioml_spatial_weight - 0.4) < 1e-12
            and abs(self.bioml_expression_weight - 0.4) < 1e-12
            and abs(self.bioml_apa_weight - 0.2) < 1e-12
        )
        if default_weights:
            params['spatial_weight'] = 0.1
            params['expression_weight'] = 0.7
            params['apa_weight'] = 0.2
        if self.highres_bioml_spatial_weight is not None:
            params['spatial_weight'] = self.highres_bioml_spatial_weight
        if self.highres_bioml_expression_weight is not None:
            params['expression_weight'] = self.highres_bioml_expression_weight
        if self.highres_bioml_apa_weight is not None:
            params['apa_weight'] = self.highres_bioml_apa_weight

        if resolved_preset == 'highres_fast' and abs(self.bioml_blend - 0.1) < 1e-12:
            params['blend'] = 0.0
        if resolved_preset == 'highres_fast' and abs(self.highres_bioml_gp_blend - 0.1) < 1e-12:
            params['highres_gp_blend'] = 0.0
        return params

    def _active_sparse_gp_defaults(self, resolved_preset: str) -> Dict:
        """Resolve run-local sparse-GP defaults without mutating the estimator.

        For ``highres_accuracy`` the number of inducing points is scaled with
        the spot count of the loaded dataset::

            n_inducing = min(500, max(100, n_spots // 100))

        so the sparse approximation stays well-conditioned across the full
        range of high-resolution tissue sizes (sub-slide to whole-slide),
        while always remaining within the [100, 500] budget supported by the
        ``SparseGPImputer``. For ``standard`` / ``highres_fast`` the sparse GP
        is not used, so the user-supplied ``self.n_inducing`` is preserved.
        """
        params = {
            'n_inducing': int(self.n_inducing),
            'inducing_method': self.sparse_gp_inducing_method,
            'length_scale': self.sparse_gp_length_scale,
            'length_scale_multiplier': self.sparse_gp_length_scale_multiplier,
            'noise_level': self.sparse_gp_noise_level,
        }
        if resolved_preset == 'highres_accuracy':
            n_spots = int(self.dataset_.n_spots) if self.dataset_ is not None else 0
            params['n_inducing'] = min(500, max(100, n_spots // 100))
            if self.sparse_gp_length_scale == 1.0:
                params['length_scale'] = 'auto'
            if abs(self.sparse_gp_noise_level - 0.1) < 1e-12:
                params['noise_level'] = 0.08
        return params

    def _resolve_run_options(
        self,
        impute: bool,
        use_bioml: Optional[bool],
    ) -> Dict:
        """Resolve preset-aware algorithm switches for the current dataset."""
        if self.dataset_ is None:
            raise ValueError("No dataset loaded")

        profile = profile_spatial_apa_matrix(
            self.dataset_.raw_counts,
            input_type=self.input_type,
        )
        resolved_preset = resolve_analysis_preset(self.analysis_preset, profile)
        self._active_bioml_params = self._active_bioml_defaults(resolved_preset)
        self._active_sparse_gp_params = self._active_sparse_gp_defaults(resolved_preset)

        use_bioml_this_run = self.use_bioml if use_bioml is None else bool(use_bioml)
        if use_bioml_this_run is None:
            use_bioml_this_run = resolved_preset in {'highres_accuracy', 'highres_fast'}

        impute_this_run = bool(impute)
        use_sparse_gp_this_run = bool(self.use_sparse_gp)
        if resolved_preset == 'highres_accuracy':
            use_sparse_gp_this_run = True
        elif resolved_preset == 'highres_fast':
            impute_this_run = False
            use_sparse_gp_this_run = False

        metadata = {
            'requested_preset': self.analysis_preset,
            'resolved_preset': resolved_preset,
            'profile': profile.to_dict(),
            'impute': bool(impute_this_run),
            'use_sparse_gp': bool(use_sparse_gp_this_run),
            'sparse_gp': {
                'n_inducing': int(self._active_sparse_gp_params['n_inducing']),
                'inducing_method': self._active_sparse_gp_params['inducing_method'],
                'length_scale': self._active_sparse_gp_params['length_scale'],
                'length_scale_multiplier': self._active_sparse_gp_params['length_scale_multiplier'],
                'noise_level': self._active_sparse_gp_params['noise_level'],
            },
            'use_bioml': bool(use_bioml_this_run),
            'bioml_weights': {
                'spatial': self._active_bioml_params['spatial_weight'],
                'expression': self._active_bioml_params['expression_weight'],
                'apa': self._active_bioml_params['apa_weight'],
            },
            'bioml_blend': self._active_bioml_params['blend'],
            'highres_bioml_gp_blend': self._active_bioml_params['highres_gp_blend'],
            'highres_bioml_apa_source': self.highres_bioml_apa_source,
            'notes': [],
        }
        if self.analysis_preset == 'auto':
            metadata['notes'].append(
                "auto resolved to highres_accuracy"
                if resolved_preset == 'highres_accuracy'
                else "auto resolved to standard"
            )
        if resolved_preset == 'highres_fast':
            metadata['notes'].append("highres_fast skips GP imputation and uncertainty")

        self.resolved_analysis_preset_ = resolved_preset
        self.dataset_profile_ = profile.to_dict()
        return metadata

    def _clip_apa_values(self, values: np.ndarray) -> np.ndarray:
        """Clip matrices according to the declared input scale."""
        if self.input_type == 'apa_index':
            return np.clip(values, 0.0, 1.0)
        return np.clip(values, 0.0, None)

    def _training_mask(self, values: np.ndarray) -> np.ndarray:
        """Return observed-value mask under the declared input semantics."""
        finite = np.isfinite(values)
        if self.input_type == 'apa_index':
            # APA index matrices encode missing cells as 0.0 (not NaN).
            # Without the >0 filter, observed_fraction=1.0 and GP uncertainty
            # is zeroed everywhere. Match qc_metrics.py's values>0 convention.
            return finite & (values > 0)
        return finite & (values > 0)

    def _resolve_sparse_gp_length_scale(
        self,
        coords: np.ndarray,
        value: Optional[Union[str, float]] = None,
        multiplier: Optional[float] = None,
    ) -> float:
        """Resolve sparse-GP length scale, optionally from local coordinate density."""
        if value is None:
            value = self.sparse_gp_length_scale
        if multiplier is None:
            multiplier = self.sparse_gp_length_scale_multiplier
        if isinstance(value, str):
            normalized = value.strip().lower()
            if normalized not in {'auto', 'none'}:
                raise ValueError("sparse_gp_length_scale must be positive, 'auto', or 'none'")
            if coords.shape[0] < 2:
                return 1.0
            k = min(max(2, int(self.bioml_n_neighbors)), coords.shape[0])
            nn = NearestNeighbors(n_neighbors=k)
            nn.fit(coords)
            distances, _ = nn.kneighbors(coords)
            positive = distances[:, -1]
            positive = positive[np.isfinite(positive) & (positive > 0)]
            if positive.size:
                base = float(np.median(positive))
            else:
                nn2 = NearestNeighbors(n_neighbors=2)
                nn2.fit(coords)
                nearest, _ = nn2.kneighbors(coords)
                vals = nearest[:, 1]
                vals = vals[np.isfinite(vals) & (vals > 0)]
                base = float(np.median(vals)) if vals.size else 1.0
            return max(1e-6, base * float(multiplier))

        length_scale = float(value)
        if length_scale <= 0:
            raise ValueError("sparse_gp_length_scale must be positive, 'auto', or 'none'")
        return length_scale

    def _load_expression_matrix(
        self,
        expression_matrix: Optional[Union[str, np.ndarray, pd.DataFrame]],
    ) -> Optional[Union[np.ndarray, pd.DataFrame]]:
        """Load an optional expression matrix while preserving table labels."""
        if expression_matrix is None:
            return None
        if isinstance(expression_matrix, (np.ndarray, pd.DataFrame)):
            return expression_matrix

        path = Path(expression_matrix)
        if path.suffix.lower() == '.npy':
            return np.load(path)

        sep = '\t' if path.suffix.lower() in {'.tsv', '.txt'} else ','
        return pd.read_csv(path, sep=sep, index_col=0)

    def _build_expression_embedding(
        self,
        expression_matrix: Optional[Union[str, np.ndarray, pd.DataFrame]],
        expression_embedding: Optional[np.ndarray],
        expression_orientation: str,
    ) -> Optional[np.ndarray]:
        """Return a spot-level expression embedding for BioML."""
        if expression_embedding is not None:
            embedding = np.asarray(expression_embedding, dtype=float)
            if embedding.ndim != 2:
                raise ValueError("expression_embedding must be 2-dimensional")
            if embedding.shape[0] != self.dataset_.n_spots:
                raise ValueError("expression_embedding must have n_spots rows")
            return embedding

        loaded = self._load_expression_matrix(expression_matrix)
        if loaded is None:
            return None

        builder = ExpressionFeatureBuilder(
            n_components=self.expression_n_components,
            orientation=expression_orientation,
        )
        embedding = builder.fit_transform(loaded)
        if embedding.shape[0] != self.dataset_.n_spots:
            raise ValueError(
                "Expression matrix produced an embedding with "
                f"{embedding.shape[0]} spots, expected {self.dataset_.n_spots}. "
                "Check expression_orientation."
            )
        return embedding

    @staticmethod
    def _confidence_from_uncertainty(
        uncertainty: Optional[np.ndarray],
        mask: np.ndarray,
    ) -> Optional[np.ndarray]:
        """Convert GP uncertainty to bounded confidence weights."""
        if uncertainty is None:
            return None
        uncertainty = np.asarray(uncertainty, dtype=float)
        finite = uncertainty[np.isfinite(uncertainty)]
        fallback = float(np.median(finite)) if finite.size else 1.0
        unc = np.nan_to_num(uncertainty, nan=fallback, posinf=fallback, neginf=fallback)
        confidence = 1.0 / (unc + 1e-6)
        finite_conf = confidence[np.isfinite(confidence) & (confidence > 0)]
        scale = float(np.median(finite_conf)) if finite_conf.size else 1.0
        confidence = confidence / max(scale, 1e-6)
        confidence = np.clip(confidence, 0.0, 10.0)
        confidence[~mask] = 0.0
        return confidence

    def _run_bioml(
        self,
        work: np.ndarray,
        coords: np.ndarray,
        uncertainty: Optional[np.ndarray],
        expression_embedding: Optional[np.ndarray],
        n_domains: int,
    ) -> Tuple[np.ndarray, Dict]:
        """Run CPU-friendly BioML refinement and domain recovery."""
        mask = self._training_mask(self.dataset_.raw_counts)
        confidence = self._confidence_from_uncertainty(uncertainty, mask)
        active_params = self._active_bioml_params or {
            'spatial_weight': self.bioml_spatial_weight,
            'expression_weight': self.bioml_expression_weight,
            'apa_weight': self.bioml_apa_weight,
            'blend': self.bioml_blend,
        }

        graph = MultiViewGraphBuilder(
            n_neighbors=self.bioml_n_neighbors,
            spatial_weight=active_params['spatial_weight'],
            expression_weight=active_params['expression_weight'],
            apa_weight=active_params['apa_weight'],
        ).build(
            coords,
            expression_embedding=expression_embedding,
            apa_matrix=work,
            uncertainty=uncertainty,
        )

        factorizer = None
        blend = 0.0 if self.bioml_domains_only else float(np.clip(active_params['blend'], 0.0, 1.0))
        refined = self._clip_apa_values(work)

        if not self.bioml_domains_only:
            factorizer = GraphRegularizedAPAFactorizer(
                rank=self.bioml_rank,
                lambda_graph=self.bioml_lambda_graph,
                lambda_l2=self.bioml_lambda_l2,
                max_iter=self.bioml_max_iter,
                random_state=42,
                preserve_observed=True,
            )
            bioml_imputed = factorizer.fit_transform(
                self.dataset_.raw_counts,
                graph_laplacian=graph.laplacian(),
                mask=mask,
                confidence=confidence,
            )
            refined = (1.0 - blend) * work + blend * bioml_imputed
            refined[mask] = self.dataset_.raw_counts[mask]
            refined = self._clip_apa_values(refined)

        if self.bioml_domain_method == 'spectral':
            labels = BioMLDomainDetector(
                method='spectral',
                n_domains=n_domains,
                random_state=42,
            ).fit_predict(graph=graph.fused)
        elif self.bioml_domain_method == 'kmeans':
            if factorizer is not None:
                spot_factors = factorizer.spot_factors_
            else:
                feature_blocks = [coords, refined.T]
                if expression_embedding is not None:
                    feature_blocks.insert(1, expression_embedding)
                spot_factors = np.column_stack(feature_blocks)
            labels = BioMLDomainDetector(
                method='kmeans',
                n_domains=n_domains,
                random_state=42,
            ).fit_predict(spot_factors=spot_factors)
        else:
            raise ValueError("bioml_domain_method must be 'spectral' or 'kmeans'")

        metadata = {
            'enabled': True,
            'mode': 'domains_only' if self.bioml_domains_only else 'value_refinement',
            'domain_method': self.bioml_domain_method,
            'n_domains': int(n_domains),
            'rank': self.bioml_rank,
            'lambda_graph': self.bioml_lambda_graph,
            'lambda_l2': self.bioml_lambda_l2,
            'max_iter': self.bioml_max_iter,
            'n_neighbors': self.bioml_n_neighbors,
            'blend': blend,
            'graph_weights_requested': {
                'spatial': active_params['spatial_weight'],
                'expression': active_params['expression_weight'],
                'apa': active_params['apa_weight'],
            },
            'graph_weights_effective': graph.weights,
            'has_expression_view': expression_embedding is not None,
            'factorization_n_iter': None if factorizer is None else int(factorizer.result_.n_iter),
            'factorization_reconstruction_error': (
                None if factorizer is None else float(factorizer.result_.reconstruction_error)
            ),
        }

        self.dataset_.set_bioml_results(
            imputed=refined,
            spot_factors=None if factorizer is None else factorizer.spot_factors_,
            gene_factors=None if factorizer is None else factorizer.gene_factors_,
            metadata=metadata,
        )

        method_name = 'spagapa_gp_bioml_domains' if self.bioml_domains_only else 'spagapa_bioml'
        return labels, {
            'labels': labels,
            'n_domains': n_domains,
            'method': method_name,
            'domain_method': self.bioml_domain_method,
            'imputed_values': refined,
            'spot_factors': None if factorizer is None else factorizer.spot_factors_,
            'gene_factors': None if factorizer is None else factorizer.gene_factors_,
            'metadata': metadata,
        }

    def _get_parent_index(self) -> Optional[np.ndarray]:
        """Return optional parent/coarse-bin labels for high-resolution graphs."""
        if self.dataset_ is None:
            return None
        obs = self.dataset_.adata.obs
        for key in ('parent_spot', 'parent_spot_id', 'parent_index', 'parent_bin'):
            if key in obs.columns:
                codes, _ = pd.factorize(obs[key], sort=True)
                return np.asarray(codes, dtype=int)
        apa_uns = self.dataset_.adata.uns.get('apa', {})
        parent = apa_uns.get('parent_index')
        if parent is not None:
            parent = np.asarray(parent, dtype=int)
            if parent.shape[0] == self.dataset_.n_spots:
                return parent
        return None

    def _run_highres_bioml(
        self,
        work: np.ndarray,
        coords: np.ndarray,
        uncertainty: Optional[np.ndarray],
        expression_embedding: Optional[np.ndarray],
        n_domains: int,
    ) -> Tuple[np.ndarray, Dict, np.ndarray]:
        """Run decoupled high-resolution BioML recovery."""
        active_params = self._active_bioml_params or {
            'spatial_weight': self.bioml_spatial_weight,
            'expression_weight': self.bioml_expression_weight,
            'apa_weight': self.bioml_apa_weight,
            'blend': self.bioml_blend,
            'highres_gp_blend': self.highres_bioml_gp_blend,
        }
        sparse_gp = work if self.dataset_.has_imputed() else None
        parent_index = self._get_parent_index()
        config = HighResBioMLConfig(
            gp_blend=active_params['highres_gp_blend'],
            spatial_weight=active_params['spatial_weight'],
            expression_weight=active_params['expression_weight'],
            apa_weight=active_params['apa_weight'],
            apa_source=self.highres_bioml_apa_source,
            expression_knn_k=self.highres_bioml_expression_knn_k,
            n_neighbors=self.bioml_n_neighbors,
            neighbor_mode=self.highres_bioml_neighbor_mode,
            adaptive_neighbor_scale=self.highres_bioml_adaptive_neighbor_scale,
            parent_weight=self.highres_bioml_parent_weight,
            parent_neighbors=self.highres_bioml_parent_neighbors,
            domain_method='spectral',
            random_state=42,
        )
        result = highres_bioml_recover(
            self.dataset_.raw_counts,
            coords,
            expression_embedding=expression_embedding,
            sparse_gp=sparse_gp,
            uncertainty=uncertainty,
            parent_index=parent_index,
            n_domains=n_domains,
            config=config,
        )

        metadata = {
            'enabled': True,
            'mode': 'highres_bioml',
            'domain_method': 'spectral',
            'n_domains': int(n_domains),
            'highres': result.metadata,
            'has_expression_view': expression_embedding is not None,
        }
        self.dataset_.set_bioml_results(
            imputed=result.recovered,
            spot_factors=None,
            gene_factors=None,
            metadata=metadata,
        )
        domain_result = {
            'labels': result.labels,
            'n_domains': n_domains,
            'method': 'spagapa_highres_bioml',
            'domain_method': 'spectral',
            'imputed_values': result.recovered,
            'metadata': metadata,
        }
        return result.labels, domain_result, result.recovered

    def _calculate_apa_indices(self, work: np.ndarray) -> Dict[str, Union[np.ndarray, bool, str]]:
        """
        Register or calculate APA indices for pipeline output.

        For ``input_type='apa_index'`` the matrix is already an APA index, so
        this step stores it explicitly instead of pretending to recompute RUD.
        For ``input_type='proximal_distal_counts'`` the dataset must contain
        AnnData layers ``proximal_counts`` and ``distal_counts`` in spots x
        genes layout.
        """
        if self.input_type == 'apa_index':
            values = self._clip_apa_values(work)
            self.dataset_.set_apa_index('APAIndex', values)
            return {
                'APAIndex': values,
                'calculated': False,
                'source': 'input_apa_index',
            }

        if self.input_type != 'proximal_distal_counts':
            raise ValueError("input_type must be 'apa_index' or 'proximal_distal_counts'")

        layers = self.dataset_.adata.layers
        if 'proximal_counts' not in layers or 'distal_counts' not in layers:
            raise ValueError(
                "input_type='proximal_distal_counts' requires AnnData layers "
                "'proximal_counts' and 'distal_counts' in spots x genes layout"
            )

        proximal = np.asarray(layers['proximal_counts']).T
        distal = np.asarray(layers['distal_counts']).T
        if proximal.shape != work.shape or distal.shape != work.shape:
            raise ValueError(
                "proximal_counts and distal_counts layers must match raw_counts "
                "shape (genes x spots) after transpose"
            )

        calculator = APAIndexCalculator()
        indices = calculator.calculate_all(proximal, distal)
        for name, values in indices.items():
            self.dataset_.set_apa_index(name, values)
        indices['calculated'] = True
        indices['source'] = 'proximal_distal_counts'
        return indices

    # ── main entry ───────────────────────────────────────────────────

    def run(
        self,
        bam_file: Optional[str] = None,
        coordinates: Optional[Union[str, np.ndarray]] = None,
        annotation: Optional[str] = None,
        apa_matrix: Optional[Union[str, np.ndarray, pd.DataFrame]] = None,
        h5ad_file: Optional[str] = None,
        matrix_orientation: str = 'genes_by_spots',
        expression_matrix: Optional[Union[str, np.ndarray, pd.DataFrame]] = None,
        expression_orientation: str = 'genes_by_spots',
        expression_embedding: Optional[np.ndarray] = None,
        dataset: Optional[APADataset] = None,
        impute: bool = True,
        quantify: bool = True,
        identify_domains: bool = True,
        use_bioml: Optional[bool] = None,
        differential_analysis: bool = False,
        detect_svapa: bool = True,
        n_domains: Optional[int] = None,
        domain_labels: Optional[np.ndarray] = None,
        fdr_threshold: float = 0.05,
        use_uncertainty_weights: bool = True,
    ) -> Dict:
        """
        Run the complete spaGAPA pipeline.

        Parameters
        ----------
        bam_file : str, optional
            Path to BAM file.
        coordinates : str or np.ndarray, optional
            Path to coordinates file or (n_spots, 2) array.
        annotation : str, optional
            Path to GTF/GFF annotation.
        apa_matrix : str or array-like, optional
            APA count/index matrix. By default rows are genes and columns are spots.
        h5ad_file : str, optional
            Existing H5AD file to load.
        matrix_orientation : str, default='genes_by_spots'
            Orientation of apa_matrix: 'genes_by_spots' or 'spots_by_genes'.
        expression_matrix : str or array-like, optional
            Optional expression matrix for BioML expression-view graph.
        expression_orientation : str, default='genes_by_spots'
            Orientation of expression_matrix: 'genes_by_spots' or 'spots_by_genes'.
        expression_embedding : np.ndarray, optional
            Precomputed spot-level expression embedding with shape
            (n_spots, n_features).
        dataset : APADataset, optional
            Pre-loaded dataset (alternative to bam_file).
        impute : bool, default=True
            Perform GP imputation.
        quantify : bool, default=True
            Calculate APA indices.
        identify_domains : bool, default=True
            Identify spatial domains.
        use_bioml : bool, optional
            Override ``self.use_bioml`` for this run.
        differential_analysis : bool, default=False
            Perform differential analysis.
        detect_svapa : bool, default=True
            Detect SVAPA genes.
        n_domains : int, optional
            Number of domains for clustering.
        domain_labels : np.ndarray, optional
            Pre-defined domain labels.
        fdr_threshold : float, default=0.05
            FDR threshold.
        use_uncertainty_weights : bool, default=True
            Use imputation uncertainty to weight downstream analyses.

        Returns
        -------
        results : dict
            Keys: dataset, imputed_values, uncertainty, apa_indices,
                  domains, differential, svapa_genes, qc_report
        """
        self._log("=" * 60)
        self._log("spaGAPA: Spatial GP-based APA Analyzer")
        self._log("=" * 60)

        # ── Step 1: Load data ──
        self._log("\n[1/7] Loading data...")
        if dataset is not None:
            self.dataset_ = dataset
        elif apa_matrix is not None or h5ad_file is not None or bam_file is not None:
            self.dataset_ = load_spatial_dataset(
                apa_matrix=apa_matrix,
                bam_file=bam_file,
                coordinates=coordinates,
                h5ad_file=h5ad_file,
                matrix_orientation=matrix_orientation,
                annotation=annotation,
            )
        else:
            raise ValueError("Provide 'dataset', 'apa_matrix', 'h5ad_file', or 'bam_file'")
        self._validate_dataset_shapes()
        self._log(f"  ✓ {self.dataset_.n_genes} genes, {self.dataset_.n_spots} spots")

        coords = self.dataset_.coords
        if coords is None:
            raise ValueError("No spatial coordinates in dataset")

        run_options = self._resolve_run_options(impute=impute, use_bioml=use_bioml)
        impute_this_run = bool(run_options['impute'])
        use_sparse_gp_this_run = bool(run_options['use_sparse_gp'])
        use_bioml_this_run = bool(run_options['use_bioml'])
        self.results_['analysis_preset'] = run_options
        self.dataset_.adata.uns['apa']['analysis_preset'] = run_options
        profile = run_options['profile']
        self._log(
            "  ✓ Analysis preset: "
            f"{run_options['requested_preset']} -> {run_options['resolved_preset']} "
            f"(spots={profile['n_spots']}, observed_fraction={profile['observed_fraction']:.3f})"
        )

        # ── Step 2: Spatial validation & filtering ──
        self._log("\n[2/7] Spatial validation + quality filtering...")
        spatial_neighbors = SpatialNeighbors(n_neighbors=self.n_neighbors)
        spatial_neighbors.fit(coords)
        self._log(f"  ✓ Built spatial graph (k={self.n_neighbors})")

        validator = SpatialValidator(
            n_neighbors=self.n_neighbors,
            support_threshold=self.min_spatial_support,
        )
        validator.fit(coords)

        qc_filter = QualityFilter(
            min_read_count=self.min_read_count,
            min_spots=self.min_spots,
            min_spatial_support=self.min_spatial_support,
        )
        self._log(f"  ✓ Validated {self.dataset_.n_genes} APA sites")

        # ── Step 3: GP imputation ──
        uncertainty = None
        if impute_this_run:
            self._log(f"\n[3/7] GP imputation (kernel={self.kernel_type})...")
            apa_matrix_values = self.dataset_.raw_counts  # (n_genes, n_spots)

            if use_sparse_gp_this_run:
                sparse_params = self._active_sparse_gp_params or {
                    'n_inducing': int(self.n_inducing),
                    'inducing_method': self.sparse_gp_inducing_method,
                    'length_scale': self.sparse_gp_length_scale,
                    'length_scale_multiplier': self.sparse_gp_length_scale_multiplier,
                    'noise_level': self.sparse_gp_noise_level,
                }
                length_scale = self._resolve_sparse_gp_length_scale(
                    coords,
                    value=sparse_params['length_scale'],
                    multiplier=sparse_params['length_scale_multiplier'],
                )
                base_imputer = SparseGPImputer(
                    n_inducing=int(sparse_params['n_inducing']),
                    inducing_method=sparse_params['inducing_method'],
                    length_scale=length_scale,
                    noise_level=sparse_params['noise_level'],
                )
                self.results_['analysis_preset']['sparse_gp']['length_scale_effective'] = float(length_scale)
            else:
                base_imputer = GPImputer(
                    kernel_type=self.kernel_type,
                    alpha=self.gp_alpha,
                    n_restarts_optimizer=self.gp_n_restarts_optimizer,
                )

            training_mask = self._training_mask(apa_matrix_values)
            self.imputer_ = base_imputer.fit_batch(
                coords,
                apa_matrix_values,
                mask=training_mask,
                n_jobs=self.gp_n_jobs,
                verbose=self.verbose,
            )
            imputed, uncertainty = self.imputer_.impute(return_uncertainty=True)
            imputed = self._clip_apa_values(imputed)

            self.dataset_.set_imputed(imputed, uncertainty)

            observed = training_mask
            if observed.sum() > 0:
                mae = float(np.mean(np.abs(
                    apa_matrix_values[observed] - imputed[observed]
                )))
                self._log(f"  ✓ Imputation done (MAE={mae:.3f})")
            else:
                self._log("  ✓ Imputation done")

            self.results_['imputed_values'] = self.dataset_.imputed
            self.results_['uncertainty'] = self.dataset_.uncertainty
        else:
            reason = "preset=highres_fast" if run_options['resolved_preset'] == 'highres_fast' else "requested"
            self._log(f"\n[3/7] Skipping imputation ({reason})")
            self.results_['imputed_values'] = None
            self.results_['uncertainty'] = None

        uncertainty = self.dataset_.uncertainty
        use_uw = use_uncertainty_weights and uncertainty is not None

        # Determine working matrix
        if self.dataset_.has_imputed():
            work = self.dataset_.imputed
        else:
            work = self.dataset_.raw_counts

        # ── Step 4: APA quantification ──
        if quantify:
            self._log("\n[4/7] APA quantification...")
            self.results_['apa_indices'] = self._calculate_apa_indices(work)
            if self.results_['apa_indices'].get('calculated'):
                self._log("  ✓ APA indices computed from proximal/distal counts")
            else:
                self._log("  ✓ Registered input APA index matrix")
        else:
            self._log("\n[4/7] Skipping quantification")
            self.results_['apa_indices'] = None

        # ── Step 5: Domain identification ──
        if identify_domains:
            self._log("\n[5/7] Identifying spatial domains...")
            if n_domains is None:
                n_domains = min(5, max(2, self.dataset_.n_spots // 20))

            if use_bioml_this_run:
                self._log(
                    "  BioML enabled "
                    f"(domain_method={self.bioml_domain_method}, "
                    f"rank={self.bioml_rank})"
                )
                expr_embedding = self._build_expression_embedding(
                    expression_matrix,
                    expression_embedding,
                    expression_orientation,
                )
                if run_options['resolved_preset'] in {'highres_accuracy', 'highres_fast'}:
                    domain_labels, domain_result, work = self._run_highres_bioml(
                        work,
                        coords,
                        uncertainty if use_uw else None,
                        expr_embedding,
                        n_domains,
                    )
                    self.results_['highres_bioml_values'] = work
                    if quantify:
                        self.results_['apa_indices'] = self._calculate_apa_indices(work)
                else:
                    domain_labels, domain_result = self._run_bioml(
                        work,
                        coords,
                        uncertainty if use_uw else None,
                        expr_embedding,
                        n_domains,
                    )
            else:
                domain_id = DomainIdentifier(
                    method='kmeans',
                    n_clusters=n_domains,
                    min_domain_size=self.min_spots,
                )
                domain_labels = domain_id.identify_domains(
                    work,
                    coords,
                    uncertainty=uncertainty if use_uw else None,
                )
                domain_result = {
                    'labels': domain_labels,
                    'n_domains': n_domains,
                    'method': 'kmeans',
                }

            # Store in dataset
            self.dataset_.set_domain_labels(domain_labels)

            self.results_['domains'] = domain_result
            self._log(f"  ✓ Identified {n_domains} spatial domains")
        else:
            self._log("\n[5/7] Skipping domain identification")
            self.results_['domains'] = None

        # ── Step 6: Differential analysis ──
        if differential_analysis:
            self._log("\n[6/7] Differential APA analysis...")
            if domain_labels is None and self.dataset_.has_domains():
                domain_labels = self.dataset_.get_domain_labels()

            if domain_labels is not None:
                analyzer = DifferentialAPAAnalyzer(method='wilcoxon')
                diff_results = analyzer.test_all_pairwise(
                    work,
                    domain_labels,
                    gene_names=self.dataset_.gene_names,
                    uncertainty=uncertainty if use_uw else None,
                )
                adjusted = {}
                for name, df in diff_results.items():
                    adjusted[name] = analyzer.adjust_pvalues(df)
                self.results_['differential'] = adjusted
                self._log("  ✓ Differential analysis complete")
            else:
                self._log("  ⚠ No domain labels; skipping")
                self.results_['differential'] = None
        else:
            self._log("\n[6/7] Skipping differential analysis")
            self.results_['differential'] = None

        # ── Step 7: SVAPA detection ──
        if detect_svapa:
            self._log(f"\n[7/7] Detecting SVAPA genes (FDR<{fdr_threshold})...")
            detector = GPTrendDetector(kernel_type=self.kernel_type)
            gene_names = self.dataset_.gene_names

            svapa_results = detector.detect_gp_trends(
                work, coords, uncertainty,
                gene_names=gene_names,
                method='likelihood_ratio',
                fdr_threshold=fdr_threshold,
            )
            self.results_['svapa_genes'] = svapa_results
            n_svapa = int(svapa_results['significant'].sum())
            self._log(f"  ✓ Detected {n_svapa} SVAPA genes")
        else:
            self._log("\n[7/7] Skipping SVAPA detection")
            self.results_['svapa_genes'] = None

        # ── QC report ──
        self._log("\nGenerating QC report...")
        qc = QCReportGenerator(dataset_name=f"spaGAPA_{self.dataset_.n_genes}genes")
        total = self.dataset_.n_spots * self.dataset_.n_genes
        observed = int((self.dataset_.adata.X > 0).sum())

        imputed_spots = observed
        if self.dataset_.has_imputed():
            imputed_spots = int((self.dataset_.imputed > 0).sum())

        qc.add_coverage_metrics(
            total_spots=total,
            observed_spots=observed,
            imputed_spots=imputed_spots,
        )

        if uncertainty is not None:
            qc.add_imputation_metrics({
                'mean_uncertainty': float(uncertainty.mean()),
                'max_uncertainty': float(uncertainty.max()),
                'min_uncertainty': float(uncertainty.min()),
            })

        self.results_['qc_report'] = qc.generate_report(format='dict')
        self.results_['dataset'] = self.dataset_

        self._log("\n" + "=" * 60)
        self._log("Pipeline complete!")
        self._log("=" * 60)

        return self.results_

    # ── accessors ────────────────────────────────────────────────────

    def get_imputed_values(self) -> Tuple[np.ndarray, np.ndarray]:
        im = self.dataset_.imputed
        un = self.dataset_.uncertainty
        if im is None:
            raise ValueError("Imputation not run")
        return im, un

    def get_svapa_genes(self, fdr_threshold: Optional[float] = None) -> pd.DataFrame:
        df = self.results_.get('svapa_genes')
        if df is None:
            raise ValueError("SVAPA detection not run")
        if fdr_threshold is not None:
            return df[df['q_value'] < fdr_threshold]
        return df[df['significant']]

    def get_domains(self) -> Dict:
        d = self.results_.get('domains')
        if d is None:
            raise ValueError("Domain identification not run")
        return d

    def save_results(self, output_dir: str):
        """Save all results."""
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        self._log(f"\nSaving results to {output_dir}...")

        if self.dataset_ is not None:
            self.dataset_.save(str(out / "dataset.h5ad"))
            self._log("  ✓ dataset.h5ad")

        imputed = self.dataset_.imputed
        if imputed is not None:
            np.save(out / "imputed_values.npy", imputed)
            self._log("  ✓ imputed_values.npy")

        uncertainty = self.dataset_.uncertainty
        if uncertainty is not None:
            np.save(out / "uncertainty.npy", uncertainty)
            self._log("  ✓ uncertainty.npy")

        domains = self.results_.get('domains')
        if domains is not None and domains.get('labels') is not None:
            labels = np.asarray(domains['labels'])
            pd.DataFrame({
                'spot': self.dataset_.spot_names,
                'domain': labels,
            }).to_csv(out / "domains.csv", index=False)
            self._log("  ✓ domains.csv")

        bioml = None
        if self.dataset_ is not None:
            bioml = self.dataset_.adata.uns.get('apa', {}).get('bioml')
        if bioml:
            bioml_imputed = self.dataset_.bioml_imputed
            if bioml_imputed is not None:
                np.save(out / "bioml_imputed_values.npy", bioml_imputed)
                self._log("  ✓ bioml_imputed_values.npy")
            if self.dataset_._OBSM_BIOML_FACTORS in self.dataset_.adata.obsm:
                np.save(
                    out / "bioml_spot_factors.npy",
                    self.dataset_.adata.obsm[self.dataset_._OBSM_BIOML_FACTORS],
                )
                self._log("  ✓ bioml_spot_factors.npy")
            if self.dataset_._VARM_BIOML_GENE_FACTORS in self.dataset_.adata.varm:
                np.save(
                    out / "bioml_gene_factors.npy",
                    self.dataset_.adata.varm[self.dataset_._VARM_BIOML_GENE_FACTORS],
                )
                self._log("  ✓ bioml_gene_factors.npy")
            (out / "bioml_metadata.json").write_text(json.dumps(bioml, indent=2))
            self._log("  ✓ bioml_metadata.json")

        preset = self.results_.get('analysis_preset')
        if preset is not None:
            (out / "analysis_preset.json").write_text(json.dumps(preset, indent=2))
            self._log("  ✓ analysis_preset.json")

        svapa = self.results_.get('svapa_genes')
        if svapa is not None:
            svapa.to_csv(out / "svapa_genes.csv", index=False)
            self._log("  ✓ svapa_genes.csv")

        qc = self.results_.get('qc_report')
        if qc is not None:
            (out / "qc_report.txt").write_text(str(qc))
            self._log("  ✓ qc_report.txt")

        self._log("All results saved.")
