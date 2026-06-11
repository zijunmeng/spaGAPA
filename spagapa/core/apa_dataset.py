"""
Core data structure for APA analysis in spatial transcriptomics.

This module provides APADataset, a tight wrapper around AnnData that stores
all APA-specific data in standard AnnData locations for ecosystem compatibility.
"""

import numpy as np
import pandas as pd
from typing import Optional, Dict, List, Union, Tuple
import logging
import warnings

try:
    import anndata
    ANNDATA_AVAILABLE = True
except ImportError:
    ANNDATA_AVAILABLE = False
    logging.warning("anndata not installed. Install with: pip install anndata")

logger = logging.getLogger(__name__)


class APADataset:
    """
    Thin wrapper around AnnData for APA spatial transcriptomics data.

    All data is stored in standard AnnData slots:
    - adata.X: (n_spots, n_genes) raw count matrix
    - adata.layers['apa_counts']: raw APA counts
    - adata.layers['apa_imputed']: GP-imputed APA values
    - adata.layers['apa_uncertainty']: GP uncertainty (std) per imputed value
    - adata.obsm['spatial']: (n_spots, 2) spatial coordinates
    - adata.obs[['spatial_x','spatial_y']]: per-spot coordinates
    - adata.obsm['RUD'] / adata.obsm['PDUI'] / adata.obsm['WUL']: APA indices
    - adata.uns['apa']: unstructured (sites, params, results)

    Delegates unknown attribute access to the underlying AnnData object,
    so you can call scanpy/squidpy functions directly:
        dataset.raw  # delegates to adata.raw
        sc.pp.normalize_total(dataset.adata)
    """

    _LAYER_RAW = 'apa_counts'
    _LAYER_IMPUTED = 'apa_imputed'
    _LAYER_UNCERTAINTY = 'apa_uncertainty'
    _LAYER_BIOML_IMPUTED = 'apa_bioml_imputed'
    _OBSM_SPATIAL = 'spatial'
    _OBSM_BIOML_FACTORS = 'bioml_factors'
    _VARM_BIOML_GENE_FACTORS = 'bioml_gene_factors'

    def __init__(self, adata: Optional['anndata.AnnData'] = None):
        if not ANNDATA_AVAILABLE:
            raise ImportError(
                "anndata is required for APADataset. "
                "Install with: pip install anndata"
            )
        if adata is None:
            adata = anndata.AnnData()
        self.adata = adata
        self._ensure_spatial_coords()
        if 'apa' not in self.adata.uns:
            self.adata.uns['apa'] = {
                'sites': None,
                'parameters': {},
                'differential': {},
                'spatial_patterns': {},
                'version': '0.1.0',
            }

    # ── delegation to adata ──────────────────────────────────────────

    def __getattr__(self, name):
        """Delegate unknown attributes to self.adata."""
        if name == 'adata':
            raise AttributeError("_adata not yet initialized")
        return getattr(self.adata, name)

    def __repr__(self) -> str:
        parts = [f"APADataset(n_spots={self.n_spots}, n_genes={self.n_genes})"]
        parts.append(f"  layers: {list(self.adata.layers.keys())}")
        parts.append(f"  obsm:   {list(self.adata.obsm.keys())}")
        parts.append(f"  imputed: {self.has_imputed()}")
        parts.append(f"  uncertainty: {self.has_uncertainty()}")
        return "\n".join(parts)

    # ── properties ───────────────────────────────────────────────────

    @property
    def n_spots(self) -> int:
        return self.adata.n_obs

    @property
    def n_genes(self) -> int:
        return self.adata.n_vars

    @property
    def gene_names(self) -> List[str]:
        return self.adata.var_names.tolist()

    @property
    def spot_names(self) -> List[str]:
        return self.adata.obs_names.tolist()

    @property
    def coords(self) -> np.ndarray:
        """Spatial coordinates, shape (n_spots, 2)."""
        return self.adata.obsm[self._OBSM_SPATIAL]

    @property
    def raw_counts(self) -> np.ndarray:
        """Raw APA counts, shape (n_genes, n_spots)."""
        layer = self._LAYER_RAW
        if layer in self.adata.layers:
            return self.adata.layers[layer].T
        return self.adata.X.T

    @property
    def imputed(self) -> Optional[np.ndarray]:
        """Imputed APA values, shape (n_genes, n_spots)."""
        layer = self._LAYER_IMPUTED
        if layer in self.adata.layers:
            return self.adata.layers[layer].T
        return None

    @property
    def uncertainty(self) -> Optional[np.ndarray]:
        """Uncertainty (std) per imputed value, shape (n_genes, n_spots)."""
        layer = self._LAYER_UNCERTAINTY
        if layer in self.adata.layers:
            return self.adata.layers[layer].T
        return None

    @property
    def bioml_imputed(self) -> Optional[np.ndarray]:
        """BioML-refined APA values, shape (n_genes, n_spots)."""
        layer = self._LAYER_BIOML_IMPUTED
        if layer in self.adata.layers:
            return self.adata.layers[layer].T
        return None

    # ── status checks ────────────────────────────────────────────────

    def has_imputed(self) -> bool:
        return self._LAYER_IMPUTED in self.adata.layers

    def has_uncertainty(self) -> bool:
        return self._LAYER_UNCERTAINTY in self.adata.layers

    def has_bioml(self) -> bool:
        return self._LAYER_BIOML_IMPUTED in self.adata.layers or self._OBSM_BIOML_FACTORS in self.adata.obsm

    def has_raw_counts(self) -> bool:
        return self._LAYER_RAW in self.adata.layers

    def has_indices(self) -> bool:
        return bool(self.adata.uns.get('apa', {}).get('indices', {}))

    def has_domains(self) -> bool:
        return 'domain' in self.adata.obs.columns

    # ── APA index storage in obsm ────────────────────────────────────

    def set_apa_index(self, name: str, values: np.ndarray):
        """
        Store an APA index (RUD, PDUI, WUL) in obsm and uns.

        Parameters
        ----------
        name : str
            Index name (e.g. 'RUD', 'PDUI', 'WUL')
        values : np.ndarray, shape (n_spots,) or (n_genes, n_spots)
            Index values. If 2D, stored as (n_spots, n_genes) in obsm.
        """
        if values.ndim == 2:
            self.adata.obsm[name] = values.T  # (n_spots, n_genes)
        else:
            self.adata.obsm[name] = values[:, np.newaxis]
        self.adata.uns['apa'].setdefault('indices', {})[name] = 'stored in obsm'

    def get_apa_index(self, name: str) -> Optional[np.ndarray]:
        """
        Get an APA index, shape (n_genes, n_spots).

        Parameters
        ----------
        name : str
            Index name (e.g. 'RUD', 'PDUI', 'WUL')
        """
        if name in self.adata.obsm:
            return self.adata.obsm[name].T
        return None

    # ── imputation ───────────────────────────────────────────────────

    def set_imputed(
        self,
        imputed: np.ndarray,
        uncertainty: Optional[np.ndarray] = None
    ):
        """
        Store GP imputation results.

        Parameters
        ----------
        imputed : np.ndarray, (n_genes, n_spots)
        uncertainty : np.ndarray, optional, (n_genes, n_spots)
        """
        self.adata.layers[self._LAYER_IMPUTED] = imputed.T
        if uncertainty is not None:
            self.adata.layers[self._LAYER_UNCERTAINTY] = uncertainty.T
        logger.info(f"Stored imputed values; uncertainty={'yes' if uncertainty is not None else 'no'}")

    def set_bioml_results(
        self,
        imputed: Optional[np.ndarray] = None,
        spot_factors: Optional[np.ndarray] = None,
        gene_factors: Optional[np.ndarray] = None,
        metadata: Optional[Dict] = None,
    ):
        """
        Store BioML outputs in AnnData-compatible locations.

        Parameters
        ----------
        imputed : np.ndarray, optional, shape (n_genes, n_spots)
            BioML-refined APA matrix.
        spot_factors : np.ndarray, optional, shape (n_spots, rank)
            Spot-level BioML embedding.
        gene_factors : np.ndarray, optional, shape (n_genes, rank)
            Gene-level BioML embedding.
        metadata : dict, optional
            Parameters and graph metadata for reproducibility.
        """
        if imputed is not None:
            imputed = np.asarray(imputed)
            if imputed.shape != (self.n_genes, self.n_spots):
                raise ValueError("imputed must have shape (n_genes, n_spots)")
            self.adata.layers[self._LAYER_BIOML_IMPUTED] = imputed.T

        if spot_factors is not None:
            spot_factors = np.asarray(spot_factors)
            if spot_factors.shape[0] != self.n_spots:
                raise ValueError("spot_factors must have n_spots rows")
            self.adata.obsm[self._OBSM_BIOML_FACTORS] = spot_factors

        if gene_factors is not None:
            gene_factors = np.asarray(gene_factors)
            if gene_factors.shape[0] != self.n_genes:
                raise ValueError("gene_factors must have n_genes rows")
            self.adata.varm[self._VARM_BIOML_GENE_FACTORS] = gene_factors

        self.adata.uns['apa']['bioml'] = metadata or {}

    # ── spatial coords ───────────────────────────────────────────────

    def _ensure_spatial_coords(self):
        """Ensure coordinates are in obsm and obs."""
        if self._OBSM_SPATIAL not in self.adata.obsm:
            return
        coords = self.adata.obsm[self._OBSM_SPATIAL]
        if 'spatial_x' not in self.adata.obs.columns:
            self.adata.obs['spatial_x'] = coords[:, 0]
        if 'spatial_y' not in self.adata.obs.columns:
            self.adata.obs['spatial_y'] = coords[:, 1]

    def set_coords(self, coords: np.ndarray):
        """
        Set spatial coordinates.

        Parameters
        ----------
        coords : np.ndarray, shape (n_spots, 2)
        """
        self.adata.obsm[self._OBSM_SPATIAL] = coords
        self.adata.obs['spatial_x'] = coords[:, 0]
        self.adata.obs['spatial_y'] = coords[:, 1]

    # ── domain labels ────────────────────────────────────────────────

    def set_domain_labels(self, labels: np.ndarray, name: str = 'domain'):
        """
        Store spatial domain labels.

        Parameters
        ----------
        labels : np.ndarray, shape (n_spots,)
        name : str
            Column name in obs (default 'domain')
        """
        self.adata.obs[name] = pd.Categorical(labels)

    def get_domain_labels(self, name: str = 'domain') -> Optional[np.ndarray]:
        """Get domain labels as integer array."""
        if name in self.adata.obs.columns:
            return self.adata.obs[name].cat.codes.values
        return None

    # ── APA sites ────────────────────────────────────────────────────

    def set_apa_sites(self, sites: pd.DataFrame):
        """Store APA site annotations."""
        self.adata.uns['apa']['sites'] = sites
        logger.info(f"Stored {len(sites)} APA sites")

    # ── factory methods ──────────────────────────────────────────────

    @classmethod
    def from_counts(
        cls,
        apa_counts: Union[np.ndarray, pd.DataFrame],
        spatial_coords: Union[np.ndarray, pd.DataFrame],
        gene_names: Optional[List[str]] = None,
        spot_names: Optional[List[str]] = None,
        apa_sites: Optional[pd.DataFrame] = None,
    ) -> 'APADataset':
        """Create from count matrix and coordinates."""
        if isinstance(apa_counts, pd.DataFrame):
            if gene_names is None:
                gene_names = apa_counts.index.tolist()
            if spot_names is None:
                spot_names = apa_counts.columns.tolist()
            apa_counts = apa_counts.values

        if isinstance(spatial_coords, pd.DataFrame):
            if spot_names is None:
                spot_names = spatial_coords.index.tolist()
            spatial_coords = spatial_coords.values

        n_genes, n_spots = apa_counts.shape
        if gene_names is None:
            gene_names = [f"Gene_{i}" for i in range(n_genes)]
        if spot_names is None:
            spot_names = [f"Spot_{i}" for i in range(n_spots)]

        adata = anndata.AnnData(
            X=apa_counts.T,
            obs=pd.DataFrame(index=spot_names),
            var=pd.DataFrame(index=gene_names),
        )
        adata.obsm[cls._OBSM_SPATIAL] = spatial_coords
        adata.obs['spatial_x'] = spatial_coords[:, 0]
        adata.obs['spatial_y'] = spatial_coords[:, 1]
        adata.layers[cls._LAYER_RAW] = apa_counts.T

        dataset = cls(adata)
        if apa_sites is not None:
            dataset.set_apa_sites(apa_sites)
        logger.info(f"Created APADataset: {n_genes} genes, {n_spots} spots")
        return dataset

    @classmethod
    def from_anndata(cls, adata: 'anndata.AnnData') -> 'APADataset':
        """Wrap an existing AnnData object."""
        return cls(adata)

    # ── I/O ──────────────────────────────────────────────────────────

    def to_anndata(self) -> 'anndata.AnnData':
        return self.adata

    def save(self, filename: str):
        self.adata.write_h5ad(filename)
        logger.info(f"Saved to {filename}")

    @classmethod
    def load(cls, filename: str) -> 'APADataset':
        logger.info(f"Loading from {filename}")
        return cls(anndata.read_h5ad(filename))

    # ── subsetting ───────────────────────────────────────────────────

    def subset_genes(self, gene_list: List[str]) -> 'APADataset':
        mask = self.adata.var_names.isin(gene_list)
        return APADataset(self.adata[:, mask].copy())

    def subset_spots(self, spot_list: List[str]) -> 'APADataset':
        mask = self.adata.obs_names.isin(spot_list)
        return APADataset(self.adata[mask, :].copy())

    # ── backward-compatible aliases ───────────────────────────────────

    def get_apa_counts(self, imputed: bool = False) -> np.ndarray:
        """(deprecated alias) Get APA counts. Use .raw_counts or .imputed."""
        return self.imputed if imputed else self.raw_counts

    def get_spatial_coords(self) -> np.ndarray:
        """(deprecated alias) Use .coords instead."""
        return self.coords

    def add_imputation(
        self, imputed_counts: np.ndarray,
        uncertainty: Optional[np.ndarray] = None,
    ):
        """(deprecated alias) Use .set_imputed() instead."""
        self.set_imputed(imputed_counts, uncertainty)

    def add_apa_sites(self, apa_sites: pd.DataFrame):
        """(deprecated alias) Use .set_apa_sites() instead."""
        self.set_apa_sites(apa_sites)

    def add_apa_indices(self, indices: Dict[str, np.ndarray]):
        """Store APA indices (RUD, PDUI, etc.)."""
        for name, arr in indices.items():
            self.set_apa_index(name, arr)
        self.adata.uns['apa'].setdefault('indices', {}).update(indices)
