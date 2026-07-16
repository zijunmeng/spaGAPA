"""
High-resolution BioML recovery for CPU-only spaGAPA workflows.

This module contains the package-level version of the decoupled highres_bioml
route that was first developed in the high-resolution benchmark runner:

* value recovery uses raw gene-mean fill plus an optional sparse-GP blend
* domain recovery uses a spatial/expression graph with an optional APA view
* adaptive graph neighborhoods are supported when parent/coarse-bin indices
  are available
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Dict, List, Optional

import numpy as np
from scipy import sparse
from sklearn.neighbors import NearestNeighbors

from .domain import BioMLDomainDetector
from .multiview_graph import MultiViewGraphBuilder


@dataclass
class HighResBioMLConfig:
    """Configuration for high-resolution decoupled BioML recovery."""

    gp_blend: float = 0.1
    spatial_weight: float = 0.1
    expression_weight: float = 0.7
    apa_weight: float = 0.2
    apa_source: str = "expression_knn"
    expression_knn_k: int = 15
    n_neighbors: int = 15
    neighbor_mode: str = "adaptive"
    adaptive_neighbor_scale: float = 10.0
    parent_weight: float = 0.0
    parent_neighbors: int = 8
    domain_method: str = "kmeans"
    random_state: int = 42

    def to_dict(self) -> Dict:
        return asdict(self)


@dataclass
class HighResBioMLResult:
    """Container for high-resolution BioML outputs."""

    recovered: np.ndarray
    labels: np.ndarray
    graph: sparse.csr_matrix
    metadata: Dict


def row_nanmean(values: np.ndarray) -> np.ndarray:
    """NaN-safe row means with zero fallback for fully missing rows."""
    arr = np.asarray(values, dtype=float)
    finite = np.isfinite(arr)
    counts = finite.sum(axis=1)
    sums = np.where(finite, arr, 0.0).sum(axis=1)
    means = np.divide(sums, np.maximum(counts, 1))
    means[counts == 0] = 0.0
    return means


def fill_missing_by_gene_mean(values: np.ndarray, clip: bool = True) -> np.ndarray:
    """Fill missing genes x spots APA values with gene-wise means."""
    filled = np.asarray(values, dtype=float).copy()
    means = row_nanmean(filled)
    missing_gene, missing_spot = np.where(~np.isfinite(filled))
    filled[missing_gene, missing_spot] = means[missing_gene]
    if clip:
        filled = np.clip(filled, 0.0, 1.0)
    return filled


def expression_knn_impute(
    values: np.ndarray,
    expression_embedding: np.ndarray,
    k: int = 15,
    max_iter: int = 10,
) -> np.ndarray:
    """
    Fill missing APA values from expression-neighbor spots.

    This is used only as a light APA proxy for high-resolution domain graphs,
    not as the final uncertainty-aware APA value model.
    """
    arr = np.asarray(values, dtype=float)
    expression = np.asarray(expression_embedding, dtype=float)
    if arr.ndim != 2:
        raise ValueError("values must have shape (n_genes, n_spots)")
    if expression.ndim != 2 or expression.shape[0] != arr.shape[1]:
        raise ValueError("expression_embedding must have n_spots rows")

    n_spots = arr.shape[1]
    if n_spots < 2:
        return fill_missing_by_gene_mean(arr)

    n_neighbors = min(max(1, int(k)), n_spots - 1)
    nn = NearestNeighbors(n_neighbors=n_neighbors + 1)
    nn.fit(np.nan_to_num(expression, nan=0.0, posinf=0.0, neginf=0.0))
    _, indices = nn.kneighbors(np.nan_to_num(expression, nan=0.0, posinf=0.0, neginf=0.0))
    neighbors = indices[:, 1:]

    filled = arr.copy()
    for _ in range(max(1, int(max_iter))):
        if np.isfinite(filled).all():
            break
        old_missing = int((~np.isfinite(filled)).sum())
        for spot_idx in range(n_spots):
            missing = ~np.isfinite(filled[:, spot_idx])
            if not missing.any():
                continue
            source = filled[:, neighbors[spot_idx]]
            means = row_nanmean(source)
            filled[missing, spot_idx] = means[missing]
        if int((~np.isfinite(filled)).sum()) == old_missing:
            break

    return fill_missing_by_gene_mean(filled)


def resolve_highres_bioml_neighbors(
    n_spots: int,
    parent_index: Optional[np.ndarray],
    config: HighResBioMLConfig,
) -> int:
    """Resolve graph KNN size for high-resolution BioML."""
    base = int(config.n_neighbors)
    if n_spots < 2:
        return 1
    if config.neighbor_mode == "fixed" or parent_index is None:
        return min(max(1, base), n_spots - 1)
    if config.neighbor_mode != "adaptive":
        raise ValueError("neighbor_mode must be 'fixed' or 'adaptive'")

    parent_index = np.asarray(parent_index)
    n_parent = max(1, int(len(np.unique(parent_index))))
    pseudo_bins_per_parent = n_spots / n_parent
    extra = max(0, int(round((pseudo_bins_per_parent - 4.0) * config.adaptive_neighbor_scale)))
    adaptive = base + extra
    return min(max(1, adaptive), n_spots - 1)


def aggregate_features_by_parent(features: np.ndarray, parent_index: np.ndarray) -> np.ndarray:
    """Average spot-level features to parent/coarse bins."""
    arr = np.asarray(features, dtype=float)
    parent_index = np.asarray(parent_index, dtype=int)
    if arr.ndim != 2:
        raise ValueError("features must be 2-dimensional")
    if arr.shape[0] != parent_index.shape[0]:
        raise ValueError("parent_index length must match number of feature rows")

    n_parent = int(parent_index.max()) + 1
    out = np.zeros((n_parent, arr.shape[1]), dtype=float)
    counts = np.bincount(parent_index, minlength=n_parent).astype(float)
    for dim in range(arr.shape[1]):
        out[:, dim] = np.bincount(parent_index, weights=arr[:, dim], minlength=n_parent)
    return out / np.maximum(counts[:, None], 1.0)


def lift_parent_graph(parent_graph: sparse.csr_matrix, parent_index: np.ndarray) -> sparse.csr_matrix:
    """Lift a sparse parent/coarse graph to high-resolution bins."""
    parent_index = np.asarray(parent_index, dtype=int)
    n_spots = len(parent_index)
    children = [
        np.where(parent_index == parent)[0]
        for parent in range(int(parent_index.max()) + 1)
    ]
    coo = parent_graph.tocoo()
    rows: List[int] = []
    cols: List[int] = []
    data: List[float] = []

    for parent_i, parent_j, value in zip(coo.row, coo.col, coo.data):
        if parent_i == parent_j or value <= 0:
            continue
        idx_i = children[int(parent_i)]
        idx_j = children[int(parent_j)]
        if idx_i.size == 0 or idx_j.size == 0:
            continue
        weight = float(value) / np.sqrt(float(idx_i.size * idx_j.size))
        block_rows = np.repeat(idx_i, idx_j.size)
        block_cols = np.tile(idx_j, idx_i.size)
        rows.extend(block_rows.tolist())
        cols.extend(block_cols.tolist())
        data.extend([weight] * len(block_rows))

    lifted = sparse.csr_matrix((data, (rows, cols)), shape=(n_spots, n_spots))
    lifted = lifted.maximum(lifted.T)
    lifted.setdiag(0.0)
    lifted.eliminate_zeros()
    return lifted


def add_parent_aware_graph(
    fused_graph: sparse.csr_matrix,
    coords: np.ndarray,
    expression_embedding: Optional[np.ndarray],
    parent_index: Optional[np.ndarray],
    config: HighResBioMLConfig,
) -> sparse.csr_matrix:
    """Fuse a high-resolution graph with an optional lifted parent graph."""
    parent_weight = float(np.clip(config.parent_weight, 0.0, 1.0))
    if parent_weight <= 0 or parent_index is None:
        return fused_graph

    parent_index = np.asarray(parent_index, dtype=int)
    if len(np.unique(parent_index)) < 2:
        return fused_graph

    parent_coords = aggregate_features_by_parent(coords, parent_index)
    parent_expr = None
    expression_weight = 0.0
    if expression_embedding is not None:
        parent_expr = aggregate_features_by_parent(expression_embedding, parent_index)
        expression_weight = 0.5

    parent_graph = MultiViewGraphBuilder(
        n_neighbors=min(max(1, config.parent_neighbors), parent_coords.shape[0] - 1),
        spatial_weight=1.0 - expression_weight,
        expression_weight=expression_weight,
        apa_weight=0.0,
    ).build(
        parent_coords,
        expression_embedding=parent_expr,
    )
    lifted = lift_parent_graph(parent_graph.fused, parent_index)
    combined = (1.0 - parent_weight) * fused_graph + parent_weight * lifted
    combined = combined.maximum(combined.T)
    combined.setdiag(0.0)
    combined.eliminate_zeros()
    return combined.tocsr()


def highres_bioml_recover(
    observed: np.ndarray,
    coords: np.ndarray,
    expression_embedding: Optional[np.ndarray],
    sparse_gp: Optional[np.ndarray],
    uncertainty: Optional[np.ndarray],
    parent_index: Optional[np.ndarray],
    n_domains: int,
    config: Optional[HighResBioMLConfig] = None,
) -> HighResBioMLResult:
    """Run decoupled high-resolution BioML value/domain recovery."""
    config = config or HighResBioMLConfig()
    observed = np.asarray(observed, dtype=float)
    coords = np.asarray(coords, dtype=float)
    if observed.ndim != 2:
        raise ValueError("observed must have shape (n_genes, n_spots)")
    if coords.ndim != 2 or coords.shape[0] != observed.shape[1]:
        raise ValueError("coords must have n_spots rows")
    if observed.shape[1] < 2:
        raise ValueError("highres_bioml_recover requires at least two spots/bins")
    if expression_embedding is not None:
        expression_embedding = np.asarray(expression_embedding, dtype=float)
        if expression_embedding.shape[0] != observed.shape[1]:
            raise ValueError("expression_embedding must have n_spots rows")

    notes: List[str] = []
    mask = np.isfinite(observed)
    raw_matrix = fill_missing_by_gene_mean(observed)

    gp_blend = float(np.clip(config.gp_blend, 0.0, 1.0))
    if gp_blend > 0 and sparse_gp is None:
        notes.append("gp_blend reset to 0 because sparse_gp is unavailable")
        gp_blend = 0.0
    if sparse_gp is not None:
        sparse_gp = np.asarray(sparse_gp, dtype=float)
        if sparse_gp.shape != observed.shape:
            raise ValueError("sparse_gp must match observed shape")

    recovered = raw_matrix.copy() if gp_blend == 0 else (1.0 - gp_blend) * raw_matrix + gp_blend * sparse_gp
    recovered[mask] = observed[mask]
    recovered = np.clip(recovered, 0.0, 1.0)

    apa_source_requested = config.apa_source
    apa_source_effective = apa_source_requested
    apa_matrix = None
    apa_uncertainty = None
    if apa_source_requested == "expression_knn":
        if expression_embedding is None:
            apa_source_effective = "raw"
            apa_matrix = raw_matrix
            notes.append("expression_knn APA source fell back to raw because expression is unavailable")
        else:
            apa_matrix = expression_knn_impute(
                observed,
                expression_embedding,
                k=config.expression_knn_k,
            )
    elif apa_source_requested == "raw":
        apa_matrix = raw_matrix
    elif apa_source_requested == "sparse_gp":
        if sparse_gp is None:
            apa_source_effective = "raw"
            apa_matrix = raw_matrix
            notes.append("sparse_gp APA source fell back to raw because sparse_gp is unavailable")
        else:
            apa_matrix = sparse_gp
            apa_uncertainty = uncertainty
    elif apa_source_requested == "none":
        apa_matrix = None
        apa_source_effective = "none"
    else:
        raise ValueError("apa_source must be 'expression_knn', 'raw', 'sparse_gp', or 'none'")

    n_neighbors = resolve_highres_bioml_neighbors(observed.shape[1], parent_index, config)
    graph = MultiViewGraphBuilder(
        n_neighbors=n_neighbors,
        spatial_weight=config.spatial_weight,
        expression_weight=config.expression_weight,
        apa_weight=config.apa_weight if apa_matrix is not None else 0.0,
    ).build(
        coords,
        expression_embedding=expression_embedding,
        apa_matrix=apa_matrix,
        uncertainty=apa_uncertainty,
    )
    fused_graph = add_parent_aware_graph(
        graph.fused,
        coords,
        expression_embedding,
        parent_index,
        config,
    )

    domain_method = config.domain_method
    if domain_method == "spectral" and fused_graph.shape[0] > 5000:
        notes.append(f"domain_method auto-switched spectral→kmeans (n_spots={fused_graph.shape[0]} > 5000)")
        domain_method = "kmeans"

    if domain_method == "spectral":
        labels = BioMLDomainDetector(
            method="spectral",
            n_domains=n_domains,
            random_state=config.random_state,
        ).fit_predict(graph=fused_graph)
    elif domain_method == "kmeans":
        from sklearn.cluster import KMeans
        recovered_T = recovered.T  # (n_spots, n_genes)
        km = KMeans(n_clusters=n_domains, random_state=config.random_state, n_init=10)
        labels = km.fit_predict(np.nan_to_num(recovered_T))
    elif domain_method == "leiden":
        try:
            import scanpy as sc
            import anndata as ad
            adata = ad.AnnData(X=np.nan_to_num(recovered.T))
            adata.obsp['connectivities'] = fused_graph
            sc.tl.leiden(adata, adjacency=fused_graph, resolution=1.0,
                          random_state=config.random_state)
            labels = adata.obs['leiden'].astype(int).values
        except ImportError:
            notes.append("Leiden requires scanpy — falling back to kmeans")
            from sklearn.cluster import KMeans
            km = KMeans(n_clusters=n_domains, random_state=config.random_state, n_init=10)
            labels = km.fit_predict(np.nan_to_num(recovered.T))
    else:
        raise ValueError(f"Unknown domain_method: {domain_method}. Use 'spectral', 'kmeans', or 'leiden'")

    metadata = {
        "mode": "highres_bioml",
        "config": config.to_dict(),
        "n_domains": int(n_domains),
        "n_neighbors_effective": int(n_neighbors),
        "gp_blend_effective": float(gp_blend),
        "apa_source_requested": apa_source_requested,
        "apa_source_effective": apa_source_effective,
        "graph_weights_requested": {
            "spatial": config.spatial_weight,
            "expression": config.expression_weight,
            "apa": config.apa_weight if apa_matrix is not None else 0.0,
        },
        "graph_weights_effective": graph.weights,
        "has_expression_view": expression_embedding is not None,
        "has_parent_index": parent_index is not None,
        "parent_weight_effective": float(np.clip(config.parent_weight, 0.0, 1.0)),
        "notes": notes,
    }
    return HighResBioMLResult(
        recovered=recovered,
        labels=labels,
        graph=fused_graph,
        metadata=metadata,
    )
