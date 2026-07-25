"""
SVAPA: Spatially Variable APA gene detection.

This module detects genes whose alternative-polyadenylation (APA) usage index
varies significantly across spatial coordinates.  It is deliberately
self-contained and lightweight: spatial neighbours are built once via a
``scipy.spatial.cKDTree`` and stored as a ``scipy.sparse`` row-normalised
weight matrix, then Moran's I (and Geary's C) are computed per gene with
permutation-based p-values and BH-FDR correction.

The implementation favours a transparent, dependency-light Moran's I:

    I = (n / S0) * (z' W z) / (z' z) ,   z = values - mean(values)

over the OBSERVED (finite) spots only, with W a symmetric kNN spatial-weights
matrix and S0 = sum of all weights.  Permutation p-values are obtained by
randomly relabelling the observed values across spots ``n_perm`` times.

Public API
----------
- :func:`morans_i`  -> (I, pvalue, n)
- :func:`gearys_c`  -> (C, pvalue, n)
- :func:`svapa`     -> DataFrame(gene, morans_i, gearys_c, pvalue, padj)
"""

from __future__ import annotations

from typing import Optional, Tuple, Union

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree

try:  # FDR correction (optional dependency)
    from statsmodels.stats.multitest import multipletests
    _HAS_SM = True
except Exception:  # pragma: no cover - statsmodels usually present
    _HAS_SM = False


# ---------------------------------------------------------------------------
# Spatial weight construction
# ---------------------------------------------------------------------------
def _build_spatial_weights(
    coords: np.ndarray,
    k: int = 8,
    symmetric: bool = True,
) -> sparse.csr_matrix:
    """Build a row-normalised symmetric kNN spatial-weights matrix.

    Parameters
    ----------
    coords : ndarray, shape (n_spots, 2)
        Spot spatial coordinates.
    k : int
        Number of nearest neighbours (excluding self).
    symmetric : bool
        Symmetrise the adjacency (W = W | W.T) so that Moran's I's S0 is
        well-defined and the statistic matches the classic formula.

    Returns
    -------
    W : scipy.sparse.csr_matrix, shape (n_spots, n_spots)
        Row-normalised spatial weights.
    """
    coords = np.asarray(coords, dtype=float)
    n = coords.shape[0]
    kk = max(1, min(int(k), n - 1))
    tree = cKDTree(coords)
    # exclude self: query k+1 then drop the first column
    _, idx = tree.query(coords, k=kk + 1)
    nbr = idx[:, 1:]  # (n, kk)
    rows = np.repeat(np.arange(n), kk)
    cols = nbr.reshape(-1)
    data = np.ones(len(rows), dtype=float)
    W = sparse.csr_matrix((data, (rows, cols)), shape=(n, n))
    if symmetric:
        W = (W + W.T)
        W.data[:] = 1.0  # binary symmetric links
    # row-normalise
    row_sums = np.asarray(W.sum(axis=1)).ravel()
    row_sums[row_sums == 0] = 1.0
    D = sparse.diags(1.0 / row_sums)
    W = D @ W
    return W.tocsr()


# ---------------------------------------------------------------------------
# Core per-gene statistics
# ---------------------------------------------------------------------------
def _morans_i_from_W(z: np.ndarray, W: sparse.csr_matrix, n: int,
                     S0: float) -> float:
    """Moran's I from centred values ``z`` and weights ``W``."""
    denom = float(z @ z)
    if denom <= 0:
        return np.nan
    Wz = W @ z
    num = float(z @ Wz)
    return (n / S0) * (num / denom)


def _gearys_c_from_W(values: np.ndarray, W: sparse.csr_matrix, n: int,
                     S0: float) -> float:
    """Geary's C from raw values ``values`` and weights ``W``.

    C = ((n-1)/S0) * (sum_ij w_ij (x_i - x_j)^2) / (sum_i (x_i - mean)^2)
    """
    mean = values.mean()
    dev = values - mean
    denom = float(dev @ dev)
    if denom <= 0:
        return np.nan
    # Geary's C = (n-1) * sum_ij w_ij (x_i-x_j)^2 / (2 * S0 * sum_i (x_i-mean)^2)
    # For symmetric W, sum_ij w_ij (x_i-x_j)^2 = 2*(sum_sq_row - sum_xWx),
    # where sum_sq_row = sum_i x_i^2 rowSum_i and sum_xWx = sum_ij w_ij x_i x_j,
    # so the factors of 2 cancel: C = (n-1)/S0 * (sum_sq_row - sum_xWx) / denom.
    Wx = W @ values
    sum_xWx = float(values @ Wx)  # sum_ij w_ij x_i x_j  (W symmetric)
    row_sum = np.asarray(W.sum(axis=1)).ravel()
    sum_sq_row = float((values ** 2) @ row_sum)  # sum_i x_i^2 rowSum_i
    C = ((n - 1) / S0) * ((sum_sq_row - sum_xWx) / denom)
    return C


def _permutation_pvalue(
    obs_stat: float,
    values: np.ndarray,
    W: sparse.csr_matrix,
    n: int,
    S0: float,
    n_perm: int,
    rng: np.random.Generator,
    geary: bool = False,
) -> float:
    """One-sided permutation p-value for Moran's I / Geary's C.

    For Moran's I we count permutations with I_perm >= I_obs (positive
    autocorrelation); for Geary's C we count C_perm <= C_obs (C small
    indicates positive autocorrelation).
    """
    if not np.isfinite(obs_stat) or n_perm <= 0:
        return np.nan
    count = 0
    perm_vals = values.copy()
    for _ in range(n_perm):
        rng.shuffle(perm_vals)
        if geary:
            s = _gearys_c_from_W(perm_vals, W, n, S0)
            if np.isfinite(s) and s <= obs_stat:
                count += 1
        else:
            z = perm_vals - perm_vals.mean()
            s = _morans_i_from_W(z, W, n, S0)
            if np.isfinite(s) and s >= obs_stat:
                count += 1
    # +1 / (n_perm+1) to avoid p == 0
    return (count + 1) / (n_perm + 1)


def morans_i(
    values_1d: np.ndarray,
    coords: np.ndarray,
    k: int = 8,
    n_perm: int = 100,
    seed: Optional[int] = None,
) -> Tuple[float, float, int]:
    """Moran's I spatial-autocorrelation for a single gene.

    Parameters
    ----------
    values_1d : ndarray, shape (n_spots,)
        Per-spot APA values for one gene.  NaNs are dropped; the statistic
        is computed over the observed (finite) spots only.
    coords : ndarray, shape (n_spots, 2)
        Spot coordinates aligned with ``values_1d``.  When ``values_1d``
        contains NaNs, both arrays are subset to the finite spots.
    k : int, default 8
        Number of nearest spatial neighbours for the weight matrix.
    n_perm : int, default 100
        Number of permutations for the p-value.  Set to 0 to skip.
    seed : int, optional
        RNG seed for permutation reproducibility.

    Returns
    -------
    I : float
        Moran's I statistic (NaN if undefined).
    pvalue : float
        One-sided permutation p-value (NaN if ``n_perm<=0``).
    n : int
        Number of observed spots used.

    Examples
    --------
    >>> import numpy as np
    >>> rng = np.random.default_rng(0)
    >>> # 10x10 grid
    >>> coords = np.mgrid[0:10, 0:10].reshape(2, -1).T.astype(float)
    >>> # spatial gradient -> high Moran's I
    >>> v = coords[:, 0] + rng.normal(0, 0.1, 100)
    >>> I, p, n = morans_i(v, coords, k=8, n_perm=50, seed=0)
    >>> I > 0.5, p < 0.05, n == 100
    (True, True, True)
    """
    values_1d = np.asarray(values_1d, dtype=float).ravel()
    coords = np.asarray(coords, dtype=float)
    if coords.ndim == 1:
        coords = coords.reshape(-1, 2)
    if values_1d.shape[0] != coords.shape[0]:
        raise ValueError(
            f"values_1d ({values_1d.shape[0]}) and coords "
            f"({coords.shape[0]}) length mismatch"
        )
    finite = np.isfinite(values_1d)
    n = int(finite.sum())
    if n < 3:
        return np.nan, np.nan, n
    v = values_1d[finite]
    c = coords[finite]
    if n - 1 < k:
        k = n - 1
    W = _build_spatial_weights(c, k=k, symmetric=True)
    S0 = float(W.sum())
    if S0 <= 0:
        return np.nan, np.nan, n
    z = v - v.mean()
    I = _morans_i_from_W(z, W, n, S0)
    if n_perm and n_perm > 0 and np.isfinite(I):
        rng = np.random.default_rng(seed)
        p = _permutation_pvalue(I, v, W, n, S0, n_perm, rng, geary=False)
    else:
        p = np.nan
    return I, p, n


def gearys_c(
    values_1d: np.ndarray,
    coords: np.ndarray,
    k: int = 8,
    n_perm: int = 100,
    seed: Optional[int] = None,
) -> Tuple[float, float, int]:
    """Geary's C spatial-autocorrelation for a single gene.

    Returns
    -------
    C : float
        Geary's C (1 = no autocorrelation; <1 positive, >1 negative).
    pvalue : float
        One-sided permutation p-value.
    n : int
        Number of observed spots.
    """
    values_1d = np.asarray(values_1d, dtype=float).ravel()
    coords = np.asarray(coords, dtype=float)
    if coords.ndim == 1:
        coords = coords.reshape(-1, 2)
    finite = np.isfinite(values_1d)
    n = int(finite.sum())
    if n < 3:
        return np.nan, np.nan, n
    v = values_1d[finite]
    c = coords[finite]
    if n - 1 < k:
        k = n - 1
    W = _build_spatial_weights(c, k=k, symmetric=True)
    S0 = float(W.sum())
    if S0 <= 0:
        return np.nan, np.nan, n
    C = _gearys_c_from_W(v, W, n, S0)
    if n_perm and n_perm > 0 and np.isfinite(C):
        rng = np.random.default_rng(seed)
        p = _permutation_pvalue(C, v, W, n, S0, n_perm, rng, geary=True)
    else:
        p = np.nan
    return C, p, n


# ---------------------------------------------------------------------------
# Matrix-level SVAPA
# ---------------------------------------------------------------------------
def svapa(
    apa_matrix: Union[np.ndarray, pd.DataFrame],
    coords: np.ndarray,
    k: int = 8,
    n_perm: int = 100,
    fdr: Optional[str] = "fdr_bh",
    seed: Optional[int] = 42,
    min_obs: int = 10,
    include_geary: bool = True,
) -> pd.DataFrame:
    """Detect spatially variable APA genes.

    Parameters
    ----------
    apa_matrix : ndarray or DataFrame, shape (n_genes, n_spots)
        Gene x spot APA index matrix.  NaN marks missing entries.  If a
        DataFrame, its index becomes the ``gene`` column.
    coords : ndarray, shape (n_spots, 2)
        Spot coordinates aligned with the matrix columns.
    k : int, default 8
        Number of nearest spatial neighbours.
    n_perm : int, default 100
        Permutations per gene for the Moran's I p-value.
    fdr : str or None, default 'fdr_bh'
        Multiple-testing correction method (passed to
        ``statsmodels.stats.multitest.multipletests``).  None skips FDR.
    seed : int, optional
        RNG seed.
    min_obs : int, default 10
        Genes with fewer than this many observed spots are dropped.
    include_geary : bool, default True
        Also compute Geary's C (no extra permutations; secondary statistic).

    Returns
    -------
    results : pd.DataFrame
        Columns: ``gene, morans_i, morans_i_pvalue, gearys_c,
        gearys_c_pvalue, n_obs, padj`` (morans-I BH-FDR).  Rows sorted by
        ``morans_i`` descending; dropped genes are excluded.

    Examples
    --------
    >>> import numpy as np
    >>> # 10x10 grid
    >>> coords = np.mgrid[0:10,0:10].reshape(2,-1).T.astype(float)
    >>> rng = np.random.default_rng(0)
    >>> # gene 0 = strong x-gradient (SV), gene 1 = noise (not SV)
    >>> M = np.vstack([coords[:,0] + rng.normal(0,0.1,100),
    ...                rng.normal(0.5,0.2,100)])
    >>> res = svapa(M, coords, k=8, n_perm=50, seed=0, min_obs=5)
    >>> res.iloc[0]['gene'] in (0, 'Gene_0')  # SV gene ranks first
    True
    """
    is_df = isinstance(apa_matrix, pd.DataFrame)
    if is_df:
        gene_names = list(apa_matrix.index)
        mat = apa_matrix.values.astype(float)
    else:
        mat = np.asarray(apa_matrix, dtype=float)
        gene_names = [f"Gene_{i}" for i in range(mat.shape[0])]
    coords = np.asarray(coords, dtype=float)
    if coords.ndim == 1:
        coords = coords.reshape(-1, 2)
    n_genes, n_spots = mat.shape
    if coords.shape[0] != n_spots:
        raise ValueError(
            f"apa_matrix has {n_spots} spots but coords has {coords.shape[0]}"
        )

    rng = np.random.default_rng(seed)
    rows = []
    for g in range(n_genes):
        v = mat[g]
        finite = np.isfinite(v)
        n = int(finite.sum())
        if n < max(min_obs, 3):
            continue
        I, pI = np.nan, np.nan
        C, pC = np.nan, np.nan
        # reuse a single weight build + permutation RNG per gene
        kk = min(k, n - 1)
        W = _build_spatial_weights(coords[finite], k=kk, symmetric=True)
        S0 = float(W.sum())
        vv = v[finite]
        z = vv - vv.mean()
        if S0 > 0:
            I = _morans_i_from_W(z, W, n, S0)
            if include_geary:
                C = _gearys_c_from_W(vv, W, n, S0)
        if n_perm and n_perm > 0:
            if np.isfinite(I):
                pI = _permutation_pvalue(I, vv, W, n, S0, n_perm,
                                         np.random.default_rng(rng.integers(2**31)),
                                         geary=False)
            if include_geary and np.isfinite(C):
                pC = _permutation_pvalue(C, vv, W, n, S0, n_perm,
                                         np.random.default_rng(rng.integers(2**31)),
                                         geary=True)
        rows.append((gene_names[g], I, pI, C, pC, n))

    out = pd.DataFrame(
        rows, columns=["gene", "morans_i", "morans_i_pvalue",
                       "gearys_c", "gearys_c_pvalue", "n_obs"]
    )

    if fdr and _HAS_SM and len(out):
        pvals = out["morans_i_pvalue"].fillna(1.0).values
        _, padj, _, _ = multipletests(pvals, method=fdr)
        out["padj"] = padj
    elif len(out):
        out["padj"] = out["morans_i_pvalue"]

    out = out.sort_values("morans_i", ascending=False, na_position="last")
    out = out.reset_index(drop=True)
    return out


__all__ = ["morans_i", "gearys_c", "svapa", "_build_spatial_weights"]
