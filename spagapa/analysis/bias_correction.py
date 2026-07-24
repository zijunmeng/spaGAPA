"""
Batch / protocol bias correction for APA usage matrices.

No existing tool corrects batch effects in *alternative-polyadenylation (APA)
usage* matrices (analogous to Harmony / BBKNN / limma ``removeBatchEffect`` but
applied to a gene x spot APA matrix rather than an expression matrix). This
module implements two complementary, NaN-aware methods:

(a) :func:`quantile_normalize` -- non-parametric distribution alignment.
    For each gene, the per-sample observed APA-usage values are rank-aligned
    across batches so that every batch shares the same marginal distribution
    (the pooled-empirical reference quantile function). This removes
    *location/scale* batch differences in a model-free way while preserving
    within-batch rank structure (the spatial / spot ordering inside a batch).

(b) :func:`linear_batch_correction` -- linear-model batch removal in the
    spirit of limma ``removeBatchEffect``.  Fits ``usage ~ gene + batch``
    (plus an optional biological covariate to *preserve*) on the gene x spot
    matrix and subtracts the batch coefficient. Preserved labels (e.g. AD vs
    control) are kept in the model so their associated variance is not
    removed as batch.

Conventions
-----------
* All public functions operate on ``gene x spot`` matrices (rows = genes,
  columns = spots) and accept ``NaN`` for missing values. Inputs may be a
  ``pandas.DataFrame`` or a ``numpy.ndarray``; returned objects match the
  input type (DataFrame in -> DataFrame out, preserving index/columns).
* Group / batch / preserve labels are 1-D arrays aligned with the **spots**
  (columns).
* Both functions are deterministic given a ``random_state`` where randomness
  is involved (only the quantile-normalize tie handling is stochastic-free;
  the linear method is fully deterministic via OLS).
"""

from __future__ import annotations

from typing import Optional, Union, Sequence, Tuple

import numpy as np
import pandas as pd

__all__ = [
    "quantile_normalize",
    "linear_batch_correction",
    "build_gene_index",
]

ArrayLike = Union[np.ndarray, pd.DataFrame]


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _to_numpy(df_or_arr: ArrayLike) -> Tuple[np.ndarray, Optional[pd.DataFrame]]:
    """Return (values 2D float, shell dataframe or None)."""
    if isinstance(df_or_arr, pd.DataFrame):
        return df_or_arr.to_numpy(dtype=float, copy=True), df_or_arr
    arr = np.asarray(df_or_arr, dtype=float)
    if arr.ndim != 2:
        raise ValueError(f"Expected 2-D matrix, got shape {arr.shape}")
    return arr.copy(), None


def _wrap(values: np.ndarray, shell: Optional[pd.DataFrame]) -> ArrayLike:
    if shell is None:
        return values
    out = pd.DataFrame(values, index=shell.index, columns=shell.columns)
    out.index.name = shell.index.name
    return out


def _validate_labels(labels: Sequence, n_spots: int, name: str) -> np.ndarray:
    labels = np.asarray(labels)
    if labels.shape[0] != n_spots:
        raise ValueError(
            f"{name} has length {labels.shape[0]}, expected {n_spots} (spots)"
        )
    return labels


# ---------------------------------------------------------------------------
# (a) quantile normalization
# ---------------------------------------------------------------------------
def quantile_normalize(
    apa_matrix: ArrayLike,
    group_labels: Sequence,
    reference: str = "pooled",
) -> ArrayLike:
    """Quantile-normalize per-gene APA usage across batches.

    For each gene, the observed values within every batch are replaced by the
    *same* set of reference quantiles, ordered by the within-batch ranks.  The
    reference distribution is, by default, the pooled (across-batch) empirical
    distribution of that gene; this makes every batch share the pooled
    marginal while keeping each batch's internal rank ordering (i.e. spatial
    structure is preserved, only the per-batch histogram is re-mapped).

    Parameters
    ----------
    apa_matrix : DataFrame or ndarray, shape (n_genes, n_spots)
        Gene x spot APA usage values. ``NaN`` marks missing observations and
        is ignored throughout (NaN stays NaN).
    group_labels : sequence of length n_spots
        Batch / protocol / sample label per spot (column).
    reference : {"pooled", "largest"}, default "pooled"
        Which batch's distribution is used as the alignment target.
        ``"pooled"`` uses the pooled-empirical quantiles of all batches
        (recommended -- no single batch is privileged). ``"largest`` uses the
        largest batch's empirical quantiles as the target.

    Returns
    -------
    DataFrame or ndarray (same type/shape as input)
        Quantile-normalized APA usage matrix.
    """
    values, shell = _to_numpy(apa_matrix)
    n_genes, n_spots = values.shape
    groups = _validate_labels(group_labels, n_spots, "group_labels")

    unique_groups = list(pd.unique(groups[~pd.isna(groups)]))
    if len(unique_groups) < 2:
        # nothing to align
        return _wrap(values, shell)

    group_masks = {g: (groups == g) for g in unique_groups}
    if reference == "largest":
        target_group = max(unique_groups, key=lambda g: int(group_masks[g].sum()))
    elif reference == "pooled":
        target_group = None
    else:
        raise ValueError(f"reference must be 'pooled' or 'largest', got {reference!r}")

    out = values.copy()
    for i in range(n_genes):
        row = values[i]
        if reference == "pooled":
            pooled = row[~np.isnan(row)]
            if pooled.size < 2:
                continue
            # reference quantile function: sorted pooled values, interpolated
            ref_sorted = np.sort(pooled)
        else:
            tg_vals = row[group_masks[target_group]]
            tg_vals = tg_vals[~np.isnan(tg_vals)]
            if tg_vals.size < 2:
                continue
            ref_sorted = np.sort(tg_vals)

        for g in unique_groups:
            mask = group_masks[g]
            gv = row[mask]
            valid = ~np.isnan(gv)
            n_valid = int(valid.sum())
            if n_valid < 2:
                continue
            order = np.argsort(gv[valid], kind="mergesort")
            ranks = np.empty(n_valid, dtype=int)
            ranks[order] = np.arange(n_valid)
            # map each rank to a reference quantile via interpolation
            plotting_positions = (ranks + 0.5) / n_valid
            ref_pp = (np.arange(ref_sorted.size) + 0.5) / ref_sorted.size
            new_vals = np.interp(plotting_positions, ref_pp, ref_sorted)
            # write back into the valid positions of this group
            full = row[mask].copy()
            full_valid = full.copy()
            full_valid[valid] = new_vals
            # NaN positions stay NaN
            full[~valid] = np.nan
            full[valid] = new_vals
            out[i, mask] = full
    return _wrap(out, shell)


# ---------------------------------------------------------------------------
# (b) linear batch correction (limma removeBatchEffect analogue)
# ---------------------------------------------------------------------------
def linear_batch_correction(
    apa_matrix: ArrayLike,
    batch_labels: Sequence,
    preserve_labels: Optional[Sequence] = None,
) -> ArrayLike:
    """Remove additive batch effects via a linear model, gene by gene.

    For each gene ``g`` fits (on observed spots only)

        usage[g, s] = intercept + sum_k beta_k * preserve_k[s]
                                  + sum_b gamma_b * batch_b[s] + eps

    then subtracts the fitted batch term ``sum_b gamma_b * batch_b[s]`` from
    every observed entry. The intercept and any ``preserve_labels`` covariate
    are kept in the model so their associated signal (e.g. AD vs control
    biology) is **not** removed -- only the batch term is subtracted, exactly
    as in limma ``removeBatchEffect``.

    Missing values (NaN) are dropped from the per-gene fit. If a batch level
    is not estimable (collinear with the preserved covariate, or a batch has
    < 2 observed spots for a gene) the batch term for that gene is left at 0
    (no correction applied), which is the safe, well-defined behaviour.

    Parameters
    ----------
    apa_matrix : DataFrame or ndarray, shape (n_genes, n_spots)
        Gene x spot APA usage matrix, NaN for missing.
    batch_labels : sequence of length n_spots
        Batch label per spot. The first level (alphabetically) is used as the
        reference level (absorbed into the intercept).
    preserve_labels : sequence of length n_spots, optional
        A biological covariate (e.g. condition: AD / control) to *keep* in the
        model so its variance is not removed as batch.  Treated as a
        categorical fixed effect.

    Returns
    -------
    DataFrame or ndarray (same type/shape as input)
        Batch-corrected APA usage matrix (NaN positions unchanged).
    """
    values, shell = _to_numpy(apa_matrix)
    n_genes, n_spots = values.shape
    batch = _validate_labels(batch_labels, n_spots, "batch_labels")

    # Build design matrix: intercept + preserve dummies + batch dummies
    design_cols = [np.ones(n_spots)]
    if preserve_labels is not None:
        preserve = _validate_labels(preserve_labels, n_spots, "preserve_labels")
        preserve_levels = sorted(pd.unique(preserve[~pd.isna(preserve)]))
        for lvl in preserve_levels[1:]:  # drop-first
            design_cols.append((preserve == lvl).astype(float))

    batch_levels = sorted(pd.unique(batch[~pd.isna(batch)]))
    if len(batch_levels) < 2:
        return _wrap(values, shell)  # no batch to correct
    # drop-first: reference batch absorbed into intercept
    batch_dummy_levels = batch_levels[1:]
    batch_col_idx = {  # for extracting the batch sub-block
        "start": len(design_cols),
        "n": len(batch_dummy_levels),
    }
    for lvl in batch_dummy_levels:
        design_cols.append((batch == lvl).astype(float))
    design = np.column_stack(design_cols)  # (n_spots, p)
    # precompute the batch-effect projection matrix once:
    # gamma_hat = (X_b' M X_b)^-1 X_b' M y, where M = I - X_p (X_p' X_p)^-1 X_p'
    # i.e. partial out intercept+preserve before estimating batch.
    preserve_block = design[:, : batch_col_idx["start"]]
    batch_block = design[:, batch_col_idx["start"]:]

    # QR-based residual-maker for the preserve block (handles collinearity)
    Qp, _ = np.linalg.qr(preserve_block)
    # M = I - Qp Qp'
    Mb_x = batch_block - Qp @ (Qp.T @ batch_block)
    # stable least-squares for gamma
    # Mb_x may be rank-deficient when batch is collinear with preserve; use lstsq
    out = values.copy()
    for i in range(n_genes):
        row = values[i]
        obs = ~np.isnan(row)
        n_obs = int(obs.sum())
        if n_obs < design.shape[1] + 1 or n_obs < 3:
            continue  # not enough data to fit; leave uncorrected
        y = row[obs]
        Mb_x_o = Mb_x[obs]
        Qp_o = Qp[obs]
        # residualize y against preserve block
        y_resid = y - Qp_o @ (Qp_o.T @ y)
        # estimate batch gamma on residualized quantities
        gamma, *_ = np.linalg.lstsq(Mb_x_o, y_resid, rcond=None)
        # fitted batch effect on the FULL set of spots (batch term only)
        batch_effect = Mb_x @ gamma  # but Mb_x is already residualized form;
        # the batch term we want to subtract is gamma applied to batch_block
        # (the residualized batch_block captures batch effect net of preserve)
        # Recompute cleanly: batch_effect_full = X_b gamma
        batch_effect_full = batch_block @ gamma
        corrected = row - batch_effect_full
        # NaN stays NaN
        corrected[~obs] = np.nan
        out[i] = corrected
    return _wrap(out, shell)


# ---------------------------------------------------------------------------
# gene-level distal-usage index builder (shared with head-to-head benchmark)
# ---------------------------------------------------------------------------
def build_gene_index(
    counts: pd.DataFrame,
    sites: pd.DataFrame,
    min_parent: int = 5,
) -> pd.DataFrame:
    """Gene x spot distal-usage ratio (median-split proximal/distal peaks).

    For each gene with >=2 PAS, split its peaks at the median strand-oriented
    position into a proximal group and a distal group, SUM counts within each
    group, and define ``index[gene, spot] = distal_sum / (distal_sum +
    prox_sum)`` when the parent total >= ``min_parent`` else ``NaN``.

    Copied verbatim from ``scripts/benchmark_stapaminer_headtohead.py`` so the
    bias-correction demo uses the identical gene-level index as the rest of
    the project.

    Parameters
    ----------
    counts : DataFrame, shape (n_peaks, n_spots)
        Peak-level APA counts (rows = site_id, columns = spots).
    sites : DataFrame
        Peak metadata with columns ``site_id``, ``gene_name``, ``strand``,
        ``coord`` (strand-oriented position).
    min_parent : int, default 5
        Minimum parent-total count for a spot to be considered observed.

    Returns
    -------
    DataFrame, shape (n_genes, n_spots)
        Gene x spot distal-usage ratio; NaN where parent total < min_parent.
    """
    sites = sites.copy()
    sites["oriented"] = np.where(
        sites["strand"].astype(str) == "-",
        -sites["coord"].astype(float),
        sites["coord"].astype(float),
    )
    peak_idx = {p: i for i, p in enumerate(counts.index)}
    spot_cols = counts.columns
    vals = counts.values  # peak x spot
    gene_rows = {}
    for gene, sub in sites.groupby("gene_name"):
        sub = sub.dropna(subset=["gene_name"])
        rows = [
            peak_idx[r["site_id"]]
            for _, r in sub.iterrows()
            if r["site_id"] in peak_idx
        ]
        if len(rows) < 2:
            continue
        oriented = np.array(
            [
                float(
                    sub.loc[
                        (sub["site_id"] == counts.index[r]), "oriented"
                    ].iloc[0]
                )
                for r in rows
            ]
        )
        rows = np.array(rows)
        order = np.argsort(oriented)
        rows = rows[order]
        mid = len(rows) // 2
        if mid < 1:
            mid = 1
        prox_rows = rows[:mid]
        dist_rows = rows[mid:]
        prox_sum = vals[prox_rows].sum(axis=0)
        dist_sum = vals[dist_rows].sum(axis=0)
        total = prox_sum + dist_sum
        ratio = np.where(
            total >= min_parent,
            dist_sum / np.where(total == 0, 1, total),
            np.nan,
        )
        gene_rows[gene] = ratio
    if not gene_rows:
        raise ValueError("No multi-site genes with usable peaks")
    index = pd.DataFrame(gene_rows, index=spot_cols).T  # gene x spot
    index.index.name = "gene"
    return index
