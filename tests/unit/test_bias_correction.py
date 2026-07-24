"""Unit tests for APA batch/protocol bias correction methods.

Both methods are exercised on synthetic data with a KNOWN injected batch
effect: a clean "biology" matrix is split into two batches, a shift is added
to one batch, and the corrector must (a) reduce the injected batch signal
and (b) recover the original biology (high corr between corrected and the
pre-batch truth). NaN handling is verified too.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from spagapa.analysis import (
    quantile_normalize,
    linear_batch_correction,
)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _make_synthetic(n_genes=60, n_spots=200, frac_missing=0.1, seed=0):
    """Return (truth matrix, batch_labels, missing mask).

    Biology and batch are constructed to be ORTHOGONAL: the per-gene spatial
    pattern is a smooth sinusoid across spots, but batch membership is
    assigned by *shuffling* (not by spot position), so removing the batch
    mean-difference does not touch the real spatial biology.  Without this,
    a method that correctly removes a batch effect can look like it is
    "hurting" recovery because the batch split happened to align with the
    biology gradient.
    """
    rng = np.random.default_rng(seed)
    # biology: per-gene spatial gradient + noise
    spots = np.linspace(0, 1, n_spots)
    base = (
        0.5
        + 0.3 * np.sin(rng.uniform(0, 2 * np.pi, n_genes)[:, None] * spots[None, :])
        + 0.05 * rng.standard_normal((n_genes, n_spots))
    )
    base = np.clip(base, 0.0, 1.0)
    mask = rng.random((n_genes, n_spots)) < frac_missing
    # half spots in batch A, half in batch B -- assigned by SHUFFLE so that
    # batch and biology are not collinear.
    spot_order = rng.permutation(n_spots)
    half = n_spots // 2
    batch = np.empty(n_spots, dtype=object)
    batch[spot_order[:half]] = "A"
    batch[spot_order[half:]] = "B"
    return base, batch, mask


def _inject_batch(truth, batch_labels, mask, shift=0.20, scale=1.30, seed=0):
    """Add a multiplicative+additive batch effect on batch B only.

    Parameters
    ----------
    truth : ndarray (n_genes, n_spots)
        Clean biology matrix.
    batch_labels : sequence of length n_spots
        Batch label per spot. The effect is injected on every spot whose
        label is the (alphabetically) last unique batch level ("B").
    mask : ndarray (n_genes, n_spots) bool
        Missing-value mask (True = missing).
    shift, scale : float
        Additive and multiplicative batch magnitude on the affected batch.

    Note: NO clipping to [0,1] is applied after injection -- clipping would
    destroy information in the affected batch (saturate values >1) and make
    recovery mechanically impossible. The truth is already in a sane [0,1]
    range and the batch-correction literature does not clip; we keep the
    injected values as-is so recovery is well-posed.
    """
    batch_labels = np.asarray(batch_labels)
    levels = sorted(np.unique(batch_labels))
    affected = levels[-1]  # "B"
    b_mask = batch_labels == affected
    contaminated = truth.copy()
    contaminated[:, b_mask] = contaminated[:, b_mask] * scale + shift
    contaminated_with_nan = contaminated.copy()
    contaminated_with_nan[mask] = np.nan
    return contaminated_with_nan, contaminated, b_mask


def _corr_robust(a, b):
    """Pearson correlation over jointly finite entries (flattened)."""
    a = np.asarray(a, dtype=float).ravel()
    b = np.asarray(b, dtype=float).ravel()
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 5:
        return np.nan
    a = a[ok] - a[ok].mean()
    b = b[ok] - b[ok].mean()
    denom = np.sqrt((a * a).sum() * (b * b).sum())
    if denom == 0:
        return np.nan
    return float((a * b).sum() / denom)


# ---------------------------------------------------------------------------
# shape / NaN / no-op tests
# ---------------------------------------------------------------------------
class TestQuantileNormalizeBasics:
    def test_preserves_shape_and_nan(self):
        truth, batch, mask = _make_synthetic()
        contaminated, _, _ = _inject_batch(truth, batch, mask)
        out = quantile_normalize(contaminated, batch)
        assert out.shape == contaminated.shape
        assert np.all(np.isnan(out[mask]))

    def test_dataframe_preserves_index(self):
        truth, batch, mask = _make_synthetic(n_genes=10, n_spots=40)
        contaminated, _, _ = _inject_batch(truth, batch, mask)
        df = pd.DataFrame(contaminated,
                          index=[f"g{i}" for i in range(contaminated.shape[0])],
                          columns=[f"s{j}" for j in range(contaminated.shape[1])])
        out = quantile_normalize(df, batch)
        assert isinstance(out, pd.DataFrame)
        assert list(out.index) == list(df.index)
        assert list(out.columns) == list(df.columns)

    def test_single_group_noop(self):
        truth, _, mask = _make_synthetic()
        truth_nan = truth.copy()
        truth_nan[mask] = np.nan
        out = quantile_normalize(truth_nan, ["X"] * truth.shape[1])
        # single batch: ranks map back to the same pooled distribution
        assert _corr_robust(out, truth_nan) > 0.999


class TestLinearBatchCorrectionBasics:
    def test_preserves_shape_and_nan(self):
        truth, batch, mask = _make_synthetic()
        contaminated, _, _ = _inject_batch(truth, batch, mask)
        out = linear_batch_correction(contaminated, batch)
        assert out.shape == contaminated.shape
        assert np.all(np.isnan(out[mask]))

    def test_single_batch_noop(self):
        truth, batch, mask = _make_synthetic()
        truth_nan = truth.copy()
        truth_nan[mask] = np.nan
        out = linear_batch_correction(truth_nan, ["A"] * truth.shape[1])
        assert _corr_robust(out, truth_nan) > 0.999


# ---------------------------------------------------------------------------
# recovery tests (the real ones)
# ---------------------------------------------------------------------------
class TestBatchRecovery:
    """Inject a known batch effect, require both methods to recover biology."""

    @pytest.mark.parametrize("shift,scale", [(0.20, 1.30), (0.10, 1.15)])
    def test_quantile_normalize_recovers_biology(self, shift, scale):
        truth, batch, mask = _make_synthetic(seed=1)
        contaminated, _, b_mask = _inject_batch(
            truth, batch, mask, shift=shift, scale=scale, seed=1)
        truth_nan = np.where(mask, np.nan, truth)
        # contaminated is far from truth
        base_corr = _corr_robust(contaminated, truth_nan)
        out = quantile_normalize(contaminated, batch)
        corr = _corr_robust(out, truth_nan)
        # recovery: corr must be substantially higher than contaminated
        assert corr > base_corr + 0.05
        assert corr > 0.85, f"QN recovery corr={corr:.3f} too low"

        # batch signal: per-gene mean difference A vs B should shrink
        def batch_mean_diff(m):
            a = np.nanmean(m[:, ~b_mask], axis=1)
            b = np.nanmean(m[:, b_mask], axis=1)
            return float(np.mean(np.abs(a - b)))
        assert batch_mean_diff(out) < batch_mean_diff(contaminated)

    @pytest.mark.parametrize("shift,scale", [(0.20, 1.30), (0.10, 1.15)])
    def test_linear_batch_correction_recovers_biology(self, shift, scale):
        truth, batch, mask = _make_synthetic(seed=2)
        contaminated, _, b_mask = _inject_batch(
            truth, batch, mask, shift=shift, scale=scale, seed=2)
        truth_nan = np.where(mask, np.nan, truth)
        base_corr = _corr_robust(contaminated, truth_nan)
        out = linear_batch_correction(contaminated, batch)
        corr = _corr_robust(out, truth_nan)
        assert corr > base_corr + 0.05
        # linear method should recover an additive+multiplicative batch well
        assert corr > 0.90, f"linear recovery corr={corr:.3f} too low"

        def batch_mean_diff(m):
            a = np.nanmean(m[:, ~b_mask], axis=1)
            b = np.nanmean(m[:, b_mask], axis=1)
            return float(np.mean(np.abs(a - b)))
        assert batch_mean_diff(out) < batch_mean_diff(contaminated)

    def test_linear_preserves_biological_covariate(self):
        """Batch correction must NOT remove a biological (condition) signal."""
        rng = np.random.default_rng(3)
        n_genes, n_spots = 40, 200
        spots = np.linspace(0, 1, n_spots)
        truth = 0.5 + 0.2 * np.sin(spots[None, :] * 3 + rng.uniform(0, 6, n_genes)[:, None])
        truth = np.clip(truth, 0, 1)
        # orthogonal design: shuffle both batch and condition assignments
        spot_order = rng.permutation(n_spots)
        batch = np.empty(n_spots, dtype=object)
        condition = np.empty(n_spots, dtype=object)
        batch[spot_order[:100]] = "A"
        batch[spot_order[100:]] = "B"
        # condition independent of batch
        cond_order = rng.permutation(n_spots)
        condition[cond_order[:100]] = "ctrl"
        condition[cond_order[100:]] = "ad"
        # add a condition effect (biology we want to KEEP)
        biology = truth.copy()
        biology[:, condition == "ad"] += 0.35
        # then add a batch effect (batch B gets a shift)
        contaminated = biology.copy()
        contaminated[:, batch == "B"] += 0.25

        out = linear_batch_correction(contaminated, batch, preserve_labels=condition)

        # biology preserved: condition separation on AD vs ctrl must survive
        def cond_diff(m):
            return float(np.nanmean(m[:, condition == "ad"])
                         - np.nanmean(m[:, condition == "ctrl"]))
        d_before = cond_diff(biology)  # the biology we want to keep
        d_after = cond_diff(out)
        # the correction should not erase the condition difference
        assert d_after > 0.5 * d_before, (
            f"linear correction erased biology: d_before={d_before:.3f} "
            f"d_after={d_after:.3f}")
