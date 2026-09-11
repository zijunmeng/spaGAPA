"""Unit tests for Mondrian (group-conditional) split-conformal calibration.

Motivation: spaGAPA's marginal conformal guarantee undercovers specific
subgroups (high-expression bins reach only ~0.82 at a 90% target).  Mondrian
conformal calibration runs a separate split-conformal quantile *within each
predefined group*, restoring per-group marginal coverage while keeping the
distribution-free, model-agnostic character.

Tests use synthetic heteroscedastic data with a KNOWN failure mode of the
global calibrator (one wide-error group) and verify the Mondrian fix.
"""
from __future__ import annotations

import numpy as np
import pytest

from spagapa.imputation.calibration import (
    MondrianConformalCalibrator,
    bin_by_quantiles,
)


# ---------------------------------------------------------------------------
# fixtures: heteroscedastic two-group data where GLOBAL conformal fails one
# group (the wide-error one) and Mondrian restores it
# ---------------------------------------------------------------------------
def _heteroscedastic_two_group(n_per_group=4000, seed=0):
    """Return (errors, stds, groups, test_errors, test_stds, test_groups).

    Group "lo": |error| ~ halfnormal(scale=0.02)
    Group "hi": |error| ~ halfnormal(scale=0.10)   (5x wider)
    GP std is weakly informative and identical across groups, so the
    pooled/global quantile is dominated by the narrow group and undercovers
    the wide group.
    """
    rng = np.random.default_rng(seed)
    err = np.concatenate([
        np.abs(rng.normal(0, 0.02, n_per_group)),
        np.abs(rng.normal(0, 0.10, n_per_group)),
    ])
    std = np.full_like(err, 0.05)
    groups = np.array(["lo"] * n_per_group + ["hi"] * n_per_group)
    perm = rng.permutation(len(err))
    err, std, groups = err[perm], std[perm], groups[perm]

    test_err = np.concatenate([
        np.abs(rng.normal(0, 0.02, n_per_group)),
        np.abs(rng.normal(0, 0.10, n_per_group)),
    ])
    test_std = np.full_like(test_err, 0.05)
    test_groups = np.array(["lo"] * n_per_group + ["hi"] * n_per_group)
    return err, std, groups, test_err, test_std, test_groups


def _coverage(pred_std, q, err):
    """Coverage of symmetric interval half-width q * std against |err|."""
    return float(np.mean(err <= q * np.maximum(pred_std, 1e-12)))


# ---------------------------------------------------------------------------
# per-group coverage (the core Mondrian property)
# ---------------------------------------------------------------------------
def test_mondrian_restores_wide_group_coverage():
    cal_err, cal_std, cal_grp, te, ts, tg = _heteroscedastic_two_group()
    mon = MondrianConformalCalibrator(alpha=0.1, mode="global")
    mon.fit(cal_err, None, cal_grp)
    lower, upper = mon.predict(np.zeros_like(te), None, tg)
    half = (upper - lower) / 2
    for g in ("lo", "hi"):
        m = tg == g
        cov = float(np.mean(te[m] <= half[m]))
        # per-group marginal coverage at (or slightly above) the 90% target
        assert cov >= 0.885, f"group {g}: coverage {cov:.4f} < 0.885"


def test_global_calibrator_undercovers_wide_group_on_same_data():
    """The motivating failure: pooled/global q undercovers the wide group."""
    from spagapa.imputation.calibration import ConformalCalibrator

    cal_err, cal_std, cal_grp, te, ts, tg = _heteroscedastic_two_group()
    glob = ConformalCalibrator(alpha=0.1, mode="global")
    glob.fit(cal_err)
    half = glob.interval_.q_hat
    cov_hi = float(np.mean(te[tg == "hi"] <= half))
    # with equal group sizes and 5x scale separation the pooled q sits near
    # the midpoint of the two scales: wide group must undercover (< 0.88)
    assert cov_hi < 0.88, f"global unexpectedly covers wide group: {cov_hi}"


# ---------------------------------------------------------------------------
# finite-sample correction per group
# ---------------------------------------------------------------------------
def test_per_group_quantile_uses_finite_sample_correction():
    # deterministic 100-point group with known quantile
    err = np.abs(np.random.default_rng(1).normal(0, 1, 100))
    groups = np.array(["a"] * 60 + ["b"] * 40)
    err = np.concatenate([np.sort(np.abs(np.arange(60) * 0.01)) + 0.001,
                          np.sort(np.abs(np.arange(40) * 0.02)) + 0.001])
    mon = MondrianConformalCalibrator(alpha=0.1, mode="global")
    mon.fit(err, None, groups)
    # group a: n=60, ceil(61*0.9)/60 -> index ceil(54.9)=55 -> "higher" of
    # the 55th order statistic (0-based: position 54 of the sorted scores)
    qa_expected = np.quantile(err[:60], min(1.0, np.ceil(61 * 0.9) / 60),
                              method="higher")
    assert mon.interval_.group_q["a"] == pytest.approx(qa_expected)


# ---------------------------------------------------------------------------
# small-group fallback + unknown-group handling
# ---------------------------------------------------------------------------
def test_small_group_falls_back_to_pooled_quantile():
    rng = np.random.default_rng(2)
    big = np.abs(rng.normal(0, 0.05, 2000))
    tiny = np.abs(rng.normal(0, 0.05, 7))  # below min_group_size
    err = np.concatenate([big, tiny])
    groups = np.array(["big"] * 2000 + ["tiny"] * 7)
    mon = MondrianConformalCalibrator(alpha=0.1, mode="global", min_group_size=50)
    mon.fit(err, None, groups)
    q_pooled = np.quantile(err, min(1.0, np.ceil((len(err) + 1) * 0.9) / len(err)),
                           method="higher")
    assert mon.interval_.group_q["tiny"] == pytest.approx(float(q_pooled))
    assert "pooled" in mon.interval_.group_q


def test_unknown_group_at_predict_uses_pooled_fallback():
    cal_err, cal_std, cal_grp, *_ = _heteroscedastic_two_group(n_per_group=500)
    mon = MondrianConformalCalibrator(alpha=0.1, mode="global")
    mon.fit(cal_err, None, cal_grp)
    pred = np.array([0.0, 0.0])
    groups = np.array(["lo", "unseen_group"])
    lower, upper = mon.predict(pred, None, groups)
    assert np.isfinite(lower).all() and np.isfinite(upper).all()
    assert (upper - lower)[1] > 0


# ---------------------------------------------------------------------------
# locally-adaptive mode: per-group q scales the GP std
# ---------------------------------------------------------------------------
def test_locally_adaptive_mode_per_group_scaling():
    cal_err, cal_std, cal_grp, te, ts, tg = _heteroscedastic_two_group(
        n_per_group=2000
    )
    mon = MondrianConformalCalibrator(alpha=0.1, mode="locally_adaptive")
    mon.fit(cal_err, cal_std, cal_grp)
    lower, upper = mon.predict(np.zeros_like(te), ts, tg)
    half = (upper - lower) / 2
    # stds are constant here, so per-group half-widths must be group-constant
    for g in ("lo", "hi"):
        m = tg == g
        assert np.allclose(half[m], half[m][0])
    # coverage per group still valid in standardized space
    for g in ("lo", "hi"):
        m = tg == g
        cov = float(np.mean(te[m] <= half[m]))
        assert cov >= 0.885


# ---------------------------------------------------------------------------
# input validation
# ---------------------------------------------------------------------------
def test_fit_before_predict_raises():
    mon = MondrianConformalCalibrator(alpha=0.1)
    with pytest.raises(RuntimeError, match="fit"):
        mon.predict(np.zeros(3), None, np.array(["a"] * 3))


def test_alpha_validation():
    with pytest.raises(ValueError, match="alpha"):
        MondrianConformalCalibrator(alpha=1.5)


def test_group_length_mismatch_raises():
    mon = MondrianConformalCalibrator(alpha=0.1)
    with pytest.raises(ValueError, match="groups"):
        mon.fit(np.abs(np.random.default_rng(3).normal(0, 1, 10)), None,
                np.array(["a"] * 9))


# ---------------------------------------------------------------------------
# quantile binning helper
# ---------------------------------------------------------------------------
def test_bin_by_quantiles_produces_balanced_bins():
    rng = np.random.default_rng(4)
    x = rng.normal(0, 1, 10000)
    bins = bin_by_quantiles(x, n_bins=5)
    assert bins.shape == x.shape
    counts = np.bincount(bins.astype(int))
    assert counts.min() >= 9000 // 5 - 200  # roughly balanced
    assert counts.max() <= 2000 + 400
    # monotone: higher x never maps to a lower bin
    order = np.argsort(x)
    assert np.all(np.diff(bins[order]) >= 0)


def test_bin_by_quantiles_handles_ties_and_nan():
    x = np.array([1.0, 1.0, 1.0, 2.0, np.nan, 3.0])
    bins = bin_by_quantiles(x, n_bins=2)
    assert np.isnan(bins[4])
    assert bins[0] == bins[1] == bins[2]
    assert bins[5] >= bins[3]
