"""Unit tests for count-native (binomial-likelihood) usage GP.

Motivation: spaGAPA models usage FRACTIONS with a Gaussian GP.  The
distal-usage index is a bounded, bimodal binomial proportion (k_distal of
n_total site counts) — the Gaussian likelihood is misspecified, which is
half the reason the per-gene MEAN baseline wins entry-wise RMSE.  The
count-native layer models the generating likelihood: empirical-logit
transform + delta-method binomial variance (heteroscedastic, per entry) +
sparse GP in logit space + sigmoid back-transform.
"""
from __future__ import annotations

import numpy as np
import pytest

from spagapa.imputation.count_gp import (
    empirical_logit,
    binomial_logit_var,
    BinomialUsageGP,
)


# ---------------------------------------------------------------------------
# pure transforms
# ---------------------------------------------------------------------------
def test_empirical_logit_known_values():
    # k=0 -> logit(0.5/(n+1)) = log(0.5 / (n+0.5)); k=n -> log((n+0.5)/0.5)
    n = 10
    assert empirical_logit(0, n) == pytest.approx(np.log(0.5 / 10.5))
    assert empirical_logit(n, n) == pytest.approx(np.log(10.5 / 0.5))
    # symmetric midpoint: p_hat = 0.5 -> logit 0
    assert empirical_logit(5, 10) == pytest.approx(0.0, abs=1e-12)


def test_binomial_logit_var_properties():
    # delta-method formula 1/(k+0.5) + 1/(n-k+0.5)
    k, n = 3, 20
    assert binomial_logit_var(k, n) == pytest.approx(1 / 3.5 + 1 / 17.5)
    # symmetry k <-> n-k
    assert binomial_logit_var(3, 20) == pytest.approx(binomial_logit_var(17, 20))
    # more data -> less variance
    assert binomial_logit_var(30, 200) < binomial_logit_var(3, 20)
    # vectorized
    v = binomial_logit_var(np.array([0, 5, 10]), np.array([10, 10, 10]))
    assert v.shape == (3,)


# ---------------------------------------------------------------------------
# synthetic recovery: smooth spatial usage field with binomial sampling
# ---------------------------------------------------------------------------
def _synthetic_field(grid=28, seed=0, mean_depth=30):
    """Smooth sigmoid field p(x); counts k ~ Binomial(n, p) with Poisson n."""
    rng = np.random.default_rng(seed)
    x = np.linspace(0, 1, grid)
    xx, yy = np.meshgrid(x, x)
    coords = np.column_stack([xx.ravel(), yy.ravel()])
    # smooth latent gradient field
    theta = 2.0 * (xx.ravel() - 0.5) + 1.0 * np.sin(3 * np.pi * yy.ravel())
    p = 1 / (1 + np.exp(-theta))
    n_tot = rng.poisson(mean_depth, size=coords.shape[0]) + 2
    k = rng.binomial(n_tot, p)
    return coords, k, n_tot, p


def test_count_gp_beats_mean_baseline_on_smooth_field():
    """THE core claim: with a binomial generative process and a spatially
    smooth usage field, the count-native GP recovers usage better than the
    per-gene mean of observed fractions."""
    rng = np.random.default_rng(7)
    n_genes = 12
    rmses = {"count_gp": [], "mean": []}
    for g in range(n_genes):
        coords, k, n_tot, p_true = _synthetic_field(seed=g)
        # mask 40% of entries (no counts observed there)
        m = rng.random(len(k)) < 0.4
        k_obs, n_obs = k.copy(), n_tot.copy()
        k_obs[m], n_obs[m] = 0, 0

        gp = BinomialUsageGP(n_inducing=80, length_scale=0.15, total_floor=2)
        gp.fit(coords, k_obs[None, :], n_obs[None, :])
        p_hat = gp.predict_usage()[0]

        frac_obs = np.where(n_obs > 0, k_obs / np.maximum(n_obs, 1), np.nan)
        mean_hat = np.nanmean(frac_obs)
        rmses["count_gp"].append(np.sqrt(np.mean((p_hat[m] - p_true[m]) ** 2)))
        rmses["mean"].append(
            np.sqrt(np.mean((mean_hat - p_true[m]) ** 2))
        )
    r_gp, r_mean = np.mean(rmses["count_gp"]), np.mean(rmses["mean"])
    assert r_gp < r_mean, (
        f"count-native GP ({r_gp:.4f}) failed to beat mean baseline "
        f"({r_mean:.4f})"
    )


def test_usage_output_bounded_and_std_positive():
    coords, k, n_tot, p_true = _synthetic_field(seed=1)
    gp = BinomialUsageGP(n_inducing=60, length_scale=0.15)
    gp.fit(coords, k[None, :], n_tot[None, :])
    p, p_std = gp.predict_usage(return_std=True)      # (1, n_spots) each
    p_hat, p_std = p[0], p_std[0]
    assert p_hat.shape == (len(k),)
    assert np.all(p_hat > 0) and np.all(p_hat < 1)
    assert np.all(p_std >= 0) and np.all(np.isfinite(p_std))


def test_low_coverage_spots_get_larger_uncertainty():
    """Depth-aware uncertainty: entries with fewer total counts should be
    assigned larger predictive std (delta method inflates through the
    binomial variance)."""
    rng = np.random.default_rng(3)
    grid = 24
    x = np.linspace(0, 1, grid)
    xx, yy = np.meshgrid(x, x)
    coords = np.column_stack([xx.ravel(), yy.ravel()])
    p = np.full(len(coords), 0.3)
    half = len(coords) // 2
    n_tot = np.concatenate([np.full(half, 4), np.full(len(coords) - half, 60)])
    k = rng.binomial(n_tot, p)
    # mask 30% of entries: uncertainty is a property of IMPUTED values
    miss = rng.random(len(k)) < 0.3
    k_m, n_m = k.copy(), n_tot.copy()
    k_m[miss], n_m[miss] = 0, 0
    gp = BinomialUsageGP(n_inducing=60, length_scale=0.3)
    gp.fit(coords, k_m[None, :], n_m[None, :])
    _, p_std = gp.predict_usage(return_std=True)
    p_std = p_std[0]
    from scipy.stats import spearmanr
    rho = spearmanr(n_tot[miss], p_std[miss]).statistic
    assert rho < -0.3, f"expected negative depth-std correlation, got {rho:.3f}"



def test_all_missing_and_zero_total_entries_are_ignored():
    coords, k, n_tot, p = _synthetic_field(seed=5, grid=12)
    k2, n2 = k.copy(), n_tot.copy()
    k2[:10], n2[:10] = 0, 0  # fully unobserved entries
    gp = BinomialUsageGP(n_inducing=30, length_scale=0.2)
    gp.fit(coords, k2[None, :], n2[None, :])
    p_hat = gp.predict_usage()[0]
    assert np.all(np.isfinite(p_hat))


def test_fit_rejects_mismatched_shapes():
    coords = np.random.default_rng(0).random((50, 2))
    with pytest.raises(ValueError, match="shape"):
        BinomialUsageGP(n_inducing=20).fit(
            coords, np.zeros((1, 50)), np.zeros((1, 40))
        )
