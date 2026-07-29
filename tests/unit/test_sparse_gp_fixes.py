"""
Unit tests for the Phase-3 GP over-smoothing / weak-uncertainty fixes:
  * length_scale_multiplier -> changes the effective length scale
  * local_noise -> changes the predictive std scaling so uncertainty tracks
    local noise rather than a flat constant.

These tests are backward-compatible: the defaults reproduce the historical
behaviour (multiplier=5, local_noise=False).
"""
import numpy as np
import pytest

from spagapa.imputation import SparseGPImputer


def _grid_coords(side=15, jitter=0.05, seed=0):
    """A jittered square grid of shape (side*side, 2)."""
    rng = np.random.default_rng(seed)
    xs, ys = np.meshgrid(np.arange(side), np.arange(side))
    coords = np.column_stack([xs.ravel(), ys.ravel()]).astype(float)
    coords += rng.normal(0, jitter, coords.shape)
    return coords


class TestLengthScaleMultiplier:
    def test_multiplier_default_is_none(self):
        """Default length_scale_multiplier=None means 'use length_scale
        directly' (backward-compatible).  local_noise defaults off."""
        imp = SparseGPImputer()
        assert imp.length_scale_multiplier is None
        assert imp.local_noise is False

    def test_multiplier_changes_effective_length_scale(self):
        """When median_nn_dist is supplied, eff_length = nn_dist * multiplier."""
        coords = _grid_coords()
        values = coords[:, 0]  # x-gradient
        nn_dist = 1.0  # arbitrary, supplied by caller
        imp2 = SparseGPImputer(length_scale_multiplier=2.0)
        imp2.fit(coords, values, median_nn_dist=nn_dist)
        imp5 = SparseGPImputer(length_scale_multiplier=5.0)
        imp5.fit(coords, values, median_nn_dist=nn_dist)
        assert imp2._effective_length_scale == pytest.approx(2.0)
        assert imp5._effective_length_scale == pytest.approx(5.0)
        assert imp5._effective_length_scale > imp2._effective_length_scale

    def test_explicit_length_scale_wins_without_median_nn(self):
        """No median_nn_dist -> eff_length falls back to self.length_scale."""
        coords = _grid_coords()
        values = coords[:, 0]
        imp = SparseGPImputer(length_scale=3.0, length_scale_multiplier=2.0)
        imp.fit(coords, values)
        # multiplier is ignored because no median_nn_dist was supplied
        assert imp._effective_length_scale == pytest.approx(3.0)

    def test_lower_multiplier_changes_predictions(self):
        """The multiplier mechanically sharpens the kernel: predictions from
        multiplier=1 and multiplier=20 on the same data must differ, and the
        high-multiplier predictions are flatter (smaller dynamic range) than
        the low-multiplier ones.  This is the lever the benchmark uses to
        combat over-smoothing (validated empirically on GSE183456/GSE220442).
        """
        rng = np.random.default_rng(2)
        coords = _grid_coords(side=20, jitter=0.0)
        # noisy spatial field: a quadrant blob + noise
        n = coords.shape[0]
        truth = np.zeros(n)
        blob = (coords[:, 0] >= 10) & (coords[:, 1] >= 10)
        truth[blob] = 1.0
        truth += rng.normal(0, 0.3, n)
        full_mask = np.ones(n, dtype=bool)

        def fit_pred(mult):
            imp = SparseGPImputer(n_inducing=80, length_scale_multiplier=mult)
            imp.fit(coords, truth, mask=full_mask, median_nn_dist=1.0)
            pred, _ = imp.predict(coords)
            return pred

        p1 = fit_pred(1.0)
        p20 = fit_pred(20.0)
        # predictions differ (multiplier changed the kernel)
        assert not np.allclose(p1, p20)
        # high-multiplier predictions are flatter (over-smoothed toward mean)
        assert (p20.max() - p20.min()) <= (p1.max() - p1.min()) + 1e-9


class TestLocalNoise:
    def test_local_noise_changes_effective_noise(self):
        """local_noise=True replaces flat 0.1 with a per-gene estimate."""
        rng = np.random.default_rng(3)
        coords = _grid_coords()
        values = coords[:, 0] + rng.normal(0, 0.5, coords.shape[0])

        imp_flat = SparseGPImputer(local_noise=False, noise_level=0.1)
        imp_flat.fit(coords, values)
        assert imp_flat._effective_noise == pytest.approx(0.1)

        imp_local = SparseGPImputer(local_noise=True, noise_level=0.1)
        imp_local.fit(coords, values)
        # Local estimate is data-driven and should NOT equal the flat 0.1.
        assert imp_local._effective_noise != pytest.approx(0.1)
        assert imp_local._effective_noise > 0

    def test_local_noise_scales_std_with_gene_noise(self):
        """Predictive std should scale with the gene's noise level.

        A noisy gene (large local variance) should get HIGHER predictive std
        than a quiet gene (tiny local variance), when local_noise=True.
        """
        coords = _grid_coords(side=15)
        rng = np.random.default_rng(4)
        # same spatial trend, different noise floors
        quiet = coords[:, 0] + rng.normal(0, 0.01, coords.shape[0])
        noisy = coords[:, 0] + rng.normal(0, 2.0, coords.shape[0])

        def fit_std(values, local_noise):
            imp = SparseGPImputer(n_inducing=30, local_noise=local_noise)
            imp.fit(coords, values)
            _, std = imp.predict(coords, return_std=True)
            return float(np.median(std))

        std_quiet_local = fit_std(quiet, True)
        std_noisy_local = fit_std(noisy, True)
        # With local noise, the noisy gene should be more uncertain.
        assert std_noisy_local > std_quiet_local

    def test_flat_noise_gives_value_independent_std(self):
        """With flat noise, predictive std depends only on coordinates.

        The std is sqrt(1 - Q_* + noise) where Q_* depends only on the
        test coordinates and inducing points (NOT on gene values), so two
        genes on the SAME coordinates AND SAME training mask must get
        identical std scaling.  (An explicit all-true mask is supplied so
        the value->mask default does not perturb the inducing points.)
        """
        coords = _grid_coords(side=12)
        rng = np.random.default_rng(5)
        quiet = coords[:, 0] + rng.normal(0, 0.01, coords.shape[0])
        noisy = coords[:, 0] + rng.normal(0, 2.0, coords.shape[0])
        full_mask = np.ones(coords.shape[0], dtype=bool)

        def fit_std(values, local_noise):
            imp = SparseGPImputer(n_inducing=25, local_noise=local_noise)
            imp.fit(coords, values, mask=full_mask)
            _, std = imp.predict(coords, return_std=True)
            return std

        std_quiet_flat = fit_std(quiet, False)
        std_noisy_flat = fit_std(noisy, False)
        np.testing.assert_allclose(std_quiet_flat, std_noisy_flat, rtol=1e-9)

    def test_local_noise_few_points_falls_back(self):
        """Very few observations should not crash the local estimate."""
        coords = np.array([[0., 0.], [1., 0.], [2., 0.], [3., 0.]])
        values = np.array([0.5, 0.6, 0.4, 0.55])
        imp = SparseGPImputer(n_inducing=3, local_noise=True)
        imp.fit(coords, values)
        assert imp._effective_noise > 0
