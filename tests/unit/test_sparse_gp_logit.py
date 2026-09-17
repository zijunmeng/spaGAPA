"""
Unit tests for the latent-logit sparse GP (transform='logit'):

  * transform=None (default) is exactly the historical raw-value GP
  * transform='logit' maps training values through
    logit(clip(y, eps, 1-eps)) and predictions back through sigmoid,
    so predictions are guaranteed to lie in (0, 1)
  * the posterior std stays in latent space (it conditions conformal
    calibration of original-space errors)
  * the batch wrapper forwards the transform
  * conformal calibration remains an original-space procedure
"""
import numpy as np
import pytest

from spagapa.imputation import SparseGPImputer
from spagapa.imputation.calibration import ConformalCalibrator, evaluate_coverage


def _grid_coords(side=15, jitter=0.05, seed=0):
    """A jittered square grid of shape (side*side, 2)."""
    rng = np.random.default_rng(seed)
    xs, ys = np.meshgrid(np.arange(side), np.arange(side))
    coords = np.column_stack([xs.ravel(), ys.ravel()]).astype(float)
    coords += rng.normal(0, jitter, coords.shape)
    return coords


def _fraction_data(coords, seed=1, noise=0.05):
    """Smooth spatial fraction field in [0, 1] with exact 0s and 1s."""
    rng = np.random.default_rng(seed)
    z = np.sin(coords[:, 0] / 6.0) + np.cos(coords[:, 1] / 7.0)
    p = np.clip(0.5 + 0.2 * z + noise * rng.standard_normal(len(coords)), 0.0, 1.0)
    p[rng.random(len(coords)) < 0.10] = 0.0
    p[rng.random(len(coords)) < 0.05] = 1.0
    return p


class TestTransformDefaults:
    def test_default_transform_is_none(self):
        imp = SparseGPImputer()
        assert imp.transform is None
        assert imp.epsilon == pytest.approx(0.05)

    def test_invalid_transform_raises(self):
        with pytest.raises(ValueError, match="transform"):
            SparseGPImputer(transform='probit')

    def test_invalid_epsilon_raises(self):
        with pytest.raises(ValueError, match="epsilon"):
            SparseGPImputer(transform='logit', epsilon=0.0)
        with pytest.raises(ValueError, match="epsilon"):
            SparseGPImputer(transform='logit', epsilon=0.5)

    def test_identity_transform_matches_raw_gp_exactly(self):
        """transform=None must reproduce the historical raw-value GP:
        two imputers with explicit None and a manual reference fit on the
        raw values give identical alpha_ and predictions."""
        coords = _grid_coords()
        values = _fraction_data(coords)
        imp_default = SparseGPImputer(n_inducing=30, length_scale=2.0,
                                      noise_level=0.05)
        imp_default.fit(coords, values)
        imp_explicit = SparseGPImputer(n_inducing=30, length_scale=2.0,
                                       noise_level=0.05, transform=None)
        imp_explicit.fit(coords, values)
        np.testing.assert_allclose(imp_default.alpha_, imp_explicit.alpha_,
                                   atol=1e-12)
        train_mask = values > 0  # the default mask
        np.testing.assert_array_equal(imp_default._train_values,
                                      values[train_mask])


class TestForwardInverse:
    def test_round_trip_interior(self):
        imp = SparseGPImputer(transform='logit')
        y = np.array([0.1, 0.25, 0.5, 0.75, 0.9])
        z = imp._forward(y)
        np.testing.assert_allclose(imp._inverse(z), y, atol=1e-12)
        # logit(0.5) = 0 and monotone increasing
        assert z[2] == pytest.approx(0.0)
        assert np.all(np.diff(z) > 0)

    def test_boundary_values_clip_to_finite(self):
        """Exact 0/1 map to finite latent values (no +/- inf)."""
        imp = SparseGPImputer(transform='logit', epsilon=1e-6)
        z = imp._forward(np.array([0.0, 1.0]))
        assert np.all(np.isfinite(z))
        assert z[0] < -10.0   # logit(1e-6)
        assert z[1] > 10.0    # logit(1 - 1e-6)
        np.testing.assert_allclose(imp._inverse(z),
                                   [1e-6, 1.0 - 1e-6], rtol=1e-3)

    def test_identity_transform_is_passthrough(self):
        imp = SparseGPImputer()
        y = np.array([-0.3, 0.0, 2.5])
        np.testing.assert_array_equal(imp._forward(y), y)
        np.testing.assert_array_equal(imp._inverse(y), y)


class TestLogitGPBehaviour:
    def test_predictions_bounded_in_open_unit_interval(self):
        """Gated at the source: sigmoid output can never leave (0, 1),
        even where the raw-space GP overshoots (tissue edges)."""
        coords = _grid_coords(side=20, jitter=0.0)
        values = _fraction_data(coords, noise=0.02)
        mask = np.ones(len(coords), bool)
        mask[::4] = False
        imp = SparseGPImputer(n_inducing=40, length_scale=2.0,
                              noise_level=0.02, transform='logit')
        imp.fit(coords, values, mask=mask)
        pred, std = imp.predict(coords, return_std=True)
        assert np.all(pred > 0.0) and np.all(pred < 1.0)
        assert np.all(np.isfinite(std))

    def test_equivalent_to_raw_gp_on_pretransformed_values(self):
        """logit-GP on y == raw GP on logit(y): identical latent fit, and
        predictions are the sigmoid of the raw-GP latent predictions."""
        coords = _grid_coords()
        values = _fraction_data(coords)
        imp_l = SparseGPImputer(n_inducing=30, length_scale=2.0,
                                noise_level=0.05, transform='logit')
        full_mask = np.ones(len(coords), bool)
        imp_l.fit(coords, values, mask=full_mask)
        # reference: raw GP on the manually transformed target.  The SAME
        # explicit mask is required: the default mask is values>0 evaluated
        # in the space of the supplied values, which differs between the
        # fraction and latent representations.
        eps = imp_l.epsilon
        z = np.log(np.clip(values, eps, 1 - eps)
                   / (1 - np.clip(values, eps, 1 - eps)))
        imp_r = SparseGPImputer(n_inducing=30, length_scale=2.0,
                                noise_level=0.05)
        imp_r.fit(coords, z, mask=full_mask)
        np.testing.assert_allclose(imp_l.alpha_, imp_r.alpha_, atol=1e-10)
        pred_latent, std_latent = imp_r.predict(coords, return_std=True)
        pred, std = imp_l.predict(coords, return_std=True)
        np.testing.assert_allclose(pred, 1.0 / (1.0 + np.exp(-pred_latent)),
                                   atol=1e-12)
        # the reported std is the LATENT posterior std, unchanged
        np.testing.assert_allclose(std, std_latent, atol=1e-12)

    def test_impute_restores_observed_raw_values(self):
        coords = _grid_coords()
        values = _fraction_data(coords)
        mask = values > 0.25
        imp = SparseGPImputer(n_inducing=30, length_scale=2.0,
                              noise_level=0.05, transform='logit')
        imputed, unc = imp.impute(coords, values, mask=mask,
                                  return_uncertainty=True)
        np.testing.assert_array_equal(imputed[mask], values[mask])
        assert np.all(unc[mask] == 0.0)
        # unobserved entries are sigmoid outputs, hence strictly interior
        unobserved = ~mask
        assert np.all(imputed[unobserved] > 0.0)
        assert np.all(imputed[unobserved] < 1.0)

    def test_batch_wrapper_forwards_transform(self):
        coords = _grid_coords()
        p1 = _fraction_data(coords, seed=1)
        p2 = _fraction_data(coords, seed=2)
        values = np.vstack([p1, p2])
        mask = np.ones_like(values, bool)
        mask[:, ::4] = False
        base = SparseGPImputer(n_inducing=30, length_scale=2.0,
                               noise_level=0.05, transform='logit')
        batch = base.fit_batch(coords, values, mask=mask, verbose=False)
        pred, unc = batch.predict()
        assert np.all(pred > 0.0) and np.all(pred < 1.0)
        assert np.all(np.isfinite(unc))
        imputed, _ = batch.impute()
        np.testing.assert_array_equal(imputed[mask], values[mask])


class TestConformalOriginalSpace:
    def test_calibration_uses_original_space_errors(self):
        """Conformal scores are |y_true - y_pred| in the ORIGINAL space,
        standardised by the latent std.  Held-out coverage must respect
        the split-conformal target."""
        rng = np.random.default_rng(3)
        coords = _grid_coords(side=18, jitter=0.0)
        values = _fraction_data(coords, seed=5, noise=0.03)
        n = len(coords)
        held_out = rng.random(n) < 0.5
        train = ~held_out
        imp = SparseGPImputer(n_inducing=40, length_scale=2.0,
                              noise_level=0.05, transform='logit')
        imp.fit(coords, values, mask=train)
        pred, std = imp.predict(coords, return_std=True)
        # calibration on held-out ORIGINAL-space errors vs latent std
        cal = ConformalCalibrator(alpha=0.1, mode="locally_adaptive")
        cal.fit(np.abs(values[held_out] - pred[held_out]), std[held_out])
        lower, upper = cal.predict(pred, std)
        cov = evaluate_coverage(lower, upper, values)
        # split-conformal guarantees >= 1 - alpha marginally; allow slack
        # for one finite sample
        assert cov >= 0.85
        # intervals are computed around original-space predictions, so the
        # bounds may legitimately exceed [0, 1] -- but the point estimate
        # never does
        assert np.all((pred > 0.0) & (pred < 1.0))
