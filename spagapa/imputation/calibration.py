"""
Post-hoc uncertainty calibration for spatial-GP imputation.

This module implements **split conformal prediction** to turn the sparse GP's
poorly-calibrated posterior standard deviation into intervals with rigorous,
reportable *marginal* coverage.  Split conformal is model-agnostic: it only
needs (i) point predictions and (ii) a held-out calibration set of true
values.  No assumptions on the GP or its hyperparameters are required.

Two interval modes are supported:

* ``mode="global"`` -- a single nonconformity quantile ``q_hat`` is applied to
  every point: ``[pred - q_hat, pred + q_hat]``.  This guarantees finite-sample
  marginal coverage at ``1 - alpha`` under exchangeability of calibration and
  test points (the standard split-conformal guarantee).

* ``mode="locally_adaptive"`` -- the per-point GP std ``s_i`` rescales the
  quantile: ``[pred - q_hat * s_i, pred + q_hat * s_i]``.  The nonconformity
  score is the standardized residual ``|y_true - y_pred| / max(s_i, s_floor)``,
  so this mode *also* has marginal coverage guarantees while producing
  tighter intervals where the GP is confident and wider intervals where it is
  not (it is "locally adaptive" in the conformal sense).  This is the more
  useful mode for spaGAPA because it preserves the GP's spatial uncertainty
  structure while correcting the overall scale.

The locally-adaptive mode benefits from a raw GP std that *tracks* error
ranking -- if the GP std is essentially flat (the current spaGAPA failure
mode, unc-error corr ~0.008), locally-adaptive degrades toward global but is
never worse.  ``ConformalCalibrator.fit`` therefore also reports the
calibration-set unc-error correlation so callers can see whether the
locally-adaptive variant is buying anything.

Reference
---------
Lei et al. (2018), "Distribution-Free Predictive Inference for Regression",
arXiv:1604.04173.  Romano et al. (2019), "Conformalized Quantile Regression",
arXiv:1905.03222 (locally-adaptive / CQR-style scaling).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Optional, Tuple

import numpy as np


@dataclass
class CalibratedInterval:
    """Result of fitting a :class:`ConformalCalibrator`.

    Attributes
    ----------
    alpha : float
        Target miscoverage (interval covers ``1 - alpha`` of points).
    mode : str
        ``"global"`` or ``"locally_adaptive"``.
    q_hat : float
        The calibrated nonconformity quantile.  For ``mode="global"`` this is
        an absolute error scale; for ``mode="locally_adaptive"`` it multiplies
        the per-point GP std.
    n_calibration : int
        Number of calibration points used.
    std_floor : float
        The GP-std floor applied to avoid division blow-ups in the
        locally-adaptive score.
    calibration_unc_error_corr : float
        Pearson correlation between GP std and |error| on the calibration
        set (NaN if mode is global or std was not provided).  Diagnostic only.
    """

    alpha: float
    mode: str
    q_hat: float
    n_calibration: int
    std_floor: float
    calibration_unc_error_corr: float = float("nan")

    def scale_for_std(self, gp_std: np.ndarray) -> np.ndarray:
        """Per-point half-width ``q_hat * max(gp_std, std_floor)`` (locally
        adaptive) or broadcast ``q_hat`` (global)."""
        if self.mode == "locally_adaptive":
            s = np.maximum(np.asarray(gp_std, dtype=float), self.std_floor)
            return self.q_hat * s
        # global: constant half-width everywhere
        return np.full(np.shape(gp_std), self.q_hat, dtype=float)


class ConformalCalibrator:
    """Split-conformal calibrator for GP posterior intervals.

    Parameters
    ----------
    alpha : float, default=0.1
        Target miscoverage.  ``alpha=0.1`` produces a 90% interval.
    mode : {"global", "locally_adaptive"}, default="locally_adaptive"
        Interval mode.  See module docstring.
    std_floor : float, optional
        Minimum GP std used in the locally-adaptive score to avoid division
        by ~0.  If ``None``, defaults to ``1e-3 * max(calibration GP std)``
        (or a tiny absolute constant if all stds are zero).

    Examples
    --------
    >>> cal = ConformalCalibrator(alpha=0.1, mode="locally_adaptive")
    >>> interval = cal.fit(cal_err, cal_std).predict(test_pred, test_std)
    >>> lower, upper = interval  # doctest: +SKIP
    """

    def __init__(
        self,
        alpha: float = 0.1,
        mode: Literal["global", "locally_adaptive"] = "locally_adaptive",
        std_floor: Optional[float] = None,
    ):
        if not (0.0 < alpha < 1.0):
            raise ValueError(f"alpha must be in (0, 1), got {alpha}")
        if mode not in ("global", "locally_adaptive"):
            raise ValueError(
                f"mode must be 'global' or 'locally_adaptive', got {mode!r}"
            )
        self.alpha = float(alpha)
        self.mode = mode
        self.std_floor = std_floor
        self.interval_: Optional[CalibratedInterval] = None

    # ------------------------------------------------------------------ fit
    def fit(
        self,
        calibration_errors: np.ndarray,
        calibration_gp_std: Optional[np.ndarray] = None,
    ) -> CalibratedInterval:
        """Compute the nonconformity quantile ``q_hat``.

        Parameters
        ----------
        calibration_errors : np.ndarray, shape (n,)
            Absolute errors ``|y_true - y_pred|`` on the held-out calibration
            set.
        calibration_gp_std : np.ndarray, shape (n,), optional
            GP posterior std on the same calibration points.  Required for
            ``mode="locally_adaptive"``; ignored for ``mode="global"``.

        Returns
        -------
        CalibratedInterval
            The fitted interval descriptor (also stored as
            ``self.interval_``).
        """
        errors = np.asarray(calibration_errors, dtype=float).ravel()
        n = errors.shape[0]
        if n == 0:
            raise ValueError("calibration_errors is empty")

        unc_corr = float("nan")
        if self.mode == "locally_adaptive":
            if calibration_gp_std is None:
                raise ValueError(
                    "mode='locally_adaptive' requires calibration_gp_std"
                )
            std = np.asarray(calibration_gp_std, dtype=float).ravel()
            if std.shape != errors.shape:
                raise ValueError(
                    "calibration_gp_std must have the same shape as "
                    "calibration_errors"
                )
            floor = self.std_floor
            if floor is None:
                smax = float(np.max(std)) if std.size else 0.0
                floor = max(1e-3 * smax, 1e-12) if smax > 0 else 1e-12
            floor = float(floor)
            s = np.maximum(std, floor)
            scores = errors / s  # standardized nonconformity scores
            if np.std(std) > 0 and np.std(errors) > 0:
                unc_corr = float(np.corrcoef(std, errors)[0, 1])
        else:  # global
            floor = float(self.std_floor) if self.std_floor is not None else 0.0
            scores = errors

        # Split-conformal quantile with the standard finite-sample correction:
        #   q_hat = ceil((n+1)(1-alpha))/n quantile, i.e. the
        #   "higher" interpolation so that coverage >= 1-alpha holds
        #   exactly under exchangeability (Lei et al. 2018, eq. 2.4).
        q_level = min(1.0, np.ceil((n + 1) * (1.0 - self.alpha)) / n)
        q_hat = float(np.quantile(scores, q_level, method="higher"))

        self.interval_ = CalibratedInterval(
            alpha=self.alpha,
            mode=self.mode,
            q_hat=q_hat,
            n_calibration=n,
            std_floor=floor,
            calibration_unc_error_corr=unc_corr,
        )
        return self.interval_

    # --------------------------------------------------------------- predict
    def predict(
        self,
        gp_pred: np.ndarray,
        gp_std: Optional[np.ndarray] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Return calibrated ``(lower, upper)`` interval bounds.

        Parameters
        ----------
        gp_pred : np.ndarray
            Point predictions.
        gp_std : np.ndarray, optional
            GP posterior std at the same points.  Required for
            ``mode="locally_adaptive"``; ignored for ``mode="global"``.
        """
        if self.interval_ is None:
            raise RuntimeError("Call fit() before predict()")
        pred = np.asarray(gp_pred, dtype=float)
        if self.mode == "locally_adaptive":
            if gp_std is None:
                raise ValueError(
                    "mode='locally_adaptive' requires gp_std in predict()"
                )
            std = np.asarray(gp_std, dtype=float)
            if std.shape != pred.shape:
                # broadcast-friendly: allow scalar -> array
                std = np.broadcast_to(std, pred.shape).astype(float)
            half = self.interval_.scale_for_std(std)
        else:
            half = self.interval_.scale_for_std(pred)
        lower = pred - half
        upper = pred + half
        return lower, upper


# ---------------------------------------------------------------------------
# Convenience: evaluate a fitted interval's empirical coverage on a test set.
# ---------------------------------------------------------------------------
def evaluate_coverage(
    lower: np.ndarray,
    upper: np.ndarray,
    y_true: np.ndarray,
) -> float:
    """Empirical marginal coverage ``mean(lower <= y_true <= upper)``."""
    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    y_true = np.asarray(y_true, dtype=float)
    if not (lower.shape == upper.shape == y_true.shape):
        raise ValueError("lower, upper, y_true must share a shape")
    inside = (lower <= y_true) & (y_true <= upper)
    if inside.size == 0:
        return float("nan")
    return float(np.mean(inside))


# ---------------------------------------------------------------------------
# Mondrian (group-conditional) split conformal.
# ---------------------------------------------------------------------------
@dataclass
class MondrianInterval:
    """Result of fitting a :class:`MondrianConformalCalibrator`.

    Per-group nonconformity quantiles.  Each listed group has its own
    finite-sample-corrected split-conformal quantile ``q_g`` computed on
    that group's calibration points, which yields a *group-conditional
    marginal* coverage guarantee: for exchangeable test points from group
    ``g``, ``P(y in interval) >= 1 - alpha``.  Groups smaller than the
    calibrator's ``min_group_size`` (and groups unseen at fit time) fall
    back to the pooled quantile ``group_q["pooled"]``, preserving the plain
    marginal guarantee.
    """

    alpha: float
    mode: str
    group_q: dict
    group_n: dict
    std_floor: float
    n_calibration: int


class MondrianConformalCalibrator:
    """Group-conditional ("Mondrian") split-conformal calibrator.

    Runs one split-conformal calibration per predefined group instead of a
    single pooled quantile.  This addresses spaGAPA's disclosed limitation
    that *marginal* coverage conceals subgroup undercoverage (high-
    expression bins reach only ~0.82 at the 90% target): Mondrian
    calibration restores the guarantee *within* each group while remaining
    distribution-free and model-agnostic (Vovk et al. 2005, "Mondrian"
    architectures; Lei et al. 2018 for the split-conformal quantile).

    Parameters
    ----------
    alpha : float, default=0.1
        Target miscoverage per group.
    mode : {"global", "locally_adaptive"}, default="locally_adaptive"
        Within each group, either absolute-error scores (global) or
        std-standardized scores (locally_adaptive), mirroring
        :class:`ConformalCalibrator`.
    min_group_size : int, default=50
        Groups with fewer calibration points use the pooled quantile
        (finite-sample correction needs enough points to be meaningful).
    std_floor : float, optional
        Same role as in :class:`ConformalCalibrator`.

    Examples
    --------
    >>> cal = MondrianConformalCalibrator(alpha=0.1, mode="global")
    >>> cal.fit(cal_err, None, groups=expr_bins)
    >>> lower, upper = cal.predict(test_pred, None, test_bins)
    """

    def __init__(
        self,
        alpha: float = 0.1,
        mode: Literal["global", "locally_adaptive"] = "locally_adaptive",
        min_group_size: int = 50,
        std_floor: Optional[float] = None,
    ):
        if not (0.0 < alpha < 1.0):
            raise ValueError(f"alpha must be in (0, 1), got {alpha}")
        if mode not in ("global", "locally_adaptive"):
            raise ValueError(
                f"mode must be 'global' or 'locally_adaptive', got {mode!r}"
            )
        self.alpha = float(alpha)
        self.mode = mode
        self.min_group_size = int(min_group_size)
        self.std_floor = std_floor
        self.interval_: Optional[MondrianInterval] = None

    @staticmethod
    def _scores(errors: np.ndarray, std: Optional[np.ndarray], mode: str,
                std_floor: Optional[float]) -> Tuple[np.ndarray, float]:
        """Nonconformity scores + the std floor actually applied."""
        if mode == "global":
            return errors, (float(std_floor) if std_floor else 0.0)
        s = np.asarray(std, dtype=float)
        floor = std_floor
        if floor is None:
            smax = float(np.max(s)) if s.size else 0.0
            floor = max(1e-3 * smax, 1e-12) if smax > 0 else 1e-12
        return errors / np.maximum(s, float(floor)), float(floor)

    @staticmethod
    def _q_level(n: int, alpha: float) -> float:
        """Finite-sample-corrected quantile level ceil((n+1)(1-a))/n."""
        return min(1.0, np.ceil((n + 1) * (1.0 - alpha)) / n)

    def fit(
        self,
        calibration_errors: np.ndarray,
        calibration_gp_std: Optional[np.ndarray],
        groups: np.ndarray,
    ) -> MondrianInterval:
        """Per-group split-conformal quantiles from a held-out set.

        ``groups`` are discrete labels (e.g. expression-quantile bins from
        :func:`bin_by_quantiles`).  Small groups pool into ``"pooled"``.
        """
        errors = np.asarray(calibration_errors, dtype=float).ravel()
        groups = np.asarray(groups).ravel()
        if groups.shape != errors.shape:
            raise ValueError(
                f"groups must have the same length as calibration_errors "
                f"({groups.shape[0]} vs {errors.shape[0]})"
            )
        if self.mode == "locally_adaptive" and calibration_gp_std is None:
            raise ValueError("mode='locally_adaptive' requires calibration_gp_std")

        all_scores, floor = self._scores(
            errors, calibration_gp_std, self.mode, self.std_floor
        )

        group_q: dict = {}
        group_n: dict = {}
        pooled_q = float(np.quantile(
            all_scores, self._q_level(errors.shape[0], self.alpha),
            method="higher",
        ))
        group_q["pooled"] = pooled_q

        for g in np.unique(groups):
            m = groups == g
            n = int(m.sum())
            group_n[str(g)] = n
            if n < self.min_group_size:
                group_q[str(g)] = pooled_q
                continue
            group_q[str(g)] = float(np.quantile(
                all_scores[m], self._q_level(n, self.alpha), method="higher"
            ))

        self.interval_ = MondrianInterval(
            alpha=self.alpha,
            mode=self.mode,
            group_q=group_q,
            group_n=group_n,
            std_floor=floor,
            n_calibration=int(errors.shape[0]),
        )
        return self.interval_

    def predict(
        self,
        gp_pred: np.ndarray,
        gp_std: Optional[np.ndarray],
        groups: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Per-point ``(lower, upper)`` using each point's group quantile.

        Unseen groups use the pooled fallback quantile.
        """
        if self.interval_ is None:
            raise RuntimeError("Call fit() before predict()")
        pred = np.asarray(gp_pred, dtype=float).ravel()
        groups = np.asarray(groups).ravel()
        if groups.shape != pred.shape:
            raise ValueError("groups must have the same length as gp_pred")

        q = np.array([
            self.interval_.group_q.get(
                str(g), self.interval_.group_q["pooled"]
            )
            for g in groups
        ])
        if self.mode == "locally_adaptive":
            if gp_std is None:
                raise ValueError("mode='locally_adaptive' requires gp_std in predict()")
            std = np.broadcast_to(
                np.asarray(gp_std, dtype=float), pred.shape
            ).astype(float)
            half = q * np.maximum(std, self.interval_.std_floor)
        else:
            half = q
        return pred - half, pred + half


def bin_by_quantiles(values: np.ndarray, n_bins: int = 5) -> np.ndarray:
    """Map a continuous conditioning variable to ``0..n_bins-1`` labels.

    Equal-frequency bins via empirical quantiles (ties collapse into the
    lower bin, NaN stays NaN).  Typical use: expression-level bins for
    Mondrian calibration of spaGAPA's disclosed high-expression-bin
    undercoverage.
    """
    x = np.asarray(values, dtype=float)
    edges = np.quantile(x[~np.isnan(x)], np.linspace(0, 1, n_bins + 1))
    edges[0], edges[-1] = -np.inf, np.inf
    out = np.digitize(x, edges[1:-1], right=True).astype(float)
    out[np.isnan(x)] = np.nan
    return out


__all__ = [
    "ConformalCalibrator",
    "CalibratedInterval",
    "MondrianConformalCalibrator",
    "MondrianInterval",
    "evaluate_coverage",
    "bin_by_quantiles",
]
