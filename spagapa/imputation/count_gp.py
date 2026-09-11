"""Count-native (binomial-likelihood) spatial usage estimation.

Motivation
----------
The distal-usage index is a **binomial proportion**: at spot ``s``, gene
``g`` has ``k`` distal-supporting site counts out of ``n`` total site
counts, and the usage index is ``k/n``.  spaGAPA's original
:class:`~spagapa.imputation.SparseGPImputer` models these fractions with a
Gaussian likelihood, which is misspecified in three ways: the support is
bounded (0, 1), the distribution is bimodal at 0/1 for low-depth entries,
and the observation noise is *depth-dependent* (a usage of 0.5 from n=4 is
not the same evidence as 0.5 from n=200).  This misspecification is half
the reason the per-gene MEAN baseline wins entry-wise RMSE on the bimodal
index (disclosed in the manuscript limitations).

Approach
--------
This module models the generating likelihood via a variance-stabilizing
GLM approximation (the "empirical logit" family):

1. ``theta = empirical_logit(k, n) = logit((k + 0.5) / (n + 1))``
2. delta-method observation variance ``Var(theta) ~= 1/(k+0.5) + 1/(n-k+0.5)``
3. a heteroscedastic sparse GP (the existing
   :class:`SparseGPImputer` with per-entry ``spot_noise``) is fit in logit
   space, where the noise floor from tiny-count entries is *exactly* the
   binomial sampling uncertainty rather than an arbitrary constant;
4. predictions are mapped back with the sigmoid, and posterior stds with
   the delta method ``std_p = sigmoid'(theta) * std_theta``.

This is one Fisher-scoring step of a full binomial-GP GLM (Laplace
approximation is the planned refinement) and is already count-aware in the
aspects that matter for the bimodal-index pathology: boundedness, depth
weighting, and heteroscedastic calibration.
"""
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

from .sparse_gp import SparseGPImputer


# ---------------------------------------------------------------------------
# pure transforms
# ---------------------------------------------------------------------------
def empirical_logit(k, n):
    """logit((k + 0.5) / (n + 1)) — variance-stabilized binomial proportion."""
    k = np.asarray(k, dtype=float)
    n = np.asarray(n, dtype=float)
    p = (k + 0.5) / (n + 1.0)
    return np.log(p / (1.0 - p))


def binomial_logit_var(k, n):
    """Delta-method variance of the empirical logit.

    ``Var(theta_hat) ~= 1/(k+0.5) + 1/(n-k+0.5)``; symmetric in k <-> n-k
    and decreasing in depth n.
    """
    k = np.asarray(k, dtype=float)
    n = np.asarray(n, dtype=float)
    return 1.0 / (k + 0.5) + 1.0 / (n - k + 0.5)


def _sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


# ---------------------------------------------------------------------------
# count-native usage GP
# ---------------------------------------------------------------------------
class BinomialUsageGP:
    """Sparse GP for binomial usage counts in empirical-logit space.

    Parameters
    ----------
    n_inducing, length_scale, noise_floor, total_floor :
        GP settings mirror :class:`SparseGPImputer`; ``noise_floor`` is the
        minimum logit-space variance used for observed entries (keeps the
        solve stable for very deep entries whose binomial variance ~ 0),
        ``total_floor`` is the minimum total count for an entry to be used
        in training (``n >= total_floor``).

    Examples
    --------
    >>> gp = BinomialUsageGP(n_inducing=100, length_scale=100.0)
    >>> gp.fit(coords, k_counts, n_totals)      # (genes x spots), (genes x spots)
    >>> usage = gp.predict_usage()              # (genes x spots) in (0, 1)
    """

    def __init__(
        self,
        n_inducing: int = 100,
        length_scale: float = 100.0,
        noise_floor: float = 0.02,
        total_floor: int = 2,
    ):
        self.n_inducing = int(n_inducing)
        self.length_scale = float(length_scale)
        self.noise_floor = float(noise_floor)
        self.total_floor = int(total_floor)
        self.batch_ = None
        self.mask_ = None
        self.coords_ = None
        self.shape_ = None

    def fit(self, coordinates: np.ndarray, k_counts: np.ndarray,
            n_totals: np.ndarray) -> "BinomialUsageGP":
        """Fit the heteroscedastic sparse GP on empirical-logit targets.

        Parameters
        ----------
        coordinates : (n_spots, 2)
        k_counts : (n_genes, n_spots) distal-supporting counts
        n_totals : (n_genes, n_spots) total site counts per entry
        """
        coords = np.asarray(coordinates, dtype=float)
        k = np.asarray(k_counts, dtype=float)
        n = np.asarray(n_totals, dtype=float)
        if k.shape != n.shape:
            raise ValueError(
                f"k_counts {k.shape} and n_totals {n.shape} must share shape"
            )
        if k.ndim != 2 or coords.ndim != 2 or coords.shape[0] != k.shape[1]:
            raise ValueError(
                f"shape mismatch: coords {coords.shape} vs counts {k.shape}"
            )
        mask = n >= self.total_floor
        theta = empirical_logit(k, n)
        var = np.maximum(binomial_logit_var(k, n), self.noise_floor)
        # unobserved entries: neutral fill (never trained on)
        theta = np.where(mask, theta, 0.0)
        var = np.where(mask, var, 1.0)

        base = SparseGPImputer(n_inducing=self.n_inducing,
                               length_scale=self.length_scale,
                               noise_level=0.05)
        self.batch_ = base.fit_batch(coords, theta, mask=mask,
                                     spot_noise=var, verbose=False)
        self.mask_, self.coords_, self.shape_ = mask, coords, k.shape
        return self

    def predict_usage(self, return_std: bool = False):
        """Back-transformed usage (and delta-method std), shape (genes, spots)."""
        if self.batch_ is None:
            raise RuntimeError("Call fit() before predict_usage()")
        if return_std:
            theta, theta_std = self.batch_.impute(return_uncertainty=True)
            p = _sigmoid(theta)
            # sigmoid'(t) = p(1-p)
            p_std = p * (1.0 - p) * theta_std
            return p, p_std
        theta = self.batch_.impute(return_uncertainty=False)
        if isinstance(theta, tuple):
            theta = theta[0]
        return _sigmoid(theta)


__all__ = ["empirical_logit", "binomial_logit_var", "BinomialUsageGP"]
