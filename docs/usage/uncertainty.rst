Uncertainty quantification
==========================

The core innovation: split-conformal prediction intervals with
distribution-free marginal coverage on top of GP posterior stds.

Concept
-------

1. Mask 20% of observed entries per gene → GP imputes them →
   ``|error| = |y_true − y_pred|`` on the mask.
2. Split masked errors 50/50 into calibration / test.
3. ConformalCalibrator computes the finite-sample-corrected quantile
   :math:`\hat{q} = \lceil (n{+}1)(1{-}\alpha) \rceil / n` of calibration
   errors; intervals = prediction ± (scaled :math:`\hat{q}`).
4. Coverage is verified on the held-out test split.

Validated performance: 11 samples, 523,174 test points — mean absolute deviation
0.21 / 0.16 / 0.10 pp at 80/90/95%, max ≤ 0.5 pp, invariant to PAS caller,
species and tissue.

Batch script
------------

.. code-block:: bash

    python scripts/calibrate_uncertainty.py \
        --apa-matrix  binned_200/apa_matrix.csv \
        --coordinates binned_200/coordinates.csv \
        --output      uncertainty_calibration/ \
        --mask-fraction 0.2 --cal-fraction 0.5 --seed 42

Python API
----------

.. code-block:: python

    from spagapa.imputation import SparseGPImputer
    from spagapa.imputation.calibration import (
        ConformalCalibrator, evaluate_coverage,
    )

    imputer = SparseGPImputer(n_inducing=200, length_scale=100.0, noise_level=0.1)
    batch = imputer.fit_batch(coords, apa_matrix, mask=observed_mask)
    predictions, uncertainty = batch.impute(return_uncertainty=True)

    calibrator = ConformalCalibrator(alpha=0.1)          # 90% intervals
    calibrator.fit(calibration_errors)                   # |y - yhat| values
    interval = calibrator.predict(predictions, uncertainty)

    print(interval.coverage(test_truth))                 # ≈ 0.90

Two interval modes: **global** (one width) and **locally adaptive**
(scaled by GP std). Locally adaptive widths track per-spot difficulty;
spatial-block splits stay at nominal.

Using uncertainty downstream
----------------------------

* **Triage**: dropping the 20% highest-uncertainty predictions removes
  ~23% of RMSE at 80% retention (paired across 5 datasets, p = 0.004).
* **Weighting**: pass ``use_uncertainty=True`` to ``SpaGAPA.run`` /
  ``spagapa run --use-uncertainty`` so differential APA and SVAPA weigh
  confident entries more.
* **Honest limits**: the guarantee is *marginal*; high-expression bins
  undercover (~0.82). Within-gene uncertainty–error correlation is modest
  (median r ≈ 0.10–0.19) — report intervals as entry-wise scores, not
  per-gene truth probabilities.
