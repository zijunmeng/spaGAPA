spaGAPA documentation
====================

**spaGAPA** — Spatial Gaussian-Process and Graph-Aware APA Analyzer — is a
Python framework for alternative polyadenylation (APA) analysis in spatial
transcriptomics, with calibrated uncertainty quantification.

.. code-block:: text

    Spatial Validation → Sparse GP Imputation (+ Conformal UQ)
        → APA Quantification → Leiden Domains → Differential APA
        → SVAPA → Bias Correction

What spaGAPA gives you
----------------------

* **Probabilistic imputation** of sparse spatial APA matrices (inducing-point
  sparse GP, :math:`O(n \cdot m^2)`; scales to 100k spots).
* **Conformal-calibrated prediction intervals** with distribution-free
  marginal coverage (validated across 11 samples / 523k held-out points;
  mean deviation ≤ 0.5 pp at 80/90/95%).
* **Caller-agnostic input**: consumes scAPAtrap *or* Sierra peak matrices.
* **Unsupervised domain detection** via a CPU-only multi-view graph
  (spatial + expression + APA) with Leiden clustering.
* **Downstream statistics**: differential APA (donor-level), SVAPA
  (Moran's I / Geary's C), APA batch correction (QN + linear removal).
* **High-resolution support**: end-to-end recipes for subcellular
  Stereo-seq data (SAW → retag → scAPAtrap → bin200).

Where to start
--------------

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   usage/installation
   usage/quickstart

.. toctree::
   :maxdepth: 2
   :caption: Practical guide

   usage/data_preparation
   usage/cli
   usage/presets
   usage/uncertainty
   usage/downstream
   usage/stereo_seq

.. toctree::
   :maxdepth: 2
   :caption: Reference

   api/modules
   usage/troubleshooting

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
