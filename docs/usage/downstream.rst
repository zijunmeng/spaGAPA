Downstream analyses
===================

Differential APA
----------------

.. code-block:: python

    from spagapa.analysis import DifferentialAPAAnalyzer

    analyzer = DifferentialAPAAnalyzer(method="t-test")
    results = analyzer.test_differential_apa(
        apa_matrix=sample_level_matrix,        # gene x n_samples
        group1_indices=control_sample_indices,
        group2_indices=ad_sample_indices,
        gene_names=gene_names,
    )
    results = analyzer.adjust_pvalues(results, method="fdr_bh")

.. important::
   The replicate unit must be the biological sample (donor/animal).
   Spot-level pooling = pseudoreplication. For n=1-per-condition designs
   (e.g. Stereo-seq pilot sections) report effect sizes only, no p-values.

Spatially variable APA (SVAPA)
------------------------------

.. code-block:: python

    from spagapa.analysis.svapa import svapa, morans_i

    score, pval, n = morans_i(apa_values, coords, k=8)
    res = svapa(apa_matrix, coords, k=8, n_perm=200, fdr="fdr_bh")

Combines GP likelihood-ratio trend detection with permutation-based
Moran's I; uncertainty weighting available via ``SpaGAPA.run``.

APA batch correction (optional)
-------------------------------

.. code-block:: python

    from spagapa.analysis.bias_correction import (
        quantile_normalize, linear_batch_correction,
    )

    corrected = quantile_normalize(apa_matrix, group_labels=batch_labels)
    corrected = linear_batch_correction(
        apa_matrix, batch_labels=batch_labels,
        preserve_labels=condition_labels,      # limma-style, keeps biology
    )

Both functions auto-detect collinear designs (batch ⊄ biology) and warn.
This module is under active evaluation (Harmony comparison pending) — use
for exploration, label as such in reports.

Domain recovery
---------------

Unsupervised, no labels required:

.. code-block:: python

    SpaGAPA(analysis_preset="highres_fast").fit_transform(dataset)
    # result.domains — Leiden labels over the fused
    # spatial+expression+APA graph

Reference anchor: MOB ST11 (260 spots, 5 annotated layers) — unsupervised
ARI = 0.60, NMI = 0.68 with the identical Leiden pipeline for both
mean- and GP-imputed inputs (fair-comparison design).

Visualization
-------------

.. code-block:: bash

    spagapa plot apa_matrix.csv -c coordinates.csv -g APOE -t spatial \
        -U uncertainty.npy -D domains.csv -o apoecombined.png
