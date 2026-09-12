CLI reference
=============

Entry point: ``spagapa`` (click group). Version 0.1.0.

``spagapa run``
---------------

Complete pipeline: validate → impute → quantify → domains → (diff) → (SVAPA).

Key options (full list via ``spagapa run --help``):

=============================================  ==========================================
Option                                         Meaning
=============================================  ==========================================
``--apa-matrix, -m PATH``                      APA matrix CSV/TSV/NPY
``--coordinates, -c PATH``                     ``spot_id,x,y`` CSV
``--dataset, -d PATH``                         Saved ``APADataset`` (``.h5ad``) instead
``--matrix-orientation``                       ``genes_by_spots`` (default) / reversed
``--analysis-preset``                          ``auto`` | ``standard`` | ``highres_accuracy`` | ``highres_fast``
``--sparse``                                    Use inducing-point sparse GP
``--sparse-gp-inducing-method``                ``kmeans`` (default) | ``random`` | ``grid``
``--sparse-gp-length-scale``                   RBF length scale, or ``auto``
``--sparse-gp-noise-level``                    Observation noise (default 0.1)
``--expression-matrix``                        Optional expression matrix for BioML
``--enable-bioml / --disable-bioml``           Override preset-driven BioML domains
``--bioml-domain-method``                      ``spectral`` (default) | ``kmeans``
``--bioml-{spatial,expression,apa}-weight``    Multi-view graph weights (0.4/0.4/0.2)
``--diff``                                     Enable differential analysis
``--no-svapa``                                 Skip SVAPA detection
``--fdr FLOAT``                                FDR threshold (default 0.05)
``--use-uncertainty``                          Weight analyses by GP uncertainty
``--output, -o DIR``                           Result directory
=============================================  ==========================================

Example — Stereo-seq binned sample, preset chosen by shape:

.. code-block:: bash

    spagapa run -m binned_200/apa_matrix.csv -c binned_200/coordinates.csv \
        --analysis-preset highres_fast \
        --expression-matrix expr_matrix.csv \
        -o out_stereo/

``spagapa impute``
------------------

GP imputation only.

.. code-block:: bash

    spagapa impute apa_matrix.csv -c coordinates.csv \
        --sparse --output-imputed imputed.npy --output-uncertainty unc.npy

``spagapa diff``
----------------

Differential APA between labeled groups (domains or conditions).

.. code-block:: bash

    spagapa diff apa_matrix.csv labels.csv -c coordinates.csv \
        --method wilcoxon --fdr 0.05 --logfc 0.5 -o diff_out/

``--method``: ``wilcoxon`` | ``t-test`` | ``permutation``.
Outputs per-domain marker tables.

.. warning::
   For *biological* condition contrasts the replicate unit must be the
   donor/animal, not the spot. Pooling spots across donors inflates n
   (pseudoreplication). Aggregate to per-gene per-sample values first when
   ``n_samples >= 3``; for ``n=1`` designs report effect sizes only.

``spagapa plot``
----------------

.. code-block:: bash

    spagapa plot apa_matrix.csv -c coordinates.csv -g APOE -t spatial
    spagapa plot apa_matrix.csv -t volcano -U unc.npy -D domains.csv

``--type``: ``spatial`` | ``volcano`` | ``heatmap`` | ``qc``.
``-U`` adds uncertainty; ``-D`` colors by domain labels.
