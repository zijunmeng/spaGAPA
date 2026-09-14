Quickstart
==========

The 5-minute path: matrix + coordinates in, domains + uncertainty out.

Input files
-----------

``apa_matrix.csv``
    Peak (or gene-level APA index) × spot matrix. Row index = peak/gene id,
    columns = spot ids. ``NaN`` = unobserved (dropout) entries.

``coordinates.csv``
    Columns ``spot_id,x,y`` matching the matrix columns.

Both files are produced by the :doc:`data_preparation` recipes.

Python API
----------

.. code-block:: python

    from spagapa import SpaGAPA

    pipeline = SpaGAPA(analysis_preset="auto")   # or highres_fast / standard
    result = pipeline.run(
        apa_matrix="data/processed/gse183456_scapatrap/apa_matrix.csv",
        coordinates="data/processed/gse183456_scapatrap/coordinates.csv",
    )

    print(result["domains"]["n_domains"])      # Leiden domains
    print(result["imputed_values"].shape)      # imputed APA matrix
    print(result["uncertainty"].shape)         # posterior std per entry

    pipeline.save_results("out_run/")

CLI
---

The same run as one command:

.. code-block:: bash

    export OPENBLAS_NUM_THREADS=8

    spagapa run \
        --apa-matrix data/processed/gse183456_scapatrap/apa_matrix.csv \
        --coordinates data/processed/gse183456_scapatrap/coordinates.csv \
        --analysis-preset auto \
        --output out_run/

Console summary prints the resolved preset, domain count, SVAPA gene count
and coverage. Results land in ``out_run/`` (imputed matrix, uncertainty,
domain labels, QC report — see ``SpaGAPA.save_results``).

Single-purpose commands
-----------------------

.. code-block:: bash

    # GP imputation only
    spagapa impute apa_matrix.csv --coordinates coordinates.csv \
        --sparse --n-inducing 200 \
        --output-imputed imputed.npy --output-uncertainty unc.npy

    # Visualise one gene's spatial APA map
    spagapa plot apa_matrix.csv --coordinates coordinates.csv \
        --gene APOE --type spatial --output apoe_map.png

    # Differential APA between annotated domains
    spagapa diff apa_matrix.csv labels.csv \
        --coordinates coordinates.csv --method wilcoxon --fdr 0.05

What happens inside ``run``
---------------------------

1. **Load & validate** — spatial support, read/spot support, neighbor
   consistency (``spagapa.calling``).
2. **GP imputation** — dense exact GP for small data; inducing-point sparse
   GP when the preset resolver flags high-resolution data
   (``spagapa.imputation``).
3. **APA quantification** — RUD / PDUI / WUL indices (``spagapa.quantification``).
4. **Domain detection** — uncertainty-weighted; BioML multi-view graph +
   Leiden when enabled by the preset (``spagapa.bioml``).
5. **Differential APA & SVAPA** — optional stages (``spagapa.analysis``).

Choose a preset once your data is larger than Visium scale — see
:doc:`presets`.
