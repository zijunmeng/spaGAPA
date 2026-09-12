Installation
============

Requirements
------------

* Python **3.10+**
* Linux (tested on RHEL 9); 64 GB RAM recommended for Visium-scale data,
  256 GB+ for large Stereo-seq chips
* ``OPENBLAS_NUM_THREADS`` **must** be bounded (8 is the validated value;
  64/32 causes segfaults in heavy linear algebra — see
  :doc:`troubleshooting`)

Install the package
-------------------

.. code-block:: bash

    conda create -n spagapa python=3.10
    conda activate spagapa
    pip install numpy scipy scikit-learn pandas matplotlib pysam statsmodels

    git clone <repo>
    cd spaGAPA/spaGAPA
    pip install -e .

This registers the ``spagapa`` console command:

.. code-block:: bash

    spagapa --version
    # spaGAPA, version 0.1.0

Optional external tools
-----------------------

spaGAPA itself is pure Python; producing an APA matrix from raw data uses
these external tools (only needed if you start from FASTQ/BAM):

====================  =====================================================
Tool                  Role
====================  =====================================================
SAW 8.2.2             Stereo-seq alignment + DNB barcode→coordinate mapping
samtools              BAM sorting/indexing/reheader
scAPAtrap (R)         PAS peak calling (primary caller)
Sierra (R)            Alternative PAS caller (validated swap)
umi_tools,            scAPAtrap dependencies (dedup/count)
featureCounts
Space Ranger          10x Visium alignment + expression (Visium inputs)
====================  =====================================================

Environment flags (important)
-----------------------------

.. code-block:: bash

    export OPENBLAS_NUM_THREADS=8      # CRITICAL: 64/32 segfaults
    export TMPDIR=/path/to/fast/local/disk   # avoid NFS for intermediates

These two flags prevent the two most common runtime failures on shared
servers (see :doc:`troubleshooting`).
