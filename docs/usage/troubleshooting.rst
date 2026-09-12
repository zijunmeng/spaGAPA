Troubleshooting
===============

Every item below was hit in production. Check here first.

Segfault / kernel kill during GP or linalg
------------------------------------------

**Cause**: ``OPENBLAS_NUM_THREADS`` at the machine default (64/32) on
many-core servers.

.. code-block:: bash

    export OPENBLAS_NUM_THREADS=8

Set it globally in the environment used for spaGAPA.

Downloads crawl at KiB/s
------------------------

A system socks proxy throttles aria2c. Disable proxy variables for the
download process:

.. code-block:: bash

    env -u http_proxy -u https_proxy -u all_proxy aria2c ...

Long S3 transfers decay to ~1 MB/s
----------------------------------

Per-transfer decay on ODP links. Restart transfers on a timer; ``aria2c
-c`` resumes from the control file (the repo ships a rotating-downloader
pattern in ``scripts/download_stereo_sra_rotating.sh``).

ENA files pass size check but corrupt
--------------------------------------

Byte-exact size with wrong MD5 happens with killed multi-connection
transfers. **Always verify MD5** against the ENA manifest; reconvert from
the local SRA when in doubt.

fasterq-dump "succeeds" instantly with no output
------------------------------------------------

You passed an accession while a local ``.sra`` was expected, or the
symlink target is relative from the wrong directory. Pass the absolute
``.sra`` **file path**.

umi_tools: "not all umis are the same length"
---------------------------------------------

umi_tools defaults to the ``RX`` tag; SAW's native RX is variable-length.
Fix = retag with fixed-length RX (``scripts/retag_scapatrap_v2.py``).

scAPAtrap: "outputDir exists" / "logf already exists"
-----------------------------------------------------

initScAPAtrap refuses pre-created output dirs and stale internal logs.
Remove both before rerunning:

.. code-block:: bash

    rm -rf scapatrap_raw/ retag_dir/scapatrap_internal.log

scAPAtrap dies on chromosome mismatch
-------------------------------------

``'chr' is not correctly specified. Valid options are NC_...`` — your
reference uses RefSeq names. Derive ``trap.params$chrs`` from the BAM
header (see :doc:`stereo_seq`).

samtools sort hangs for hours on big BAMs
-----------------------------------------

Sort temp files default to the output directory; on NFS this thrashes.
Point ``-T`` at local disk. Symptom: 11% CPU, ~1 MB/s writes.

R dies silently in findTails / export
-------------------------------------

Memory ceiling (0.5–0.7 TB per 1e9 reads). Serialize samples; if
``scAPAtrapData.rda`` was saved, re-export the sparse counts with a
streaming writer instead of a giant ``data.frame``.

GP imputation slow on big chips
-------------------------------

You are on the exact GP. Let ``--analysis-preset auto`` resolve to
``highres_fast`` (or set it explicitly), which switches to the sparse GP
with inducing points.

Differential APA p-values look too good
---------------------------------------

Pseudoreplication: spots pooled across donors. Aggregate to sample level
first (see :doc:`downstream`).
