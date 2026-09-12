Stereo-seq high-resolution guide
================================

End-to-end from raw reads to spaGAPA on a subcellular Stereo-seq chip.
(Command details live in :doc:`data_preparation`; this page is the
operational checklist and the lessons that cost us days to learn.)

Pipeline at a glance
--------------------

.. code-block:: text

    FASTQ + chip mask ─► SAW 8.2.2 count ─► retag (CB/UB/RX/GX/GN + chr)
        ─► scAPAtrap (TenX=TRUE) ─► bin200 ─► spaGAPA highres_fast

Resource planning
-----------------

===========================  ===========================================
Stage                        Practical guidance
===========================  ===========================================
SAW count                    ``--threads=32 --memory=200``; ~2 h per
                             ~600 Gb sample
retag (pysam)                ~3.5 h per 1e9 reads; 55–60 GB BAM out
scAPAtrap                    dedup ~5 h (1e9 reads); findTails is the
                             long pole. Peak RSS 400–700 GB per sample —
                             **serialize large samples**, do not run
                             several findTails concurrently
bin200                       ~1 h per 300M-row counts file
spaGAPA                      ~10 min (15k spots, 20k+ PAS)
===========================  ===========================================

Memory rule of thumb: one scAPAtrap findTails/generatescExpMa pass per
1e9 reads wants 0.5–0.7 TB RAM. Two concurrent samples = OOM territory on
a 1 TB node. If you must overlap, suspend one process (``kill -STOP`` /
``kill -CONT``) instead of losing hours of work.

Reading SAW output layout
-------------------------

SAW minor builds differ; the annotated BAM lives in one of:

.. code-block:: text

    <id>/outs/bam/annotated_bam/<SN>....target.bam          # newer
    <id>/STEREO_ANALYSIS_WORKFLOW_PROCESSING/ANNOTATION/... # older

Probe both (the repo's pipeline scripts do exactly this).

Binning sensitivity
-------------------

Cross-bin Pearson r between bin sizes is part of routine QC (50↔100,
50↔200 over representative peaks); Moran's I decays gently with bin size.
bin200 is the default production resolution (≈ Visium spot scale).

What "good" looks like
----------------------

* PAS yield: 21k–23k curated peaks per large chip section (mouse brain,
  human organoid, rat thymus all land in this range).
* 3′-end enrichment: >50% of gene-annotated reads within 500 bp of TES
  (polyA-capture signature; ~73% within 2 kb).
* Post-bin observed fraction: ~4–10% — this sparsity is exactly what the
  sparse GP + conformal layer is built for.

Known pitfalls (all field-tested)
---------------------------------

* **Chip masks expire.** STOmics OSS retention is finite; prefer GEO-
  deposited ``barcodeToPos.h5``. Probe availability before committing to
  a dataset.
* **fasterq-dump needs the file path**, not the accession (remote fetch
  is ~100× slower and fails silently).
* **umi_tools reads tag RX**, not UB; SAW-native RX lengths can be
  heterogeneous → assert. The v2 retag writes fixed-length ``RX=UB``.
* **scAPAtrap must run with ``TenX=TRUE``** on retagged BAMs (tag-mode
  umi_tools; otherwise it parses read names and crashes).
* **samtools sort on NFS self-throttles** for >40 GB BAMs; always pass
  ``-T /local/disk/tmp``.
* **Fasterq/sort intermediates** double disk usage; budget ≥ 3× the SRA
  size per sample.
