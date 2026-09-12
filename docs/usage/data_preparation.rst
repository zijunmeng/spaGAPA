Data preparation
================

spaGAPA consumes a **PAS × spot matrix + coordinates**. This page shows the
two validated production chains that produce it from raw data. All commands
are the exact invocations used on the project's datasets.

Route A — 10x Visium (Space Ranger → scAPAtrap)
-----------------------------------------------

.. code-block:: bash

    # 1. Align + expression (standard 10x workflow)
    spaceranger count --id=SAMPLE --transcriptome=refdata-gex-GRCh38-2024-A \
        --fastqs=fastqs/ --sample=SAMPLE --create-bam true

    # 2. Call PAS peaks with scAPAtrap (R)
    #    tools = samtools/umi_tools/featureCounts/STAR paths
    #    trap.params$chrs = chr1..chr22,chrX,chrY ; readlength = 91
    scAPAtrap(tools=tools, trap.params=trap.params,
              inputBam="SAMPLE/outs/possorted_genome_bam.bam",
              outputDir="scapatrap_raw/")

    # 3. Convert to spaGAPA format (peak x spot usage matrix + coordinates)
    #    repo scripts: run_scapatrap_spaceranger.py + conversion utilities

Sierra-called BAMs work identically (caller-agnostic by design — validated
by the scAPAtrap↔Sierra swap; coverage guarantees unchanged, max deviation
0.5 pp across 3 datasets).

Route B — Stereo-seq, subcellular (SAW → retag → scAPAtrap → bin200)
---------------------------------------------------------------------

Step 1 — SAW count (needs raw FASTQ + chip mask)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    saw count \
        --id=GSM_SAMPLE_saw --sn=<CHIP_SN> \
        --omics=transcriptomics \
        --kit-version="Stereo-seq T FF V1.3" \
        --sequencing-type="PE75_50+100" \
        --chip-mask=masks/<CHIP_SN>.barcodeToPos.h5 \
        --organism=human --tissue=brain \
        --fastqs=fastq_saw/<CHIP_SN>/ \
        --reference=/path/to/SAW_refs/Homo_sapiens_index \
        --threads-num=32 --memory=200

Reference indexes: official ``reference-data-{human,mouse,rat}.tar.gz``
from the SAW download centre (md5-verified on arrival).

**Chip mask acquisition** — the hard part; two sources:

* **GEO-deposited masks** (best): 2025+ studies increasingly deposit
  ``barcodeToPos.h5`` in sample supplementary files — download directly:

  .. code-block:: bash

      aria2c -c -x4 -s4 -o D02266B4.barcodeToPos.h5 \
        "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSMnnnGSMXXXX/suppl/GSMXXXX_SN.barcodeToPos.h5"

* **STOmics chip-mask portal** (login required): masks exist only during a
  limited cloud retention window; older chips are purged
  (github.com/STOmics/SAW/issues/268). The official ``download_mask`` CLI
  retrieves by chip SN once configured with portal credentials.

Step 2 — retag the SAW BAM for scAPAtrap/umi_tools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The SAW annotation BAM stores barcode as ``Cx``/``Cy``, UMI as ``UR``.
scAPAtrap's TenX path needs 10x-style tags, and umi_tools requires a
**fixed-length** UMI tag:

.. code-block:: python

    # scripts/retag_scapatrap_v2.py — per read:
    #   CB = "{Cx}_{Cy}"                spatial barcode
    #   UB = RX = (UR + 'NNNNN')[:5]    UMI padded/truncated to fixed length
    #   GX = GI ; GN = GS               gene id / symbol

.. code-block:: bash

    python scripts/retag_scapatrap_v2.py \
        SAW/outs/bam/annotated_bam/<SN>....target.bam retag.bam

    samtools reheader chr_header.sam retag.bam > retag.chr.bam  # SN: → chrSN
    samtools sort -@16 -T /local/disk/tmp -o final.bam retag.chr.bam
    samtools index final.bam

Step 3 — scAPAtrap with TenX tag mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: r

    trap.params <- setTrapParams(print=FALSE)
    trap.params$TenX <- TRUE     # REQUIRED: reads UB/CB tags and passes
                                 # --extract-umi-method=tag to umi_tools
    trap.params$chrs <- ...      # your genome's names (see note below)
    trap.params$readlength <- 100

**Chromosome naming**: the STOmics *rat* reference uses RefSeq names
(``NC_051336.1``…), *not* ``chr1..chr20``. Derive the list from the BAM
header instead of hard-coding:

.. code-block:: r

    hdr <- Rsamtools::scanBamHeader(bam)
    trap.params$chrs <- grep("^NC_", names(hdr$targets), value=TRUE)

Step 4 — bin aggregation (DNB → bin200)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    python scripts/stereo_bin200.py \
        --raw-dir scapatrap_raw/ \
        --out-dir binned_200/ \
        --dataset gse293464 --sample GSM8882884 --bin 200

Outputs ``apa_matrix.csv`` (site × bin raw UMI), ``coordinates.csv``,
``apa_sites.csv``, ``qc_summary.json`` — the exact quickstart inputs.

Reference-scale outputs from this chain: mouse AD brain 21,455 PAS ×
15,235 bins (~10.3% observed); human retina organoid 22,762 PAS × 14,905
bins; rat thymus 23,138 PAS × 1.26e8 DNBs before binning.

Downloading SRA runs (fast, resumable)
--------------------------------------

.. code-block:: bash

    # proven pattern; proxy must be DISABLED (socks proxy throttles to KiB/s)
    env -u http_proxy -u https_proxy -u all_proxy \
      aria2c -c -x4 -s4 https://sra-pub-run-odp.s3.amazonaws.com/sra/SRRXXXX/SRRXXXX

    # convert — pass the FILE PATH, never the accession (accession silently
    # goes to slow NCBI remote):
    fasterq-dump -e 16 --split-files --temp $TMPDIR -O out/ SRRXXXX.sra
