#!/usr/bin/env python
"""Cache read 3'-end-to-TES signed distances for Supplementary Figure S12
Panel B. Streams the Stereo-seq BAM, samples up to TARGET gene-annotated
reads, computes strand-aware distance, caches a histogram + raw sample to CSV.

Mirrors the approach in gse263789_stereo_pilot/figures/fig1_3prime_enrichment.py
but only persists the data needed for the supp-figure histogram.
"""
import os, sys
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import numpy as np

GTF = "/s1/SHARE/01_software/SAW_refs/Mus_musculus_index/genes/Mus_musculus.GRCm38.93.saw.gtf"
BAM = ("/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output/"
       "gse263789_stereo_pilot/gsm8199179_pilot/STEREO_ANALYSIS_WORKFLOW_PROCESSING/"
       "ANNOTATION/SS200000769BR.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam")
# Fallback to the retag forward BAM if the SAW BAM is missing.
BAM_FALLBACK = ("/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output/"
                "gse263789_stereo_pilot/scapatrap_input.retag.chr.UniqSorted.dedup.forward.bam")
OUT_DIR = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output/supplementary_figures/_cache"
os.makedirs(OUT_DIR, exist_ok=True)
OUT_NPY = os.path.join(OUT_DIR, "s12_tes_distances.npy")
TARGET_READS = 500_000  # enough for a smooth histogram
HIST_BINS = np.arange(-3000, 3001, 50)  # signed dist in bp


def pick_bam():
    for p in [BAM, BAM_FALLBACK]:
        if os.path.exists(p):
            return p
    raise FileNotFoundError("no BAM found")


def load_genes():
    genes = {}
    with open(GTF) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            chrom, start, end, strand = f[0], int(f[3]), int(f[4]), f[6]
            attr = f[8]
            i = attr.find("gene_id")
            if i < 0:
                continue
            j = attr.find('"', i)
            k = attr.find('"', j + 1)
            gid = attr[j + 1:k]
            genes[gid] = (chrom, start, end, strand)
    return genes


def main():
    import pysam
    if os.path.exists(OUT_NPY):
        print(f"[skip] cache exists -> {OUT_NPY}", flush=True)
        return
    bam_path = pick_bam()
    print(f"[bam] using {bam_path}", flush=True)
    genes = load_genes()
    print(f"[gtf] {len(genes):,} genes loaded", flush=True)

    bam = pysam.AlignmentFile(bam_path, "rb")
    filled = 0
    n_total = 0
    n_annotated = 0
    reservoir = []  # reservoir sample of signed distances
    rng = np.random.RandomState(42)
    for r in bam:
        n_total += 1
        if filled >= TARGET_READS and n_total % 500_000 == 0:
            # Checkpoint: stop early once we have enough.
            break
        gi = r.get_tag("GI") if r.has_tag("GI") else None
        if gi is None or gi not in genes:
            continue
        n_annotated += 1
        chrom, gstart, gend, gstrand = genes[gi]
        if gstrand == "+":
            tes = gend
            r3 = r.reference_end  # forward read 3'
        else:
            tes = gstart
            r3 = r.reference_start  # reverse read 3'
        if r3 is None:
            continue
        dist = int(r3) - int(tes)
        # Reservoir-sample to cap memory while keeping uniform coverage.
        if filled < TARGET_READS:
            reservoir.append(dist)
            filled += 1
        else:
            j = rng.randint(0, filled + 1)
            if j < TARGET_READS:
                reservoir[j] = dist
        if n_total % 500_000 == 0:
            print(f"  ...scanned {n_total:,} reads, annotated kept {filled:,}", flush=True)

    bam.close()
    arr = np.asarray(reservoir, dtype=np.int64)
    np.save(OUT_NPY, arr)
    print(f"[done] cached {arr.size:,} distances -> {OUT_NPY}", flush=True)
    print(f"       median={np.median(arr):.0f}bp mean={arr.mean():.0f}bp "
          f"frac within ±200bp of TES={np.mean(np.abs(arr) <= 200):.3f}", flush=True)


if __name__ == "__main__":
    main()
