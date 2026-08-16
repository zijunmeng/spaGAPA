#!/usr/bin/env python
"""Extract splice junctions from a Space Ranger possorted_genome_bam as a
regtools-style bed12 file for Sierra::FindPeaks.

Usage: python extract_junctions.py --dataset gse220442_gsm6801751 [--nproc 16]

Counts every N CIGAR gap of primary alignments (skips secondary/supplementary
and unmapped). Chromosomes are split into parts for parallelism; each junction
(intron start, end) is assigned to exactly the part containing its start, so
boundary reads are never double-counted. Only junctions with >= MIN_KEEP
supporting reads are kept (Sierra's default junction mask cutoff is
min.jcutoff=50, so lower-count junctions can never be used anyway).
"""
import os, pickle, argparse
from collections import Counter
import multiprocessing as mp
import pysam

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from datasets import get as get_dataset

MIN_KEEP = 25
ANCHOR = 8
PART = 400  # total work units

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--dataset", required=True,
                    help="dataset key from datasets.py (bam/outdir resolved from it)")
parser.add_argument("--bam", default=None, help="override BAM path")
parser.add_argument("--out-dir", default=None, help="override output dir")
parser.add_argument("--nproc", type=int, default=16)
args = parser.parse_args()
ds = get_dataset(args.dataset)
BAM = args.bam or ds["bam"]
OUT_DIR = args.out_dir or ds["outdir"]
OUT_BED = os.path.join(OUT_DIR, "junctions.bed")
TMP_DIR = os.path.join(OUT_DIR, "junc_tmp")
NPROC = args.nproc


def worker(job):
    chrom, part, nparts, start, end = job
    counts = Counter()
    n_reads = 0
    bam = pysam.AlignmentFile(BAM, "rb", threads=1)
    for read in bam.fetch(chrom, start, end):
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        cigar = read.cigartuples
        if not cigar:
            continue
        has_n = False
        for op, _ in cigar:
            if op == 3:
                has_n = True
                break
        if not has_n:
            continue
        n_reads += 1
        pos = read.reference_start
        for op, length in cigar:
            if op in (0, 2, 7, 8):      # M, D, =, X
                pos += length
            elif op == 3:              # N -> intron [pos, pos+length)
                if start <= pos < end:
                    counts[(pos, pos + length)] += 1
                pos += length
    bam.close()
    out = os.path.join(TMP_DIR, f"{chrom}.part{part}.pkl")
    with open(out, "wb") as f:
        pickle.dump((chrom, dict(counts), n_reads), f)
    return f"{chrom}:{start + 1}-{end} reads={n_reads} juncs={len(counts)}"


def main():
    os.makedirs(TMP_DIR, exist_ok=True)
    done = set(os.listdir(TMP_DIR))
    bam = pysam.AlignmentFile(BAM, "rb")
    chrs, total_len = [], 0
    for c in bam.references:
        if not c.startswith("chr") or c == "chrM" or c == "chrEBV" or "_" in c or "." in c:
            continue
        chrs.append(c)
        total_len += bam.get_reference_length(c)
    jobs = []
    for c in chrs:
        ln = bam.get_reference_length(c)
        nparts = max(1, round(PART * ln / total_len))
        for p in range(nparts):
            step = ln // nparts + 1
            start, end = p * step, min(ln, (p + 1) * step)
            if f"{c}.part{p}.pkl" in done:
                continue
            jobs.append((c, p, nparts, start, end))
    bam.close()
    print(f"[junc] {len(chrs)} chromosomes, {len(jobs)} jobs to run", flush=True)
    if jobs:
        with mp.Pool(NPROC) as pool:
            for i, msg in enumerate(pool.imap_unordered(worker, jobs, chunksize=1)):
                print(f"[junc] ({i + 1}/{len(jobs)}) {msg}", flush=True)

    by_chrom = {}
    for fn in sorted(os.listdir(TMP_DIR)):
        if not fn.endswith(".pkl"):
            continue
        with open(os.path.join(TMP_DIR, fn), "rb") as f:
            chrom, cnt, _ = pickle.load(f)
        d = by_chrom.setdefault(chrom, Counter())
        for (s, e), v in cnt.items():
            d[(s, e)] += v
    n_lines = 0
    with open(OUT_BED, "w") as out:
        for chrom in sorted(by_chrom):
            for (s, e), v in sorted(by_chrom[chrom].items()):
                if v < MIN_KEEP:
                    continue
                bs = e - s - ANCHOR  # second block start relative to chromStart
                out.write("\t".join(map(str, [
                    chrom, s - ANCHOR, e + ANCHOR, "junc", v, ".",
                    s - ANCHOR, e + ANCHOR, "0", 2, f"{ANCHOR},{ANCHOR}", f"0,{bs}",
                ])) + "\n")
                n_lines += 1
    total = sum(sum(d.values()) for d in by_chrom.values())
    kept = sum(v for d in by_chrom.values() for v in d.values() if v >= MIN_KEEP)
    print(f"[junc] raw junction-spanning reads total={total}; wrote {n_lines} "
          f"junctions (count>={MIN_KEEP}, covering {kept} read-introns) -> {OUT_BED}", flush=True)


if __name__ == "__main__":
    main()
