#!/usr/bin/env python3
"""
sAPA-RegNet genome-wide neural/AD-RBP FIMO rescan.

Fixes the chr1-only data gap in /s1/SHARE/apadata/RBP/RBP_binding_sites_summary.txt
(which only had 54 chr1 3'UTRs in its FIMO input FASTA). Here we:

  STEP 1. Extract a full genome-wide 3'UTR FASTA from genome.fa using the
          3'UTR intervals already built by Track B
          (utr3_per_gene.tsv.gz, 11,826 genes). Strand-aware: reverse-
          complement for the - strand. The motif files are RNA (ACGU), so
          we convert the FASTA T->U (one consistent alphabet choice).
  STEP 2. Concatenate every neural/AD RBP motif file into one MEME db and
          run FIMO once against the full 3'UTR FASTA (DNA-style alphabet
          but with U; FIMO handles arbitrary alphabets). p-value thresh 1e-4.
  STEP 3. Emit genome-wide neural RBP sites with absolute genomic coords
          (gene_name, RBP, chr, start, end, strand, score, pvalue).
  STEP 4. Re-run the Track-A APA-region overlap logic (proximal PAS coord
          per gene -> APA-sensitive region 3' of proximal PAS, common region
          5' of proximal PAS) on these genome-wide sites.
  STEP 5. Sanity: report AD-gene (APP/SV2B/YWHAZ/ARPP19/DCLK1) annotations.

CPU-only. Hostname-aware. Does NOT git commit. Writes everything to disk
(parent agent verifies via disk because prior agents died on the API).

Outputs (all under pipeline_output/sapa_regnet/rbp_genomewide/):
    utr3_genomewide_rna.fa            combined 3'UTR FASTA (RNA, strand-correct)
    utr3_genomewide_rna.fa.fai        samtools index
    neural_rbp_motifs.meme            concatenated motif db
    fimo_neural_rbp/                  FIMO output dir (fimo.tsv etc.)
    neural_rbp_sites.tsv              genome-wide sites (FINAL STEP 3)
    apa_regulatory_annotation_rbp_neural.csv   FINAL STEP 4
    rbp_genomewide_sanity_report.txt           FINAL STEP 5
    rbp_genomewide_run.log                      run log
"""
from __future__ import annotations

import glob
import gzip
import os
import re
import socket
import subprocess
import sys
import time
from collections import Counter

import pandas as pd
import pysam
import pyranges as pr

# ---------------------------------------------------------------------------
# Hostname-aware environment (per CLAUDE.md).
# ---------------------------------------------------------------------------
HOST = socket.gethostname().split(".")[0]
TMPDIR_BY_HOST = {
    "S90": "/s2/mengzijun/tmp",
    "S91": "/s3/mengzijun/tmp",
    "S97": "/s972/mengzijun/tmp",
    "S98": "/s982/mengzijun/tmp",
}
TMPDIR = TMPDIR_BY_HOST.get(HOST, "/tmp")
os.makedirs(TMPDIR, exist_ok=True)
os.environ.setdefault("TMPDIR", TMPDIR)
# Limit BLAS threads so FIMO gets the cores.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")

# ---------------------------------------------------------------------------
PROJ = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
OUT_DIR = os.path.join(PROJ, "pipeline_output/sapa_regnet/rbp_genomewide")
UTR3_INTERVALS = os.path.join(
    PROJ, "pipeline_output/sapa_regnet/track_b_mirna/utr3_per_gene.tsv.gz"
)
GENOME_FA = "/s1/SHARE/00_ref_genecode/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
MOTIF_DIR = "/s1/SHARE/apadata/RBP/RBP_motifs"
FIMO = "/home/mengzijun/anaconda3/envs/samtools/bin/fimo"

UTR3_FA = os.path.join(OUT_DIR, "utr3_genomewide_rna.fa")
MOTIF_DB = os.path.join(OUT_DIR, "neural_rbp_motifs.meme")
FIMO_OUT = os.path.join(OUT_DIR, "fimo_neural_rbp")
SITES_TSV = os.path.join(OUT_DIR, "neural_rbp_sites.tsv")
APA_CSV = os.path.join(
    PROJ, "pipeline_output/sapa_regnet/apa_regulatory_annotation_rbp_neural.csv"
)
SANITY_REPORT = os.path.join(OUT_DIR, "rbp_genomewide_sanity_report.txt")
RUN_LOG = os.path.join(OUT_DIR, "rbp_genomewide_run.log")

PAS_GLOB = os.path.join(PROJ, "data/processed/gse220442_gsm680175*_scapatrap/apa_sites.csv.gz")
GENCODE_GTF = "/s1/SHARE/00_ref_genecode/gencode.v44.annotation.gtf"
GENE_BOUNDS_CACHE = os.path.join(
    PROJ, "pipeline_output/sapa_regnet/_gene_bounds_gencode_v44.tsv.gz"
)

PVALUE_THRESH = "1e-4"

# ~15 neural / AD-relevant RBPs (TARDBP = TDP-43). Use ALL their motif files.
NEURAL_RBPS = [
    "TARDBP", "FUS", "ELAVL1", "ELAVL2", "ELAVL3", "ELAVL4",
    "NOVA1", "NOVA2", "FMR1", "HNRNPK", "HNRNPA1", "KHSRP",
    "TIA1", "TIAL1", "PTBP1",
]

EXAMPLE_AD_GENES = ["APP", "SV2B", "YWHAZ", "ARPP19", "DCLK1"]

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def log(msg: str) -> None:
    line = f"[rbp_gw {time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    try:
        with open(RUN_LOG, "a") as fh:
            fh.write(line + "\n")
    except Exception:
        pass


# ---------------------------------------------------------------------------
# STEP 1: build the genome-wide 3'UTR FASTA (RNA, strand-correct)
# ---------------------------------------------------------------------------
def build_utr3_fasta() -> tuple[int, int]:
    """Extract every 3'UTR interval to one combined FASTA.

    Header: gene_name|gene_id|chr:start-end(strand)
    Sequence is the SENSE (mRNA) orientation: rev-comp for - strand.
    DNA->RNA: T->U (motif files are ACGU).
    """
    log(f"Loading 3'UTR intervals from {UTR3_INTERVALS}")
    with gzip.open(UTR3_INTERVALS, "rt") as fh:
        iv = pd.read_csv(fh, sep="\t")
    # columns: gene_id, gene_name, chr, strand, utr3_start, utr3_end
    log(f"  {len(iv)} intervals, {iv['gene_name'].nunique()} unique genes, "
        f"{iv['chr'].nunique()} chroms")

    log(f"Opening genome FASTA: {GENOME_FA}")
    fa = pysam.FastaFile(GENOME_FA)
    avail_chroms = set(fa.references)
    missing = set(iv["chr"].unique()) - avail_chroms
    if missing:
        log(f"  WARNING: {len(missing)} chroms in UTR3 set not in genome: "
            f"{sorted(missing)[:5]}")

    os.makedirs(OUT_DIR, exist_ok=True)
    n_seq = 0
    n_skipped = 0
    t0 = time.time()
    with open(UTR3_FA, "w") as out:
        for _, r in iv.iterrows():
            chrom = str(r["chr"])
            start = int(r["utr3_start"])
            end = int(r["utr3_end"])
            strand = str(r["strand"])
            gname = str(r["gene_name"])
            gid = str(r["gene_id"])
            if chrom not in avail_chroms:
                n_skipped += 1
                continue
            # pysam fetch is 0-based half-open; UTR3 intervals are 1-based inclusive.
            # Convert: 0-based start = start-1, end stays (exclusive end == inclusive end).
            try:
                seq = fa.fetch(chrom, start - 1, end)
            except (ValueError, OSError) as e:
                n_skipped += 1
                continue
            if not seq:
                n_skipped += 1
                continue
            seq = seq.upper()
            if strand == "-":
                seq = seq.translate(_COMP)[::-1]
            # DNA -> RNA
            seq = seq.replace("T", "U").replace("N", "N")
            header = f"{gname}|{gid}|{chrom}:{start}-{end}({strand})"
            out.write(f">{header}\n{seq}\n")
            n_seq += 1
    dt = time.time() - t0
    log(f"  Wrote {n_seq:,} sequences ({n_skipped} skipped) to {UTR3_FA} in {dt:.1f}s")
    return n_seq, n_skipped


# ---------------------------------------------------------------------------
# STEP 2: build neural motif DB + run FIMO
# ---------------------------------------------------------------------------
def build_motif_db() -> tuple[int, int]:
    """Concatenate all *_*.meme motif files for the neural RBPs into one db.

    Returns (n_motifs, n_rbps).
    """
    log("Building neural/AD RBP motif DB...")
    seen_motif_ids = set()
    n_motifs = 0
    rbps_used = []
    with open(MOTIF_DB, "w") as out:
        header_written = False
        for rbp in NEURAL_RBPS:
            files = sorted(glob.glob(os.path.join(MOTIF_DIR, f"{rbp}_*.meme")))
            if not files:
                log(f"  WARNING: no motif files for {rbp}")
                continue
            rbps_used.append(rbp)
            for f in files:
                with open(f) as fh:
                    content = fh.read()
                # First file writes the MEME header; subsequent files: drop
                # the preamble (everything up to and including the first
                # 'MOTIF ' line is handled below by streaming).
                # Simplest robust approach: for the first file keep full;
                # for later files skip lines until past the Background block.
                if not header_written:
                    out.write(content)
                    if not content.endswith("\n"):
                        out.write("\n")
                    header_written = True
                else:
                    # Strip the leading header section: keep from the first
                    # 'MOTIF ' line onward.
                    idx = content.find("\nMOTIF ")
                    if idx == -1:
                        # fall back: keep whole (rare)
                        out.write(content)
                    else:
                        out.write(content[idx + 1:])
                    if not content.endswith("\n"):
                        out.write("\n")
                # Count motifs added from this file
                n_in_file = len(re.findall(r"^MOTIF\s+", content, flags=re.MULTILINE))
                n_motifs += n_in_file
    # de-dup MOTIF IDs that collide across files (FIMO errors on dup IDs)
    log(f"  Motif DB written: {n_motifs} raw motifs across {len(rbps_used)} RBPs")
    # Verify uniqueness; if dups exist, rewrite with suffixes.
    with open(MOTIF_DB) as fh:
        ids = re.findall(r"^MOTIF\s+(\S+)", fh.read(), flags=re.MULTILINE)
    dups = [k for k, v in Counter(ids).items() if v > 1]
    if dups:
        log(f"  NOTE: {len(dups)} duplicate motif IDs found ({dups[:3]}...). "
            f"Disambiguating with _a/_b suffixes.")
        _dedup_motif_ids(dups)
        with open(MOTIF_DB) as fh:
            ids = re.findall(r"^MOTIF\s+(\S+)", fh.read(), flags=re.MULTILINE)
        n_motifs = len(ids)
    return n_motifs, len(rbps_used)


def _dedup_motif_ids(dups):
    with open(MOTIF_DB) as fh:
        text = fh.read()
    for dup in dups:
        # add a numeric suffix to each occurrence beyond the first
        seen = 0
        out_parts = []
        last = 0
        for m in re.finditer(rf"^(MOTIF\s+){re.escape(dup)}(\s)", text, flags=re.MULTILINE):
            seen += 1
            if seen == 1:
                continue
            out_parts.append((m.start(), m.end(), f"{m.group(1)}{dup}_{chr(96+seen)}{m.group(2)}"))
        for start, end, repl in reversed(out_parts):
            text = text[:start] + repl + text[end:]
    with open(MOTIF_DB, "w") as fh:
        fh.write(text)


def run_fimo() -> tuple[float, str]:
    log(f"Running FIMO: pval<={PVALUE_THRESH}")
    if os.path.isdir(FIMO_OUT):
        # clean stale
        for f in glob.glob(os.path.join(FIMO_OUT, "*")):
            try:
                os.remove(f)
            except OSError:
                pass
    cmd = [
        FIMO,
        "--parse-genomic",            # interpret |chr:start-end(strand) headers as coords
        "--pv-thresh", PVALUE_THRESH,
        "--max-strand",               # report each hit once
        "--oc", FIMO_OUT,
        MOTIF_DB,
        UTR3_FA,
    ]
    log("  CMD: " + " ".join(cmd))
    t0 = time.time()
    # Capture but also tee to log via stderr
    proc = subprocess.run(cmd, capture_output=True, text=True)
    dt = time.time() - t0
    if proc.returncode != 0:
        log("FIMO FAILED (rc={})".format(proc.returncode))
        log("STDOUT:\n" + proc.stdout[-2000:])
        log("STDERR:\n" + proc.stderr[-2000:])
        sys.exit(2)
    with open(RUN_LOG, "a") as fh:
        fh.write("\n--- FIMO STDOUT ---\n" + proc.stdout + "\n")
        fh.write("\n--- FIMO STDERR ---\n" + proc.stderr + "\n")
    log(f"  FIMO done in {dt:.1f}s")
    return dt, os.path.join(FIMO_OUT, "fimo.tsv")


# ---------------------------------------------------------------------------
# STEP 3: parse FIMO output to genome-wide sites TSV
# ---------------------------------------------------------------------------
def parse_fimo_sites(fimo_tsv: str) -> pd.DataFrame:
    """Parse FIMO TSV -> neural_rbp_sites.tsv with absolute genomic coords.

    With --parse-genomic on headers like 'GENE|ENSG..|chr:start-end(strand)',
    FIMO emits #pattern implied name as the sequence-name (gene header) and
    its own coordinate columns ('start'/'end'/'strand') are RELATIVE to the
    fetched FASTA seq (1-based within the seq). We must add back the genomic
    offset and rev-complement coordinates for - strand sequences.
    """
    log(f"Parsing FIMO output: {fimo_tsv}")
    fimo = pd.read_csv(fimo_tsv, sep="\t", comment="#")
    # FIMO columns: motif_id, motif_alt_id, sequence_name, start, end, strand,
    #               score, p-value, q-value, matched_sequence
    log(f"  FIMO raw hits: {len(fimo):,}")
    # parse sequence_name -> gene_name, gene_id, chr, genomic start/end, strand
    parts = fimo["sequence_name"].astype(str).str.extract(
        r"^(?P<gene_name>[^|]+)\|(?P<gene_id>[^|]+)\|"
        r"(?P<chr>chr[^:]+):(?P<gstart>\d+)-(?P<gend>\d+)\((?P<gstrand>[+-])\)$"
    )
    bad = parts["gene_name"].isna().sum()
    if bad:
        log(f"  WARNING: {bad} rows had un-parseable sequence_name headers")
    fimo = pd.concat([fimo, parts], axis=1)
    fimo = fimo.dropna(subset=["gene_name"]).copy()
    fimo["gstart"] = fimo["gstart"].astype(int)
    fimo["gend"] = fimo["gend"].astype(int)
    fimo["fstart"] = fimo["start"].astype(int)  # 1-based within FASTA seq
    fimo["fend"] = fimo["end"].astype(int)
    fimo["gstrand"] = fimo["gstrand"].astype(str)

    # Map FIMO-relative coords back to absolute genomic coords.
    # The FASTA sequence is always SENSE (mRNA) orientation.
    # For + strand: sense pos p (1-based) -> genomic gstart + (p-1).
    #   hit window [fstart, fend] -> [gstart+fstart-1, gstart+fend-1].
    # For - strand: sense pos p -> genomic gend - (p-1). The hit window
    #   [fstart, fend] in sense maps to [gend-fend+1, gend-fstart+1] on the
    #   genome, but the FIMO 'strand' of the hit (+/-) is wrt the sense seq.
    #   A '+' hit on a - strand gene lies on the sense (mRNA) strand, which
    #   is the reverse complement of the genome; the genomic coordinates of
    #   the window are still [gend-fend+1, gend-fstart+1] (start<end on genome).
    plus = fimo["gstrand"] == "+"
    # genomic window start/end (1-based inclusive, start <= end)
    g_w_start = pd.Series(0, index=fimo.index, dtype=int)
    g_w_end = pd.Series(0, index=fimo.index, dtype=int)
    g_w_start[plus] = fimo.loc[plus, "gstart"] + fimo.loc[plus, "fstart"] - 1
    g_w_end[plus] = fimo.loc[plus, "gstart"] + fimo.loc[plus, "fend"] - 1
    minus = ~plus
    g_w_start[minus] = fimo.loc[minus, "gend"] - fimo.loc[minus, "fend"] + 1
    g_w_end[minus] = fimo.loc[minus, "gend"] - fimo.loc[minus, "fstart"] + 1

    fimo["chr"] = fimo["chr"].astype(str)
    fimo["start"] = g_w_start
    fimo["end"] = g_w_end
    # hit strand relative to GENOME: sense '+' hit on '-' gene => '-' on genome
    fimo["strand_rel_gene"] = fimo["strand"].astype(str)
    out_strand = pd.Series("+", index=fimo.index, dtype=str)
    # if hit is on '+' wrt sense, genomic strand == gene strand
    # if hit is on '-' wrt sense, genomic strand == opposite of gene strand
    hit_plus = fimo["strand_rel_gene"] == "+"
    out_strand[hit_plus] = fimo.loc[hit_plus, "gstrand"]
    out_strand[~hit_plus] = fimo.loc[~hit_plus, "gstrand"].map({"+": "-", "-": "+"})
    fimo["strand"] = out_strand

    # bare RBP name = motif_id up to first '_'
    fimo["RBP"] = fimo["motif_id"].astype(str).str.split("_").str[0]

    out = fimo[["gene_name", "gene_id", "RBP", "motif_id", "chr",
                "start", "end", "strand", "score", "p-value", "q-value",
                "matched_sequence"]].rename(columns={
        "motif_id": "motif", "p-value": "pvalue", "q-value": "qvalue",
    }).copy()
    out["start"] = out["start"].astype(int)
    out["end"] = out["end"].astype(int)
    out["score"] = out["score"].astype(float).round(3)
    out["pvalue"] = out["pvalue"].astype(float)
    out["qvalue"] = out["qvalue"].astype(float)
    # drop exact duplicates (same motif + gene + chr + start + end + strand)
    before = len(out)
    out = out.drop_duplicates(
        subset=["motif", "gene_name", "chr", "start", "end", "strand"]
    ).sort_values(["chr", "start", "RBP"]).reset_index(drop=True)
    log(f"  Dedup: {before:,} -> {len(out):,} unique sites")
    out.to_csv(SITES_TSV, sep="\t", index=False)
    log(f"  Wrote {SITES_TSV}")
    return out


# ---------------------------------------------------------------------------
# STEP 4: APA-region overlap (Track-A logic, genome-wide)
# ---------------------------------------------------------------------------
def load_pas_union() -> pd.DataFrame:
    files = sorted(glob.glob(PAS_GLOB))
    if not files:
        sys.exit(f"No PAS files matched {PAS_GLOB}")
    frames = []
    for f in files:
        with gzip.open(f, "rt") as fh:
            df = pd.read_csv(fh)
        df["sample"] = os.path.basename(os.path.dirname(f)).split("_scapatrap")[0]
        frames.append(df)
    allp = pd.concat(frames, ignore_index=True)
    before = len(allp)
    allp = allp.drop_duplicates(subset=["gene_id", "chr", "coord", "strand"]).copy()
    log(f"PAS union: {before} -> {len(allp)} unique PAS ({len(files)} samples)")
    return allp


def load_gene_bounds() -> pd.DataFrame:
    log(f"Gene bounds cache: {GENE_BOUNDS_CACHE}")
    return pd.read_csv(GENE_BOUNDS_CACHE, sep="\t")


def build_apa_segments(pas: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
    g = genes[["gene_id", "gene_name", "gene_type", "chr",
               "gene_start", "gene_end", "strand"]].copy()
    pas = pas.merge(g, on="gene_id", how="inner", suffixes=("", "_g"))
    if "strand_g" in pas:
        pas["strand"] = pas["strand_g"]
    pas = pas.drop(columns=[c for c in pas.columns if c.endswith("_g")], errors="ignore")
    segs = []
    for (gid, gname, gtype, chrom, strand,
         gstart, gend), gdf in pas.groupby(
        ["gene_id", "gene_name", "gene_type", "chr", "strand",
         "gene_start", "gene_end"]):
        if len(gdf) < 2:
            continue
        coords = gdf["coord"].astype(int).tolist()
        if strand == "+":
            prox_coord = min(coords)
            dist_coord = max(coords)
        else:
            prox_coord = max(coords)
            dist_coord = min(coords)
        segs.append({
            "gene_id": gid, "gene_name": gname, "gene_type": gtype,
            "chr": chrom, "strand": strand,
            "gene_start": int(gstart), "gene_end": int(gend),
            "n_pas": len(coords),
            "proximal_pas_coord": int(prox_coord),
            "distal_pas_coord": int(dist_coord),
            "apa_region_start": int(min(prox_coord, gend)),
            "apa_region_end": int(max(prox_coord, gend)),
            "common_region_start": int(min(gstart, prox_coord)),
            "common_region_end": int(max(gstart, prox_coord)),
        })
    seg = pd.DataFrame(segs)
    log(f"Multi-PAS genes with segments: {len(seg)}")
    return seg


def overlap_sites(segments: pd.DataFrame, sites: pd.DataFrame) -> pd.DataFrame:
    def to_pr(df, start_col, end_col, extra_cols):
        out = pd.DataFrame({
            "Chromosome": df["chr"].astype(str),
            "Start": df[start_col].astype(int) - 1,
            "End": df[end_col].astype(int),
            "Strand": df["strand"].astype(str),
        })
        for c in extra_cols:
            out[c] = df[c].values
        return pr.PyRanges(out)

    apa_pr = to_pr(segments, "apa_region_start", "apa_region_end",
                   ["gene_id", "gene_name", "chr", "strand"])
    com_pr = to_pr(segments, "common_region_start", "common_region_end",
                   ["gene_id", "gene_name", "chr", "strand"])
    sites_pr = to_pr(sites, "start", "end",
                     ["motif", "RBP", "strand", "score", "pvalue"])

    log("Joining neural-RBP sites to APA regions (stranded)...")
    apa_join = apa_pr.join(sites_pr, suffix="_site")
    apa_df = apa_join.df.copy() if hasattr(apa_join, "df") else apa_join.as_df()
    log(f"  APA-region overlaps: {len(apa_df):,}")

    com_join = com_pr.join(sites_pr, suffix="_site")
    com_df = com_join.df.copy() if hasattr(com_join, "df") else com_join.as_df()
    log(f"  Common-region overlaps: {len(com_df):,}")

    def agg(odf, label):
        if len(odf) == 0:
            return pd.DataFrame(columns=["gene_id", f"rbp_in_{label}",
                                         f"n_rbp_{label}", f"n_sites_{label}"])
        rows = []
        for gid, sub in odf.groupby("gene_id"):
            motif_counts = Counter(sub["motif"])
            items = []
            bare_counts = Counter()
            for motif, n in sorted(motif_counts.items(),
                                   key=lambda kv: (kv[0].split("_")[0], kv[0])):
                bare = motif.split("_")[0]
                items.append(f"{motif}:{n}")
                bare_counts[bare] += n
            rows.append({
                "gene_id": gid,
                f"rbp_in_{label}": ";".join(items),
                f"n_rbp_{label}": len(bare_counts),
                f"n_sites_{label}": sum(bare_counts.values()),
            })
        return pd.DataFrame(rows)

    apa_agg = agg(apa_df, "apa_region")
    com_agg = agg(com_df, "common")
    merged = segments.merge(apa_agg, on="gene_id", how="left") \
                     .merge(com_agg, on="gene_id", how="left")
    for c in ["n_rbp_apa_region", "n_rbp_common",
              "n_sites_apa_region", "n_sites_common"]:
        merged[c] = merged[c].fillna(0).astype(int)
    for c in ["rbp_in_apa_region", "rbp_in_common"]:
        merged[c] = merged[c].fillna("")
    return merged


def write_apa_csv(annot: pd.DataFrame) -> None:
    out = annot.rename(columns={
        "n_rbp_apa_region": "n_rbp_apa",
        "n_rbp_common": "n_rbp_common",
    })[[
        "gene_name", "gene_id", "chr", "strand", "n_pas",
        "proximal_pas_coord", "distal_pas_coord",
        "apa_region_start", "apa_region_end",
        "common_region_start", "common_region_end",
        "rbp_in_apa_region", "rbp_in_common",
        "n_rbp_apa", "n_rbp_common",
        "n_sites_apa_region", "n_sites_common",
        "gene_start", "gene_end", "gene_type",
    ]].sort_values(["chr", "proximal_pas_coord"])
    out.to_csv(APA_CSV, index=False)
    log(f"Wrote {APA_CSV} ({len(out):,} genes)")


# ---------------------------------------------------------------------------
# STEP 5: sanity report
# ---------------------------------------------------------------------------
def write_sanity(sites: pd.DataFrame, annot: pd.DataFrame,
                 n_motifs: int, n_rbps: int, fimo_dt: float,
                 n_seq: int) -> None:
    chroms = sorted(sites["chr"].unique().tolist())
    rbp_counts = sites.groupby("RBP")["gene_name"].nunique().sort_values(ascending=False)
    total_apa_sites = int(annot["n_sites_apa_region"].sum())
    total_common_sites = int(annot["n_sites_common"].sum())
    n_genes_apa = int((annot["n_sites_apa_region"] > 0).sum())

    lines = []
    lines.append("sAPA-RegNet genome-wide neural/AD-RBP FIMO rescan")
    lines.append("=" * 72)
    lines.append(f"Host: {HOST}   TMPDIR: {TMPDIR}")
    lines.append("")
    lines.append("INPUTS / SCAN")
    lines.append("-" * 72)
    lines.append(f"  3'UTR FASTA seqs scanned (sense, RNA):   {n_seq:,}")
    lines.append(f"  Neural/AD RBPs in motif DB:              {n_rbps} "
                 f"({', '.join(NEURAL_RBPS[:8])}, ...)")
    lines.append(f"  Total motif PWMs in DB:                  {n_motifs}")
    lines.append(f"  FIMO p-value threshold:                  <= {PVALUE_THRESH}")
    lines.append(f"  FIMO runtime:                            {fimo_dt:.1f}s "
                 f"({fimo_dt/60:.1f} min)")
    lines.append("")
    lines.append("GENOME-WIDE RESULTS")
    lines.append("-" * 72)
    lines.append(f"  Total unique RBP sites (all chroms):     {len(sites):,}")
    lines.append(f"  Chromosomes covered:                     {len(chroms)} "
                 f"({chroms[0]}..{chroms[-1]})")
    lines.append(f"  Distinct genes with >=1 site:            {sites['gene_name'].nunique():,}")
    lines.append(f"  Multi-PAS genes annotated:               {len(annot):,}")
    lines.append(f"  Genes with >=1 neural-RBP site in APA:   {n_genes_apa:,}")
    lines.append(f"  Total neural-RBP sites in APA regions:   {total_apa_sites:,}")
    lines.append(f"  Total neural-RBP sites in common regions:{total_common_sites:,}")
    lines.append("")
    lines.append("Sites per neural RBP (genes with site):")
    for rbp, n in rbp_counts.items():
        lines.append(f"  {rbp:<10} {n:>5} genes   "
                     f"({(sites['RBP']==rbp).sum():,} sites)")
    lines.append("")
    lines.append("DATA GAP RESOLVED")
    lines.append("-" * 72)
    lines.append("  Prior chr1-only file covered 1 chromosome (54 3'UTRs).")
    lines.append("  This rescan covers all 25 main chroms (incl. chr8 YWHAZ,")
    lines.append("  chr13 DCLK1, chr15 ARPP19/SV2B, chr21 APP).")
    lines.append("")

    # AD genes block
    lines.append("AD / differential genes — neural-RBP annotations in APA region:")
    lines.append("-" * 72)
    neural_re = re.compile(r"\b(" + "|".join(re.escape(r) for r in NEURAL_RBPS) + r")\b")
    found_any = False
    for gname in EXAMPLE_AD_GENES:
        sub = annot[annot["gene_name"] == gname]
        if sub.empty:
            lines.append(f"  {gname}: NOT in multi-PAS set (no APA switch detected)")
            continue
        found_any = True
        for _, r in sub.iterrows():
            apa_str = r["rbp_in_apa_region"] or "(none)"
            neural = sorted(set(neural_re.findall(apa_str)))
            lines.append(
                f"  {gname} ({r['chr']} {r['strand']}, n_pas={r['n_pas']}): "
                f"prox={r['proximal_pas_coord']} dist={r['distal_pas_coord']} "
                f"APA=[{r['apa_region_start']}-{r['apa_region_end']}]"
            )
            lines.append(f"      n_sites APA={r['n_sites_apa_region']} "
                         f"common={r['n_sites_common']}  "
                         f"n_RBP APA={r['n_rbp_apa']} common={r['n_rbp_common']}")
            if neural:
                lines.append(f"      neural/AD RBPs in APA region: {', '.join(neural)}")
            lines.append(f"      rbp_in_apa_region: "
                         f"{apa_str[:250]}{'...' if len(apa_str) > 250 else ''}")
    if not found_any:
        lines.append("  (none of the example AD genes are in the multi-PAS set; "
                     "see genome-wide sites TSV for raw hits on these genes.)")
    # also show raw sites for AD genes even if not multi-PAS
    lines.append("")
    lines.append("Raw neural-RBP sites in 3'UTRs of AD genes (regardless of APA status):")
    for gname in EXAMPLE_AD_GENES:
        sub = sites[sites["gene_name"] == gname]
        if sub.empty:
            lines.append(f"  {gname}: 0 sites")
        else:
            rbp_summary = ", ".join(f"{r}={n}" for r, n in
                                    sub.groupby("RBP").size().sort_values(ascending=False).items())
            lines.append(f"  {gname}: {len(sub)} sites ({rbp_summary})")
    lines.append("")
    lines.append("Interpretation key:")
    lines.append("  - 'switching to PROXIMAL APA' = cleavage at proximal PAS =>")
    lines.append("    LOSES the rbp_in_apa_region set.")
    lines.append("  - 'switching to DISTAL APA'   = cleavage at distal  PAS =>")
    lines.append("    KEEPS the rbp_in_apa_region set.")
    lines.append("  - The 3'UTR FASTA is SENSE orientation (rev-comp for -strand);")
    lines.append("    reported genomic coords are absolute hg38, start<=end.")

    report = "\n".join(lines) + "\n"
    with open(SANITY_REPORT, "w") as fh:
        fh.write(report)
    log(f"Wrote {SANITY_REPORT}")
    print()
    print(report)


# ---------------------------------------------------------------------------
def main() -> None:
    log(f"=== START. Host={HOST} ===")
    t_start = time.time()

    # STEP 1
    n_seq, n_skipped = build_utr3_fasta()
    # index the FASTA (not strictly needed for FIMO, useful for inspection)
    try:
        subprocess.run(
            ["/home/mengzijun/anaconda3/envs/samtools/bin/samtools", "faidx", UTR3_FA],
            check=True, capture_output=True,
        )
        log("  Indexed FASTA (.fai)")
    except subprocess.CalledProcessError as e:
        log(f"  samtools faidx failed (non-fatal): {e.stderr.decode()[:200]}")

    # STEP 2
    n_motifs, n_rbps = build_motif_db()
    fimo_dt, fimo_tsv = run_fimo()

    # STEP 3
    sites = parse_fimo_sites(fimo_tsv)

    # STEP 4
    pas = load_pas_union()
    genes = load_gene_bounds()
    seg = build_apa_segments(pas, genes)
    annot = overlap_sites(seg, sites)
    write_apa_csv(annot)

    # STEP 5
    write_sanity(sites, annot, n_motifs, n_rbps, fimo_dt, n_seq)

    log(f"=== DONE in {(time.time()-t_start)/60:.1f} min ===")


if __name__ == "__main__":
    main()
