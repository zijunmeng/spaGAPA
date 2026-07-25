#!/usr/bin/env python3
"""
sAPA-RegNet Component 1 / Track B -- site-level miRNA APA regulatory annotation.

For each in-scope multi-PAS gene, determine which miRNA binding sites (predicted
by LOCAL miRanda on mature.fa) fall in the "APA-sensitive region" (3'UTR
sequence 3' of the proximal PAS, lost on a proximal APA switch) vs the "common
region" (3'UTR sequence 5' of the proximal PAS, always present).

The miRNA scan is run on the gene's 3'UTR sequence (gencode v44 3'UTR exons,
strand-correct), partitioned at the proximal PAS coordinate. This keeps miRanda
runtime bounded (3'UTRs are ~1-5 kb, not the full 100 kb gene body) and matches
where miRNAs biologically bind.

Outputs:
  pipeline_output/sapa_regnet/apa_regulatory_annotation_mirna.csv
  pipeline_output/sapa_regnet/track_b_mirna/  (per-gene miranda raw + parsed, logs)

Switching to PROXIMAL APA = cleavage at proximal PAS => LOSES mirna_in_apa_region.
"""
from __future__ import annotations

import gzip
import os
import re
import subprocess
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
from pyfaidx import Fasta

# ---------------------------------------------------------------------------
# Config (hostname-aware paths per CLAUDE.md; we are on S91 here)
# ---------------------------------------------------------------------------
HOSTNAME = os.popen("hostname -s").read().strip()
TMPDIR = {
    "S90": "/s2/mengzijun/tmp",
    "S91": "/s3/mengzijun/tmp",
    "S97": "/s972/mengzijun/tmp",
    "S98": "/s982/mengzijun/tmp",
}.get(HOSTNAME, "/tmp")
os.environ["TMPDIR"] = TMPDIR
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
Path(TMPDIR).mkdir(parents=True, exist_ok=True)

PROJ = Path("/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA")
OUT_ROOT = PROJ / "pipeline_output" / "sapa_regnet"
OUT_DIR = OUT_ROOT / "track_b_mirna"
OUT_DIR.mkdir(parents=True, exist_ok=True)

APA_SITES = PROJ / "data" / "processed" / "gse220442_gsm6801751_scapatrap" / "apa_sites.csv.gz"
DIFF_APA = PROJ / "pipeline_output" / "gse220442_differential_apa_unified" / "results.csv"
SVAPA = PROJ / "pipeline_output" / "svapa" / "gse220442_svapa_genelevel.csv"
GENE_BOUNDS = OUT_ROOT / "_gene_bounds_gencode_v44.tsv.gz"  # reuse from Track A
UTR3 = OUT_DIR / "utr3_per_gene.tsv.gz"  # built by build_utr3_gencode.py

GENOME_FA = Path("/s1/SHARE/00_ref_genecode/refdata-gex-GRCh38-2024-A/fasta/genome.fa")
MATURE_FA = Path("/s1/SHARE/apadata/MicroRNA/mature.fa")
MATURE_HSA = OUT_DIR / "mature_hsa.fa"  # filtered, built once
MIRANDA = Path.home() / "anaconda3/envs/samtools/bin/miranda"

MIRANDA_SCORE = 160.0   # high-confidence site cutoff (miRanda "good" hits).
# sc 140 (repo default) is too permissive on 3'UTRs (~2200 unique miRNAs bind a
# single 16 kb 3'UTR = near-degenerate); sc 160 yields confident, interpretable
# site sets (e.g. APP=19, MAPT=71 unique miRNAs).
MIRANDA_ENERGY = 1.0    # kcal/mol (miranda default; no extra energy cutoff)

# Scope sizes (CPU budget). With 16-way parallel miranda the full sig diff-APA
# set runs in minutes, so we include ALL significant diff-APA genes (padj<0.05)
# plus the top SVAPA genes, rather than only the top-50 pooled subset.
N_DIFF_REPLICATED = None        # all direction_agree==1 & sig (subset of sig)
N_DIFF_POOLED_TOP = None        # None => use ALL sig diff-APA genes
N_SVAPA_TOP = 50
N_WORKERS = min(16, max(2, (os.cpu_count() or 8) // 8))  # parallel miranda procs


# ---------------------------------------------------------------------------
# 1. Load PAS + build per-gene regions (mirror Track A logic exactly)
# ---------------------------------------------------------------------------
def load_apa_sites() -> pd.DataFrame:
    with gzip.open(APA_SITES, "rt") as fh:
        df = pd.read_csv(fh)
    df["coord"] = df["coord"].astype(int)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return df


def load_gene_bounds() -> pd.DataFrame:
    """gencode v44 gene bounds (gene_start/gene_end in genome coords)."""
    with gzip.open(GENE_BOUNDS, "rt") as fh:
        gb = pd.read_csv(fh, sep="\t")
    gb = gb.rename(columns={"gene_start": "gencode_start", "gene_end": "gencode_end"})
    return gb[["gene_id", "gencode_start", "gencode_end"]]


def compute_regions(sites: pd.DataFrame, gene_bounds: pd.DataFrame) -> pd.DataFrame:
    """For each multi-PAS gene compute proximal/distal PAS and the two regions,
    using the GENCODE gene bounds for the TSS-side common-region boundary.

    Strand-aware (gene body runs TSS -> TES):
      + strand: TSS at low coord, TES at high coord.
         proximal PAS (near TSS) = MIN coord ; distal PAS = MAX coord
         APA-sensitive region = [proximal, distal]   (3' of proximal, toward TES)
         common region        = [gene_start, proximal]
      - strand: TSS at high coord, TES at low coord.
         proximal PAS (near TSS) = MAX coord ; distal PAS = MIN coord
         APA-sensitive region = [distal, proximal]   (3' of proximal, toward TES)
         common region        = [proximal, gene_end]
    Region coordinates are reported as (start<=end) in genome forward space;
    strand correction happens only at FASTA extraction.
    """
    gb = gene_bounds.set_index("gene_id")
    rows = []
    for gene_id, g in sites.groupby("gene_id"):
        if len(g) < 2:
            continue
        if gene_id not in gb.index:
            continue
        gene_name = g["gene_name"].iloc[0]
        chrom = g["chr"].iloc[0]
        strand = g["strand"].iloc[0]
        coords = g["coord"].tolist()
        gstart = int(gb.loc[gene_id, "gencode_start"])
        gend = int(gb.loc[gene_id, "gencode_end"])

        if strand == "+":
            prox = min(coords)
            dist = max(coords)
            apa_start, apa_end = prox, dist
            com_start, com_end = gstart, prox
        else:
            prox = max(coords)
            dist = min(coords)
            apa_start, apa_end = dist, prox
            com_start, com_end = prox, gend

        rows.append({
            "gene_id": gene_id,
            "gene_name": gene_name,
            "chr": chrom,
            "strand": strand,
            "n_pas": len(coords),
            "proximal_pas_coord": prox,
            "distal_pas_coord": dist,
            "apa_region_start": apa_start,
            "apa_region_end": apa_end,
            "common_region_start": com_start,
            "common_region_end": com_end,
            "gene_start": gstart,
            "gene_end": gend,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 2. Scope: which genes to run
# ---------------------------------------------------------------------------
def scope_gene_names() -> set[str]:
    scope = set()
    # Tier 1 + Tier 2: differential APA -- include ALL significant genes
    diff = pd.read_csv(DIFF_APA)
    sig = diff[(diff["padj_ttest"].fillna(1) < 0.05)]
    rep = sig[sig["direction_agree"] == 1]["gene"].tolist()
    scope.update(rep)
    # all significant pooled genes (CPU-feasible with parallel miranda)
    if N_DIFF_POOLED_TOP is None:
        scope.update(sig["gene"].tolist())
    else:
        pooled = sig.assign(a=sig["log2fc"].abs()).sort_values("a", ascending=False)
        scope.update(pooled.head(N_DIFF_POOLED_TOP)["gene"].tolist())
    # Tier 3: SVAPA (by morans_i, nominal p<0.05)
    try:
        sv = pd.read_csv(SVAPA)
        sv_sig = sv[sv["pvalue"].fillna(1) < 0.05].sort_values("morans_i", ascending=False)
        scope.update(sv_sig.head(N_SVAPA_TOP)["gene_name"].tolist())
    except Exception as e:
        print(f"[warn] SVAPA load failed: {e}", file=sys.stderr)
    return scope


# ---------------------------------------------------------------------------
# 3. FASTA extraction (strand-correct via pyfaidx)
# ---------------------------------------------------------------------------
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def fetch_seq(genome: Fasta, chrom: str, start: int, end: int, strand: str) -> str:
    """1-based inclusive coords; returns 5'->3' on the gene strand.
    pyfaidx Sequence[start:end] is 0-based half-open; with as_raw=True it
    returns a plain str."""
    chrom = chrom if chrom in genome else chrom.replace("chr", "")
    seq = genome[chrom][start - 1:end]
    if hasattr(seq, "seq"):
        seq = seq.seq
    seq = str(seq).upper()
    if strand == "-":
        seq = revcomp(seq)
    return seq


def load_utr3() -> pd.DataFrame:
    """Collapsed 3'UTR intervals per gene (built by build_utr3_gencode.py)."""
    with gzip.open(UTR3, "rt") as fh:
        return pd.read_csv(fh, sep="\t")


def build_region_seq(genome: Fasta, chrom: str, strand: str,
                     intervals: list, proximal_pas: int):
    """Given a gene's merged 3'UTR intervals (genome forward coords), split them
    at the proximal PAS into APA-sensitive (3' of proximal) and common (5' of
    proximal) pieces, fetch each piece's strand-correct sequence, concatenate in
    transcription order (5'->3'), and return:
        apa_seq, apa_map (list of (seq_pos0, genome_start, genome_end))
        com_seq, com_map

    seq positions are 1-based within the concatenated region (miranda reports
    1-based target positions), so maps store genome coords per 1-based offset.
    """
    apa_pieces = []  # (genome_start, genome_end)
    com_pieces = []
    for (s, e) in sorted(intervals):
        if strand == "+":
            # transcription 5'->3' = low->high coord.
            # 3' of proximal PAS (toward TES) = coords > proximal_pas  => APA-sensitive
            if s > proximal_pas:
                apa_pieces.append((s, e))
            elif e <= proximal_pas:
                com_pieces.append((s, e))
            else:  # interval straddles proximal PAS
                com_pieces.append((s, proximal_pas))
                apa_pieces.append((proximal_pas + 1, e))
        else:  # '-': transcription 5'->3' = high->low coord.
            # 3' of proximal PAS (toward TES) = coords < proximal_pas => APA-sensitive
            if e < proximal_pas:
                apa_pieces.append((s, e))
            elif s >= proximal_pas:
                com_pieces.append((s, e))
            else:
                apa_pieces.append((s, proximal_pas - 1))
                com_pieces.append((proximal_pas, e))

    def assemble(pieces):
        # order in transcription direction so concatenation is 5'->3'
        if strand == "+":
            pieces = sorted(pieces)
        else:
            pieces = sorted(pieces, reverse=True)
        seq_parts = []
        posmap = {}  # 1-based offset -> genome coord
        offset = 1
        for (s, e) in pieces:
            sub = fetch_seq(genome, chrom, s, e, strand)
            for i, base in enumerate(sub):
                # genome coord of this base on the forward strand
                if strand == "+":
                    gcoord = s + i
                else:
                    gcoord = e - i
                posmap[offset] = gcoord
                offset += 1
            seq_parts.append(sub)
        return "".join(seq_parts), posmap

    apa_seq, apa_map = assemble(apa_pieces)
    com_seq, com_map = assemble(com_pieces)
    return apa_seq, apa_map, com_seq, com_map, apa_pieces, com_pieces


# ---------------------------------------------------------------------------
# 4. miRanda: build hsa mature + run per region + parse
# ---------------------------------------------------------------------------
def build_hsa_mature():
    if MATURE_HSA.exists() and MATURE_HSA.stat().st_size > 0:
        return
    n = 0
    with open(MATURE_FA) as fh, open(MATURE_HSA, "w") as out:
        keep = False
        for line in fh:
            if line.startswith(">"):
                keep = line.startswith(">hsa-")
                if keep:
                    n += 1
            if keep:
                out.write(line)
    print(f"[info] wrote {n} hsa mature miRNAs -> {MATURE_HSA}", file=sys.stderr)


_HIT_RE = re.compile(r"^>(\S+)\t(\S+)\t([0-9.]+)\t(-?[0-9.]+)\t(\d+)\s+(\d+)\t(\d+)\s+(\d+)\t(\d+)\t")


def run_miranda(target_fa: Path, out_txt: Path) -> bool:
    """Run miranda hsa-mature vs a single target FASTA; return True if ran."""
    if out_txt.exists() and out_txt.stat().st_size > 0:
        return True
    cmd = [
        str(MIRANDA), str(MATURE_HSA), str(target_fa),
        "-sc", str(MIRANDA_SCORE),
        "-en", str(MIRANDA_ENERGY),
        "-out", str(out_txt),
    ]
    try:
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       check=False, timeout=1800)
        return True
    except subprocess.TimeoutExpired:
        print(f"[warn] miranda timeout on {target_fa.name}", file=sys.stderr)
        return False


def _miranda_worker(args):
    """Module-level worker for ProcessPoolExecutor (must be picklable)."""
    fa, out = args
    ok = run_miranda(Path(fa), Path(out))
    return (fa, ok)


def parse_miranda(out_txt: Path) -> list[dict]:
    """Parse per-hit lines: >mirna<TAB>target<TAB>score<TAB>energy<TAB>
    q_start q_end<TAB>r_start r_end<TAB>align_len<TAB>%<TAB>%"""
    hits = []
    if not out_txt.exists():
        return hits
    with open(out_txt, errors="replace") as fh:
        for line in fh:
            if line.startswith(">") and not line.startswith(">>"):
                m = _HIT_RE.match(line.rstrip("\n"))
                if not m:
                    continue
                mirna, target, score, energy = m.group(1), m.group(2), float(m.group(3)), float(m.group(4))
                q_start, q_end = int(m.group(5)), int(m.group(6))
                r_start, r_end = int(m.group(7)), int(m.group(8))
                align_len = int(m.group(9))
                hits.append({
                    "mirna": mirna,
                    "target_region_id": target,
                    "score": score,
                    "energy_kcal": energy,
                    "target_start": r_start,
                    "target_end": r_end,
                    "mirna_start": q_start,
                    "mirna_end": q_end,
                    "align_len": align_len,
                })
    return hits


# ---------------------------------------------------------------------------
# 5. Main
# ---------------------------------------------------------------------------
def main():
    t0 = time.time()
    print(f"[info] host={HOSTNAME} tmpdir={TMPDIR}", file=sys.stderr)
    if not MIRANDA.exists():
        print(f"[FATAL] miranda not found at {MIRANDA}", file=sys.stderr)
        sys.exit(1)

    build_hsa_mature()
    if not UTR3.exists():
        print(f"[FATAL] 3'UTR table not found: {UTR3}. Run build_utr3_gencode.py first.", file=sys.stderr)
        sys.exit(1)
    genome = Fasta(str(GENOME_FA), sequence_always_upper=True, as_raw=True)

    sites = load_apa_sites()
    gene_bounds = load_gene_bounds()
    regions = compute_regions(sites, gene_bounds)
    utr3 = load_utr3()
    # gene_id -> list of merged 3'UTR intervals (genome forward coords)
    utr_by_gene = {gid: list(zip(g["utr3_start"].tolist(), g["utr3_end"].tolist()))
                   for gid, g in utr3.groupby("gene_id")}
    print(f"[info] {len(regions)} multi-PAS genes total; "
          f"{len(utr_by_gene)} genes with 3'UTR annotation", file=sys.stderr)

    scope = scope_gene_names()
    print(f"[info] scope gene-name set size: {len(scope)}", file=sys.stderr)

    # Scope = multi-PAS genes whose gene_name is in the scope set
    in_scope = regions[regions["gene_name"].isin(scope)].copy()
    # require a 3'UTR annotation
    in_scope = in_scope[in_scope["gene_id"].isin(utr_by_gene)].copy()
    in_scope = in_scope.sort_values(["gene_name"]).reset_index(drop=True)
    print(f"[info] in-scope multi-PAS genes WITH 3'UTR to annotate: {len(in_scope)}",
          file=sys.stderr)
    in_scope.to_csv(OUT_DIR / "_scope_genes.tsv", sep="\t", index=False)

    # ---- Phase 1: build all per-gene region FASTAs + position maps (fast) ----
    gene_tasks = []   # list of dicts with everything needed for parse/assembly
    for _, r in in_scope.iterrows():
        gn = r["gene_name"]
        gid = r["gene_id"]
        intervals = utr_by_gene.get(gid, [])
        apa_seq, apa_map, com_seq, com_map, apa_pieces, com_pieces = build_region_seq(
            genome, r["chr"], r["strand"], intervals, r["proximal_pas_coord"])

        apa_fa = OUT_DIR / f"{gn}.apa.fa"
        com_fa = OUT_DIR / f"{gn}.common.fa"
        with open(apa_fa, "w") as fh:
            fh.write(f">{gn}__apa\n{apa_seq}\n")
        with open(com_fa, "w") as fh:
            fh.write(f">{gn}__common\n{com_seq}\n")
        gene_tasks.append({
            "row": r, "apa_map": apa_map, "com_map": com_map,
            "apa_seq_len": len(apa_seq), "com_seq_len": len(com_seq),
            "apa_fa": apa_fa, "com_fa": com_fa,
            "apa_out": OUT_DIR / f"{gn}.apa.miranda.txt",
            "com_out": OUT_DIR / f"{gn}.common.miranda.txt",
        })
    print(f"[info] Phase 1 done: built {len(gene_tasks)} gene FASTA pairs "
          f"({time.time()-t0:.0f}s)", file=sys.stderr)

    # ---- Phase 2: run miranda in parallel across all (gene, region) jobs ----
    jobs = []  # (fa_path_str, out_path_str)
    for t in gene_tasks:
        if t["apa_seq_len"] >= 8:
            jobs.append((str(t["apa_fa"]), str(t["apa_out"])))
        if t["com_seq_len"] >= 8:
            jobs.append((str(t["com_fa"]), str(t["com_out"])))
    print(f"[info] Phase 2: running {len(jobs)} miranda jobs with "
          f"{N_WORKERS} workers...", file=sys.stderr)
    t_mir = time.time()
    done = 0
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        futs = {ex.submit(_miranda_worker, j): j for j in jobs}
        for fut in as_completed(futs):
            fa, ok = fut.result()
            done += 1
            if done % 10 == 0:
                print(f"[info]   miranda {done}/{len(jobs)} "
                      f"({time.time()-t_mir:.0f}s)", file=sys.stderr)
    print(f"[info] Phase 2 done: miranda runtime {time.time()-t_mir:.0f}s "
          f"({(time.time()-t_mir)/60:.1f} min) for {len(jobs)} jobs",
          file=sys.stderr)

    # ---- Phase 3: parse + assemble per-gene output ----
    all_rows = []        # final per-gene summary
    all_sites_rows = []  # per-site detail (long)
    n_apa_sites = 0
    n_common_sites = 0
    n_genes_with_apa_mir = 0
    n_no_apa_utr = 0

    def genome_coords_for(h, posmap):
        gs = posmap.get(h["target_start"])
        ge = posmap.get(h["target_end"])
        if gs is not None and ge is not None:
            return (min(gs, ge), max(gs, ge))
        return (None, None)

    for t in gene_tasks:
        r = t["row"]
        gn = r["gene_name"]
        gid = r["gene_id"]
        apa_hits = parse_miranda(t["apa_out"]) if t["apa_seq_len"] >= 8 else []
        com_hits = parse_miranda(t["com_out"]) if t["com_seq_len"] >= 8 else []

        # dedup miRNA within a region (keep best score)
        apa_best = {}
        for h in apa_hits:
            if h["mirna"] not in apa_best or h["score"] > apa_best[h["mirna"]]["score"]:
                apa_best[h["mirna"]] = h
        com_best = {}
        for h in com_hits:
            if h["mirna"] not in com_best or h["score"] > com_best[h["mirna"]]["score"]:
                com_best[h["mirna"]] = h

        apa_mirnas = sorted(apa_best.keys())
        com_mirnas = sorted(m for m in com_best.keys() if m not in apa_best)

        for m, h in apa_best.items():
            gs, ge = genome_coords_for(h, t["apa_map"])
            all_sites_rows.append({
                "gene_name": gn, "gene_id": gid, "region": "apa_sensitive",
                "mirna": m, "score": h["score"], "energy_kcal": h["energy_kcal"],
                "site_start_genome": gs, "site_end_genome": ge,
                "chr": r["chr"], "strand": r["strand"],
                "proximal_pas_coord": r["proximal_pas_coord"],
            })
        for m, h in com_best.items():
            gs, ge = genome_coords_for(h, t["com_map"])
            all_sites_rows.append({
                "gene_name": gn, "gene_id": gid, "region": "common",
                "mirna": m, "score": h["score"], "energy_kcal": h["energy_kcal"],
                "site_start_genome": gs, "site_end_genome": ge,
                "chr": r["chr"], "strand": r["strand"],
                "proximal_pas_coord": r["proximal_pas_coord"],
            })

        n_apa_sites += len(apa_best)
        n_common_sites += len(com_best)
        if apa_best:
            n_genes_with_apa_mir += 1
        if t["apa_seq_len"] == 0:
            n_no_apa_utr += 1

        all_rows.append({
            "gene_name": gn,
            "gene_id": gid,
            "chr": r["chr"],
            "strand": r["strand"],
            "n_pas": r["n_pas"],
            "proximal_pas_coord": r["proximal_pas_coord"],
            "distal_pas_coord": r["distal_pas_coord"],
            "apa_region_start": r["apa_region_start"],
            "apa_region_end": r["apa_region_end"],
            "common_region_start": r["common_region_start"],
            "common_region_end": r["common_region_end"],
            "mirna_in_apa_region": ";".join(apa_mirnas),
            "mirna_in_common": ";".join(com_mirnas),
            "n_mirna_apa": len(apa_mirnas),
            "n_mirna_common": len(com_mirnas),
            "utr3_apa_len_bp": t["apa_seq_len"],
            "utr3_common_len_bp": t["com_seq_len"],
            "gene_start": r["gene_start"],
            "gene_end": r["gene_end"],
        })
    n_done = len(gene_tasks)

    summary = pd.DataFrame(all_rows)
    sites_long = pd.DataFrame(all_sites_rows)

    out_csv = OUT_ROOT / "apa_regulatory_annotation_mirna.csv"
    summary.to_csv(out_csv, index=False)
    sites_long.to_csv(OUT_DIR / "mirna_sites_long.csv", index=False)

    elapsed = time.time() - t0
    # ---- sanity / report ----
    rep = []
    rep.append("sAPA-RegNet Component 1 / Track B -- miRNA site-level APA annotation")
    rep.append("=" * 70)
    rep.append(f"Host: {HOSTNAME}   TMPDIR: {TMPDIR}")
    rep.append(f"Genome: {GENOME_FA.name}   miRanda: {MIRANDA} (v3.3a)")
    rep.append(f"miRanda score cutoff: {MIRANDA_SCORE}   energy cutoff: {MIRANDA_ENERGY}")
    rep.append(f"Parallel miranda workers: {N_WORKERS}")
    rep.append(f"hsa mature miRNAs used: {sum(1 for _ in open(MATURE_HSA) if _.startswith('>'))}")
    rep.append("")
    rep.append(f"Multi-PAS genes total: {len(regions)}")
    rep.append(f"Scope (all sig diff-APA padj<0.05, SVAPA top{N_SVAPA_TOP}): "
               f"{len(scope)} gene names")
    rep.append(f"Genes annotated (multi-PAS in scope WITH 3'UTR): {n_done}")
    rep.append(f"  (genes whose 3'UTR lies entirely 5' of proximal PAS: {n_no_apa_utr})")
    rep.append(f"Total miRNA sites in APA-sensitive regions (3'UTR 3' of prox PAS): {n_apa_sites}")
    rep.append(f"Total miRNA sites in common regions (3'UTR 5' of prox PAS):        {n_common_sites}")
    rep.append(f"Genes with >=1 miRNA in APA region: {n_genes_with_apa_mir} / {n_done}")
    rep.append(f"miRanda runtime: {elapsed:.0f}s ({elapsed/60:.1f} min)")
    rep.append("")
    rep.append("Note: 'switching to PROXIMAL APA' = cleavage at proximal PAS =>")
    rep.append("      LOSES the mirna_in_apa_region set. Switching to DISTAL keeps them.")
    rep.append("")
    rep.append("Example AD differential genes:")
    diff = pd.read_csv(DIFF_APA)
    ad_examples = ["APP", "MAPT", "AAK1", "SV2B", "YWHAZ", "APLP2", "APBA2", "ANK2"]
    for gn in ad_examples:
        sub = summary[summary["gene_name"] == gn]
        if sub.empty:
            continue
        s = sub.iloc[0]
        rep.append(f"  {gn} ({s['chr']} {s['strand']}, n_pas={s['n_pas']}): "
                   f"prox={s['proximal_pas_coord']} dist={s['distal_pas_coord']} "
                   f"APA-region=[{s['apa_region_start']}-{s['apa_region_end']}]")
        rep.append(f"      n_mirna APA={s['n_mirna_apa']} common={s['n_mirna_common']}")
        apa_list = s['mirna_in_apa_region']
        if apa_list:
            rep.append(f"      mirna_in_apa_region: {apa_list[:200]}")
        else:
            rep.append(f"      mirna_in_apa_region: (none)")
    rep.append("")
    rep.append("Files written:")
    rep.append(f"  {out_csv}")
    rep.append(f"  {OUT_DIR / 'mirna_sites_long.csv'}  (per-site detail, genome coords)")
    rep.append(f"  {UTR3}  (gencode v44 per-gene 3'UTR intervals)")
    rep.append(f"  {OUT_DIR / '_scope_genes.tsv'}")
    rep.append(f"  {OUT_DIR}/<gene>.(apa|common).fa + .miranda.txt  (raw per-gene)")

    report_txt = OUT_DIR / "track_b_mirna_report.txt"
    report_txt.write_text("\n".join(rep))
    print("\n" + "\n".join(rep))
    print(f"\n[done] report -> {report_txt}", file=sys.stderr)


if __name__ == "__main__":
    main()
