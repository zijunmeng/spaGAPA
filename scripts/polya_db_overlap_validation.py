#!/usr/bin/env python3
"""
Experiment 2a: orthogonal validation of scAPAtrap PAS calls against the
PolyA_DB v4.1 catalog (Yu et al., NAR 2026, PMID 41316728).

Question
--------
Do the PAS coordinates called de novo by scAPAtrap on our 9-GSE + Stereo-seq
corpus land on independently curated polyadenylation sites?

Protocol (per sample)
---------------------
1. Load peaks_meta.csv.gz (peakID, chr, start, end, strand, coord).  The
   representative cleavage position is `coord` (= peak end on '+', peak start
   on '-'; verified against every input file).
2. Load the species-matched PolyA_DB v4.1 `*.PAS.max.tsv` catalog and index
   PAS positions per (chromosome, strand).
3. Observed metrics (same chr + same strand only):
   - nearest-distance distribution from each scAPAtrap PAS to the closest
     catalog PAS;
   - hit rate within +/-50 bp of a catalog PAS.
4. Random-background control: same number of sites per chromosome, drawn
   uniformly over the catalog's per-chromosome coordinate range, strands
   assigned from the sample's own per-chromosome strand proportions.
   Primary draw seed=42; 100 replicates (seeds 42..141) give mean +/- sd and
   an empirical p-value for the observed hit rate.
* Mouse Visium samples GSE169749 / GSE263303 were aligned to GRCm39; their
  PAS coordinates are lifted to mm10 with UCSC liftOver
  (mm39ToMm10.over.chain) before comparison.  Unlifted sites are counted and
  excluded from the hit-rate denominator.

Outputs (pipeline_output/polya_db_overlap/)
-------------------------------------------
overlap_by_sample.csv     per-sample hit rates, distances, background, enrichment
distance_distribution.png histogram + CDF, observed vs background, human/mouse
distance_distributions.npz raw pooled distance arrays behind the figure
overlap_summary.json      machine-readable summary incl. provenance + strata

Notes
-----
* PolyA_DB v4.1 coordinates are hg38 (human) / mm10 (mouse); our scAPAtrap
  runs use the same builds (chr-prefixed names match; chrM/ scaffolds absent
  from the catalog are dropped and reported).
* Run with ~/anaconda3/envs/spagapa/bin/python.
"""
from __future__ import annotations

import hashlib
import json
import os
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
OUT_DIR = os.path.join(REPO, "pipeline_output", "polya_db_overlap")
DB_DIR = os.path.join(REPO, "data", "external", "polya_db_v4")

SPECIES_DB = {
    "human": os.path.join(DB_DIR, "human", "HumanPas", "hg38.PAS.max.tsv"),
    "mouse": os.path.join(DB_DIR, "mouse", "MousePas", "mm10.PAS.max.tsv"),
}
DB_URL = "https://exon.apps.wistar.org/polya_db/v4/download/4.1/"

FLANK = 50            # +/- bp window for the primary hit-rate metric
N_REPS = 100          # random-background replicates (seeds 42..42+N_REPS-1)
BASE_SEED = 42
HIST_MAX = 1000       # bp; distance range shown in the figure
FIG_BINS = 50

# Curated representative corpus: 12 Visium samples spanning the 8 Visium GSE
# accessions (same list as scripts/calibrate_uncertainty_all_datasets.py),
# + 4 GSE293464 retinal-organoid Stereo-seq samples (9th GSE),
# + 2 GSE263789 Stereo-seq mouse-brain samples.
# (sample, gse, species, tissue, platform, genome build, peaks_meta path)
DB_BUILD = "mm10"
_S = "raw_scapatrap/peaks_meta.csv.gz"
SAMPLES = [
    ("gse237183_gsm7596587", "GSE237183", "human", "glioma", "visium",
    "hg38", f"pipeline_output/gse237183_GSM7596587_scapatrap/{_S}"),
    ("gse237183_gsm7596595", "GSE237183", "human", "glioma", "visium",
    "hg38", f"pipeline_output/gse237183_GSM7596595_scapatrap/{_S}"),
    ("gse237183_gsm7596604", "GSE237183", "human", "glioma", "visium",
    "hg38", f"pipeline_output/gse237183_GSM7596604_scapatrap/{_S}"),
    ("gse183456_gsm6047774", "GSE183456", "human", "kidney", "visium",
    "hg38", f"pipeline_output/gse183456_GSM6047774_scapatrap/{_S}"),
    ("gse179572_gsm5420751", "GSE179572", "human", "brain_metastasis", "visium",
    "hg38", f"pipeline_output/gse179572_GSM5420751_scapatrap/{_S}"),
    ("gse220442_gsm6801751", "GSE220442", "human", "AD_brain_PFC", "visium",
    "hg38", f"pipeline_output/gse220442_GSM6801751_scapatrap/{_S}"),
    ("gse220442_gsm6801753", "GSE220442", "human", "AD_brain_PFC", "visium",
    "hg38", f"pipeline_output/gse220442_GSM6801753_scapatrap/{_S}"),
    ("gse206391_gsm6252925", "GSE206391", "human", "skin_psoriasis", "visium",
    "hg38", f"pipeline_output/gse206391_GSM6252925_scapatrap/{_S}"),
    ("gse338525_gsm9876373", "GSE338525", "human", "liver", "visium",
    "hg38", f"pipeline_output/gse338525_GSM9876373_scapatrap/{_S}"),
    ("gse338525_gsm9876374", "GSE338525", "human", "liver", "visium",
    "hg38", f"pipeline_output/gse338525_GSM9876374_scapatrap/{_S}"),
    ("gse169749_gsm5213483", "GSE169749", "mouse", "colon_DSS", "visium",
    "GRCm39", f"pipeline_output/gse169749_gsm5213483_scapatrap/{_S}"),
    ("gse263303_gsm8189356", "GSE263303", "mouse", "brain_Nf1", "visium",
    "GRCm39", f"pipeline_output/gse263303_GSM8189356_scapatrap/{_S}"),
    ("gse293464_gsm8882884", "GSE293464", "human", "retinal_organoid", "stereo",
    "hg38", "pipeline_output/gse293464_retina/GSM8882884_scapatrap_raw/peaks_meta.csv.gz"),
    ("gse293464_gsm8882885", "GSE293464", "human", "retinal_organoid", "stereo",
    "hg38", "pipeline_output/gse293464_retina/GSM8882885_scapatrap_raw/peaks_meta.csv.gz"),
    ("gse293464_gsm8882886", "GSE293464", "human", "retinal_organoid", "stereo",
    "hg38", "pipeline_output/gse293464_retina/GSM8882886_scapatrap_raw/peaks_meta.csv.gz"),
    ("gse293464_gsm8882887", "GSE293464", "human", "retinal_organoid", "stereo",
    "hg38", "pipeline_output/gse293464_retina/GSM8882887_scapatrap_raw/peaks_meta.csv.gz"),
    ("gse263789_stereo_pilot", "GSE263789", "mouse", "AD_brain", "stereo",
    DB_BUILD, "pipeline_output/gse263789_stereo_pilot/scapatrap_raw/peaks_meta.csv.gz"),
    ("gse263789_wt_control", "GSE263789", "mouse", "AD_brain_WT", "stereo",
    DB_BUILD, "pipeline_output/gse263789_wt_control/scapatrap_raw/peaks_meta.csv.gz"),
]


LIFTOVER_BIN = "/s1/mengzijun/pkgs/kent/liftOver"
LIFTOVER_CHAIN = os.path.join(DB_DIR, "mm39ToMm10.over.chain.gz")


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_polya_db(path: str):
    """Parse PAS_ID = 'chr:strand:pos' into per-(chr,strand) sorted arrays."""
    df = pd.read_csv(path, sep="\t", usecols=["PAS_ID"])
    parts = df["PAS_ID"].str.split(":", expand=True)
    df = df.assign(chr=parts[0], strand=parts[1], pos=pd.to_numeric(parts[2]))
    db = {}
    for (chrom, strand), grp in df.groupby(["chr", "strand"]):
        db[(chrom, strand)] = np.sort(grp["pos"].to_numpy(dtype=np.int64))
    # per-chromosome coordinate range (any strand) for background sampling
    bounds = df.groupby("chr")["pos"].agg(["min", "max"])
    log(f"  {os.path.basename(path)}: {len(df):,} PAS, "
        f"{df['chr'].nunique()} chromosomes")
    return db, bounds, len(df)


def nearest_distances(positions: np.ndarray, db_sorted: np.ndarray) -> np.ndarray:
    """Distance from each query position to the nearest catalog position."""
    idx = np.searchsorted(db_sorted, positions)
    left = db_sorted[np.clip(idx - 1, 0, len(db_sorted) - 1)]
    right = db_sorted[np.clip(idx, 0, len(db_sorted) - 1)]
    return np.minimum(np.abs(positions - left), np.abs(positions - right))


def load_sample(path: str, db):
    df = pd.read_csv(path)
    n_raw = len(df)
    chrs_db = {c for (c, _s) in db}
    df = df[df["chr"].isin(chrs_db)]
    dropped = sorted(set(pd.read_csv(path)["chr"]) - chrs_db)
    return df[["chr", "strand", "coord"]], n_raw, dropped

def liftover_mm39_mm10(sample: str, df: pd.DataFrame):
    """UCSC-liftOver the sample's (chr, strand, coord) sites GRCm39 -> mm10."""
    import subprocess
    import tempfile

    tmp = tempfile.mkdtemp(prefix=f"lift_{sample}_", dir=OUT_DIR)
    bed_in = os.path.join(tmp, "in.bed")
    bed_out = os.path.join(tmp, "out.bed")
    bed_un = os.path.join(tmp, "unmapped.bed")
    with open(bed_in, "w") as fh:
        for row in df.itertuples(index=False):
            fh.write(f"{row.chr}\t{int(row.coord) - 1}\t{int(row.coord)}\t"
                     f"{row.chr}:{row.strand}:{int(row.coord)}\t0\t"
                     f"{row.strand}\n")
    subprocess.run([LIFTOVER_BIN, bed_in, LIFTOVER_CHAIN, bed_out, bed_un],
                   check=True, capture_output=True)
    out = pd.read_csv(bed_out, sep="\t", header=None,
                      names=["chr", "start", "end", "name", "score", "strand"])
    out["coord"] = out["end"]  # 1-based position restored from BED end
    out = out[["chr", "strand", "coord"]]
    n_unmapped = len(df) - len(out)
    for f in (bed_in, bed_out, bed_un):
        os.remove(f)
    os.rmdir(tmp)
    return out.reset_index(drop=True), n_unmapped


def random_background(sample: pd.DataFrame, bounds: pd.DataFrame, seed: int):
    """Same #sites/chromosome, uniform over the catalog chr range,
    strands drawn from the sample's per-chromosome strand proportions."""
    rng = np.random.default_rng(seed)
    rows = []
    for chrom, grp in sample.groupby("chr"):
        n = len(grp)
        lo, hi = bounds.loc[chrom, "min"], bounds.loc[chrom, "max"]
        pos = rng.integers(lo, hi + 1, size=n)
        p_plus = float((grp["strand"] == "+").mean())
        strand = np.where(rng.random(n) < p_plus, "+", "-")
        rows.append(pd.DataFrame({"chr": chrom, "strand": strand, "coord": pos}))
    return pd.concat(rows, ignore_index=True)


def distances_vs_db(query: pd.DataFrame, db) -> np.ndarray:
    out = []
    for (chrom, strand), grp in query.groupby(["chr", "strand"]):
        arr = db.get((chrom, strand))
        if arr is None or arr.size == 0:
            continue
        out.append(nearest_distances(grp["coord"].to_numpy(np.int64), arr))
    return np.concatenate(out) if out else np.array([], dtype=np.int64)


def make_figure(pooled: dict, out_path: str) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5))
    colors = {"human": "#1f77b4", "mouse": "#d62728"}
    bins = np.linspace(0, HIST_MAX, FIG_BINS + 1)
    centers = 0.5 * (bins[:-1] + bins[1:])

    ax = axes[0]
    for sp in ("human", "mouse"):
        d = pooled[("observed", sp)].astype(float)
        d = d[d <= HIST_MAX]
        w = np.histogram(d, bins=bins, density=True)[0]
        ax.plot(centers, w, color=colors[sp], lw=2,
                label=f"{sp} observed (n={len(d):,})")
        b = pooled[("background", sp)].astype(float)
        b = b[b <= HIST_MAX]
        wb = np.histogram(b, bins=bins, density=True)[0]
        ax.plot(centers, wb, color=colors[sp], lw=1.4, ls="--", alpha=0.8,
                label=f"{sp} random background (n={len(b):,})")
    ax.set_yscale("log")
    ax.set_xlabel("Distance to nearest PolyA_DB v4.1 PAS (bp)")
    ax.set_ylabel("Density (log scale)")
    ax.set_title("A  Nearest-distance distribution (0-1000 bp)")
    ax.axvline(FLANK, color="grey", lw=0.8, ls=":")
    ax.legend(fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    for sp in ("human", "mouse"):
        d = np.sort(pooled[("observed", sp)].astype(float))
        cdf = np.arange(1, len(d) + 1) / len(d)
        ax.plot(d, cdf, color=colors[sp], lw=2, label=f"{sp} observed")
        b = np.sort(pooled[("background", sp)].astype(float))
        cdfb = np.arange(1, len(b) + 1) / len(b)
        ax.plot(b, cdfb, color=colors[sp], lw=1.4, ls="--", alpha=0.8,
                label=f"{sp} random background")
    ax.axvline(FLANK, color="grey", lw=0.8, ls=":")
    hit = {sp: float((pooled[("observed", sp)] <= FLANK).mean())
           for sp in ("human", "mouse")}
    ax.set_xlim(0, HIST_MAX)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Distance to nearest PolyA_DB v4.1 PAS (bp)")
    ax.set_ylabel("Cumulative fraction")
    ax.set_title("B  Empirical CDF (50 bp hit: "
                 f"human {hit['human']:.1%}, mouse {hit['mouse']:.1%})",
                 fontsize=10)

    fig.suptitle("scAPAtrap PAS vs PolyA_DB v4.1 catalog - orthogonal overlap "
                 f"validation ({len(SAMPLES)} samples, 9 GSE + Stereo-seq)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    log("Loading PolyA_DB v4.1 catalogs ...")
    catalogs, catalog_meta = {}, {}
    for sp, path in SPECIES_DB.items():
        db, bounds, n = load_polya_db(path)
        catalogs[sp] = (db, bounds)
        catalog_meta[sp] = {
            "file": os.path.relpath(path, REPO),
            "n_pas": int(n),
            "sha256": sha256(path),
            "source_url": DB_URL + ("HumanPas.zip" if sp == "human"
                                    else "MousePas.zip"),
        }

    pooled: dict = {}
    per_sample_rows, dropped_chrs, lift_stats = [], {}, {}
    for sample, gse, sp, tissue, platform, build, rel in SAMPLES:
        path = os.path.join(REPO, rel)
        db, bounds = catalogs[sp]
        peaks, n_raw, dropped = load_sample(path, db)
        n_unlifted = 0
        if build == "GRCm39":
            peaks, n_unlifted = liftover_mm39_mm10(sample, peaks)
            lift_stats[sample] = {
                "from": "GRCm39", "to": "mm10",
                "n_input": int(n_raw), "n_lifted": int(len(peaks)),
                "n_unlifted": int(n_unlifted),
            }
            # some sites land on mm10 non-primary contigs absent from the
            # catalog (e.g. chr4_GL456216_random) -- drop and report them
            chrs_db = {c for (c, _s) in db}
            n_alt = int((~peaks["chr"].isin(chrs_db)).sum())
            if n_alt:
                lift_stats[sample]["n_mapped_to_nonprimary_contig"] = n_alt
                peaks = peaks[peaks["chr"].isin(chrs_db)].reset_index(drop=True)

        if dropped:
            dropped_chrs[sample] = dropped
        obs = distances_vs_db(peaks, db)

        # random background: primary seed-42 draw + N_REPS replicates
        bg_dists, bg_hit = [], []
        for i in range(N_REPS):
            bg = random_background(peaks, bounds, BASE_SEED + i)
            d = distances_vs_db(bg, db)
            bg_dists.append(d)
            bg_hit.append(float((d <= FLANK).mean()))
        bg_primary = bg_dists[0]  # seed == BASE_SEED == 42
        bg_hit = np.array(bg_hit)
        obs_hit = float((obs <= FLANK).mean())
        p_emp = float((1 + np.sum(bg_hit >= obs_hit)) / (N_REPS + 1))

        row = {
            "sample": sample, "gse": gse, "species": sp, "tissue": tissue,
            "platform": platform, "genome_build": build,
            "n_peaks_raw": n_raw, "n_peaks_used": int(len(peaks)),
            "n_unlifted": int(n_unlifted),
            "db_pas_total": catalog_meta[sp]["n_pas"],
            "hit_rate_50bp": obs_hit,
            "bg_hit_rate_50bp_seed42": float(bg_hit[0]),
            "bg_hit_rate_50bp_mean": float(bg_hit.mean()),
            "bg_hit_rate_50bp_sd": float(bg_hit.std()),
            "enrichment_fold": obs_hit / bg_hit.mean(),
            "empirical_p_100reps": p_emp,
            "median_dist_bp": float(np.median(obs)),
            "mean_dist_bp": float(obs.mean()),
            "p90_dist_bp": float(np.percentile(obs, 90)),
            "frac_dist_le_100bp": float((obs <= 100).mean()),
            "frac_dist_le_500bp": float((obs <= 500).mean()),
            "frac_dist_le_1000bp": float((obs <= 1000).mean()),
        }
        per_sample_rows.append(row)
        key_o, key_b = ("observed", sp), ("background", sp)
        pooled.setdefault(key_o, []).append(obs)
        pooled.setdefault(key_b, []).append(bg_primary)
        log(f"{sample} ({build}): n={len(peaks):,} hit@50bp={obs_hit:.3%} "
            f"bg={bg_hit.mean():.3%}±{bg_hit.std():.3%} "
            f"enrich={row['enrichment_fold']:.1f}x p={p_emp:.4f} "
            f"median_d={row['median_dist_bp']:.0f}bp")

    res = pd.DataFrame(per_sample_rows)
    res.to_csv(os.path.join(OUT_DIR, "overlap_by_sample.csv"), index=False)

    pooled = {k: np.concatenate(v) for k, v in pooled.items()}
    np.savez_compressed(
        os.path.join(OUT_DIR, "distance_distributions.npz"),
        **{"_".join(k): v for k, v in pooled.items()})
    make_figure(pooled, os.path.join(OUT_DIR, "distance_distribution.png"))

    def strat(df: pd.DataFrame, col: str) -> dict:
        out = {}
        for key, grp in df.groupby(col):
            out[str(key)] = {
                "n_samples": int(len(grp)),
                "n_pas_total": int(grp["n_peaks_used"].sum()),
                "mean_hit_rate_50bp": float(grp["hit_rate_50bp"].mean()),
                "min_hit_rate_50bp": float(grp["hit_rate_50bp"].min()),
                "max_hit_rate_50bp": float(grp["hit_rate_50bp"].max()),
                "mean_bg_hit_rate_50bp": float(
                    grp["bg_hit_rate_50bp_mean"].mean()),
                "mean_enrichment_fold": float(
                    grp["enrichment_fold"].mean()),
                "median_of_median_dist_bp": float(
                    grp["median_dist_bp"].median()),
            }
        return out

    obs_all = np.concatenate([pooled[("observed", s)] for s in ("human", "mouse")])
    bg_all = np.concatenate([pooled[("background", s)]
                             for s in ("human", "mouse")])
    summary = {
        "experiment": "2a_polya_db_overlap",
        "question": ("Do scAPAtrap PAS calls overlap independently curated "
                     "PolyA_DB v4.1 polyadenylation sites?"),
        "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
        "polya_db": {
            "version": "v4.1 (NAR 2026 database issue, PMID 41316728)",
            "site": "https://exon.apps.wistar.org/polya_db/v4/",
            "catalogs": catalog_meta,
            "genome_builds": {"human": "hg38", "mouse": "mm10"},
        },
        "samples": {
            "n_total": len(SAMPLES),
            "n_gse_accessions": len({s[1] for s in SAMPLES}),
            "gse_accessions": sorted({s[1] for s in SAMPLES}),
            "platforms": sorted({s[4] for s in SAMPLES}),
            "genome_builds": sorted({s[5] for s in SAMPLES}),
            "liftover": {
                "tool": "UCSC liftOver (kent)",
                "chain": "mm39ToMm10.over.chain.gz "
                         "(hgdownload.soe.ucsc.edu/goldenPath/mm39/liftOver/)",
                "note": "mouse Visium samples aligned to GRCm39 lifted to "
                        "mm10 before comparison; unlifted sites excluded",
                "per_sample": lift_stats,
            },
            "dropped_chromosomes_not_in_db": dropped_chrs,
        },
        "method": {
            "pas_position": "peaks_meta.coord (peak end on '+', start on '-')",
            "match_rule": "same chromosome + same strand, nearest catalog PAS",
            "hit_window_bp": FLANK,
            "background": ("per chromosome, same site count, uniform over the "
                           "catalog coordinate range, strands from the "
                           "sample's per-chr strand proportions"),
            "n_background_replicates": N_REPS,
            "seed_primary_draw": BASE_SEED,
            "seed_replicates": f"{BASE_SEED}..{BASE_SEED + N_REPS - 1}",
        },
        "headline": {
            "pooled_hit_rate_50bp": float((obs_all <= FLANK).mean()),
            "pooled_bg_hit_rate_50bp": float((bg_all <= FLANK).mean()),
            "pooled_enrichment_fold": float(
                (obs_all <= FLANK).mean() / (bg_all <= FLANK).mean()),
            "pooled_median_dist_bp": float(np.median(obs_all)),
            "pooled_frac_le_500bp": float((obs_all <= 500).mean()),
            "by_species": {
                sp: {
                    "hit_rate_50bp": float(
                        (pooled[("observed", sp)] <= FLANK).mean()),
                    "bg_hit_rate_50bp": float(
                        (pooled[("background", sp)] <= FLANK).mean()),
                    "median_dist_bp": float(
                        np.median(pooled[("observed", sp)])),
                    "frac_dist_le_500bp": float(
                        (pooled[("observed", sp)] <= 500).mean()),
                } for sp in ("human", "mouse")
            },
        },
        "stratified": {
            "by_species": strat(res, "species"),
            "by_tissue": strat(res, "tissue"),
            "by_platform": strat(res, "platform"),
        },
        "per_sample": per_sample_rows,
        "outputs": [
            "pipeline_output/polya_db_overlap/overlap_by_sample.csv",
            "pipeline_output/polya_db_overlap/distance_distribution.png",
            "pipeline_output/polya_db_overlap/distance_distributions.npz",
            "pipeline_output/polya_db_overlap/overlap_summary.json",
        ],
    }
    with open(os.path.join(OUT_DIR, "overlap_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    log("Done. Outputs in pipeline_output/polya_db_overlap/")
    print(res[["sample", "hit_rate_50bp", "bg_hit_rate_50bp_mean",
               "enrichment_fold", "median_dist_bp"]].to_string(index=False))


if __name__ == "__main__":
    main()
