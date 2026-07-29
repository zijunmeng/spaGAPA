#!/usr/bin/env python
"""Build marker-defined weak biological labels for a prepared spaGAPA dataset.

The labels produced here are not pathology ground truth. They are an
interpretable marker proxy for biological consistency benchmarks before manual
ROI annotation is available.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_MARKERS = {
    "tumor_epithelial": [
        "EPCAM",
        "KRT7",
        "KRT8",
        "KRT18",
        "KRT19",
        "MUC1",
        "CEACAM5",
        "SLC34A2",
        "NAPSA",
        "TACSTD2",
        "ALDH1A1",
    ],
    "stromal_fibroblast": [
        "COL1A1",
        "COL1A2",
        "COL3A1",
        "DCN",
        "LUM",
        "FN1",
        "ACTA2",
        "TAGLN",
    ],
    "immune_myeloid": [
        "PTPRC",
        "LST1",
        "C1QA",
        "C1QB",
        "C1QC",
        "CD68",
        "LYZ",
        "FCER1G",
    ],
    "lymphoid": [
        "CD3D",
        "CD3E",
        "TRAC",
        "CD2",
        "MS4A1",
        "CD79A",
        "NKG7",
        "GNLY",
    ],
    "endothelial": [
        "PECAM1",
        "VWF",
        "KDR",
        "CLDN5",
        "RAMP2",
        "EMCN",
    ],
    "brain_glial": [
        "GFAP",
        "MBP",
        "PLP1",
        "AQP4",
        "SLC1A2",
        "SNAP25",
        "RBFOX3",
        "OLIG1",
        "OLIG2",
    ],
    "proliferative": [
        "MKI67",
        "TOP2A",
        "PCNA",
        "TYMS",
        "STMN1",
    ],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-dir", required=True)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--score-threshold", type=float, default=0.25)
    parser.add_argument("--margin-threshold", type=float, default=0.10)
    parser.add_argument(
        "--top-genes",
        type=int,
        default=0,
        help="Optional top variable genes to keep before scoring. Default keeps all markers only.",
    )
    parser.add_argument(
        "--metadata-layer-alias",
        action="store_true",
        help="Also write metadata layer/layer_source columns for legacy benchmark runners.",
    )
    return parser.parse_args()


def read_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    return df


def load_expression(dataset_dir: Path) -> pd.DataFrame:
    expr = read_table(dataset_dir / "expression_matrix.csv")
    genes = pd.read_csv(dataset_dir / "expression_genes.csv")
    if "gene" in genes.columns and "gene_symbol" in genes.columns:
        gene_to_symbol = dict(zip(genes["gene"].astype(str), genes["gene_symbol"].astype(str)))
        expr["__gene_symbol__"] = [gene_to_symbol.get(g, g).upper() for g in expr.index.astype(str)]
        expr = expr.set_index("__gene_symbol__", drop=True)
        expr.index.name = None
        expr = expr.groupby(expr.index).sum()
    else:
        expr.index = expr.index.astype(str).str.upper()
        expr = expr.groupby(expr.index).sum()
    return expr


def log_normalize(expr: pd.DataFrame) -> pd.DataFrame:
    counts = expr.astype(float)
    library = counts.sum(axis=0).replace(0, np.nan)
    cpm = counts.divide(library, axis=1) * 1e4
    return np.log1p(cpm.fillna(0.0))


def zscore_rows(expr: pd.DataFrame) -> pd.DataFrame:
    values = expr.values.astype(float)
    mean = values.mean(axis=1, keepdims=True)
    std = values.std(axis=1, keepdims=True)
    std[std == 0] = 1.0
    z = (values - mean) / std
    z = np.clip(z, -5.0, 5.0)
    return pd.DataFrame(z, index=expr.index, columns=expr.columns)


def score_markers(z_expr: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    scores = {}
    used = {}
    for label, markers in DEFAULT_MARKERS.items():
        present = [marker for marker in markers if marker in z_expr.index]
        used[label] = present
        if present:
            scores[label] = z_expr.loc[present].mean(axis=0)
        else:
            scores[label] = pd.Series(0.0, index=z_expr.columns)
    return pd.DataFrame(scores), used


def assign_labels(scores: pd.DataFrame, score_threshold: float, margin_threshold: float) -> pd.DataFrame:
    ordered = np.sort(scores.values, axis=1)
    max_score = ordered[:, -1]
    second_score = ordered[:, -2] if scores.shape[1] > 1 else np.zeros_like(max_score)
    margin = max_score - second_score
    best = scores.idxmax(axis=1).astype(str).values
    weak = np.where((max_score >= score_threshold) & (margin >= margin_threshold), best, "ambiguous")
    confidence = 1.0 / (1.0 + np.exp(-margin))
    out = pd.DataFrame(
        {
            "spot_id": scores.index.astype(str),
            "marker_weak_label": weak,
            "marker_weak_confidence": confidence,
            "marker_weak_margin": margin,
            "marker_weak_max_score": max_score,
            "marker_weak_second_score": second_score,
        }
    )
    for col in scores.columns:
        out[f"score_{col}"] = scores[col].values
    return out


def plot_labels(coords: pd.DataFrame, labels: pd.DataFrame, figure_dir: Path) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    merged = coords.merge(labels, on="spot_id", how="inner")
    label_order = [x for x in DEFAULT_MARKERS if x in set(merged["marker_weak_label"])]
    if "ambiguous" in set(merged["marker_weak_label"]):
        label_order.append("ambiguous")
    colors = plt.get_cmap("tab10")(np.linspace(0, 1, max(10, len(label_order))))
    color_map = {label: colors[i] for i, label in enumerate(label_order)}

    fig, ax = plt.subplots(figsize=(6.5, 5.8))
    for label in label_order:
        sub = merged[merged["marker_weak_label"] == label]
        ax.scatter(sub["x"], sub["y"], s=12, color=color_map[label], label=f"{label} ({len(sub)})", linewidths=0)
    ax.set_aspect("equal", adjustable="box")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.legend(frameon=False, fontsize=7, loc="center left", bbox_to_anchor=(1.02, 0.5))
    ax.set_title("GSE179572 marker-defined weak labels")
    fig.tight_layout()
    fig.savefig(figure_dir / "marker_weak_label_map.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    score_cols = [c for c in labels.columns if c.startswith("score_")]
    n_cols = 3
    n_rows = int(np.ceil(len(score_cols) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.0 * n_cols, 3.6 * n_rows))
    axes = np.asarray(axes).reshape(-1)
    for ax, col in zip(axes, score_cols):
        vals = merged[col].values
        sc = ax.scatter(merged["x"], merged["y"], c=vals, cmap="viridis", s=10, linewidths=0)
        ax.set_aspect("equal", adjustable="box")
        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(col.replace("score_", ""))
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    for ax in axes[len(score_cols) :]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(figure_dir / "marker_score_maps.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    counts = labels["marker_weak_label"].value_counts().reindex(label_order).fillna(0)
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    ax.bar(counts.index.astype(str), counts.values, color=[color_map[x] for x in counts.index])
    ax.set_ylabel("Spots")
    ax.set_xticklabels(counts.index.astype(str), rotation=35, ha="right")
    ax.set_title("Marker weak-label spot counts")
    fig.tight_layout()
    fig.savefig(figure_dir / "marker_weak_label_counts.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def update_metadata(dataset_dir: Path, labels: pd.DataFrame, layer_alias: bool) -> None:
    metadata_path = dataset_dir / "metadata.csv"
    metadata = pd.read_csv(metadata_path)
    keep = [c for c in metadata.columns if not c.startswith("marker_weak_")]
    keep = [c for c in keep if c not in {"layer", "layer_source"}]
    metadata = metadata[keep]
    merged = metadata.merge(
        labels[
            [
                "spot_id",
                "marker_weak_label",
                "marker_weak_confidence",
                "marker_weak_margin",
                "marker_weak_max_score",
            ]
        ],
        on="spot_id",
        how="left",
    )
    if layer_alias:
        merged["layer"] = merged["marker_weak_label"].fillna("ambiguous").astype(str)
        merged["layer_source"] = "marker_weak_label_not_pathology_gold_standard"
    merged.to_csv(metadata_path, index=False)


def update_qc(dataset_dir: Path, labels: pd.DataFrame, used_markers: dict[str, list[str]]) -> None:
    path = dataset_dir / "qc_summary.json"
    qc = json.loads(path.read_text()) if path.exists() else {}
    counts = labels["marker_weak_label"].value_counts().to_dict()
    qc["weak_biological_labels_ready"] = True
    qc["weak_biological_label_source"] = "marker_defined_expression_proxy"
    qc["weak_biological_label_warning"] = "Not pathology ROI or gold-standard manual annotation."
    qc["n_marker_weak_label_spots"] = int(labels.shape[0])
    qc["marker_weak_label_counts"] = {str(k): int(v) for k, v in counts.items()}
    qc["marker_sets_used"] = {key: value for key, value in used_markers.items()}
    path.write_text(json.dumps(qc, indent=2))


def main() -> None:
    args = parse_args()
    dataset_dir = Path(args.dataset_dir)
    output_dir = Path(args.output_dir) if args.output_dir else dataset_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    coords = pd.read_csv(dataset_dir / "coordinates.csv")
    coords["spot_id"] = coords["spot_id"].astype(str)
    expr = load_expression(dataset_dir)
    shared = [spot for spot in coords["spot_id"].astype(str) if spot in expr.columns]
    if len(shared) < 10:
        raise ValueError("Too few shared spots between expression and coordinates")
    coords = coords[coords["spot_id"].isin(shared)].copy()
    expr = expr.loc[:, shared]

    marker_genes = sorted({m for markers in DEFAULT_MARKERS.values() for m in markers})
    marker_present = [gene for gene in marker_genes if gene in expr.index]
    expr = expr.loc[marker_present]
    norm = log_normalize(expr)
    z_expr = zscore_rows(norm)
    scores, used_markers = score_markers(z_expr)
    labels = assign_labels(scores, args.score_threshold, args.margin_threshold)
    labels.to_csv(output_dir / "marker_weak_labels.csv", index=False)

    plot_labels(coords, labels, output_dir / "figures")
    if output_dir.resolve() == dataset_dir.resolve():
        update_metadata(dataset_dir, labels, args.metadata_layer_alias)
        update_qc(dataset_dir, labels, used_markers)

    summary = {
        "dataset_dir": str(dataset_dir),
        "output_dir": str(output_dir),
        "score_threshold": float(args.score_threshold),
        "margin_threshold": float(args.margin_threshold),
        "n_spots": int(labels.shape[0]),
        "label_counts": {str(k): int(v) for k, v in labels["marker_weak_label"].value_counts().items()},
        "used_markers": used_markers,
        "figures": [
            str(output_dir / "figures" / "marker_weak_label_map.png"),
            str(output_dir / "figures" / "marker_score_maps.png"),
            str(output_dir / "figures" / "marker_weak_label_counts.png"),
        ],
    }
    (output_dir / "marker_weak_label_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
