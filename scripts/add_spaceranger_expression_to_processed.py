#!/usr/bin/env python
"""Add Space Ranger expression counts to a prepared spaGAPA APA dataset."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATASET_DIR = PACKAGE_ROOT / "data/processed/gse179572_gsm5420751_scapatrap"
DEFAULT_H5 = PACKAGE_ROOT / "pipeline_output/gse179572_GSM5420751_sr/outs/filtered_feature_bc_matrix.h5"


def make_unique_names(names: list[str], feature_ids: list[str]) -> list[str]:
    """Return stable unique gene labels while preserving common symbols."""
    seen: dict[str, int] = {}
    unique: list[str] = []
    for name, feature_id in zip(names, feature_ids):
        if name not in seen:
            seen[name] = 1
            unique.append(name)
        else:
            seen[name] += 1
            unique.append(f"{name}|{feature_id}")
    return unique


def read_visium_h5(path: Path) -> pd.DataFrame:
    import h5py
    from scipy import sparse

    with h5py.File(path, "r") as handle:
        matrix_group = handle["matrix"]
        shape = tuple(int(x) for x in matrix_group["shape"][:])
        matrix = sparse.csc_matrix(
            (
                matrix_group["data"][:],
                matrix_group["indices"][:],
                matrix_group["indptr"][:],
            ),
            shape=shape,
        )
        barcodes = [
            x.decode() if isinstance(x, bytes) else str(x)
            for x in matrix_group["barcodes"][:]
        ]
        genes = [
            x.decode() if isinstance(x, bytes) else str(x)
            for x in matrix_group["features/name"][:]
        ]
        feature_ids = [
            x.decode() if isinstance(x, bytes) else str(x)
            for x in matrix_group["features/id"][:]
        ]

    unique_genes = make_unique_names(genes, feature_ids)
    expression = pd.DataFrame.sparse.from_spmatrix(
        matrix,
        index=pd.Index(unique_genes, name="gene"),
        columns=pd.Index(barcodes, dtype=str),
    )
    expression.attrs["gene_symbols"] = genes
    expression.attrs["feature_ids"] = feature_ids
    return expression


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-dir", default=str(DEFAULT_DATASET_DIR))
    parser.add_argument("--expression-h5", default=str(DEFAULT_H5))
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    dataset_dir = Path(args.dataset_dir)
    expression_h5 = Path(args.expression_h5)
    expression_path = dataset_dir / "expression_matrix.csv"
    gene_meta_path = dataset_dir / "expression_genes.csv"

    if expression_path.exists() and not args.force:
        raise FileExistsError(f"{expression_path} exists; use --force to overwrite")
    if not dataset_dir.exists():
        raise FileNotFoundError(dataset_dir)
    if not expression_h5.exists():
        raise FileNotFoundError(expression_h5)

    apa_cols = pd.read_csv(dataset_dir / "apa_matrix.csv", index_col=0, nrows=0).columns.astype(str)
    expression = read_visium_h5(expression_h5)
    shared = [spot for spot in apa_cols if spot in expression.columns]
    if len(shared) != len(apa_cols):
        missing = len(apa_cols) - len(shared)
        raise ValueError(f"Expression matrix is missing {missing} APA spots")

    expression = expression.loc[:, shared]
    expression.to_csv(expression_path)

    gene_meta = pd.DataFrame(
        {
            "gene": expression.index.astype(str),
            "gene_symbol": expression.attrs.get("gene_symbols", expression.index.astype(str)),
            "feature_id": expression.attrs.get("feature_ids", [""] * expression.shape[0]),
        }
    )
    gene_meta.to_csv(gene_meta_path, index=False)

    qc_path = dataset_dir / "qc_summary.json"
    qc = json.loads(qc_path.read_text()) if qc_path.exists() else {}
    qc.update(
        {
            "expression_ready": True,
            "n_expression_genes": int(expression.shape[0]),
            "n_expression_spots": int(expression.shape[1]),
            "n_expression_duplicate_gene_symbols": int(
                pd.Series(expression.attrs.get("gene_symbols", [])).duplicated().sum()
            ),
        }
    )
    qc.setdefault("files", {})["expression_matrix"] = "expression_matrix.csv"
    qc["files"]["expression_genes"] = "expression_genes.csv"
    notes = qc.setdefault("notes", [])
    notes.append(
        "Expression matrix was exported from Space Ranger filtered_feature_bc_matrix.h5 "
        "and aligned to apa_matrix.csv spot order."
    )
    notes.append(
        "Duplicate Space Ranger gene symbols were made unique as gene|feature_id; "
        "the original symbols are stored in expression_genes.csv."
    )
    qc_path.write_text(json.dumps(qc, indent=2))

    metadata_path = dataset_dir / "metadata.csv"
    metadata = pd.read_csv(metadata_path)
    metadata["expression_ready"] = True
    metadata["label_status"] = metadata.get("label_status", "missing_biological_labels")
    metadata.to_csv(metadata_path, index=False)

    print(f"Wrote {expression_path}")
    print(f"Expression genes: {expression.shape[0]}, spots: {expression.shape[1]}")


if __name__ == "__main__":
    main()
