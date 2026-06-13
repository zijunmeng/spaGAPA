#!/usr/bin/env python
"""Prepare a Visium expression-only candidate dataset for spaGAPA.

This script intentionally does not create ``apa_matrix.csv``. Expression-only
datasets are useful candidates for future APA calling from BAM/FASTQ, but they
must not be counted as true APA benchmarks until real APA/PAS evidence is added.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def read_expression_csv(path: Path, orientation: str) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    if orientation == "spots_by_genes":
        df = df.T
    return df


def read_visium_h5(path: Path) -> pd.DataFrame:
    try:
        import h5py
        from scipy import sparse
    except ImportError as exc:  # pragma: no cover - dependency guard
        raise ImportError("h5py and scipy are required to read 10x H5 files") from exc

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
        barcodes = [x.decode() if isinstance(x, bytes) else str(x) for x in matrix_group["barcodes"][:]]
        names = matrix_group["features/name"][:]
        genes = [x.decode() if isinstance(x, bytes) else str(x) for x in names]

    return pd.DataFrame(
        matrix.toarray(),
        index=pd.Index(genes, dtype=str),
        columns=pd.Index(barcodes, dtype=str),
    )


def read_visium_coordinates(spatial_dir: Path) -> pd.DataFrame:
    candidates = [
        spatial_dir / "tissue_positions.csv",
        spatial_dir / "tissue_positions_list.csv",
    ]
    path = next((candidate for candidate in candidates if candidate.exists()), None)
    if path is None:
        raise FileNotFoundError(f"No tissue_positions file found in {spatial_dir}")

    first = pd.read_csv(path, nrows=1, header=None)
    has_header = str(first.iloc[0, 0]).lower() in {"barcode", "spot_id"}
    if first.shape[1] == 6 and not has_header:
        df = pd.read_csv(
            path,
            header=None,
            names=["spot_id", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"],
        )
    else:
        df = pd.read_csv(path)
        rename = {}
        if "barcode" in df.columns and "spot_id" not in df.columns:
            rename["barcode"] = "spot_id"
        if "pxl_row_in_fullres" in df.columns:
            rename["pxl_row_in_fullres"] = "pxl_row"
        if "pxl_col_in_fullres" in df.columns:
            rename["pxl_col_in_fullres"] = "pxl_col"
        df = df.rename(columns=rename)

    required = {"spot_id", "in_tissue", "pxl_row", "pxl_col"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Missing required coordinate columns: {missing}")

    df = df[df["in_tissue"].astype(int) == 1].copy()
    return pd.DataFrame(
        {
            "spot_id": df["spot_id"].astype(str),
            "x": df["pxl_col"].astype(float),
            "y": df["pxl_row"].astype(float),
        }
    )


def select_top_variable_genes(expression: pd.DataFrame, n_genes: int | None) -> pd.DataFrame:
    if n_genes is None or n_genes <= 0 or expression.shape[0] <= n_genes:
        return expression
    variances = expression.var(axis=1)
    genes = variances.sort_values(ascending=False).head(n_genes).index
    return expression.loc[genes]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--visium-dir", default=None, help="Directory containing filtered_feature_bc_matrix.h5 and spatial/")
    parser.add_argument("--expression-h5", default=None)
    parser.add_argument("--expression-csv", default=None)
    parser.add_argument("--expression-orientation", default="genes_by_spots", choices=["genes_by_spots", "spots_by_genes"])
    parser.add_argument("--spatial-dir", default=None)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--source", default="unknown")
    parser.add_argument("--platform", default="10x Visium")
    parser.add_argument("--species", default="unknown")
    parser.add_argument("--tissue", default="unknown")
    parser.add_argument("--top-variable-genes", type=int, default=3000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    visium_dir = Path(args.visium_dir) if args.visium_dir else None
    expression_h5 = Path(args.expression_h5) if args.expression_h5 else None
    expression_csv = Path(args.expression_csv) if args.expression_csv else None
    spatial_dir = Path(args.spatial_dir) if args.spatial_dir else None

    if visium_dir is not None:
        expression_h5 = expression_h5 or (visium_dir / "filtered_feature_bc_matrix.h5")
        spatial_dir = spatial_dir or (visium_dir / "spatial")
    if expression_h5 is None and expression_csv is None:
        raise ValueError("Provide --visium-dir, --expression-h5, or --expression-csv")
    if spatial_dir is None:
        raise ValueError("Provide --visium-dir or --spatial-dir")

    if expression_csv is not None:
        expression = read_expression_csv(expression_csv, args.expression_orientation)
    else:
        expression = read_visium_h5(expression_h5)
    coords = read_visium_coordinates(spatial_dir)

    shared_spots = [spot for spot in coords["spot_id"].astype(str) if spot in expression.columns]
    if not shared_spots:
        raise ValueError("No shared spots between expression matrix and spatial coordinates")
    expression = expression.loc[:, shared_spots]
    expression = select_top_variable_genes(expression, args.top_variable_genes)
    coords = coords.set_index("spot_id").loc[shared_spots].reset_index()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    expression.to_csv(output_dir / "expression_matrix.csv")
    coords.to_csv(output_dir / "coordinates.csv", index=False)

    metadata = pd.DataFrame(
        {
            "spot_id": shared_spots,
            "dataset": args.dataset_id,
            "source": args.source,
            "platform": args.platform,
            "species": args.species,
            "tissue": args.tissue,
            "apa_ready": False,
            "label_status": "missing_biological_labels",
        }
    )
    metadata.to_csv(output_dir / "metadata.csv", index=False)

    qc = {
        "dataset": args.dataset_id,
        "source": args.source,
        "platform": args.platform,
        "species": args.species,
        "tissue": args.tissue,
        "apa_source": "expression_only_candidate",
        "apa_ready": False,
        "expression_ready": True,
        "biological_labels_ready": False,
        "site_level_ready": False,
        "n_spots": int(expression.shape[1]),
        "n_expression_genes": int(expression.shape[0]),
        "files": {
            "expression_matrix": "expression_matrix.csv",
            "coordinates": "coordinates.csv",
            "metadata": "metadata.csv",
        },
        "notes": [
            "No apa_matrix.csv was generated.",
            "Do not count this dataset as a true APA benchmark until APA/PAS calls are added.",
        ],
    }
    (output_dir / "qc_summary.json").write_text(json.dumps(qc, indent=2))

    print(f"Prepared expression-only candidate dataset: {output_dir}")
    print(f"Expression genes: {expression.shape[0]}, spots: {expression.shape[1]}")
    print("APA ready: false")


if __name__ == "__main__":
    main()
