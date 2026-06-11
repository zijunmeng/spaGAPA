#!/usr/bin/env python
"""Run a small spaGAPA smoke test on the prepared stAPAminer MOB dataset."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from spagapa import SpaGAPA
from spagapa.io import load_spatial_dataset


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        default="spaGAPA/data/processed/stapaminer_mob",
        help="Prepared stAPAminer MOB directory.",
    )
    parser.add_argument(
        "--output-dir",
        default="spaGAPA/benchmark_results/real/stapaminer_mob_smoke",
        help="Output directory for smoke-test results.",
    )
    parser.add_argument(
        "--n-genes",
        type=int,
        default=30,
        help="Number of sufficiently observed genes to test.",
    )
    parser.add_argument(
        "--min-observed-spots",
        type=int,
        default=80,
        help="Minimum non-NA observed spots per selected gene.",
    )
    parser.add_argument(
        "--sparse",
        action="store_true",
        help="Use sparse GP instead of dense GP.",
    )
    parser.add_argument(
        "--n-domains",
        type=int,
        default=3,
        help="Number of spatial domains for the smoke-test pipeline.",
    )
    return parser.parse_args()


def select_gene_subset(
    matrix_file: Path,
    n_genes: int,
    min_observed_spots: int,
) -> pd.DataFrame:
    matrix = pd.read_csv(matrix_file, index_col=0)
    observed = matrix.notna().sum(axis=1)
    has_missing = observed < matrix.shape[1]
    selected = observed[(observed >= min_observed_spots) & has_missing]
    selected = selected.sort_values(ascending=False)
    if selected.empty:
        selected = observed[observed >= min_observed_spots].sort_values(ascending=False)
    if selected.empty:
        raise ValueError(
            f"No genes have at least {min_observed_spots} observed spots"
        )
    return matrix.loc[selected.head(n_genes).index]


def main() -> None:
    args = parse_args()
    data_dir = Path(args.data_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    subset = select_gene_subset(
        data_dir / "apa_matrix.csv",
        n_genes=args.n_genes,
        min_observed_spots=args.min_observed_spots,
    )
    subset_file = output_dir / "apa_matrix_subset.csv"
    subset.to_csv(subset_file)

    dataset = load_spatial_dataset(
        apa_matrix=subset_file,
        coordinates=data_dir / "coordinates.csv",
        matrix_orientation="genes_by_spots",
    )

    pipeline = SpaGAPA(
        input_type="apa_index",
        kernel_type="matern",
        use_sparse_gp=args.sparse,
        n_inducing=40,
        n_neighbors=6,
        min_spots=10,
        verbose=True,
    )
    results = pipeline.run(
        dataset=dataset,
        impute=True,
        quantify=True,
        identify_domains=True,
        differential_analysis=True,
        detect_svapa=True,
        n_domains=args.n_domains,
        fdr_threshold=0.1,
    )

    imputed = results["imputed_values"]
    uncertainty = results["uncertainty"]
    svapa = results["svapa_genes"]
    domains = results["domains"]["labels"]

    np.save(output_dir / "spagapa_imputed.npy", imputed)
    np.save(output_dir / "spagapa_uncertainty.npy", uncertainty)
    svapa.to_csv(output_dir / "spagapa_svapa.csv", index=False)
    pd.DataFrame({
        "spot_id": dataset.spot_names,
        "domain": domains,
    }).to_csv(output_dir / "spagapa_domains.csv", index=False)

    differential_summary = {}
    for name, df in results["differential"].items():
        out = output_dir / f"spagapa_differential_{name}.csv"
        df.to_csv(out, index=False)
        differential_summary[name] = {
            "n_genes": int(len(df)),
            "min_padj": float(df["padj"].min()) if len(df) else None,
        }

    summary = {
        "dataset": "stAPAminer_MOB",
        "n_genes": int(dataset.n_genes),
        "n_spots": int(dataset.n_spots),
        "raw_missing_rate": float(np.isnan(dataset.raw_counts).mean()),
        "imputed_shape": list(imputed.shape),
        "uncertainty_shape": list(uncertainty.shape),
        "mean_uncertainty": float(np.nanmean(uncertainty)),
        "n_domains": int(len(np.unique(domains))),
        "n_svapa_significant": int(svapa["significant"].sum()),
        "use_sparse_gp": bool(args.sparse),
        "differential": differential_summary,
    }
    (output_dir / "smoke_summary.json").write_text(
        json.dumps(summary, indent=2)
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
