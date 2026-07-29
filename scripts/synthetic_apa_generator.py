#!/usr/bin/env python3
"""Synthetic spatial APA data generator for the runtime/scalability benchmark.

Generates, for a given number of spots N, a realistic spatial-transcriptomics
APA imputation problem with three artifacts that all four methods consume:

  1. ``coordinates.csv``  - spot x,y on a near-square grid (column ``x``,``y``,
     row index = spot id ``S0..S{N-1}``)
  2. ``apa_index.csv``     - gene x spot APA index matrix, ~80% NaN, with a
     subset of spatially-structured genes (radial gradient) so spatial methods
     have real signal to recover.
  3. ``expression.csv``    - gene x spot integer COUNT matrix with lognormal
     per-gene means and 4-quadrant domain structure (needed by stAPAminer's
     expression-KNN and spvAPA's WNN / SCTransform pipeline).

Determinism: a fixed ``seed`` plus the spot count N makes the dataset at a
given scale identical across runs and across methods (fairness: everyone sees
the same masked matrix).

Design notes
------------
* The grid is laid out as ``side = ceil(sqrt(N))`` and the first N cells are
  used, so coordinates are integer lattice points (Visium-like spacing = 1).
* ``n_genes=2000`` matches the real Visium head-to-head scale (1801 genes on
  GSE183456) and is large enough that spvAPA's SCTransform+PCA+WNN pipeline is
  well-conditioned (it errors below ~1000 informative genes).
* APA sparsity target ~0.20 observed (i.e. ~80% NaN), matching the real
  gene-level distal-usage index observed fraction (~0.1-0.3 in processed data).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def generate_synthetic_apa(
    n_spots: int,
    n_genes: int = 2000,
    observed_fraction: float = 0.20,
    n_structured_frac: float = 0.10,
    seed: int = 42,
) -> dict[str, np.ndarray]:
    """Generate a synthetic spatial APA dataset.

    Returns a dict with keys: ``coordinates`` (n_spots,2),
    ``apa_index`` (n_genes,n_spots) float with NaN for missing,
    ``expression`` (n_genes,n_spots) int >= 0.
    """
    rng = np.random.default_rng(seed)

    # ---- near-square integer grid ----
    side = int(np.ceil(np.sqrt(n_spots)))
    spots = np.arange(n_spots)
    rows, cols = np.unravel_index(spots, (side, side))
    coordinates = np.column_stack([cols, rows]).astype(float)

    # 4-quadrant domain labels (drives expression block structure)
    domain = (cols < side // 2).astype(int) + (rows < side // 2).astype(int)

    # ---- expression: lognormal gene means x domain effects ----
    gene_means = np.exp(rng.normal(2.5, 1.0, n_genes))  # ~mean count ~12-25
    domain_effect = rng.normal(0.0, 0.6, (n_genes, 4))
    lam = gene_means[:, None] * np.exp(domain_effect[:, domain])
    expression = np.random.poisson(lam).astype(np.int64)
    # guarantee no all-zero genes (SCTransform / PCA dislike them)
    zero_genes = np.where(expression.sum(axis=1) == 0)[0]
    if len(zero_genes):
        expression[zero_genes] = np.maximum(
            1, np.random.poisson(gene_means[zero_genes], size=len(zero_genes))
        ).astype(np.int64)

    # ---- APA index: ~80% NaN, structured genes have radial gradient ----
    apa_index = np.full((n_genes, n_spots), np.nan)
    observed = rng.random((n_genes, n_spots)) < observed_fraction
    n_structured = max(1, int(n_genes * n_structured_frac))
    structured = set(rng.choice(n_genes, n_structured, replace=False))
    base_value = rng.random(n_genes)
    cdist2 = (cols - side / 2.0) ** 2 + (rows - side / 2.0) ** 2
    cdist2 = cdist2 / (cdist2.max() + 1e-9)
    for gi in range(n_genes):
        if gi in structured:
            # radial gradient: distal-usage rises from center outward
            truth = np.clip(0.2 + 0.6 * cdist2 + rng.normal(0, 0.05, n_spots), 0, 1)
        else:
            truth = np.clip(
                base_value[gi] + rng.normal(0, 0.1, n_spots), 0, 1
            )
        apa_index[gi, observed[gi]] = truth[observed[gi]]

    return {
        "coordinates": coordinates,
        "apa_index": apa_index,
        "expression": expression,
    }


def write_dataset(out_dir: Path, data: dict[str, np.ndarray]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    coords = pd.DataFrame(
        {"spot_id": [f"S{i}" for i in range(len(data["coordinates"]))],
         "x": data["coordinates"][:, 0],
         "y": data["coordinates"][:, 1]},
    )
    spots = coords["spot_id"].tolist()
    genes = [f"G{i}" for i in range(data["apa_index"].shape[0])]
    coords.to_csv(out_dir / "coordinates.csv", index=False)
    pd.DataFrame(data["apa_index"], index=genes, columns=spots).to_csv(
        out_dir / "apa_index.csv"
    )
    pd.DataFrame(data["expression"], index=genes, columns=spots).to_csv(
        out_dir / "expression.csv"
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-spots", type=int, required=True)
    ap.add_argument("--n-genes", type=int, default=2000)
    ap.add_argument("--observed-fraction", type=float, default=0.20)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    data = generate_synthetic_apa(
        n_spots=args.n_spots,
        n_genes=args.n_genes,
        observed_fraction=args.observed_fraction,
        seed=args.seed,
    )
    write_dataset(Path(args.out_dir), data)
    obs_frac = np.isfinite(data["apa_index"]).mean()
    print(
        f"generated {args.n_spots} spots x {args.n_genes} genes | "
        f"observed {obs_frac:.3f} | mean expr {data['expression'].mean():.1f} "
        f"-> {args.out_dir}"
    )


if __name__ == "__main__":
    main()
