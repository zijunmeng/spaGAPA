#!/usr/bin/env python
"""Profile spaGAPA GP imputation runtime on a prepared dataset.

This script separates fit and impute/predict time, and compares exact GP
parallel settings against sparse GP. It is intended for engineering decisions,
not publication figures.
"""

from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

os.makedirs("/tmp/matplotlib-spagapa", exist_ok=True)
os.environ["MPLCONFIGDIR"] = "/tmp/matplotlib-spagapa"

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import psutil

from spagapa.imputation import GPImputer, SparseGPImputer


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--n-genes", type=int, default=40)
    parser.add_argument("--min-observed-spots", type=int, default=120)
    parser.add_argument("--mask-fraction", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--kernel", default="matern", choices=["matern", "rbf", "auto"])
    parser.add_argument("--alpha", type=float, default=1e-10)
    parser.add_argument("--n-restarts", type=int, default=1)
    parser.add_argument("--n-jobs", default="1,2,4,8")
    parser.add_argument("--profile-sparse", action="store_true")
    parser.add_argument("--sparse-n-inducing", type=int, default=100)
    parser.add_argument("--sparse-length-scale", type=float, default=1.0)
    parser.add_argument("--sparse-noise-level", type=float, default=0.1)
    return parser.parse_args()


def read_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    return df


def load_inputs(data_dir: Path) -> tuple[pd.DataFrame, np.ndarray]:
    apa = read_table(data_dir / "apa_matrix.csv")
    coords = pd.read_csv(data_dir / "coordinates.csv")
    coords = coords.set_index("spot_id" if "spot_id" in coords.columns else coords.columns[0])
    coords.index = coords.index.astype(str)
    common = [spot for spot in apa.columns if spot in coords.index]
    apa = apa.loc[:, common]
    coords_arr = coords.loc[common][["x", "y"]].values.astype(float)
    return apa, coords_arr


def select_genes(apa: pd.DataFrame, n_genes: int, min_observed_spots: int) -> pd.DataFrame:
    observed = apa.notna().sum(axis=1)
    eligible = observed[observed >= min_observed_spots]
    if eligible.empty:
        eligible = observed[observed >= max(3, min_observed_spots // 2)]
    return apa.loc[eligible.sort_values(ascending=False).head(n_genes).index]


def random_holdout(values: np.ndarray, mask_fraction: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    observed = np.argwhere(np.isfinite(values))
    n_mask = max(1, int(len(observed) * mask_fraction))
    chosen = observed[rng.choice(len(observed), size=n_mask, replace=False)]
    holdout = np.zeros(values.shape, dtype=bool)
    holdout[chosen[:, 0], chosen[:, 1]] = True
    return holdout


def compute_rmse(true_values: np.ndarray, pred: np.ndarray, holdout: np.ndarray) -> float:
    diff = true_values[holdout] - pred[holdout]
    return float(np.sqrt(np.mean(diff * diff)))


def profile_exact(
    coords: np.ndarray,
    train: np.ndarray,
    mask: np.ndarray,
    true_values: np.ndarray,
    holdout: np.ndarray,
    kernel: str,
    alpha: float,
    n_restarts: int,
    n_jobs: int,
) -> dict:
    process = psutil.Process(os.getpid())
    start_mem = process.memory_info().rss / 1024 / 1024
    imputer = GPImputer(
        kernel_type=kernel,
        alpha=alpha,
        n_restarts_optimizer=n_restarts,
    )
    t0 = time.perf_counter()
    batch = imputer.fit_batch(
        coords,
        np.nan_to_num(train, nan=0.0),
        mask=mask,
        n_jobs=n_jobs,
        verbose=False,
    )
    t_fit = time.perf_counter() - t0
    t1 = time.perf_counter()
    pred, unc = batch.impute(return_uncertainty=True)
    t_impute = time.perf_counter() - t1
    peak_mem = process.memory_info().rss / 1024 / 1024
    return {
        "method": "exact_gp",
        "n_jobs": int(n_jobs),
        "n_restarts": int(n_restarts),
        "fit_s": float(t_fit),
        "impute_s": float(t_impute),
        "total_s": float(t_fit + t_impute),
        "rss_delta_mb": float(max(0.0, peak_mem - start_mem)),
        "rmse_holdout": compute_rmse(true_values, np.clip(pred, 0.0, 1.0), holdout),
        "uncertainty_mean_holdout": float(np.nanmean(unc[holdout])) if unc is not None else np.nan,
    }


def profile_sparse(
    coords: np.ndarray,
    train: np.ndarray,
    mask: np.ndarray,
    true_values: np.ndarray,
    holdout: np.ndarray,
    n_inducing: int,
    length_scale: float,
    noise_level: float,
) -> dict:
    process = psutil.Process(os.getpid())
    start_mem = process.memory_info().rss / 1024 / 1024
    imputer = SparseGPImputer(
        n_inducing=n_inducing,
        inducing_method="kmeans",
        length_scale=length_scale,
        noise_level=noise_level,
    )
    t0 = time.perf_counter()
    batch = imputer.fit_batch(
        coords,
        np.nan_to_num(train, nan=0.0),
        mask=mask,
        n_jobs=1,
        verbose=False,
    )
    t_fit = time.perf_counter() - t0
    t1 = time.perf_counter()
    pred, unc = batch.impute(return_uncertainty=True)
    t_impute = time.perf_counter() - t1
    peak_mem = process.memory_info().rss / 1024 / 1024
    return {
        "method": "sparse_gp",
        "n_jobs": 1,
        "n_restarts": 0,
        "fit_s": float(t_fit),
        "impute_s": float(t_impute),
        "total_s": float(t_fit + t_impute),
        "rss_delta_mb": float(max(0.0, peak_mem - start_mem)),
        "rmse_holdout": compute_rmse(true_values, np.clip(pred, 0.0, 1.0), holdout),
        "uncertainty_mean_holdout": float(np.nanmean(unc[holdout])) if unc is not None else np.nan,
        "n_inducing": int(n_inducing),
        "length_scale": float(length_scale),
        "noise_level": float(noise_level),
    }


def plot(rows: pd.DataFrame, output_dir: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
    labels = rows["method"] + "_j" + rows["n_jobs"].astype(str)
    x = np.arange(len(rows))
    for ax, metric, title in zip(
        axes,
        ["total_s", "rmse_holdout", "rss_delta_mb"],
        ["Total Runtime", "Holdout RMSE", "RSS Delta MB"],
    ):
        ax.bar(x, rows[metric].values, color="#4c78a8")
        ax.set_title(title)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=35, ha="right")
        ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_dir / "gp_runtime_profile.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    data_dir = Path(args.data_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    apa, coords = load_inputs(data_dir)
    selected = select_genes(apa, args.n_genes, args.min_observed_spots)
    values = selected.values.astype(float)
    holdout = random_holdout(values, args.mask_fraction, args.seed)
    train = values.copy()
    train[holdout] = np.nan
    mask = np.isfinite(train)

    rows = []
    for n_jobs in [int(x) for x in args.n_jobs.split(",") if x.strip()]:
        rows.append(
            profile_exact(
                coords,
                train,
                mask,
                values,
                holdout,
                args.kernel,
                args.alpha,
                args.n_restarts,
                n_jobs,
            )
        )

    if args.profile_sparse:
        rows.append(
            profile_sparse(
                coords,
                train,
                mask,
                values,
                holdout,
                args.sparse_n_inducing,
                args.sparse_length_scale,
                args.sparse_noise_level,
            )
        )

    df = pd.DataFrame(rows)
    df.to_csv(output_dir / "gp_runtime_profile.csv", index=False)
    (output_dir / "gp_runtime_profile.json").write_text(
        json.dumps(
            {
                "data_dir": str(data_dir),
                "n_genes": int(values.shape[0]),
                "n_spots": int(values.shape[1]),
                "mask_fraction": float(args.mask_fraction),
                "seed": int(args.seed),
                "results": df.to_dict(orient="records"),
            },
            indent=2,
        )
    )
    plot(df, output_dir)
    print(df.to_string(index=False))
    print(f"Saved profile to: {output_dir}")


if __name__ == "__main__":
    main()
