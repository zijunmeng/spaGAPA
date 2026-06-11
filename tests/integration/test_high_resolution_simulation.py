"""Integration test for high-resolution pseudo-bin benchmark."""

from __future__ import annotations

import subprocess
import sys

import numpy as np
import pandas as pd


def test_high_resolution_simulation_runner_toy(tmp_path):
    dataset_dir = tmp_path / "toy_prepared"
    output_dir = tmp_path / "highres_out"
    dataset_dir.mkdir()

    rng = np.random.default_rng(7)
    grid = 4
    x = np.repeat(np.arange(grid), grid).astype(float)
    y = np.tile(np.arange(grid), grid).astype(float)
    spots = [f"Spot_{i}" for i in range(len(x))]
    layers = np.where(x < 2, "inner", "outer")

    apa = np.vstack(
        [
            np.clip(x / max(x.max(), 1.0) + rng.normal(scale=0.02, size=len(x)), 0, 1),
            np.clip(y / max(y.max(), 1.0) + rng.normal(scale=0.02, size=len(x)), 0, 1),
            np.clip((x + y) / max((x + y).max(), 1.0), 0, 1),
            np.full(len(x), 0.5),
        ]
    )
    apa[:, ::5] = np.nan
    expression = np.vstack([x, y, x + y, np.ones(len(x))])

    pd.DataFrame(
        apa,
        index=[f"Gene_{i}" for i in range(apa.shape[0])],
        columns=spots,
    ).to_csv(dataset_dir / "apa_matrix.csv")
    pd.DataFrame(
        expression,
        index=[f"Expr_{i}" for i in range(expression.shape[0])],
        columns=spots,
    ).to_csv(dataset_dir / "expression_matrix.csv")
    pd.DataFrame({"spot_id": spots, "x": x, "y": y}).to_csv(
        dataset_dir / "coordinates.csv",
        index=False,
    )
    pd.DataFrame(
        {
            "spot_id": spots,
            "layer": layers,
            "dataset": "toy_highres",
            "source": "pytest",
        }
    ).to_csv(dataset_dir / "metadata.csv", index=False)

    cmd = [
        sys.executable,
        "scripts/run_high_resolution_simulation.py",
        "--data-dir",
        str(dataset_dir),
        "--output-dir",
        str(output_dir),
        "--dataset-name",
        "toy_highres",
        "--n-genes",
        "3",
        "--min-observed-spots",
        "4",
        "--subbins-per-spot",
        "2",
        "--methods",
        "raw,expression_knn,highres_bioml",
        "--capture-rate",
        "0.6",
        "--dropout-rate",
        "0.1",
        "--sparse-n-inducing",
        "8",
        "--bioml-n-neighbors",
        "4",
    ]
    result = subprocess.run(cmd, cwd=".", capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr + result.stdout

    summary = pd.read_csv(output_dir / "highres_results_summary.csv")
    assert set(summary["method"]) == {"raw", "expression_knn", "highres_bioml"}
    assert int(summary["n_bins"].iloc[0]) == 32
    highres_row = summary[summary["method"] == "highres_bioml"].iloc[0]
    assert highres_row["domain_source"] == "model_domain"
    assert np.isfinite(highres_row["rmse_holdout"])
    assert (output_dir / "figures" / "highres_method_comparison.png").exists()
    assert (output_dir / "figures" / "highres_scaling.png").exists()
    assert (output_dir / "figures" / "highres_domain_maps_last.png").exists()
