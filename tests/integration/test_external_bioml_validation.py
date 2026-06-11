"""Integration tests for external BioML validation runner."""

from __future__ import annotations

import subprocess
import sys

import numpy as np
import pandas as pd


def test_external_bioml_validation_runner_toy(tmp_path):
    dataset_dir = tmp_path / "toy_dataset"
    output_dir = tmp_path / "external_validation"
    dataset_dir.mkdir()

    rng = np.random.default_rng(42)
    grid = 5
    x = np.repeat(np.arange(grid), grid).astype(float)
    y = np.tile(np.arange(grid), grid).astype(float)
    spot_names = [f"Spot_{i}" for i in range(grid * grid)]
    layers = np.where(x < 2, "inner", np.where(x > 2, "outer", "middle"))

    values = []
    for gene_idx in range(6):
        base = (x / max(x.max(), 1.0)) if gene_idx % 2 == 0 else (y / max(y.max(), 1.0))
        noise = rng.normal(scale=0.03, size=len(spot_names))
        values.append(np.clip(base + noise, 0.0, 1.0))
    apa = np.vstack(values)
    apa[:, ::6] = np.nan

    expression = np.vstack(
        [
            x,
            y,
            x + y,
            (x > 2).astype(float),
            rng.normal(size=len(spot_names)),
            np.ones(len(spot_names)),
        ]
    )

    pd.DataFrame(
        apa,
        index=[f"Gene_{i}" for i in range(apa.shape[0])],
        columns=spot_names,
    ).to_csv(dataset_dir / "apa_matrix.csv")
    pd.DataFrame(
        expression,
        index=[f"Expr_{i}" for i in range(expression.shape[0])],
        columns=spot_names,
    ).to_csv(dataset_dir / "expression_matrix.csv")
    pd.DataFrame({"spot_id": spot_names, "x": x, "y": y}).to_csv(
        dataset_dir / "coordinates.csv",
        index=False,
    )
    pd.DataFrame(
        {
            "spot_id": spot_names,
            "layer": layers,
            "dataset": "toy_external",
            "source": "pytest",
        }
    ).to_csv(dataset_dir / "metadata.csv", index=False)

    cmd = [
        sys.executable,
        "scripts/run_external_bioml_validation.py",
        "--dataset-dir",
        str(dataset_dir),
        "--output-dir",
        str(output_dir),
        "--dataset-name",
        "toy_external",
        "--n-genes",
        "4",
        "--min-observed-spots",
        "5",
        "--methods",
        "raw,stapaminer_knn_expression",
        "--skip-downstream",
    ]
    result = subprocess.run(cmd, cwd=".", capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr + result.stdout

    summary = pd.read_csv(output_dir / "external_validation_summary.csv")
    assert set(summary["method"]) == {"raw", "stapaminer_knn_expression"}
    assert (output_dir / "figures" / "layer_separation_metrics.png").exists()
    assert (output_dir / "figures" / "ref_package_experiment_coverage.png").exists()
    assert (output_dir / "toy_external" / "site_level_validation_checklist.json").exists()
