"""Integration test for the pipeline-level high-resolution formal suite."""

from __future__ import annotations

import json
import subprocess
import sys

import numpy as np
import pandas as pd


def test_pipeline_highres_suite_runner_toy(tmp_path):
    dataset_dir = tmp_path / "toy_prepared"
    output_dir = tmp_path / "pipeline_highres_suite"
    dataset_dir.mkdir()

    rng = np.random.default_rng(17)
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
        "scripts/run_pipeline_highres_suite.py",
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
        "--max-parent-spots",
        "12",
        "--seeds",
        "11,12",
        "--formal-subbins",
        "2",
        "--scaling-subbins",
        "2",
        "--dropout-rates",
        "0.10,0.20",
        "--capture-rate",
        "0.6",
        "--dropout-rate",
        "0.1",
        "--sparse-n-inducing",
        "8",
        "--bioml-n-neighbors",
        "4",
        "--bioml-rank",
        "3",
        "--bioml-max-iter",
        "4",
        "--skip-graph-sweep",
        "--skip-blend-sweep",
    ]
    result = subprocess.run(cmd, cwd=".", capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr + result.stdout

    long = pd.read_csv(output_dir / "pipeline_highres_suite_results_long.csv")
    assert set(long["run_type"]) == {"formal", "scaling", "dropout"}
    assert set(long["method"]) == {
        "runner_highres_bioml",
        "pipeline_highres_accuracy",
        "pipeline_highres_fast",
    }
    assert set(long.loc[long["run_type"] == "formal", "seed"]) == {11, 12}
    assert np.isfinite(long["rmse_holdout"]).all()
    assert np.isfinite(long["layer_ari"]).all()

    decision = json.loads((output_dir / "decision_summary.json").read_text())
    assert "pipeline_highres_accuracy" in decision["comparisons"]
    assert "pipeline_highres_fast" in decision["comparisons"]

    assert (output_dir / "pipeline_highres_suite_overall_summary.csv").exists()
    assert (output_dir / "figures" / "pipeline_highres_suite_formal.png").exists()
    assert (output_dir / "figures" / "pipeline_highres_suite_scaling.png").exists()
    assert (output_dir / "figures" / "pipeline_highres_suite_dropout.png").exists()
