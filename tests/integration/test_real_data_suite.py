"""Integration test for real-data suite readiness orchestrator."""

from __future__ import annotations

import json
import subprocess
import sys

import numpy as np
import pandas as pd


def test_real_data_suite_readiness_toy(tmp_path):
    processed_root = tmp_path / "processed"
    dataset_dir = processed_root / "toy_real"
    output_dir = tmp_path / "real_data_suite"
    dataset_dir.mkdir(parents=True)

    genes = [f"Gene_{i}" for i in range(4)]
    spots = [f"Spot_{i}" for i in range(6)]
    x = np.arange(len(spots), dtype=float)
    y = x % 2
    apa = np.vstack(
        [
            np.linspace(0.1, 0.8, len(spots)),
            np.linspace(0.8, 0.1, len(spots)),
            np.full(len(spots), 0.5),
            x / max(x.max(), 1.0),
        ]
    )
    apa[:, ::5] = np.nan

    pd.DataFrame(apa, index=genes, columns=spots).to_csv(dataset_dir / "apa_matrix.csv")
    pd.DataFrame({"spot_id": spots, "x": x, "y": y}).to_csv(
        dataset_dir / "coordinates.csv",
        index=False,
    )
    pd.DataFrame(
        {
            "spot_id": spots,
            "layer": ["L1", "L1", "L2", "L2", "L3", "L3"],
            "dataset": "toy_real",
        }
    ).to_csv(dataset_dir / "metadata.csv", index=False)
    pd.DataFrame(apa, index=genes, columns=spots).to_csv(dataset_dir / "expression_matrix.csv")

    cmd = [
        sys.executable,
        "scripts/run_real_data_suite.py",
        "--processed-root",
        str(processed_root),
        "--output-dir",
        str(output_dir),
    ]
    result = subprocess.run(cmd, cwd=".", text=True, capture_output=True, check=False)
    assert result.returncode == 0, result.stderr + result.stdout

    readiness = pd.read_csv(output_dir / "real_data_readiness.csv")
    assert readiness.loc[0, "dataset"] == "toy_real"
    assert bool(readiness.loc[0, "external_validation_ready"]) is True
    assert bool(readiness.loc[0, "highres_validation_ready"]) is True

    summary = json.loads((output_dir / "real_data_readiness_summary.json").read_text())
    assert summary["n_datasets"] == 1
    assert summary["n_external_validation_ready"] == 1
    assert summary["bib_ready"] is False

    decision = json.loads((output_dir / "decision_summary.json").read_text())
    assert decision["ready_datasets"]["external_validation"] == ["toy_real"]
    assert decision["ready_datasets"]["highres_validation"] == ["toy_real"]
    assert (output_dir / "real_data_suite_run_manifest.csv").exists()
