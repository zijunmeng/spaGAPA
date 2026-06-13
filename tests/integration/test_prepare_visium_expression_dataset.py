"""Integration test for expression-only Visium candidate preparation."""

from __future__ import annotations

import json
import subprocess
import sys

import numpy as np
import pandas as pd


def test_prepare_visium_expression_dataset_from_csv(tmp_path):
    spatial_dir = tmp_path / "spatial"
    output_dir = tmp_path / "processed_brain_candidate"
    spatial_dir.mkdir()

    spots = [f"Spot_{i}" for i in range(5)]
    genes = [f"Gene_{i}" for i in range(4)]
    expr = pd.DataFrame(
        np.arange(len(genes) * len(spots)).reshape(len(genes), len(spots)),
        index=genes,
        columns=spots,
    )
    expr_path = tmp_path / "expression.csv"
    expr.to_csv(expr_path)

    positions = pd.DataFrame(
        {
            "spot_id": spots,
            "in_tissue": [1, 1, 1, 0, 1],
            "array_row": [0, 0, 1, 1, 2],
            "array_col": [0, 1, 0, 1, 0],
            "pxl_row": [10, 20, 30, 40, 50],
            "pxl_col": [11, 21, 31, 41, 51],
        }
    )
    positions.to_csv(spatial_dir / "tissue_positions.csv", index=False)

    cmd = [
        sys.executable,
        "scripts/prepare_visium_expression_dataset.py",
        "--expression-csv",
        str(expr_path),
        "--spatial-dir",
        str(spatial_dir),
        "--output-dir",
        str(output_dir),
        "--dataset-id",
        "toy_brain_candidate",
        "--source",
        "pytest",
        "--species",
        "human",
        "--tissue",
        "brain",
        "--top-variable-genes",
        "3",
    ]
    result = subprocess.run(cmd, cwd=".", text=True, capture_output=True, check=False)
    assert result.returncode == 0, result.stderr + result.stdout

    expression = pd.read_csv(output_dir / "expression_matrix.csv", index_col=0)
    coords = pd.read_csv(output_dir / "coordinates.csv")
    metadata = pd.read_csv(output_dir / "metadata.csv")
    qc = json.loads((output_dir / "qc_summary.json").read_text())

    assert expression.shape == (3, 4)
    assert coords.shape[0] == 4
    assert set(coords["spot_id"]) == {"Spot_0", "Spot_1", "Spot_2", "Spot_4"}
    assert bool(metadata["apa_ready"].iloc[0]) is False
    assert qc["apa_ready"] is False
    assert qc["apa_source"] == "expression_only_candidate"
    assert not (output_dir / "apa_matrix.csv").exists()
