"""Tests for real-data benchmark readiness registry."""

from __future__ import annotations

import numpy as np
import pandas as pd

from spagapa.benchmark import (
    check_prepared_dataset,
    discover_prepared_datasets,
    statuses_to_dataframe,
    summarize_bib_readiness,
)


def write_prepared_dataset(path, *, with_expression=True, with_sites=True):
    path.mkdir(parents=True)
    genes = [f"Gene_{i}" for i in range(3)]
    spots = [f"Spot_{i}" for i in range(5)]
    values = np.linspace(0.1, 0.9, len(genes) * len(spots)).reshape(len(genes), len(spots))
    values[:, ::4] = np.nan

    pd.DataFrame(values, index=genes, columns=spots).to_csv(path / "apa_matrix.csv")
    pd.DataFrame(
        {"spot_id": spots, "x": np.arange(len(spots)), "y": np.arange(len(spots)) % 2}
    ).to_csv(path / "coordinates.csv", index=False)
    pd.DataFrame(
        {
            "spot_id": spots,
            "layer": ["L1", "L1", "L2", "L2", "L2"],
            "dataset": path.name,
        }
    ).to_csv(path / "metadata.csv", index=False)
    if with_expression:
        pd.DataFrame(values, index=genes, columns=spots).to_csv(path / "expression_matrix.csv")
    if with_sites:
        pd.DataFrame({"site": ["PAS1"], "gene": ["Gene_0"]}).to_csv(path / "apa_sites.csv", index=False)
        pd.DataFrame({"site": ["PAS1"], "Spot_0": [3]}).to_csv(path / "apa_site_counts.csv", index=False)


def test_check_prepared_dataset_ready(tmp_path):
    dataset_dir = tmp_path / "dataset_a"
    write_prepared_dataset(dataset_dir)

    status = check_prepared_dataset(dataset_dir)

    assert status.required_ready is True
    assert status.external_validation_ready is True
    assert status.highres_validation_ready is True
    assert status.metaapa_readiness == "partial_site_level"
    assert status.n_genes == 3
    assert status.n_spots == 5
    assert 0.0 < status.observed_fraction < 1.0


def test_discover_and_summarize_bib_readiness(tmp_path):
    write_prepared_dataset(tmp_path / "dataset_a")
    write_prepared_dataset(tmp_path / "dataset_b", with_expression=False)

    statuses = discover_prepared_datasets(tmp_path)
    table = statuses_to_dataframe(statuses)
    summary = summarize_bib_readiness(statuses)

    assert len(statuses) == 2
    assert set(table["dataset"]) == {"dataset_a", "dataset_b"}
    assert summary["n_datasets"] == 2
    assert summary["n_external_validation_ready"] == 2
    assert summary["n_highres_validation_ready"] == 1
    assert summary["bib_ready"] is False


def test_missing_required_dataset_not_ready(tmp_path):
    dataset_dir = tmp_path / "incomplete"
    dataset_dir.mkdir()
    pd.DataFrame({"spot_id": ["Spot_0"], "x": [0.0], "y": [0.0]}).to_csv(
        dataset_dir / "coordinates.csv",
        index=False,
    )

    status = check_prepared_dataset(dataset_dir)

    assert status.required_ready is False
    assert "apa_matrix.csv" in status.missing_required
    assert status.external_validation_ready is False
