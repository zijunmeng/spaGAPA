"""Real-data benchmark registry and readiness checks.

This module keeps BIB-oriented real-data benchmarks from turning into a pile of
manual shell commands. It validates prepared dataset directories, summarizes
available evidence, and decides which benchmark suites can run.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable

import pandas as pd


REQUIRED_FILES = ("apa_matrix.csv", "coordinates.csv")
OPTIONAL_FILES = (
    "metadata.csv",
    "expression_matrix.csv",
    "stapaminer_rud_imputed.csv",
    "stapaminer_rud_raw.csv",
    "apa_sites.csv",
    "apa_sites.csv.gz",
    "apa_site_counts.csv",
    "apa_site_counts.csv.gz",
    "qc_summary.json",
)


@dataclass
class PreparedDatasetStatus:
    """Readiness status for one prepared real spatial APA dataset."""

    dataset: str
    data_dir: str
    exists: bool
    required_ready: bool
    external_validation_ready: bool
    highres_validation_ready: bool
    metaapa_readiness: str
    n_genes: int | None = None
    n_spots: int | None = None
    observed_fraction: float | None = None
    has_metadata: bool = False
    has_layer_labels: bool = False
    has_expression: bool = False
    has_stapaminer_original: bool = False
    has_site_table: bool = False
    has_site_counts: bool = False
    missing_required: list[str] = field(default_factory=list)
    missing_recommended: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _read_csv_shape(path: Path) -> tuple[int | None, int | None, float | None, str | None]:
    """Read a prepared matrix shape and observed fraction."""
    try:
        df = pd.read_csv(path, index_col=0)
    except Exception as exc:  # pragma: no cover - defensive metadata path
        return None, None, None, str(exc)
    values = df.values
    return int(df.shape[0]), int(df.shape[1]), float(pd.notna(values).mean()), None


def _read_metadata_columns(path: Path) -> set[str]:
    if not path.exists():
        return set()
    try:
        return set(pd.read_csv(path, nrows=5).columns.astype(str))
    except Exception:
        return set()


def check_prepared_dataset(
    data_dir: str | Path,
    *,
    dataset_name: str | None = None,
    layer_column: str = "layer",
) -> PreparedDatasetStatus:
    """Check whether a prepared dataset can enter BIB benchmark suites."""
    data_dir = Path(data_dir)
    dataset = dataset_name or data_dir.name
    exists = data_dir.exists() and data_dir.is_dir()
    missing_required = [name for name in REQUIRED_FILES if not (data_dir / name).exists()]

    has_metadata = (data_dir / "metadata.csv").exists()
    metadata_cols = _read_metadata_columns(data_dir / "metadata.csv")
    has_layer_labels = has_metadata and layer_column in metadata_cols
    has_expression = (data_dir / "expression_matrix.csv").exists()
    has_stapaminer_original = (data_dir / "stapaminer_rud_imputed.csv").exists()
    has_site_table = (data_dir / "apa_sites.csv").exists() or (data_dir / "apa_sites.csv.gz").exists()
    has_site_counts = (data_dir / "apa_site_counts.csv").exists() or (data_dir / "apa_site_counts.csv.gz").exists()

    n_genes = None
    n_spots = None
    observed_fraction = None
    notes: list[str] = []
    if exists and not missing_required:
        n_genes, n_spots, observed_fraction, error = _read_csv_shape(data_dir / "apa_matrix.csv")
        if error is not None:
            notes.append(f"Failed to read apa_matrix.csv: {error}")

    missing_recommended = []
    if not has_metadata:
        missing_recommended.append("metadata.csv")
    if not has_layer_labels:
        missing_recommended.append(f"metadata.{layer_column}")
    if not has_expression:
        missing_recommended.append("expression_matrix.csv")
    if not has_stapaminer_original:
        missing_recommended.append("stapaminer_rud_imputed.csv")
    if not has_site_table:
        missing_recommended.append("apa_sites.csv(.gz)")
    if not has_site_counts:
        missing_recommended.append("apa_site_counts.csv(.gz)")

    required_ready = exists and not missing_required and n_genes is not None and n_spots is not None
    external_validation_ready = required_ready and has_metadata and has_layer_labels
    highres_validation_ready = external_validation_ready and has_expression

    if has_site_table and has_site_counts:
        metaapa_readiness = "partial_site_level"
    elif has_site_table or has_site_counts:
        metaapa_readiness = "incomplete_site_level"
    else:
        metaapa_readiness = "missing_site_level"

    if not exists:
        notes.append("Dataset directory does not exist.")
    if missing_required:
        notes.append("Missing required files: " + ", ".join(missing_required))
    if external_validation_ready and not has_expression:
        notes.append("External layer validation can run, but expression-aware baselines will be skipped.")
    if metaapa_readiness != "partial_site_level":
        notes.append("metaAPA-style upstream site-level validation is not yet complete.")

    return PreparedDatasetStatus(
        dataset=dataset,
        data_dir=str(data_dir),
        exists=exists,
        required_ready=required_ready,
        external_validation_ready=external_validation_ready,
        highres_validation_ready=highres_validation_ready,
        metaapa_readiness=metaapa_readiness,
        n_genes=n_genes,
        n_spots=n_spots,
        observed_fraction=observed_fraction,
        has_metadata=has_metadata,
        has_layer_labels=has_layer_labels,
        has_expression=has_expression,
        has_stapaminer_original=has_stapaminer_original,
        has_site_table=has_site_table,
        has_site_counts=has_site_counts,
        missing_required=missing_required,
        missing_recommended=missing_recommended,
        notes=notes,
    )


def discover_prepared_datasets(
    processed_root: str | Path,
    *,
    layer_column: str = "layer",
) -> list[PreparedDatasetStatus]:
    """Discover and check all prepared dataset directories under a root."""
    processed_root = Path(processed_root)
    if not processed_root.exists():
        return []
    statuses = []
    for data_dir in sorted(path for path in processed_root.iterdir() if path.is_dir()):
        statuses.append(check_prepared_dataset(data_dir, layer_column=layer_column))
    return statuses


def statuses_to_dataframe(statuses: Iterable[PreparedDatasetStatus]) -> pd.DataFrame:
    """Convert readiness statuses to a tabular report."""
    rows = [status.to_dict() for status in statuses]
    if not rows:
        return pd.DataFrame(
            columns=[
                "dataset",
                "data_dir",
                "required_ready",
                "external_validation_ready",
                "highres_validation_ready",
                "metaapa_readiness",
            ]
        )
    return pd.DataFrame(rows)


def summarize_bib_readiness(statuses: Iterable[PreparedDatasetStatus]) -> dict[str, Any]:
    """Summarize whether the current real-data collection is BIB-ready."""
    statuses = list(statuses)
    external_ready = [status for status in statuses if status.external_validation_ready]
    highres_ready = [status for status in statuses if status.highres_validation_ready]
    site_ready = [status for status in statuses if status.metaapa_readiness == "partial_site_level"]
    bib_ready = len(external_ready) >= 3 and len(highres_ready) >= 1 and len(site_ready) >= 2
    return {
        "n_datasets": len(statuses),
        "n_external_validation_ready": len(external_ready),
        "n_highres_validation_ready": len(highres_ready),
        "n_site_level_partial_ready": len(site_ready),
        "bib_ready": bool(bib_ready),
        "next_required": (
            "Ready for BIB-scale real-data benchmark."
            if bib_ready
            else "Add and standardize more external real datasets before claiming BIB-scale real-data evidence."
        ),
    }
