"""
Analysis preset resolution for user-facing spaGAPA modes.

The implementation keeps internal algorithm choices explicit while exposing a
small public surface:

``auto`` | ``standard`` | ``highres_accuracy`` | ``highres_fast``.
"""

from dataclasses import asdict, dataclass
from typing import Dict, List

import numpy as np


VALID_ANALYSIS_PRESETS = ("auto", "standard", "highres_accuracy", "highres_fast")


@dataclass
class DatasetProfile:
    """Compact data-shape summary used by the preset resolver."""

    n_genes: int
    n_spots: int
    finite_fraction: float
    observed_fraction: float
    positive_fraction: float
    median_observed_per_spot: float
    median_positive_per_spot: float
    highres_like: bool
    highres_reasons: List[str]

    def to_dict(self) -> Dict:
        return asdict(self)


def profile_spatial_apa_matrix(
    values: np.ndarray,
    input_type: str = "apa_index",
    highres_spot_threshold: int = 1000,
    sparse_observed_fraction_threshold: float = 0.35,
) -> DatasetProfile:
    """
    Summarize a genes x spots APA matrix for automatic preset selection.

    High-resolution-like status is intentionally based on data shape rather
    than platform name, because pseudo-bins, beads, cells, and segmentation
    bins can all expose the same sparse-APA failure mode.
    """
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 2:
        raise ValueError("values must have shape (n_genes, n_spots)")

    n_genes, n_spots = arr.shape
    finite = np.isfinite(arr)
    positive = finite & (arr > 0)

    if input_type == "apa_index":
        observed = finite
    else:
        observed = positive

    total = max(1, n_genes * n_spots)
    finite_fraction = float(finite.sum() / total)
    observed_fraction = float(observed.sum() / total)
    positive_fraction = float(positive.sum() / total)
    median_observed_per_spot = float(np.median(observed.sum(axis=0))) if n_spots else 0.0
    median_positive_per_spot = float(np.median(positive.sum(axis=0))) if n_spots else 0.0

    reasons: List[str] = []
    if n_spots >= highres_spot_threshold:
        reasons.append(f"n_spots>={highres_spot_threshold}")
    if observed_fraction < sparse_observed_fraction_threshold:
        reasons.append(f"observed_fraction<{sparse_observed_fraction_threshold}")

    # Some APA-index matrices encode missing/low-coverage bins as zeros rather
    # than NaN. Treat this as a high-resolution hint only when the data also has
    # enough bins to make dense exact-GP defaults risky.
    low_positive_gene_count = max(5.0, 0.10 * max(float(n_genes), 1.0))
    if (
        n_spots >= 300
        and positive_fraction < 0.15
        and median_positive_per_spot <= low_positive_gene_count
    ):
        reasons.append("sparse_positive_APA_signal")

    return DatasetProfile(
        n_genes=int(n_genes),
        n_spots=int(n_spots),
        finite_fraction=finite_fraction,
        observed_fraction=observed_fraction,
        positive_fraction=positive_fraction,
        median_observed_per_spot=median_observed_per_spot,
        median_positive_per_spot=median_positive_per_spot,
        highres_like=bool(reasons),
        highres_reasons=reasons,
    )


def resolve_analysis_preset(requested_preset: str, profile: DatasetProfile) -> str:
    """Resolve ``auto`` to a concrete analysis preset."""
    if requested_preset not in VALID_ANALYSIS_PRESETS:
        valid = ", ".join(VALID_ANALYSIS_PRESETS)
        raise ValueError(f"analysis_preset must be one of: {valid}")
    if requested_preset == "auto":
        return "highres_accuracy" if profile.highres_like else "standard"
    return requested_preset
