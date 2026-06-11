"""Unit tests for user-facing analysis preset resolution."""

import numpy as np

from spagapa.presets import profile_spatial_apa_matrix, resolve_analysis_preset


def test_auto_preset_resolves_standard_for_dense_small_matrix():
    values = np.ones((10, 50), dtype=float)
    profile = profile_spatial_apa_matrix(values, input_type="apa_index")

    assert profile.highres_like is False
    assert resolve_analysis_preset("auto", profile) == "standard"


def test_auto_preset_resolves_highres_for_many_bins():
    values = np.ones((10, 1200), dtype=float)
    profile = profile_spatial_apa_matrix(values, input_type="apa_index")

    assert profile.highres_like is True
    assert "n_spots>=1000" in profile.highres_reasons
    assert resolve_analysis_preset("auto", profile) == "highres_accuracy"


def test_auto_preset_resolves_highres_for_sparse_observations():
    values = np.full((10, 80), np.nan, dtype=float)
    values[:, :20] = 0.5
    profile = profile_spatial_apa_matrix(values, input_type="apa_index")

    assert profile.observed_fraction == 0.25
    assert profile.highres_like is True
    assert "observed_fraction<0.35" in profile.highres_reasons
    assert resolve_analysis_preset("auto", profile) == "highres_accuracy"


def test_explicit_preset_is_not_changed_by_profile():
    values = np.ones((10, 1200), dtype=float)
    profile = profile_spatial_apa_matrix(values, input_type="apa_index")

    assert resolve_analysis_preset("standard", profile) == "standard"
    assert resolve_analysis_preset("highres_fast", profile) == "highres_fast"
