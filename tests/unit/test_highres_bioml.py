"""Unit tests for high-resolution decoupled BioML recovery."""

import numpy as np

from spagapa.bioml import (
    HighResBioMLConfig,
    expression_knn_impute,
    highres_bioml_recover,
    resolve_highres_bioml_neighbors,
)


def make_highres_inputs():
    rng = np.random.default_rng(12)
    parent_coords = np.array(
        [
            [0.0, 0.0],
            [0.0, 1.0],
            [3.0, 0.0],
            [3.0, 1.0],
        ]
    )
    subbins = 3
    parent_index = np.repeat(np.arange(parent_coords.shape[0]), subbins)
    coords = np.repeat(parent_coords, subbins, axis=0)
    coords = coords + rng.normal(scale=0.03, size=coords.shape)

    expression = np.column_stack(
        [
            coords[:, 0],
            coords[:, 1],
            (parent_index >= 2).astype(float),
        ]
    )
    n_spots = coords.shape[0]
    apa = np.zeros((5, n_spots), dtype=float)
    apa[:, parent_index < 2] = 0.25
    apa[:, parent_index >= 2] = 0.75
    apa += rng.normal(scale=0.02, size=apa.shape)
    apa = np.clip(apa, 0.0, 1.0)
    apa[:, ::4] = np.nan
    sparse_gp = np.nan_to_num(apa, nan=np.nanmean(apa))
    uncertainty = np.full_like(apa, 0.1)
    return apa, coords, expression, sparse_gp, uncertainty, parent_index


def test_expression_knn_impute_fills_missing_values():
    apa, _, expression, _, _, _ = make_highres_inputs()

    filled = expression_knn_impute(apa, expression, k=3)

    assert filled.shape == apa.shape
    assert np.isfinite(filled).all()
    assert filled.min() >= 0.0
    assert filled.max() <= 1.0


def test_resolve_highres_neighbors_adaptive_uses_parent_density():
    parent_index = np.repeat(np.arange(8), 8)
    config = HighResBioMLConfig(n_neighbors=6, adaptive_neighbor_scale=10.0)

    fixed = resolve_highres_bioml_neighbors(64, parent_index, HighResBioMLConfig(
        n_neighbors=6,
        neighbor_mode="fixed",
    ))
    adaptive = resolve_highres_bioml_neighbors(64, parent_index, config)

    assert fixed == 6
    assert adaptive > fixed


def test_highres_bioml_recover_returns_domains_and_metadata():
    apa, coords, expression, sparse_gp, uncertainty, parent_index = make_highres_inputs()
    config = HighResBioMLConfig(
        gp_blend=0.3,
        n_neighbors=3,
        expression_knn_k=3,
        adaptive_neighbor_scale=3.0,
    )

    result = highres_bioml_recover(
        apa,
        coords,
        expression,
        sparse_gp=sparse_gp,
        uncertainty=uncertainty,
        parent_index=parent_index,
        n_domains=2,
        config=config,
    )

    assert result.recovered.shape == apa.shape
    assert np.isfinite(result.recovered).all()
    assert result.labels.shape == (coords.shape[0],)
    assert len(np.unique(result.labels)) == 2
    assert result.graph.shape == (coords.shape[0], coords.shape[0])
    assert result.metadata["mode"] == "highres_bioml"
    assert result.metadata["apa_source_effective"] == "expression_knn"
