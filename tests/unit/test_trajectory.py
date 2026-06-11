"""Unit tests for spatial trajectory analysis (Tasks 3.6–3.9)."""

import pytest
import numpy as np
import pandas as pd
from spagapa.analysis.trajectory import (
    TrajectoryBuilder,
    TrajectoryAnalyzer,
    build_spatial_trajectory,
    analyse_apa_trajectory,
)


@pytest.fixture
def grid_data():
    """Simple 10×10 grid with a linear APA gradient."""
    np.random.seed(42)
    side = 10
    x = np.repeat(np.arange(side), side).astype(float)
    y = np.tile(np.arange(side), side).astype(float)
    coords = np.column_stack([x, y])
    n_spots = len(coords)
    n_genes = 8
    apa = np.random.rand(n_genes, n_spots) * 0.3 + 0.35
    # Gene 0: increasing gradient along x
    apa[0] = x / (side - 1)
    # Gene 1: decreasing gradient
    apa[1] = 1.0 - x / (side - 1)
    return coords, apa


class TestTrajectoryBuilder:

    def test_build_trajectory_returns_path(self, grid_data):
        coords, _ = grid_data
        builder = TrajectoryBuilder(n_neighbors=4)
        path = builder.build_trajectory(coords, coords[0], coords[-1])
        assert len(path) >= 2
        assert path[0] == 0 or np.allclose(coords[path[0]], coords[0])

    def test_order_spots_pseudospace_range(self, grid_data):
        coords, _ = grid_data
        builder = TrajectoryBuilder(n_neighbors=4)
        path = builder.build_trajectory(coords, coords[0], coords[-1])
        _, ps = builder.order_spots_along_path(path, coords)
        assert ps[0] == pytest.approx(0.0)
        assert ps[-1] == pytest.approx(1.0)
        assert np.all(np.diff(ps) >= 0)

    def test_compute_trajectory_distance(self, grid_data):
        coords, _ = grid_data
        builder = TrajectoryBuilder(n_neighbors=4)
        path = builder.build_trajectory(coords, coords[0], coords[-1])
        dist = builder.compute_trajectory_distance(path, coords)
        assert dist > 0

    def test_validate_trajectory_ok(self, grid_data):
        coords, _ = grid_data
        builder = TrajectoryBuilder(n_neighbors=4)
        path = builder.build_trajectory(coords, coords[0], coords[-1])
        assert builder.validate_trajectory(path, coords, min_length=3)

    def test_validate_trajectory_too_short(self, grid_data):
        coords, _ = grid_data
        builder = TrajectoryBuilder()
        assert not builder.validate_trajectory(np.array([0, 1]), coords, min_length=5)

    def test_same_start_end(self, grid_data):
        coords, _ = grid_data
        builder = TrajectoryBuilder(n_neighbors=4)
        path = builder.build_trajectory(coords, coords[5], coords[5])
        assert len(path) == 1


class TestTrajectoryAnalyzer:

    def test_fit_apa_curves_shape(self, grid_data):
        coords, apa = grid_data
        path, ps = build_spatial_trajectory(coords, coords[0], coords[-1])
        analyzer = TrajectoryAnalyzer(n_bins=10)
        curves = analyzer.fit_apa_curves(apa, path, ps)
        assert curves.shape[0] == 10          # bins
        assert curves.shape[1] == apa.shape[0]  # genes

    def test_fit_apa_curves_gene_names(self, grid_data):
        coords, apa = grid_data
        path, ps = build_spatial_trajectory(coords, coords[0], coords[-1])
        names = [f'G{i}' for i in range(apa.shape[0])]
        analyzer = TrajectoryAnalyzer(n_bins=10)
        curves = analyzer.fit_apa_curves(apa, path, ps, gene_names=names)
        assert list(curves.columns) == names

    def test_detect_switch_points_returns_dict(self, grid_data):
        coords, apa = grid_data
        path, ps = build_spatial_trajectory(coords, coords[0], coords[-1])
        analyzer = TrajectoryAnalyzer(n_bins=10)
        curves = analyzer.fit_apa_curves(apa, path, ps)
        switches = analyzer.detect_switch_points(curves)
        assert isinstance(switches, dict)
        assert set(switches.keys()) == set(curves.columns)

    def test_rank_dynamic_genes_variance(self, grid_data):
        coords, apa = grid_data
        path, ps = build_spatial_trajectory(coords, coords[0], coords[-1])
        analyzer = TrajectoryAnalyzer(n_bins=10)
        curves = analyzer.fit_apa_curves(apa, path, ps)
        ranked = analyzer.rank_dynamic_genes(curves, metric='variance')
        assert 'gene' in ranked.columns
        assert 'score' in ranked.columns
        assert ranked['score'].iloc[0] >= ranked['score'].iloc[-1]

    def test_rank_dynamic_genes_range(self, grid_data):
        coords, apa = grid_data
        path, ps = build_spatial_trajectory(coords, coords[0], coords[-1])
        analyzer = TrajectoryAnalyzer(n_bins=10)
        curves = analyzer.fit_apa_curves(apa, path, ps)
        ranked = analyzer.rank_dynamic_genes(curves, metric='range')
        assert len(ranked) == apa.shape[0]

    def test_classify_patterns(self, grid_data):
        coords, apa = grid_data
        path, ps = build_spatial_trajectory(coords, coords[0], coords[-1])
        analyzer = TrajectoryAnalyzer(n_bins=10)
        curves = analyzer.fit_apa_curves(apa, path, ps)
        patterns = analyzer.classify_trajectory_patterns(curves)
        assert 'gene' in patterns.columns
        assert 'pattern' in patterns.columns
        valid_patterns = {'increasing', 'decreasing', 'switch', 'stable'}
        assert set(patterns['pattern']).issubset(valid_patterns)

    def test_gradient_gene_classified_correctly(self, grid_data):
        """Gene 0 (x-gradient) should be 'increasing'."""
        coords, apa = grid_data
        path, ps = build_spatial_trajectory(coords, coords[0], coords[-1])
        names = [f'G{i}' for i in range(apa.shape[0])]
        analyzer = TrajectoryAnalyzer(n_bins=10)
        curves = analyzer.fit_apa_curves(apa, path, ps, gene_names=names)
        patterns = analyzer.classify_trajectory_patterns(curves)
        g0_pattern = patterns.loc[patterns['gene'] == 'G0', 'pattern'].values[0]
        assert g0_pattern in ('increasing', 'switch')  # gradient → increasing


class TestConvenienceFunctions:

    def test_build_spatial_trajectory(self, grid_data):
        coords, _ = grid_data
        path, ps = build_spatial_trajectory(coords, coords[0], coords[-1])
        assert len(path) >= 2
        assert len(ps) == len(path)
        assert ps[0] == pytest.approx(0.0)
        assert ps[-1] == pytest.approx(1.0)

    def test_analyse_apa_trajectory_keys(self, grid_data):
        coords, apa = grid_data
        result = analyse_apa_trajectory(apa, coords, coords[0], coords[-1])
        for key in ('path_indices', 'pseudospace', 'curves',
                    'switches', 'ranked_genes', 'patterns'):
            assert key in result

    def test_analyse_apa_trajectory_with_names(self, grid_data):
        coords, apa = grid_data
        names = [f'GENE{i}' for i in range(apa.shape[0])]
        result = analyse_apa_trajectory(apa, coords, coords[0], coords[-1],
                                        gene_names=names)
        assert list(result['curves'].columns) == names
