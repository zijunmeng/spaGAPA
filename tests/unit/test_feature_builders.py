"""
Unit tests for expression and geometry feature builders.
"""

import numpy as np
import pandas as pd

from spagapa.imputation import ExpressionFeatureBuilder, GeometryFeatureBuilder


class TestExpressionFeatureBuilder:
    def test_pca_embedding_shape(self):
        rng = np.random.default_rng(42)
        expression = pd.DataFrame(
            rng.normal(size=(50, 12)),
            index=[f"g{i}" for i in range(50)],
            columns=[f"s{i}" for i in range(12)],
        )

        builder = ExpressionFeatureBuilder(n_components=5)
        embedding = builder.fit_transform(expression)

        assert embedding.shape == (12, 5)

    def test_hvg_selection(self):
        rng = np.random.default_rng(42)
        expression = pd.DataFrame(
            rng.normal(size=(100, 10)),
            index=[f"g{i}" for i in range(100)],
            columns=[f"s{i}" for i in range(10)],
        )

        builder = ExpressionFeatureBuilder(
            n_components=4,
            use_hvg=True,
            n_top_genes=20,
        )
        embedding = builder.fit_transform(expression)

        assert embedding.shape == (10, 4)

    def test_spots_by_genes_orientation(self):
        rng = np.random.default_rng(42)
        expression = rng.normal(size=(8, 30))

        builder = ExpressionFeatureBuilder(
            n_components=3,
            orientation="spots_by_genes",
        )
        embedding = builder.fit_transform(expression)

        assert embedding.shape == (8, 3)


class TestGeometryFeatureBuilder:
    def test_build_radius_and_theta(self):
        coordinates = np.array(
            [
                [0.0, 0.0],
                [1.0, 0.0],
                [0.0, 1.0],
                [1.0, 1.0],
            ]
        )
        builder = GeometryFeatureBuilder(use_theta=True)
        features = builder.build(coordinates)

        assert features["xy"].shape == (4, 2)
        assert features["radius"].shape == (4, 1)
        assert features["theta"].shape == (4, 1)

    def test_build_pseudolayer(self):
        coordinates = np.array(
            [
                [0.0, 0.0],
                [0.2, 0.1],
                [1.0, 1.0],
                [1.2, 1.1],
            ]
        )
        labels = np.array(["inner", "inner", "outer", "outer"])
        builder = GeometryFeatureBuilder(use_theta=False)
        features = builder.build(coordinates, layer_labels=labels)

        assert "pseudolayer" in features
        assert features["pseudolayer"].shape == (4, 1)
        assert set(features["layer_order"].tolist()) == {"inner", "outer"}
