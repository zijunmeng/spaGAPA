"""
Unit tests for CPU-friendly BioML components.
"""

import numpy as np
import pytest

from spagapa.bioml import (
    BioMLDomainDetector,
    GraphRegularizedAPAFactorizer,
    MultiViewGraphBuilder,
)


@pytest.fixture
def bioml_data():
    rng = np.random.default_rng(42)
    n_spots = 18
    n_genes = 6

    coords_a = rng.normal(loc=0.0, scale=0.2, size=(9, 2))
    coords_b = rng.normal(loc=3.0, scale=0.2, size=(9, 2))
    coordinates = np.vstack([coords_a, coords_b])

    expression = np.vstack(
        [
            rng.normal(loc=0.0, scale=0.2, size=(9, 4)),
            rng.normal(loc=2.0, scale=0.2, size=(9, 4)),
        ]
    )

    apa = np.zeros((n_genes, n_spots), dtype=float)
    apa[:, :9] = rng.normal(loc=0.25, scale=0.03, size=(n_genes, 9))
    apa[:, 9:] = rng.normal(loc=0.75, scale=0.03, size=(n_genes, 9))
    apa[0, 1] = np.nan
    apa[3, 12] = np.nan
    mask = np.isfinite(apa)
    uncertainty = np.full_like(apa, 0.1)
    uncertainty[:, 9:] = 0.2
    return coordinates, expression, apa, mask, uncertainty


class TestMultiViewGraphBuilder:
    def test_build_multiview_graph(self, bioml_data):
        coordinates, expression, apa, _, uncertainty = bioml_data
        builder = MultiViewGraphBuilder(n_neighbors=4)
        graph = builder.build(
            coordinates,
            expression_embedding=expression,
            apa_matrix=apa,
            uncertainty=uncertainty,
        )

        assert graph.spatial.shape == (18, 18)
        assert graph.expression.shape == (18, 18)
        assert graph.apa.shape == (18, 18)
        assert graph.fused.shape == (18, 18)
        assert graph.fused.nnz > 0
        np.testing.assert_allclose(graph.fused.toarray(), graph.fused.toarray().T)
        assert pytest.approx(sum(graph.weights.values())) == 1.0

    def test_laplacian_row_sums_are_zero(self, bioml_data):
        coordinates, expression, apa, _, uncertainty = bioml_data
        graph = MultiViewGraphBuilder(n_neighbors=4).build(
            coordinates,
            expression_embedding=expression,
            apa_matrix=apa,
            uncertainty=uncertainty,
        )
        lap = graph.laplacian()

        assert lap.shape == (18, 18)
        np.testing.assert_allclose(np.asarray(lap.sum(axis=1)).ravel(), 0.0, atol=1e-8)

    def test_invalid_neighbors_raises(self):
        with pytest.raises(ValueError):
            MultiViewGraphBuilder(n_neighbors=0)


class TestGraphRegularizedAPAFactorizer:
    def test_fit_transform_preserves_observed(self, bioml_data):
        coordinates, expression, apa, mask, uncertainty = bioml_data
        graph = MultiViewGraphBuilder(n_neighbors=4).build(
            coordinates,
            expression_embedding=expression,
            apa_matrix=apa,
            uncertainty=uncertainty,
        )
        confidence = 1.0 / (uncertainty + 1e-6)

        factorizer = GraphRegularizedAPAFactorizer(
            rank=3,
            lambda_graph=0.2,
            max_iter=8,
            random_state=42,
        )
        imputed = factorizer.fit_transform(
            apa,
            graph_laplacian=graph.laplacian(),
            mask=mask,
            confidence=confidence,
        )

        assert imputed.shape == apa.shape
        assert factorizer.gene_factors_.shape == (6, 3)
        assert factorizer.spot_factors_.shape == (18, 3)
        assert factorizer.result_.n_iter >= 1
        np.testing.assert_allclose(imputed[mask], apa[mask])
        assert np.isfinite(imputed).all()

    def test_invalid_rank_raises(self):
        with pytest.raises(ValueError):
            GraphRegularizedAPAFactorizer(rank=0)


class TestBioMLDomainDetector:
    def test_kmeans_domain_detection(self, bioml_data):
        coordinates, expression, apa, mask, uncertainty = bioml_data
        graph = MultiViewGraphBuilder(n_neighbors=4).build(
            coordinates,
            expression_embedding=expression,
            apa_matrix=apa,
            uncertainty=uncertainty,
        )
        factorizer = GraphRegularizedAPAFactorizer(rank=3, max_iter=6).fit(
            apa,
            graph_laplacian=graph.laplacian(),
            mask=mask,
            confidence=1.0 / (uncertainty + 1e-6),
        )

        labels = BioMLDomainDetector(method="kmeans", n_domains=2).fit_predict(
            factorizer.spot_factors_
        )

        assert labels.shape == (18,)
        assert len(np.unique(labels)) == 2

    def test_spectral_domain_detection(self, bioml_data):
        coordinates, expression, apa, _, uncertainty = bioml_data
        graph = MultiViewGraphBuilder(n_neighbors=4).build(
            coordinates,
            expression_embedding=expression,
            apa_matrix=apa,
            uncertainty=uncertainty,
        )
        labels = BioMLDomainDetector(method="spectral", n_domains=2).fit_predict(
            graph=graph.fused
        )

        assert labels.shape == (18,)
        assert len(np.unique(labels)) == 2

    def test_invalid_method_raises(self):
        with pytest.raises(ValueError):
            BioMLDomainDetector(method="invalid")
