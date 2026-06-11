"""
Unit tests for expression-informed GP imputation.
"""

import numpy as np
import pytest

from spagapa.imputation import (
    ExpressionGPImputer,
    ExpressionGPImputerBatch,
    GeometryFeatureBuilder,
)


class TestExpressionGPImputer:
    @pytest.fixture
    def data(self):
        rng = np.random.default_rng(42)
        x = np.repeat(np.arange(5), 5)
        y = np.tile(np.arange(5), 5)
        coordinates = np.column_stack([x, y]).astype(float)

        radius = np.linalg.norm(coordinates - coordinates.mean(axis=0), axis=1)
        values = np.sin(radius * 0.5) + rng.normal(0, 0.05, size=len(radius))
        mask = np.ones(len(values), dtype=bool)
        mask[rng.choice(len(values), size=8, replace=False)] = False

        expression_embedding = np.column_stack(
            [
                np.sin(x * 0.3),
                np.cos(y * 0.3),
                radius,
            ]
        )
        return coordinates, values, mask, expression_embedding

    def test_spatial_radial_predict(self, data):
        coordinates, values, mask, _ = data
        geometry = GeometryFeatureBuilder(use_theta=True).build(coordinates)

        imputer = ExpressionGPImputer(
            variant="spatial_radial",
            use_theta=True,
            alpha=1e-3,
        )
        imputer.fit(coordinates, values, mask=mask, geometry_features=geometry)
        pred, std = imputer.predict(coordinates, geometry_features=geometry, return_std=True)

        assert pred.shape == values.shape
        assert std.shape == values.shape
        assert np.all(std >= 0)

    def test_additive_impute_preserves_observed(self, data):
        coordinates, values, mask, expression_embedding = data

        imputer = ExpressionGPImputer(
            variant="additive",
            alpha=1e-3,
            lambda_expr=0.5,
        )
        imputed, uncertainty = imputer.impute(
            coordinates,
            values,
            mask=mask,
            expression_embedding=expression_embedding,
            return_uncertainty=True,
        )

        assert imputed.shape == values.shape
        assert uncertainty.shape == values.shape
        np.testing.assert_allclose(imputed[mask], values[mask])
        np.testing.assert_allclose(uncertainty[mask], 0.0)

    def test_additive_requires_expression_embedding(self, data):
        coordinates, values, mask, _ = data
        imputer = ExpressionGPImputer(variant="additive")

        with pytest.raises(ValueError):
            imputer.fit(coordinates, values, mask=mask)

    def test_product_impute_preserves_observed(self, data):
        coordinates, values, mask, expression_embedding = data

        imputer = ExpressionGPImputer(
            variant="product",
            alpha=1e-3,
            lambda_expr=0.5,
            product_offset=1.0,
        )
        imputed, uncertainty = imputer.impute(
            coordinates,
            values,
            mask=mask,
            expression_embedding=expression_embedding,
            return_uncertainty=True,
        )

        assert imputed.shape == values.shape
        assert uncertainty.shape == values.shape
        np.testing.assert_allclose(imputed[mask], values[mask])
        np.testing.assert_allclose(uncertainty[mask], 0.0)
        assert np.all(uncertainty[~mask] >= 0)

    def test_product_requires_expression_embedding(self, data):
        coordinates, values, mask, _ = data
        imputer = ExpressionGPImputer(variant="product")

        with pytest.raises(ValueError):
            imputer.fit(coordinates, values, mask=mask)

    def test_adaptive_additive_impute_preserves_observed(self, data):
        coordinates, values, mask, expression_embedding = data

        imputer = ExpressionGPImputer(
            variant="adaptive_additive",
            alpha=1e-3,
            lambda_expr=0.5,
            gate_tau=0.8,
        )
        imputed, uncertainty = imputer.impute(
            coordinates,
            values,
            mask=mask,
            expression_embedding=expression_embedding,
            return_uncertainty=True,
        )

        assert imputed.shape == values.shape
        assert uncertainty.shape == values.shape
        np.testing.assert_allclose(imputed[mask], values[mask])
        np.testing.assert_allclose(uncertainty[mask], 0.0)
        assert np.all(uncertainty[~mask] >= 0)

    def test_adaptive_additive_requires_expression_embedding(self, data):
        coordinates, values, mask, _ = data
        imputer = ExpressionGPImputer(variant="adaptive_additive")

        with pytest.raises(ValueError):
            imputer.fit(coordinates, values, mask=mask)

    def test_layer_local_impute_preserves_observed(self, data):
        coordinates, values, mask, expression_embedding = data

        imputer = ExpressionGPImputer(
            variant="layer_local",
            alpha=1e-3,
            lambda_expr=0.5,
            product_offset=0.5,
            local_k=5,
            layer_gate_mode="radius",
        )
        imputed, uncertainty = imputer.impute(
            coordinates,
            values,
            mask=mask,
            expression_embedding=expression_embedding,
            return_uncertainty=True,
        )

        assert imputed.shape == values.shape
        assert uncertainty.shape == values.shape
        np.testing.assert_allclose(imputed[mask], values[mask])
        np.testing.assert_allclose(uncertainty[mask], 0.0)
        assert np.all(uncertainty[~mask] >= 0)

    def test_layer_local_requires_expression_embedding(self, data):
        coordinates, values, mask, _ = data
        imputer = ExpressionGPImputer(variant="layer_local")

        with pytest.raises(ValueError):
            imputer.fit(coordinates, values, mask=mask)

    def test_layer_local_pseudolayer_fallback_runs_without_labels(self, data):
        coordinates, values, mask, expression_embedding = data

        imputer = ExpressionGPImputer(
            variant="layer_local",
            alpha=1e-3,
            lambda_expr=0.5,
            product_offset=0.5,
            local_k=5,
            layer_gate_mode="pseudolayer",
        )
        pred, std = imputer.impute(
            coordinates,
            values,
            mask=mask,
            expression_embedding=expression_embedding,
            return_uncertainty=True,
        )

        assert pred.shape == values.shape
        assert std.shape == values.shape
        assert imputer._pseudolayer_edges_fit is not None

    def test_layer_local_gate_is_symmetric(self, data):
        coordinates, values, mask, expression_embedding = data

        imputer = ExpressionGPImputer(
            variant="layer_local",
            alpha=1e-3,
            local_k=5,
            layer_gate_mode="radius",
        )
        imputer.fit(
            coordinates,
            values,
            mask=mask,
            expression_embedding=expression_embedding,
        )
        gate = imputer._local_gate_kernel(
            imputer._train_local,
            imputer._train_local,
            imputer._local_gate_tau_fit,
        )

        assert gate.shape == (int(mask.sum()), int(mask.sum()))
        np.testing.assert_allclose(gate, gate.T)
        assert np.all(gate >= 0)
        assert np.all(gate <= 1)

    def test_adaptive_and_additive_outputs_differ(self, data):
        coordinates, values, mask, expression_embedding = data

        additive = ExpressionGPImputer(
            variant="additive",
            alpha=1e-3,
            lambda_expr=0.5,
        )
        adaptive = ExpressionGPImputer(
            variant="adaptive_additive",
            alpha=1e-3,
            lambda_expr=0.5,
            gate_tau=0.8,
        )

        pred_add, _ = additive.impute(
            coordinates,
            values,
            mask=mask,
            expression_embedding=expression_embedding,
        )
        pred_adaptive, _ = adaptive.impute(
            coordinates,
            values,
            mask=mask,
            expression_embedding=expression_embedding,
        )

        assert not np.allclose(pred_add[~mask], pred_adaptive[~mask])

    def test_adaptive_different_gate_tau_outputs_differ(self, data):
        coordinates, values, mask, expression_embedding = data

        adaptive_small_tau = ExpressionGPImputer(
            variant="adaptive_additive",
            alpha=1e-3,
            lambda_expr=0.5,
            gate_tau=0.3,
        )
        adaptive_large_tau = ExpressionGPImputer(
            variant="adaptive_additive",
            alpha=1e-3,
            lambda_expr=0.5,
            gate_tau=2.0,
        )

        pred_small, _ = adaptive_small_tau.impute(
            coordinates,
            values,
            mask=mask,
            expression_embedding=expression_embedding,
        )
        pred_large, _ = adaptive_large_tau.impute(
            coordinates,
            values,
            mask=mask,
            expression_embedding=expression_embedding,
        )

        assert not np.allclose(pred_small[~mask], pred_large[~mask])

    def test_product_and_additive_outputs_differ(self, data):
        coordinates, values, mask, expression_embedding = data

        additive = ExpressionGPImputer(
            variant="additive",
            alpha=1e-3,
            lambda_expr=0.5,
        )
        product = ExpressionGPImputer(
            variant="product",
            alpha=1e-3,
            lambda_expr=0.5,
            product_offset=1.0,
        )

        pred_add, _ = additive.impute(
            coordinates,
            values,
            mask=mask,
            expression_embedding=expression_embedding,
        )
        pred_product, _ = product.impute(
            coordinates,
            values,
            mask=mask,
            expression_embedding=expression_embedding,
        )

        assert not np.allclose(pred_add[~mask], pred_product[~mask])

    def test_negative_product_offset_raises(self):
        with pytest.raises(ValueError):
            ExpressionGPImputer(variant="product", product_offset=-0.1)

    def test_invalid_gate_mode_raises(self):
        with pytest.raises(ValueError):
            ExpressionGPImputer(variant="adaptive_additive", gate_mode="layer")

    def test_non_positive_gate_tau_raises(self):
        with pytest.raises(ValueError):
            ExpressionGPImputer(variant="adaptive_additive", gate_tau=0.0)

    def test_invalid_layer_gate_mode_raises(self):
        with pytest.raises(ValueError):
            ExpressionGPImputer(variant="layer_local", layer_gate_mode="metadata")

    def test_non_positive_local_k_raises(self):
        with pytest.raises(ValueError):
            ExpressionGPImputer(variant="layer_local", local_k=0)

    def test_predict_before_fit_raises(self, data):
        coordinates, _, _, expression_embedding = data
        imputer = ExpressionGPImputer(variant="additive")

        with pytest.raises(ValueError):
            imputer.predict(coordinates, expression_embedding=expression_embedding)


class TestExpressionGPImputerBatch:
    def test_batch_impute(self):
        rng = np.random.default_rng(42)
        coordinates = rng.random((20, 2)) * 5
        expression_embedding = rng.normal(size=(20, 4))

        values = np.zeros((3, 20), dtype=float)
        mask = np.zeros_like(values, dtype=bool)
        for gene_idx in range(values.shape[0]):
            observed = rng.choice(20, size=8, replace=False)
            mask[gene_idx, observed] = True
            values[gene_idx, observed] = rng.normal(loc=gene_idx, scale=0.2, size=len(observed))

        base = ExpressionGPImputer(
            variant="product",
            alpha=1e-3,
            lambda_expr=1.0,
            product_offset=1.0,
        )
        batch = base.fit_batch(
            coordinates,
            values,
            mask=mask,
            expression_embedding=expression_embedding,
            n_jobs=1,
            verbose=False,
        )

        assert isinstance(batch, ExpressionGPImputerBatch)
        assert len(batch.imputers_) == 3

        imputed, uncertainty = batch.impute(return_uncertainty=True)
        assert imputed.shape == values.shape
        assert uncertainty.shape == values.shape
        np.testing.assert_allclose(imputed[mask], values[mask])
        np.testing.assert_allclose(uncertainty[mask], 0.0)
