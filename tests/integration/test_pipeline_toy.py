"""Toy end-to-end tests for the BIB-ready pipeline contract."""

import numpy as np
import pandas as pd
import pytest

from spagapa import APADataset, SpaGAPA
from spagapa.io import load_spatial_dataset


def make_toy_dataset(n_genes=4, grid_size=5):
    """Create a small genes x spots APA-index dataset."""
    x = np.repeat(np.arange(grid_size), grid_size)
    y = np.tile(np.arange(grid_size), grid_size)
    coords = np.column_stack([x, y]).astype(float)
    n_spots = len(coords)

    values = np.zeros((n_genes, n_spots), dtype=float)
    values[0] = (x - x.min()) / (x.max() - x.min())
    if n_genes > 1:
        values[1] = (y - y.min()) / (y.max() - y.min())
    center_dist = np.sqrt((x - x.mean()) ** 2 + (y - y.mean()) ** 2)
    if n_genes > 2:
        values[2] = 1.0 - center_dist / center_dist.max()
    if n_genes > 3:
        values[3] = 0.45

    dropout = np.zeros_like(values, dtype=bool)
    dropout[:, ::4] = True
    values[dropout] = 0.0

    return APADataset.from_counts(
        values,
        coords,
        gene_names=[f"Gene_{i}" for i in range(n_genes)],
        spot_names=[f"Spot_{i}" for i in range(n_spots)],
    )


@pytest.mark.parametrize("use_sparse_gp", [False, True])
def test_pipeline_toy_end_to_end(use_sparse_gp):
    dataset = make_toy_dataset()
    pipeline = SpaGAPA(
        kernel_type='matern',
        use_sparse_gp=use_sparse_gp,
        n_inducing=6,
        n_neighbors=4,
        min_spots=3,
        verbose=False,
    )

    results = pipeline.run(
        dataset=dataset,
        impute=True,
        quantify=True,
        identify_domains=True,
        differential_analysis=True,
        detect_svapa=True,
        n_domains=2,
        fdr_threshold=0.5,
    )

    assert results['imputed_values'].shape == (dataset.n_genes, dataset.n_spots)
    assert results['uncertainty'].shape == (dataset.n_genes, dataset.n_spots)
    assert results['apa_indices']['APAIndex'].shape == (dataset.n_genes, dataset.n_spots)
    assert results['apa_indices']['calculated'] is False
    assert results['domains']['labels'].shape == (dataset.n_spots,)
    assert set(results['differential'].keys()) == {'Domain0_vs_Domain1'}
    assert {'gene', 'pvalue', 'padj', 'log2fc'}.issubset(
        results['differential']['Domain0_vs_Domain1'].columns
    )
    assert len(results['svapa_genes']) == dataset.n_genes
    assert results['dataset'].raw_counts.shape == (dataset.n_genes, dataset.n_spots)


def test_load_spatial_dataset_from_tables(tmp_path):
    dataset = make_toy_dataset(n_genes=3, grid_size=4)
    matrix = pd.DataFrame(
        dataset.raw_counts,
        index=dataset.gene_names,
        columns=dataset.spot_names,
    )
    coords = pd.DataFrame(
        dataset.coords,
        index=dataset.spot_names,
        columns=['x', 'y'],
    )

    matrix_file = tmp_path / "apa_matrix.csv"
    coords_file = tmp_path / "coordinates.csv"
    matrix.to_csv(matrix_file)
    coords.assign(barcode=coords.index).to_csv(coords_file, index=False)

    loaded = load_spatial_dataset(
        apa_matrix=matrix_file,
        coordinates=coords_file,
        matrix_orientation='genes_by_spots',
    )

    assert loaded.raw_counts.shape == dataset.raw_counts.shape
    np.testing.assert_allclose(loaded.raw_counts, dataset.raw_counts)
    np.testing.assert_allclose(loaded.coords, dataset.coords)
