"""Toy end-to-end tests for the BIB-ready pipeline contract."""

import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

from spagapa import APADataset, SpaGAPA
from spagapa.cli import main as cli_main
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


def test_pipeline_toy_bioml_mainline():
    dataset = make_toy_dataset(n_genes=5, grid_size=5)
    rng = np.random.default_rng(42)
    expression_embedding = np.column_stack(
        [
            dataset.coords[:, 0],
            dataset.coords[:, 1],
            rng.normal(scale=0.01, size=dataset.n_spots),
        ]
    )

    pipeline = SpaGAPA(
        kernel_type='matern',
        n_neighbors=4,
        min_spots=3,
        use_bioml=True,
        bioml_rank=3,
        bioml_max_iter=4,
        bioml_n_neighbors=4,
        expression_n_components=3,
        verbose=False,
    )

    results = pipeline.run(
        dataset=dataset,
        expression_embedding=expression_embedding,
        impute=True,
        quantify=True,
        identify_domains=True,
        use_bioml=True,
        differential_analysis=False,
        detect_svapa=False,
        n_domains=2,
    )

    domains = results['domains']
    assert domains['method'] == 'spagapa_bioml'
    assert domains['labels'].shape == (dataset.n_spots,)
    assert domains['imputed_values'].shape == (dataset.n_genes, dataset.n_spots)
    assert domains['spot_factors'].shape == (dataset.n_spots, 3)
    assert results['dataset'].has_bioml()
    assert results['dataset'].bioml_imputed.shape == (dataset.n_genes, dataset.n_spots)
    assert results['dataset'].adata.uns['apa']['bioml']['has_expression_view'] is True


def test_analysis_preset_auto_resolves_standard_on_toy_data():
    dataset = make_toy_dataset(n_genes=4, grid_size=4)
    pipeline = SpaGAPA(
        analysis_preset='auto',
        n_neighbors=4,
        min_spots=3,
        verbose=False,
    )

    results = pipeline.run(
        dataset=dataset,
        impute=False,
        quantify=True,
        identify_domains=False,
        differential_analysis=False,
        detect_svapa=False,
    )

    preset = results['analysis_preset']
    assert preset['requested_preset'] == 'auto'
    assert preset['resolved_preset'] == 'standard'
    assert preset['use_bioml'] is False
    assert preset['use_sparse_gp'] is False
    assert dataset.adata.uns['apa']['analysis_preset']['resolved_preset'] == 'standard'


def test_analysis_preset_highres_accuracy_enables_sparse_bioml():
    dataset = make_toy_dataset(n_genes=5, grid_size=5)
    expression_embedding = np.column_stack(
        [
            dataset.coords[:, 0],
            dataset.coords[:, 1],
            dataset.coords[:, 0] + dataset.coords[:, 1],
        ]
    )
    pipeline = SpaGAPA(
        analysis_preset='highres_accuracy',
        n_neighbors=4,
        n_inducing=8,
        min_spots=3,
        bioml_rank=3,
        bioml_max_iter=3,
        bioml_n_neighbors=4,
        verbose=False,
    )

    results = pipeline.run(
        dataset=dataset,
        expression_embedding=expression_embedding,
        impute=True,
        quantify=True,
        identify_domains=True,
        differential_analysis=False,
        detect_svapa=False,
        n_domains=2,
    )

    preset = results['analysis_preset']
    assert preset['resolved_preset'] == 'highres_accuracy'
    assert preset['use_sparse_gp'] is True
    assert preset['use_bioml'] is True
    assert preset['bioml_weights'] == {'spatial': 0.2, 'expression': 0.6, 'apa': 0.2}
    assert results['imputed_values'].shape == (dataset.n_genes, dataset.n_spots)
    assert results['domains']['method'] == 'spagapa_bioml'


def test_analysis_preset_highres_fast_skips_gp_but_runs_bioml():
    dataset = make_toy_dataset(n_genes=5, grid_size=5)
    expression_embedding = np.column_stack(
        [
            dataset.coords[:, 0],
            dataset.coords[:, 1],
            dataset.coords[:, 0] + dataset.coords[:, 1],
        ]
    )
    pipeline = SpaGAPA(
        analysis_preset='highres_fast',
        n_neighbors=4,
        min_spots=3,
        bioml_rank=3,
        bioml_max_iter=3,
        bioml_n_neighbors=4,
        verbose=False,
    )

    results = pipeline.run(
        dataset=dataset,
        expression_embedding=expression_embedding,
        impute=True,
        quantify=True,
        identify_domains=True,
        differential_analysis=False,
        detect_svapa=False,
        n_domains=2,
    )

    preset = results['analysis_preset']
    assert preset['resolved_preset'] == 'highres_fast'
    assert preset['impute'] is False
    assert preset['use_bioml'] is True
    assert results['imputed_values'] is None
    assert results['uncertainty'] is None
    assert results['domains']['method'] == 'spagapa_bioml'


def test_cli_run_with_bioml_from_tables(tmp_path):
    dataset = make_toy_dataset(n_genes=4, grid_size=4)
    matrix = pd.DataFrame(
        dataset.raw_counts,
        index=dataset.gene_names,
        columns=dataset.spot_names,
    )
    coords = pd.DataFrame(
        {
            'barcode': dataset.spot_names,
            'x': dataset.coords[:, 0],
            'y': dataset.coords[:, 1],
        }
    )
    expression = pd.DataFrame(
        np.vstack(
            [
                dataset.coords[:, 0],
                dataset.coords[:, 1],
                dataset.coords[:, 0] + dataset.coords[:, 1],
                np.ones(dataset.n_spots),
            ]
        ),
        index=[f"Expr_{i}" for i in range(4)],
        columns=dataset.spot_names,
    )

    matrix_file = tmp_path / "apa_matrix.csv"
    coords_file = tmp_path / "coordinates.csv"
    expression_file = tmp_path / "expression.csv"
    output_dir = tmp_path / "cli_results"
    matrix.to_csv(matrix_file)
    coords.to_csv(coords_file, index=False)
    expression.to_csv(expression_file)

    runner = CliRunner()
    result = runner.invoke(
        cli_main,
        [
            'run',
            '--apa-matrix',
            str(matrix_file),
            '--coordinates',
            str(coords_file),
            '--expression-matrix',
            str(expression_file),
            '--enable-bioml',
            '--analysis-preset',
            'standard',
            '--no-svapa',
            '--n-domains',
            '2',
            '--bioml-rank',
            '3',
            '--bioml-max-iter',
            '3',
            '--bioml-n-neighbors',
            '4',
            '--n-neighbors',
            '4',
            '--output',
            str(output_dir),
            '--quiet',
        ],
    )

    assert result.exit_code == 0, result.output
    assert (output_dir / "dataset.h5ad").exists()
    assert (output_dir / "domains.csv").exists()
    assert (output_dir / "analysis_preset.json").exists()
    assert (output_dir / "bioml_metadata.json").exists()
    assert (output_dir / "bioml_spot_factors.npy").exists()
