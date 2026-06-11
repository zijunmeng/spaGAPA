"""
Command-line interface for spaGAPA.

Usage:
    spagapa run       Run the complete pipeline
    spagapa impute    Run GP imputation only
    spagapa diff      Differential APA analysis
    spagapa plot      Generate visualization
"""

import sys
from pathlib import Path
from typing import Optional

import click
import numpy as np
import pandas as pd


@click.group()
@click.version_option(version='0.1.0', prog_name='spaGAPA')
def main():
    """spaGAPA: Spatial GP-based APA Analyzer.

    A comprehensive toolkit for analyzing alternative polyadenylation
    in spatial transcriptomics data using Gaussian process imputation
    with uncertainty quantification.
    """
    pass


# ── shared options ───────────────────────────────────────────────────

_opt_neighbors = click.option(
    '--n-neighbors', '-k', default=6, show_default=True,
    help='Number of spatial neighbors.'
)
_opt_kernel = click.option(
    '--kernel', '-K', default='matern', show_default=True,
    type=click.Choice(['rbf', 'matern', 'auto']),
    help='GP kernel type.'
)
_opt_sparse = click.option(
    '--sparse/--no-sparse', default=False, show_default=True,
    help='Use sparse GP approximation.'
)
_opt_fdr = click.option(
    '--fdr', '-F', default=0.05, show_default=True,
    help='FDR threshold.'
)
_opt_verbose = click.option(
    '--verbose/--quiet', '-v/-q', default=True,
    help='Print progress messages.'
)
_opt_output = click.option(
    '--output', '-o', default='./spagapa_results', show_default=True,
    help='Output directory.'
)
_opt_coords = click.option(
    '--coordinates', '-c', required=True,
    help='Path to spatial coordinates CSV/TSV.'
)
_opt_uncertainty = click.option(
    '--use-uncertainty/--no-uncertainty', default=True, show_default=True,
    help='Use imputation uncertainty to weight downstream analyses.'
)


# ── run ──────────────────────────────────────────────────────────────

@main.command()
@click.option('--bam', '-b', default=None, help='Path to BAM file.')
@click.option('--annotation', '-a', default=None, help='Path to GTF/GFF annotation.')
@click.option('--apa-matrix', '-m', default=None, help='Path to APA matrix CSV/TSV/NPY.')
@click.option('--matrix-orientation', default='genes_by_spots', show_default=True,
              type=click.Choice(['genes_by_spots', 'spots_by_genes']),
              help='Orientation of --apa-matrix.')
@_opt_coords
@click.option('--dataset', '-d', default=None,
              help='Path to saved APADataset (.h5ad).')
@click.option('--gp-alpha', default=1e-10, show_default=True,
              help='Exact GP noise/regularization parameter.')
@click.option('--gp-n-restarts', default=1, show_default=True,
              help='Number of GP hyperparameter optimizer restarts.')
@click.option('--expression-matrix', default=None,
              help='Optional expression matrix CSV/TSV/NPY for BioML.')
@click.option('--expression-orientation', default='genes_by_spots', show_default=True,
              type=click.Choice(['genes_by_spots', 'spots_by_genes']),
              help='Orientation of --expression-matrix.')
@click.option('--no-impute', is_flag=True, help='Skip GP imputation.')
@click.option('--no-quantify', is_flag=True, help='Skip APA quantification.')
@click.option('--no-domains', is_flag=True, help='Skip domain identification.')
@click.option('--enable-bioml/--disable-bioml', default=False, show_default=True,
              help='Use BioML multi-view graph domain recovery.')
@click.option('--bioml-domain-method', default='spectral', show_default=True,
              type=click.Choice(['spectral', 'kmeans']),
              help='BioML domain detector.')
@click.option('--bioml-rank', default=8, show_default=True,
              help='BioML low-rank dimension.')
@click.option('--bioml-lambda-graph', default=0.5, show_default=True,
              help='BioML graph regularization strength.')
@click.option('--bioml-lambda-l2', default=1e-2, show_default=True,
              help='BioML ridge regularization strength.')
@click.option('--bioml-max-iter', default=20, show_default=True,
              help='BioML maximum ALS iterations.')
@click.option('--bioml-n-neighbors', default=15, show_default=True,
              help='BioML graph KNN size.')
@click.option('--bioml-blend', default=0.1, show_default=True,
              help='Blend BioML factorized matrix into GP imputed values.')
@click.option('--bioml-spatial-weight', default=0.4, show_default=True,
              help='BioML spatial graph weight.')
@click.option('--bioml-expression-weight', default=0.4, show_default=True,
              help='BioML expression graph weight.')
@click.option('--bioml-apa-weight', default=0.2, show_default=True,
              help='BioML APA graph weight.')
@click.option('--expression-n-components', default=10, show_default=True,
              help='PCA dimensions for --expression-matrix.')
@click.option('--diff', is_flag=True, help='Enable differential analysis.')
@click.option('--no-svapa', is_flag=True, help='Skip SVAPA detection.')
@click.option('--n-domains', default=None, type=int, help='Number of spatial domains.')
@_opt_neighbors
@_opt_kernel
@_opt_sparse
@_opt_fdr
@_opt_uncertainty
@_opt_verbose
@_opt_output
def run(**kwargs):
    """Run the complete spaGAPA pipeline."""
    from spagapa.pipeline import SpaGAPA
    from spagapa.core import APADataset

    # Load dataset
    dataset = None
    if kwargs['dataset']:
        click.echo(f"Loading dataset from {kwargs['dataset']}...")
        dataset = APADataset.load(kwargs['dataset'])

    spa = SpaGAPA(
        n_neighbors=kwargs['n_neighbors'],
        kernel_type=kwargs['kernel'],
        gp_alpha=kwargs['gp_alpha'],
        gp_n_restarts_optimizer=kwargs['gp_n_restarts'],
        use_sparse_gp=kwargs['sparse'],
        use_bioml=kwargs['enable_bioml'],
        bioml_rank=kwargs['bioml_rank'],
        bioml_lambda_graph=kwargs['bioml_lambda_graph'],
        bioml_lambda_l2=kwargs['bioml_lambda_l2'],
        bioml_max_iter=kwargs['bioml_max_iter'],
        bioml_n_neighbors=kwargs['bioml_n_neighbors'],
        bioml_blend=kwargs['bioml_blend'],
        bioml_domain_method=kwargs['bioml_domain_method'],
        bioml_spatial_weight=kwargs['bioml_spatial_weight'],
        bioml_expression_weight=kwargs['bioml_expression_weight'],
        bioml_apa_weight=kwargs['bioml_apa_weight'],
        expression_n_components=kwargs['expression_n_components'],
        verbose=kwargs['verbose'],
    )

    results = spa.run(
        bam_file=kwargs['bam'],
        apa_matrix=kwargs['apa_matrix'],
        coordinates=kwargs['coordinates'],
        annotation=kwargs['annotation'],
        matrix_orientation=kwargs['matrix_orientation'],
        expression_matrix=kwargs['expression_matrix'],
        expression_orientation=kwargs['expression_orientation'],
        dataset=dataset,
        impute=not kwargs['no_impute'],
        quantify=not kwargs['no_quantify'],
        identify_domains=not kwargs['no_domains'],
        use_bioml=kwargs['enable_bioml'],
        differential_analysis=kwargs['diff'],
        detect_svapa=not kwargs['no_svapa'],
        n_domains=kwargs['n_domains'],
        fdr_threshold=kwargs['fdr'],
        use_uncertainty_weights=kwargs['use_uncertainty'],
    )

    spa.save_results(kwargs['output'])
    click.echo(f"\n✅ Results saved to {kwargs['output']}/")

    # Summary
    svapa = results.get('svapa_genes')
    if svapa is not None:
        n_sig = int(svapa['significant'].sum())
        click.echo(f"   SVAPA genes detected: {n_sig}")

    domains = results.get('domains')
    if domains is not None:
        click.echo(f"   Spatial domains: {domains['n_domains']}")
        if domains.get('method') == 'spagapa_bioml':
            click.echo(f"   BioML domain method: {domains.get('domain_method')}")

    qc = results.get('qc_report')
    if qc and isinstance(qc, dict):
        cov = qc.get('coverage_metrics', {})
        if cov:
            click.echo(f"   Coverage: {cov.get('observed_spots', '?')} spots")


# ── impute ───────────────────────────────────────────────────────────

@main.command()
@click.argument('input_file', type=click.Path(exists=True))
@_opt_coords
@_opt_kernel
@_opt_sparse
@_opt_neighbors
@click.option('--output-imputed', '-I', default='imputed.npy',
              help='Output path for imputed values.')
@click.option('--output-uncertainty', '-U', default='uncertainty.npy',
              help='Output path for uncertainty.')
@_opt_verbose
def impute(input_file, coordinates, kernel, sparse, n_neighbors,
           output_imputed, output_uncertainty, verbose):
    """Run GP imputation on an APA matrix.

    INPUT_FILE: Path to APA count matrix (.npy or .csv).
    """
    from spagapa.imputation import GPImputer, SparseGPImputer

    # Load data
    if input_file.endswith('.csv'):
        df = pd.read_csv(input_file, index_col=0)
        apa = df.values
    else:
        apa = np.load(input_file)

    coords = pd.read_csv(coordinates, index_col=0).values

    click.echo(f"Loaded: {apa.shape[0]} genes × {apa.shape[1]} spots")
    click.echo(f"Coordinates: {coords.shape[0]} spots")

    # Impute
    if sparse:
        imputer = SparseGPImputer(kernel_type=kernel)
    else:
        imputer = GPImputer(kernel_type=kernel)

    click.echo(f"Fitting GP ({kernel}{' sparse' if sparse else ''})...")
    imputer.fit_batch(coords, apa.T, n_jobs=-1, verbose=verbose)

    imputed_T, uncertainty_T = imputer.predict(coords, return_std=True)
    imputed = imputed_T.T  # (genes, spots)
    uncertainty = uncertainty_T.T

    np.save(output_imputed, imputed)
    np.save(output_uncertainty, uncertainty)

    mae = np.nanmean(np.abs(apa - imputed))
    click.echo(f"\n✅ Imputation complete (MAE={mae:.3f})")
    click.echo(f"   Imputed → {output_imputed}")
    click.echo(f"   Uncertainty → {output_uncertainty}")


# ── diff ─────────────────────────────────────────────────────────────

@main.command()
@click.argument('apa_file', type=click.Path(exists=True))
@click.argument('labels_file', type=click.Path(exists=True))
@_opt_coords
@click.option('--uncertainty-file', '-U', default=None,
              help='Path to uncertainty matrix (.npy).')
@click.option('--method', '-m', default='wilcoxon',
              type=click.Choice(['wilcoxon', 't-test', 'permutation']),
              help='Statistical test method.')
@_opt_fdr
@click.option('--logfc', default=0.5, show_default=True,
              help='Log2 fold change threshold.')
@_opt_output
@_opt_verbose
def diff(apa_file, labels_file, coordinates, uncertainty_file, method,
         fdr, logfc, output, verbose):
    """Differential APA analysis between spatial domains.

    APA_FILE: Path to APA index matrix (.npy or .csv).
    LABELS_FILE: Path to domain labels (.csv or .npy).
    """
    from spagapa.analysis import find_domain_markers

    apa = _load_matrix(apa_file)
    labels = _load_vector(labels_file)
    coords = pd.read_csv(coordinates, index_col=0).values
    uncertainty = np.load(uncertainty_file) if uncertainty_file else None

    click.echo(f"APA: {apa.shape[0]} genes × {apa.shape[1]} spots")
    click.echo(f"Domains: {len(np.unique(labels))}")
    if uncertainty is not None:
        click.echo("Using uncertainty weights ✓")

    markers = find_domain_markers(
        apa, labels, method=method, padj_threshold=fdr,
        logfc_threshold=logfc, uncertainty=uncertainty,
    )

    Path(output).mkdir(parents=True, exist_ok=True)
    for domain, df in markers.items():
        path = Path(output) / f"domain_{domain}_markers.csv"
        df.to_csv(path, index=False)
        click.echo(f"  Domain {domain}: {len(df)} markers → {path.name}")


# ── plot ─────────────────────────────────────────────────────────────

@main.command()
@click.argument('apa_file', type=click.Path(exists=True))
@_opt_coords
@click.option('--gene', '-g', default=None, help='Gene name to plot.')
@click.option('--gene-index', '-i', default=0, help='Gene row index.')
@click.option('--type', '-t', default='spatial',
              type=click.Choice(['spatial', 'volcano', 'heatmap', 'qc']),
              help='Plot type.')
@click.option('--uncertainty-file', '-U', default=None,
              help='Path to uncertainty matrix (.npy).')
@click.option('--domains-file', '-D', default=None,
              help='Path to domain labels for coloring.')
@_opt_output
@click.option('--dpi', default=300, show_default=True)
def plot(apa_file, coordinates, gene, gene_index, type, uncertainty_file,
         domains_file, output, dpi):
    """Generate visualizations.

    APA_FILE: Path to APA matrix (.npy or .csv).
    """
    from spagapa.visualization import SpatialPlotter, StatisticalPlotter, QCPlotter
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    apa = _load_matrix(apa_file)
    coords = pd.read_csv(coordinates, index_col=0).values
    uncertainty = np.load(uncertainty_file) if uncertainty_file else None
    domains = _load_vector(domains_file) if domains_file else None

    Path(output).mkdir(parents=True, exist_ok=True)

    if type == 'spatial':
        plotter = SpatialPlotter()
        idx = list(apa.index).index(gene) if gene else gene_index
        vals = apa.values[idx] if hasattr(apa, 'values') else apa[idx]
        fig = plotter.plot_spatial_apa(
            coords[:, 0], coords[:, 1], vals,
            title=gene or f"Gene_{idx}",
            domain_labels=domains,
        )
    elif type == 'qc' and uncertainty is not None:
        plotter = QCPlotter()
        fig = plotter.plot_imputation_quality(
            apa.values if hasattr(apa, 'values') else apa,
            uncertainty, coords,
        )
    else:
        click.echo(f"Plot type '{type}' requires additional configuration.")
        return

    path = Path(output) / f"{type}_{gene or gene_index}.png"
    fig.savefig(path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    click.echo(f"✅ Saved → {path}")


# ── helpers ──────────────────────────────────────────────────────────

def _load_matrix(path: str) -> np.ndarray:
    if path.endswith('.csv'):
        return pd.read_csv(path, index_col=0).values
    return np.load(path)


def _load_vector(path: str) -> np.ndarray:
    if path.endswith('.csv'):
        return pd.read_csv(path, index_col=0).values.flatten()
    return np.load(path).flatten()


if __name__ == '__main__':
    main()
