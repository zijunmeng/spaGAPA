"""
Prepare real benchmark datasets for spaGAPA.

This script:
1. Loads 10x Visium H5 files
2. Extracts spatial coordinates
3. Runs scAPAtrap (via R) to get APA matrices
4. Saves standardised CSV files for run_real_benchmark.py

Usage:
    conda activate spagapa
    python scripts/prepare_benchmark_data.py --data_dir data/real --output_dir data/processed
"""

import argparse
import os
import numpy as np
import pandas as pd


def load_visium_coordinates(spatial_dir: str) -> pd.DataFrame:
    """Load tissue positions from 10x Visium spatial folder."""
    pos_file = os.path.join(spatial_dir, 'tissue_positions_list.csv')
    if not os.path.exists(pos_file):
        pos_file = os.path.join(spatial_dir, 'tissue_positions.csv')

    if not os.path.exists(pos_file):
        raise FileNotFoundError(f"tissue_positions file not found in {spatial_dir}")

    df = pd.read_csv(pos_file, header=None,
                     names=['barcode', 'in_tissue', 'row', 'col', 'y', 'x'])
    # Keep only spots in tissue
    df = df[df['in_tissue'] == 1].copy()
    df = df[['barcode', 'x', 'y']].set_index('barcode')
    return df


def load_visium_h5(h5_file: str) -> pd.DataFrame:
    """Load filtered feature-barcode matrix from H5 file."""
    try:
        import scanpy as sc
        adata = sc.read_10x_h5(h5_file)
        return pd.DataFrame(
            adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X,
            index=adata.obs_names,
            columns=adata.var_names
        )
    except ImportError:
        raise ImportError("scanpy required: pip install scanpy")


def simulate_apa_from_expression(expr_df: pd.DataFrame,
                                  dropout_rate: float = 0.3,
                                  seed: int = 42) -> pd.DataFrame:
    """
    Simulate APA indices from gene expression (placeholder until scAPAtrap runs).

    In real usage, replace this with actual scAPAtrap output.
    APA index ≈ normalised expression (proxy for distal site usage).
    """
    np.random.seed(seed)
    # Normalise each gene to [0,1]
    expr = expr_df.values.astype(float)
    col_min = expr.min(axis=0)
    col_max = expr.max(axis=0)
    denom = col_max - col_min
    denom[denom == 0] = 1
    apa = (expr - col_min) / denom

    # Add dropout
    mask = np.random.rand(*apa.shape) < dropout_rate
    apa[mask] = np.nan

    return pd.DataFrame(apa, index=expr_df.index, columns=expr_df.columns)


def prepare_dataset(data_dir: str, output_dir: str, name: str,
                    dropout_rate: float = 0.3):
    """Prepare one dataset."""
    print(f"\nPreparing {name}...")
    os.makedirs(output_dir, exist_ok=True)

    h5_file      = os.path.join(data_dir, 'filtered_feature_bc_matrix.h5')
    spatial_dir  = os.path.join(data_dir, 'spatial')

    if not os.path.exists(h5_file):
        print(f"  {h5_file} not found — skipping {name}")
        return False

    # Load expression
    print("  Loading expression matrix...")
    expr_df = load_visium_h5(h5_file)

    # Load coordinates
    print("  Loading spatial coordinates...")
    coord_df = load_visium_coordinates(spatial_dir)

    # Align
    common = expr_df.index.intersection(coord_df.index)
    expr_df  = expr_df.loc[common]
    coord_df = coord_df.loc[common]

    # Select highly variable genes (top 500 by variance)
    variances = expr_df.var(axis=0)
    top_genes = variances.nlargest(500).index
    expr_df = expr_df[top_genes]

    # Simulate APA (replace with scAPAtrap output in real usage)
    print("  Generating APA matrix (simulated from expression)...")
    apa_df = simulate_apa_from_expression(expr_df, dropout_rate=dropout_rate)

    # Save
    apa_path   = os.path.join(output_dir, 'apa_matrix.csv')
    coord_path = os.path.join(output_dir, 'coordinates.csv')

    apa_df.to_csv(apa_path)
    coord_df.to_csv(coord_path)

    print(f"  Saved: {apa_path}  ({apa_df.shape[0]} spots × {apa_df.shape[1]} genes)")
    print(f"  Saved: {coord_path}")
    return True


def main():
    parser = argparse.ArgumentParser(description='Prepare spaGAPA benchmark data')
    parser.add_argument('--data_dir',   default='data/real',      help='Raw data directory')
    parser.add_argument('--output_dir', default='data/processed',  help='Output directory')
    parser.add_argument('--dropout',    type=float, default=0.3,   help='Dropout rate')
    args = parser.parse_args()

    datasets = {
        'mob':   os.path.join(args.data_dir, 'mob'),
        'brain': os.path.join(args.data_dir, 'brain'),
        'embryo': os.path.join(args.data_dir, 'embryo'),
    }

    for name, data_dir in datasets.items():
        out = os.path.join(args.output_dir, name)
        prepare_dataset(data_dir, out, name, args.dropout)

    print("\nDone. Run: python scripts/run_real_benchmark.py")


if __name__ == '__main__':
    main()
