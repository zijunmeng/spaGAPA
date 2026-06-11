"""
Quality control plots for APA analysis.

This module provides functions for visualizing quality control metrics.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Optional, Tuple, List
import warnings


class QCPlotter:
    """Create quality control plots."""
    
    def __init__(self, figsize=(8, 6), dpi=300):
        self.figsize = figsize
        self.dpi = dpi
        
    def plot_imputation_quality(
        self,
        observed: np.ndarray,
        imputed: np.ndarray,
        uncertainties: Optional[np.ndarray] = None,
        title: str = "Imputation Quality",
        save: Optional[str] = None
    ) -> plt.Figure:
        """Plot imputation quality metrics."""
        if uncertainties is not None:
            fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        else:
            fig, axes = plt.subplots(1, 2, figsize=(10, 4))
            axes = list(axes) + [None]
            
        # 1. Observed vs Imputed
        ax = axes[0]
        ax.scatter(observed, imputed, alpha=0.5, s=10)
        lims = [min(observed.min(), imputed.min()), 
                max(observed.max(), imputed.max())]
        ax.plot(lims, lims, 'r--', alpha=0.5, label='Perfect')
        ax.set_xlabel('Observed')
        ax.set_ylabel('Imputed')
        ax.set_title('Observed vs Imputed')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # 2. Residuals
        ax = axes[1]
        residuals = imputed - observed
        ax.hist(residuals, bins=50, alpha=0.7, edgecolor='black')
        ax.axvline(0, color='red', linestyle='--', label='Zero')
        ax.set_xlabel('Residual (Imputed - Observed)')
        ax.set_ylabel('Frequency')
        ax.set_title(f'Residuals (MAE={np.abs(residuals).mean():.3f})')
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        
        # 3. Uncertainty distribution
        if uncertainties is not None and axes[2] is not None:
            ax = axes[2]
            ax.hist(uncertainties, bins=50, alpha=0.7, 
                   edgecolor='black', color='orange')
            ax.set_xlabel('Uncertainty (Std Dev)')
            ax.set_ylabel('Frequency')
            ax.set_title(f'Uncertainty (Mean={uncertainties.mean():.3f})')
            ax.grid(True, alpha=0.3, axis='y')
            
        plt.suptitle(title, fontsize=14, y=1.02)
        plt.tight_layout()
        
        if save:
            fig.savefig(save, dpi=self.dpi, bbox_inches='tight')
            
        return fig
        
    def plot_spatial_support(
        self,
        coords: np.ndarray,
        support_scores: np.ndarray,
        threshold: float = 0.3,
        title: str = "Spatial Support",
        save: Optional[str] = None
    ) -> plt.Axes:
        """Plot spatial support scores."""
        fig, ax = plt.subplots(figsize=self.figsize)
        
        # Color by support score
        scatter = ax.scatter(
            coords[:, 0], coords[:, 1],
            c=support_scores,
            cmap='RdYlGn',
            s=20, alpha=0.7,
            vmin=0, vmax=1
        )
        
        # Threshold line in colorbar
        cbar = plt.colorbar(scatter, ax=ax, label='Support Score')
        
        ax.set_xlabel('X coordinate')
        ax.set_ylabel('Y coordinate')
        ax.set_title(title)
        ax.set_aspect('equal', adjustable='box')
        
        if save:
            fig.savefig(save, dpi=self.dpi, bbox_inches='tight')
            
        return ax
        
    def plot_dropout_stats(
        self,
        apa_matrix: np.ndarray,
        coords: np.ndarray,
        title: str = "Dropout Statistics",
        save: Optional[str] = None
    ) -> plt.Figure:
        """Plot dropout statistics."""
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        
        # 1. Dropout rate per spot
        ax = axes[0]
        dropout_per_spot = np.isnan(apa_matrix).mean(axis=0)
        scatter = ax.scatter(
            coords[:, 0], coords[:, 1],
            c=dropout_per_spot * 100,
            cmap='Reds',
            s=20, alpha=0.7
        )
        plt.colorbar(scatter, ax=ax, label='Dropout Rate (%)')
        ax.set_xlabel('X coordinate')
        ax.set_ylabel('Y coordinate')
        ax.set_title('Dropout Rate per Spot')
        ax.set_aspect('equal', adjustable='box')
        
        # 2. Dropout rate per gene
        ax = axes[1]
        dropout_per_gene = np.isnan(apa_matrix).mean(axis=1)
        ax.hist(dropout_per_gene * 100, bins=50, alpha=0.7, edgecolor='black')
        ax.axvline(dropout_per_gene.mean() * 100, color='red', 
                  linestyle='--', label=f'Mean={dropout_per_gene.mean()*100:.1f}%')
        ax.set_xlabel('Dropout Rate (%)')
        ax.set_ylabel('Number of Genes')
        ax.set_title('Dropout Rate Distribution')
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.suptitle(title, fontsize=14, y=1.02)
        plt.tight_layout()
        
        if save:
            fig.savefig(save, dpi=self.dpi, bbox_inches='tight')
            
        return fig
        
    def generate_qc_report(
        self,
        dataset,
        output_dir: str = "."
    ):
        """Generate comprehensive QC report."""
        print("Generating QC report...")
        
        # Plot 1: Imputation quality
        if hasattr(dataset, 'imputed') and dataset.imputed is not None:
            observed = dataset.apa_matrix[~np.isnan(dataset.apa_matrix)]
            imputed_vals = dataset.imputed[~np.isnan(dataset.apa_matrix)]
            uncertainties = dataset.uncertainties[~np.isnan(dataset.apa_matrix)] if hasattr(dataset, 'uncertainties') else None
            
            self.plot_imputation_quality(
                observed, imputed_vals, uncertainties,
                save=f"{output_dir}/qc_imputation.png"
            )
            
        # Plot 2: Dropout statistics
        self.plot_dropout_stats(
            dataset.apa_matrix,
            dataset.spatial_coords,
            save=f"{output_dir}/qc_dropout.png"
        )
        
        print(f"QC report saved to {output_dir}/")


def plot_imputation_quality(
    observed: np.ndarray,
    imputed: np.ndarray,
    uncertainties: Optional[np.ndarray] = None,
    figsize: Tuple[float, float] = (15, 4),
    save: Optional[str] = None
) -> plt.Figure:
    """Plot imputation quality (convenience function)."""
    plotter = QCPlotter(figsize=figsize)
    return plotter.plot_imputation_quality(
        observed, imputed, uncertainties, save=save
    )
