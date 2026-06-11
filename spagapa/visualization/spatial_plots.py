"""
Spatial visualization for APA data.

This module provides functions for creating publication-quality spatial plots
of APA usage patterns.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize, ListedColormap
from typing import Optional, Union, List, Tuple, Dict
import warnings


class SpatialPlotter:
    """
    Create spatial plots for APA data.
    
    This class provides methods to visualize APA usage patterns in spatial
    transcriptomics data.
    
    Parameters
    ----------
    figsize : tuple, default=(8, 6)
        Figure size (width, height) in inches
    dpi : int, default=300
        Resolution for saved figures
    style : str, default='default'
        Matplotlib style to use
        
    Attributes
    ----------
    figsize : tuple
        Default figure size
    dpi : int
        Default DPI for saving
    """
    
    def __init__(
        self,
        figsize: Tuple[float, float] = (8, 6),
        dpi: int = 300,
        style: str = 'default'
    ):
        self.figsize = figsize
        self.dpi = dpi
        self.style = style
        
        # Set style
        if style != 'default':
            plt.style.use(style)
            
    def plot_spatial_apa(
        self,
        coords: np.ndarray,
        apa_values: np.ndarray,
        gene_name: str = "Gene",
        cmap: str = 'viridis',
        size: float = 10,
        alpha: float = 0.8,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        title: Optional[str] = None,
        xlabel: str = "X coordinate",
        ylabel: str = "Y coordinate",
        colorbar_label: str = "APA Index",
        ax: Optional[plt.Axes] = None,
        save: Optional[str] = None
    ) -> plt.Axes:
        """
        Plot spatial distribution of APA usage for a single gene.
        
        Parameters
        ----------
        coords : np.ndarray, shape (n_spots, 2)
            Spatial coordinates (x, y)
        apa_values : np.ndarray, shape (n_spots,)
            APA values for each spot
        gene_name : str, default="Gene"
            Gene name for title
        cmap : str, default='viridis'
            Colormap name
        size : float, default=10
            Point size
        alpha : float, default=0.8
            Point transparency
        vmin, vmax : float, optional
            Color scale limits
        title : str, optional
            Custom title
        xlabel, ylabel : str
            Axis labels
        colorbar_label : str
            Colorbar label
        ax : plt.Axes, optional
            Existing axes to plot on
        save : str, optional
            Path to save figure
            
        Returns
        -------
        ax : plt.Axes
            Matplotlib axes object
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=self.figsize)
        else:
            fig = ax.figure
            
        # Handle NaN values
        valid_mask = ~np.isnan(apa_values)
        coords_valid = coords[valid_mask]
        values_valid = apa_values[valid_mask]
        
        # Plot
        scatter = ax.scatter(
            coords_valid[:, 0],
            coords_valid[:, 1],
            c=values_valid,
            cmap=cmap,
            s=size,
            alpha=alpha,
            vmin=vmin,
            vmax=vmax,
            edgecolors='none'
        )
        
        # Colorbar
        cbar = plt.colorbar(scatter, ax=ax, label=colorbar_label)
        
        # Labels and title
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if title is None:
            title = f"{gene_name} APA Usage"
        ax.set_title(title)
        
        # Equal aspect ratio
        ax.set_aspect('equal', adjustable='box')
        
        # Save if requested
        if save:
            fig.savefig(save, dpi=self.dpi, bbox_inches='tight')
            
        return ax
        
    def plot_spatial_domains(
        self,
        coords: np.ndarray,
        domain_labels: np.ndarray,
        domain_names: Optional[Dict[int, str]] = None,
        cmap: str = 'tab10',
        size: float = 10,
        alpha: float = 0.8,
        show_boundaries: bool = False,
        title: str = "Spatial Domains",
        xlabel: str = "X coordinate",
        ylabel: str = "Y coordinate",
        ax: Optional[plt.Axes] = None,
        save: Optional[str] = None
    ) -> plt.Axes:
        """
        Plot spatial domains with different colors.
        
        Parameters
        ----------
        coords : np.ndarray, shape (n_spots, 2)
            Spatial coordinates
        domain_labels : np.ndarray, shape (n_spots,)
            Domain label for each spot
        domain_names : dict, optional
            Mapping from domain ID to name
        cmap : str, default='tab10'
            Colormap for domains
        size : float, default=10
            Point size
        alpha : float, default=0.8
            Point transparency
        show_boundaries : bool, default=False
            Whether to draw domain boundaries
        title : str
            Plot title
        xlabel, ylabel : str
            Axis labels
        ax : plt.Axes, optional
            Existing axes
        save : str, optional
            Path to save figure
            
        Returns
        -------
        ax : plt.Axes
            Matplotlib axes object
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=self.figsize)
        else:
            fig = ax.figure
            
        # Get unique domains
        unique_domains = np.unique(domain_labels)
        n_domains = len(unique_domains)
        
        # Get colormap
        if isinstance(cmap, str):
            cmap_obj = plt.get_cmap(cmap, n_domains)
            colors = [cmap_obj(i) for i in range(n_domains)]
        else:
            colors = cmap
            
        # Plot each domain
        for i, domain in enumerate(unique_domains):
            mask = domain_labels == domain
            ax.scatter(
                coords[mask, 0],
                coords[mask, 1],
                c=[colors[i]],
                s=size,
                alpha=alpha,
                label=domain_names.get(domain, f"Domain {domain}") if domain_names else f"Domain {domain}",
                edgecolors='none'
            )
            
        # Draw boundaries if requested
        if show_boundaries:
            self._draw_domain_boundaries(ax, coords, domain_labels)
            
        # Legend
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
        
        # Labels
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_aspect('equal', adjustable='box')
        
        # Save if requested
        if save:
            fig.savefig(save, dpi=self.dpi, bbox_inches='tight')
            
        return ax
        
    def plot_apa_comparison(
        self,
        coords: np.ndarray,
        apa_matrix: np.ndarray,
        genes: List[str],
        gene_names: Optional[List[str]] = None,
        ncols: int = 3,
        cmap: str = 'viridis',
        size: float = 5,
        alpha: float = 0.8,
        share_colorbar: bool = False,
        suptitle: Optional[str] = None,
        save: Optional[str] = None
    ) -> plt.Figure:
        """
        Plot APA usage for multiple genes in a grid.
        
        Parameters
        ----------
        coords : np.ndarray, shape (n_spots, 2)
            Spatial coordinates
        apa_matrix : np.ndarray, shape (n_genes, n_spots)
            APA matrix
        genes : list of int or str
            Gene indices or names to plot
        gene_names : list of str, optional
            Gene names for titles
        ncols : int, default=3
            Number of columns in grid
        cmap : str, default='viridis'
            Colormap
        size : float, default=5
            Point size
        alpha : float, default=0.8
            Point transparency
        share_colorbar : bool, default=False
            Use same color scale for all plots
        suptitle : str, optional
            Super title for figure
        save : str, optional
            Path to save figure
            
        Returns
        -------
        fig : plt.Figure
            Matplotlib figure object
        """
        n_genes = len(genes)
        nrows = int(np.ceil(n_genes / ncols))
        
        fig, axes = plt.subplots(
            nrows, ncols,
            figsize=(ncols * 4, nrows * 3.5)
        )
        axes = np.atleast_2d(axes).flatten()
        
        # Determine color scale
        if share_colorbar:
            all_values = []
            for gene_idx in genes:
                if isinstance(gene_idx, str) and gene_names:
                    gene_idx = gene_names.index(gene_idx)
                values = apa_matrix[gene_idx, :]
                all_values.extend(values[~np.isnan(values)])
            vmin, vmax = np.min(all_values), np.max(all_values)
        else:
            vmin, vmax = None, None
            
        # Plot each gene
        for i, gene_idx in enumerate(genes):
            if isinstance(gene_idx, str) and gene_names:
                gene_idx = gene_names.index(gene_idx)
                
            gene_name = gene_names[gene_idx] if gene_names else f"Gene {gene_idx}"
            
            self.plot_spatial_apa(
                coords,
                apa_matrix[gene_idx, :],
                gene_name=gene_name,
                cmap=cmap,
                size=size,
                alpha=alpha,
                vmin=vmin,
                vmax=vmax,
                ax=axes[i]
            )
            
        # Hide unused axes
        for i in range(n_genes, len(axes)):
            axes[i].axis('off')
            
        # Super title
        if suptitle:
            fig.suptitle(suptitle, fontsize=14, y=1.02)
            
        plt.tight_layout()
        
        # Save if requested
        if save:
            fig.savefig(save, dpi=self.dpi, bbox_inches='tight')
            
        return fig
        
    def plot_differential_spatial(
        self,
        coords: np.ndarray,
        apa_values: np.ndarray,
        pvalues: np.ndarray,
        gene_name: str = "Gene",
        threshold: float = 0.05,
        cmap: str = 'RdBu_r',
        size: float = 10,
        alpha: float = 0.8,
        title: Optional[str] = None,
        save: Optional[str] = None
    ) -> plt.Axes:
        """
        Plot spatial distribution with significance highlighting.
        
        Parameters
        ----------
        coords : np.ndarray, shape (n_spots, 2)
            Spatial coordinates
        apa_values : np.ndarray, shape (n_spots,)
            APA values
        pvalues : np.ndarray, shape (n_spots,)
            P-values for each spot
        gene_name : str
            Gene name
        threshold : float, default=0.05
            Significance threshold
        cmap : str, default='RdBu_r'
            Colormap
        size : float, default=10
            Point size
        alpha : float, default=0.8
            Point transparency
        title : str, optional
            Custom title
        save : str, optional
            Path to save figure
            
        Returns
        -------
        ax : plt.Axes
            Matplotlib axes object
        """
        fig, ax = plt.subplots(figsize=self.figsize)
        
        # Separate significant and non-significant
        sig_mask = pvalues < threshold
        
        # Plot non-significant (gray)
        ax.scatter(
            coords[~sig_mask, 0],
            coords[~sig_mask, 1],
            c='lightgray',
            s=size,
            alpha=alpha * 0.5,
            label='Not significant',
            edgecolors='none'
        )
        
        # Plot significant (colored by APA)
        if sig_mask.sum() > 0:
            scatter = ax.scatter(
                coords[sig_mask, 0],
                coords[sig_mask, 1],
                c=apa_values[sig_mask],
                cmap=cmap,
                s=size * 1.5,
                alpha=alpha,
                label=f'Significant (p < {threshold})',
                edgecolors='black',
                linewidths=0.5
            )
            plt.colorbar(scatter, ax=ax, label='APA Index')
            
        ax.legend()
        ax.set_xlabel('X coordinate')
        ax.set_ylabel('Y coordinate')
        if title is None:
            title = f"{gene_name} - Differential APA"
        ax.set_title(title)
        ax.set_aspect('equal', adjustable='box')
        
        if save:
            fig.savefig(save, dpi=self.dpi, bbox_inches='tight')
            
        return ax
        
    def _draw_domain_boundaries(
        self,
        ax: plt.Axes,
        coords: np.ndarray,
        labels: np.ndarray
    ):
        """Draw boundaries between domains (simplified)."""
        # This is a simplified version
        # For production, consider using alpha shapes or contours
        from scipy.spatial import Delaunay
        
        try:
            tri = Delaunay(coords)
            
            # Find edges between different domains
            for simplex in tri.simplices:
                for i in range(3):
                    p1, p2 = simplex[i], simplex[(i+1) % 3]
                    if labels[p1] != labels[p2]:
                        ax.plot(
                            [coords[p1, 0], coords[p2, 0]],
                            [coords[p1, 1], coords[p2, 1]],
                            'k-', linewidth=0.5, alpha=0.3
                        )
        except:
            warnings.warn("Could not draw domain boundaries")


def plot_spatial_apa(
    coords: np.ndarray,
    apa_values: np.ndarray,
    gene_name: str = "Gene",
    cmap: str = 'viridis',
    size: float = 10,
    figsize: Tuple[float, float] = (8, 6),
    save: Optional[str] = None,
    **kwargs
) -> plt.Axes:
    """
    Plot spatial APA distribution (convenience function).
    
    Parameters
    ----------
    coords : np.ndarray, shape (n_spots, 2)
        Spatial coordinates
    apa_values : np.ndarray, shape (n_spots,)
        APA values
    gene_name : str
        Gene name
    cmap : str, default='viridis'
        Colormap
    size : float, default=10
        Point size
    figsize : tuple, default=(8, 6)
        Figure size
    save : str, optional
        Path to save figure
    **kwargs
        Additional arguments passed to SpatialPlotter.plot_spatial_apa
        
    Returns
    -------
    ax : plt.Axes
        Matplotlib axes object
        
    Examples
    --------
    >>> import numpy as np
    >>> from spagapa.visualization import plot_spatial_apa
    >>> 
    >>> coords = np.random.rand(100, 2) * 100
    >>> apa_values = np.random.rand(100)
    >>> 
    >>> ax = plot_spatial_apa(coords, apa_values, gene_name="GENE1")
    >>> plt.show()
    """
    plotter = SpatialPlotter(figsize=figsize)
    return plotter.plot_spatial_apa(
        coords, apa_values, gene_name, cmap, size, save=save, **kwargs
    )


def plot_spatial_domains(
    coords: np.ndarray,
    domain_labels: np.ndarray,
    domain_names: Optional[Dict[int, str]] = None,
    cmap: str = 'tab10',
    size: float = 10,
    figsize: Tuple[float, float] = (8, 6),
    save: Optional[str] = None,
    **kwargs
) -> plt.Axes:
    """
    Plot spatial domains (convenience function).
    
    Parameters
    ----------
    coords : np.ndarray, shape (n_spots, 2)
        Spatial coordinates
    domain_labels : np.ndarray, shape (n_spots,)
        Domain labels
    domain_names : dict, optional
        Domain names
    cmap : str, default='tab10'
        Colormap
    size : float, default=10
        Point size
    figsize : tuple, default=(8, 6)
        Figure size
    save : str, optional
        Path to save figure
    **kwargs
        Additional arguments
        
    Returns
    -------
    ax : plt.Axes
        Matplotlib axes object
        
    Examples
    --------
    >>> import numpy as np
    >>> from spagapa.visualization import plot_spatial_domains
    >>> 
    >>> coords = np.random.rand(100, 2) * 100
    >>> domains = np.repeat([0, 1, 2, 3], 25)
    >>> 
    >>> ax = plot_spatial_domains(coords, domains)
    >>> plt.show()
    """
    plotter = SpatialPlotter(figsize=figsize)
    return plotter.plot_spatial_domains(
        coords, domain_labels, domain_names, cmap, size, save=save, **kwargs
    )
