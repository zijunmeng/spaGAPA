"""
Statistical plots for APA analysis.

This module provides functions for creating statistical visualizations
such as volcano plots, heatmaps, and box plots.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, Union, List, Tuple, Dict
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist


class StatisticalPlotter:
    """
    Create statistical plots for APA analysis.
    
    Parameters
    ----------
    figsize : tuple, default=(8, 6)
        Default figure size
    dpi : int, default=300
        Resolution for saved figures
    style : str, default='default'
        Matplotlib style
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
        
        if style != 'default':
            plt.style.use(style)
            
    def plot_volcano(
        self,
        logfc: np.ndarray,
        pvalues: np.ndarray,
        gene_names: Optional[List[str]] = None,
        threshold_fc: float = 0.5,
        threshold_p: float = 0.05,
        label_top: int = 10,
        title: str = "Volcano Plot",
        xlabel: str = "log2 Fold Change",
        ylabel: str = "-log10(p-value)",
        ax: Optional[plt.Axes] = None,
        save: Optional[str] = None
    ) -> plt.Axes:
        """
        Create volcano plot for differential APA analysis.
        
        Parameters
        ----------
        logfc : np.ndarray
            Log2 fold changes
        pvalues : np.ndarray
            P-values
        gene_names : list of str, optional
            Gene names for labeling
        threshold_fc : float, default=0.5
            Log2FC threshold for significance
        threshold_p : float, default=0.05
            P-value threshold
        label_top : int, default=10
            Number of top genes to label
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
            
        # Calculate -log10(p-value)
        neg_log_p = -np.log10(pvalues + 1e-300)  # Add small value to avoid log(0)
        
        # Classify genes
        sig_up = (logfc > threshold_fc) & (pvalues < threshold_p)
        sig_down = (logfc < -threshold_fc) & (pvalues < threshold_p)
        not_sig = ~(sig_up | sig_down)
        
        # Plot
        ax.scatter(
            logfc[not_sig], neg_log_p[not_sig],
            c='gray', s=10, alpha=0.5, label='Not significant'
        )
        ax.scatter(
            logfc[sig_down], neg_log_p[sig_down],
            c='blue', s=20, alpha=0.7, label='Down-regulated'
        )
        ax.scatter(
            logfc[sig_up], neg_log_p[sig_up],
            c='red', s=20, alpha=0.7, label='Up-regulated'
        )
        
        # Threshold lines
        ax.axhline(-np.log10(threshold_p), color='black', linestyle='--', 
                  linewidth=1, alpha=0.5)
        ax.axvline(threshold_fc, color='black', linestyle='--', 
                  linewidth=1, alpha=0.5)
        ax.axvline(-threshold_fc, color='black', linestyle='--', 
                  linewidth=1, alpha=0.5)
        
        # Label top genes
        if gene_names and label_top > 0:
            # Get top genes by p-value
            top_indices = np.argsort(pvalues)[:label_top]
            for idx in top_indices:
                if sig_up[idx] or sig_down[idx]:
                    ax.annotate(
                        gene_names[idx],
                        (logfc[idx], neg_log_p[idx]),
                        xytext=(5, 5),
                        textcoords='offset points',
                        fontsize=8,
                        alpha=0.7
                    )
        
        # Labels and legend
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(frameon=False)
        ax.grid(True, alpha=0.3)
        
        # Add counts
        n_up = sig_up.sum()
        n_down = sig_down.sum()
        ax.text(
            0.02, 0.98,
            f"Up: {n_up}\nDown: {n_down}",
            transform=ax.transAxes,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        )
        
        if save:
            fig.savefig(save, dpi=self.dpi, bbox_inches='tight')
            
        return ax
        
    def plot_heatmap(
        self,
        matrix: np.ndarray,
        row_labels: Optional[List[str]] = None,
        col_labels: Optional[List[str]] = None,
        cluster_rows: bool = True,
        cluster_cols: bool = True,
        cmap: str = 'RdBu_r',
        center: Optional[float] = None,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        title: str = "Heatmap",
        xlabel: str = "",
        ylabel: str = "",
        cbar_label: str = "Value",
        figsize: Optional[Tuple[float, float]] = None,
        save: Optional[str] = None
    ) -> Tuple[plt.Figure, plt.Axes]:
        """
        Create clustered heatmap.
        
        Parameters
        ----------
        matrix : np.ndarray, shape (n_rows, n_cols)
            Data matrix
        row_labels : list of str, optional
            Row labels
        col_labels : list of str, optional
            Column labels
        cluster_rows : bool, default=True
            Cluster rows
        cluster_cols : bool, default=True
            Cluster columns
        cmap : str, default='RdBu_r'
            Colormap
        center : float, optional
            Value to center colormap
        vmin, vmax : float, optional
            Color scale limits
        title : str
            Plot title
        xlabel, ylabel : str
            Axis labels
        cbar_label : str
            Colorbar label
        figsize : tuple, optional
            Figure size
        save : str, optional
            Path to save figure
            
        Returns
        -------
        fig : plt.Figure
            Matplotlib figure
        ax : plt.Axes
            Matplotlib axes
        """
        if figsize is None:
            figsize = (max(8, matrix.shape[1] * 0.3), 
                      max(6, matrix.shape[0] * 0.2))
            
        # Clustering
        row_order = np.arange(matrix.shape[0])
        col_order = np.arange(matrix.shape[1])
        
        if cluster_rows and matrix.shape[0] > 1:
            row_linkage = hierarchy.linkage(
                pdist(matrix, metric='euclidean'),
                method='average'
            )
            row_order = hierarchy.leaves_list(row_linkage)
            
        if cluster_cols and matrix.shape[1] > 1:
            col_linkage = hierarchy.linkage(
                pdist(matrix.T, metric='euclidean'),
                method='average'
            )
            col_order = hierarchy.leaves_list(col_linkage)
            
        # Reorder matrix
        matrix_ordered = matrix[row_order, :][:, col_order]
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Plot heatmap
        im = ax.imshow(
            matrix_ordered,
            cmap=cmap,
            aspect='auto',
            vmin=vmin,
            vmax=vmax,
            interpolation='nearest'
        )
        
        # Colorbar
        cbar = plt.colorbar(im, ax=ax, label=cbar_label)
        
        # Labels
        if row_labels:
            row_labels_ordered = [row_labels[i] for i in row_order]
            ax.set_yticks(np.arange(len(row_labels_ordered)))
            ax.set_yticklabels(row_labels_ordered, fontsize=8)
        else:
            ax.set_yticks([])
            
        if col_labels:
            col_labels_ordered = [col_labels[i] for i in col_order]
            ax.set_xticks(np.arange(len(col_labels_ordered)))
            ax.set_xticklabels(col_labels_ordered, rotation=90, 
                             ha='right', fontsize=8)
        else:
            ax.set_xticks([])
            
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        
        plt.tight_layout()
        
        if save:
            fig.savefig(save, dpi=self.dpi, bbox_inches='tight')
            
        return fig, ax
        
    def plot_boxplot(
        self,
        apa_values: Union[np.ndarray, List[np.ndarray]],
        groups: Union[np.ndarray, List[str]],
        gene_name: str = "Gene",
        ylabel: str = "APA Index",
        title: Optional[str] = None,
        colors: Optional[List[str]] = None,
        show_points: bool = True,
        ax: Optional[plt.Axes] = None,
        save: Optional[str] = None
    ) -> plt.Axes:
        """
        Create box plot of APA values across groups.
        
        Parameters
        ----------
        apa_values : np.ndarray or list of arrays
            APA values (either single array with groups, or list of arrays)
        groups : np.ndarray or list of str
            Group labels
        gene_name : str
            Gene name for title
        ylabel : str
            Y-axis label
        title : str, optional
            Custom title
        colors : list of str, optional
            Colors for each group
        show_points : bool, default=True
            Show individual points
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
            
        # Prepare data
        if isinstance(apa_values, np.ndarray):
            # Single array with group labels
            unique_groups = np.unique(groups)
            data = [apa_values[groups == g] for g in unique_groups]
            labels = [str(g) for g in unique_groups]
        else:
            # List of arrays
            data = apa_values
            labels = groups if isinstance(groups, list) else [str(g) for g in groups]
            
        # Box plot
        bp = ax.boxplot(
            data,
            labels=labels,
            patch_artist=True,
            showfliers=False
        )
        
        # Colors
        if colors is None:
            colors = plt.cm.Set3(np.linspace(0, 1, len(data)))
            
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
            
        # Show points
        if show_points:
            for i, d in enumerate(data):
                x = np.random.normal(i + 1, 0.04, size=len(d))
                ax.scatter(x, d, alpha=0.3, s=10, color='black')
                
        # Labels
        ax.set_ylabel(ylabel)
        ax.set_xlabel("Group")
        if title is None:
            title = f"{gene_name} APA Distribution"
        ax.set_title(title)
        ax.grid(True, alpha=0.3, axis='y')
        
        if save:
            fig.savefig(save, dpi=self.dpi, bbox_inches='tight')
            
        return ax
        
    def plot_violin(
        self,
        apa_values: Union[np.ndarray, List[np.ndarray]],
        groups: Union[np.ndarray, List[str]],
        gene_name: str = "Gene",
        ylabel: str = "APA Index",
        title: Optional[str] = None,
        colors: Optional[List[str]] = None,
        ax: Optional[plt.Axes] = None,
        save: Optional[str] = None
    ) -> plt.Axes:
        """
        Create violin plot of APA values across groups.
        
        Parameters
        ----------
        apa_values : np.ndarray or list of arrays
            APA values
        groups : np.ndarray or list of str
            Group labels
        gene_name : str
            Gene name
        ylabel : str
            Y-axis label
        title : str, optional
            Custom title
        colors : list of str, optional
            Colors for each group
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
            
        # Prepare data for seaborn
        if isinstance(apa_values, np.ndarray):
            df = pd.DataFrame({
                'APA': apa_values,
                'Group': groups
            })
        else:
            # List of arrays
            df_list = []
            for i, (values, group) in enumerate(zip(apa_values, groups)):
                df_list.append(pd.DataFrame({
                    'APA': values,
                    'Group': group
                }))
            df = pd.concat(df_list, ignore_index=True)
            
        # Violin plot
        sns.violinplot(
            data=df,
            x='Group',
            y='APA',
            palette=colors,
            ax=ax
        )
        
        # Labels
        ax.set_ylabel(ylabel)
        ax.set_xlabel("Group")
        if title is None:
            title = f"{gene_name} APA Distribution"
        ax.set_title(title)
        ax.grid(True, alpha=0.3, axis='y')
        
        if save:
            fig.savefig(save, dpi=self.dpi, bbox_inches='tight')
            
        return ax


def plot_volcano(
    logfc: np.ndarray,
    pvalues: np.ndarray,
    gene_names: Optional[List[str]] = None,
    threshold_fc: float = 0.5,
    threshold_p: float = 0.05,
    figsize: Tuple[float, float] = (8, 6),
    save: Optional[str] = None,
    **kwargs
) -> plt.Axes:
    """
    Create volcano plot (convenience function).
    
    Parameters
    ----------
    logfc : np.ndarray
        Log2 fold changes
    pvalues : np.ndarray
        P-values
    gene_names : list of str, optional
        Gene names
    threshold_fc : float, default=0.5
        Log2FC threshold
    threshold_p : float, default=0.05
        P-value threshold
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
    >>> from spagapa.visualization import plot_volcano
    >>> 
    >>> logfc = np.random.randn(100)
    >>> pvalues = np.random.rand(100)
    >>> 
    >>> ax = plot_volcano(logfc, pvalues)
    >>> plt.show()
    """
    plotter = StatisticalPlotter(figsize=figsize)
    return plotter.plot_volcano(
        logfc, pvalues, gene_names, threshold_fc, threshold_p, 
        save=save, **kwargs
    )


def plot_heatmap(
    matrix: np.ndarray,
    row_labels: Optional[List[str]] = None,
    col_labels: Optional[List[str]] = None,
    cluster_rows: bool = True,
    cluster_cols: bool = True,
    figsize: Optional[Tuple[float, float]] = None,
    save: Optional[str] = None,
    **kwargs
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Create clustered heatmap (convenience function).
    
    Parameters
    ----------
    matrix : np.ndarray
        Data matrix
    row_labels : list of str, optional
        Row labels
    col_labels : list of str, optional
        Column labels
    cluster_rows : bool, default=True
        Cluster rows
    cluster_cols : bool, default=True
        Cluster columns
    figsize : tuple, optional
        Figure size
    save : str, optional
        Path to save figure
    **kwargs
        Additional arguments
        
    Returns
    -------
    fig : plt.Figure
        Matplotlib figure
    ax : plt.Axes
        Matplotlib axes
        
    Examples
    --------
    >>> import numpy as np
    >>> from spagapa.visualization import plot_heatmap
    >>> 
    >>> matrix = np.random.rand(20, 10)
    >>> fig, ax = plot_heatmap(matrix)
    >>> plt.show()
    """
    plotter = StatisticalPlotter(figsize=figsize or (8, 6))
    return plotter.plot_heatmap(
        matrix, row_labels, col_labels, cluster_rows, cluster_cols,
        figsize=figsize, save=save, **kwargs
    )
