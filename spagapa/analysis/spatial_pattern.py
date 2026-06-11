"""
Spatial pattern discovery for APA analysis.

This module provides methods to identify spatially variable APA genes and
discover spatial expression patterns.
"""

import numpy as np
import pandas as pd
from typing import Optional, Union, Tuple, List, Dict
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import AgglomerativeClustering
from sklearn.preprocessing import StandardScaler


class SpatialPatternAnalyzer:
    """
    Analyze spatial patterns in APA data.
    
    This class provides methods to:
    1. Calculate Moran's I for spatial autocorrelation
    2. Identify spatially variable APA genes (SVAPA)
    3. Cluster genes by spatial patterns
    
    Parameters
    ----------
    n_permutations : int, default=1000
        Number of permutations for significance testing
    alpha : float, default=0.05
        Significance level
    random_state : int, optional
        Random state for reproducibility
        
    Attributes
    ----------
    morans_i_ : pd.DataFrame
        Moran's I statistics for each gene
    svapa_genes_ : list
        List of spatially variable APA genes
    """
    
    def __init__(
        self,
        n_permutations: int = 1000,
        alpha: float = 0.05,
        random_state: Optional[int] = None
    ):
        self.n_permutations = n_permutations
        self.alpha = alpha
        self.random_state = random_state
        
        self.morans_i_ = None
        self.svapa_genes_ = None
        
        if random_state is not None:
            np.random.seed(random_state)
            
    def compute_morans_i(
        self,
        apa_values: np.ndarray,
        spatial_weights: np.ndarray
    ) -> Tuple[float, float, float]:
        """
        Compute Moran's I statistic for spatial autocorrelation.
        
        Moran's I measures spatial autocorrelation:
        - I > 0: positive spatial autocorrelation (similar values cluster)
        - I ≈ 0: random spatial pattern
        - I < 0: negative spatial autocorrelation (dissimilar values cluster)
        
        Parameters
        ----------
        apa_values : np.ndarray, shape (n_spots,)
            APA values for one gene
        spatial_weights : np.ndarray, shape (n_spots, n_spots)
            Spatial weight matrix (row-normalized)
            
        Returns
        -------
        morans_i : float
            Moran's I statistic
        expected_i : float
            Expected value under null hypothesis
        variance_i : float
            Variance under null hypothesis
        """
        # Remove NaN values
        valid_mask = ~np.isnan(apa_values)
        if valid_mask.sum() < 3:
            return np.nan, np.nan, np.nan
            
        values = apa_values[valid_mask]
        weights = spatial_weights[np.ix_(valid_mask, valid_mask)]
        
        n = len(values)
        
        # Standardize values
        mean_val = np.mean(values)
        values_centered = values - mean_val
        
        # Compute Moran's I
        numerator = np.sum(weights * np.outer(values_centered, values_centered))
        denominator = np.sum(values_centered ** 2)
        
        W = np.sum(weights)
        morans_i = (n / W) * (numerator / denominator)
        
        # Expected value and variance
        expected_i = -1.0 / (n - 1)
        
        # Simplified variance (assuming normality)
        S1 = 0.5 * np.sum((weights + weights.T) ** 2)
        S2 = np.sum((np.sum(weights, axis=1) + np.sum(weights, axis=0)) ** 2)
        S3 = (np.sum(values_centered ** 4) / n) / (denominator / n) ** 2
        S4 = (n ** 2 - 3 * n + 3) * S1 - n * S2 + 3 * W ** 2
        S5 = (n ** 2 - n) * S1 - 2 * n * S2 + 6 * W ** 2
        
        variance_i = (n * S4 - S3 * S5) / ((n - 1) * (n - 2) * (n - 3) * W ** 2)
        variance_i -= expected_i ** 2
        
        return morans_i, expected_i, variance_i
        
    def test_spatial_autocorrelation(
        self,
        apa_matrix: np.ndarray,
        spatial_coords: np.ndarray,
        gene_names: Optional[List[str]] = None,
        n_neighbors: int = 6
    ) -> pd.DataFrame:
        """
        Test spatial autocorrelation for all genes.
        
        Parameters
        ----------
        apa_matrix : np.ndarray, shape (n_genes, n_spots)
            APA index matrix
        spatial_coords : np.ndarray, shape (n_spots, 2)
            Spatial coordinates
        gene_names : list of str, optional
            Gene names
        n_neighbors : int, default=6
            Number of neighbors for spatial weights
            
        Returns
        -------
        results : pd.DataFrame
            Results with columns:
            - gene: gene name
            - morans_i: Moran's I statistic
            - expected_i: expected value
            - z_score: standardized statistic
            - pvalue: p-value
        """
        # Build spatial weight matrix
        from ..spatial import build_knn_graph
        distances, indices = build_knn_graph(spatial_coords, k=n_neighbors)
        
        # Convert to sparse adjacency matrix
        from scipy.sparse import csr_matrix
        n_spots = len(spatial_coords)
        row_ind = np.repeat(np.arange(n_spots), indices.shape[1])
        col_ind = indices.flatten()
        data = np.ones(len(row_ind))
        spatial_graph = csr_matrix(
            (data, (row_ind, col_ind)),
            shape=(n_spots, n_spots)
        )
        
        # Row-normalize
        row_sums = np.array(spatial_graph.sum(axis=1)).flatten()
        row_sums[row_sums == 0] = 1
        spatial_weights = spatial_graph.multiply(1.0 / row_sums[:, np.newaxis])
        spatial_weights = spatial_weights.toarray()
        
        n_genes = apa_matrix.shape[0]
        
        if gene_names is None:
            gene_names = [f"Gene_{i}" for i in range(n_genes)]
            
        results_list = []
        
        for i, gene in enumerate(gene_names):
            apa_values = apa_matrix[i, :]
            
            morans_i, expected_i, variance_i = self.compute_morans_i(
                apa_values, spatial_weights
            )
            
            if np.isnan(morans_i):
                continue
                
            # Z-score and p-value
            if variance_i > 0:
                z_score = (morans_i - expected_i) / np.sqrt(variance_i)
                pvalue = 2 * (1 - stats.norm.cdf(np.abs(z_score)))
            else:
                z_score = np.nan
                pvalue = 1.0
                
            results_list.append({
                'gene': gene,
                'morans_i': morans_i,
                'expected_i': expected_i,
                'z_score': z_score,
                'pvalue': pvalue
            })
            
        results = pd.DataFrame(results_list)
        self.morans_i_ = results
        
        return results
        
    def identify_svapa_genes(
        self,
        apa_matrix: np.ndarray,
        spatial_coords: np.ndarray,
        gene_names: Optional[List[str]] = None,
        fdr_threshold: float = 0.05,
        n_neighbors: int = 6
    ) -> List[str]:
        """
        Identify spatially variable APA genes (SVAPA).
        
        SVAPA genes show significant spatial autocorrelation in their
        APA usage patterns.
        
        Parameters
        ----------
        apa_matrix : np.ndarray, shape (n_genes, n_spots)
            APA index matrix
        spatial_coords : np.ndarray, shape (n_spots, 2)
            Spatial coordinates
        gene_names : list of str, optional
            Gene names
        fdr_threshold : float, default=0.05
            FDR threshold for significance
        n_neighbors : int, default=6
            Number of neighbors
            
        Returns
        -------
        svapa_genes : list of str
            List of SVAPA gene names
        """
        # Test spatial autocorrelation
        results = self.test_spatial_autocorrelation(
            apa_matrix, spatial_coords, gene_names, n_neighbors
        )
        
        # FDR correction
        from statsmodels.stats.multitest import multipletests
        _, padj, _, _ = multipletests(
            results['pvalue'].values,
            alpha=fdr_threshold,
            method='fdr_bh'
        )
        results['padj'] = padj
        
        # Filter significant genes
        svapa_mask = (results['padj'] < fdr_threshold) & (results['morans_i'] > 0)
        svapa_genes = results.loc[svapa_mask, 'gene'].tolist()
        
        self.svapa_genes_ = svapa_genes
        
        return svapa_genes
        
    def compute_weighted_morans_i(
        self,
        apa_values: np.ndarray,
        spatial_weights: np.ndarray,
        uncertainty: np.ndarray,
    ) -> Tuple[float, float, float]:
        """
        Compute uncertainty-weighted Moran's I.

        Higher-uncertainty observations receive lower weight, so noisy
        spots contribute less to spatial autocorrelation detection.

        Parameters
        ----------
        apa_values : np.ndarray, shape (n_spots,)
        spatial_weights : np.ndarray, shape (n_spots, n_spots)
        uncertainty : np.ndarray, shape (n_spots,)
            Per-spot uncertainty std.
        """
        valid = ~np.isnan(apa_values) & ~np.isnan(uncertainty) & (uncertainty > 0)
        if valid.sum() < 3:
            return np.nan, np.nan, np.nan

        values = apa_values[valid]
        w_mat = spatial_weights[np.ix_(valid, valid)]
        unc = uncertainty[valid]
        alpha = 1.0 / (unc + 1e-8)
        alpha = alpha / alpha.sum()  # normalized

        n = len(values)
        w_mean = np.average(values, weights=alpha)
        values_c = values - w_mean

        # Weighted Moran's I
        num = np.sum(
            w_mat * np.outer(values_c, values_c)
            * np.outer(alpha, alpha)
        )
        denom = np.sum(alpha * values_c ** 2)
        W = np.sum(w_mat)

        if denom < 1e-12:
            return np.nan, np.nan, np.nan

        morans_i = (n / W) * (num / denom)
        expected_i = -1.0 / (n - 1)

        # Approximate variance using standard formula on alpha-scaled values
        S1 = 0.5 * np.sum((w_mat + w_mat.T) ** 2)
        S2 = np.sum((np.sum(w_mat, axis=1) + np.sum(w_mat, axis=0)) ** 2)
        S3 = (np.sum(values_c ** 4 * alpha) / n) / (denom / n) ** 2
        S4 = (n ** 2 - 3 * n + 3) * S1 - n * S2 + 3 * W ** 2
        S5 = (n ** 2 - n) * S1 - 2 * n * S2 + 6 * W ** 2
        denom_v = (n - 1) * (n - 2) * (n - 3) * W ** 2
        if denom_v < 1e-12:
            return morans_i, expected_i, 0.0
        variance_i = (n * S4 - S3 * S5) / denom_v - expected_i ** 2

        return morans_i, expected_i, max(variance_i, 1e-12)

    def test_spatial_autocorrelation_weighted(
        self,
        apa_matrix: np.ndarray,
        spatial_coords: np.ndarray,
        uncertainty: np.ndarray,
        gene_names: Optional[List[str]] = None,
        n_neighbors: int = 6,
    ) -> pd.DataFrame:
        """
        Test spatial autocorrelation with uncertainty weighting.

        Parameters
        ----------
        apa_matrix : np.ndarray, (n_genes, n_spots)
        spatial_coords : np.ndarray, (n_spots, 2)
        uncertainty : np.ndarray, (n_genes, n_spots)
            Per-gene-per-spot uncertainty.
        gene_names : list of str, optional
        n_neighbors : int, default=6

        Returns
        -------
        results : pd.DataFrame
            Columns: gene, morans_i, morans_i_weighted, z_score, pvalue
        """
        from ..spatial import build_knn_graph
        from scipy.sparse import csr_matrix

        distances, indices = build_knn_graph(spatial_coords, k=n_neighbors)
        n_spots = len(spatial_coords)
        row_ind = np.repeat(np.arange(n_spots), indices.shape[1])
        col_ind = indices.flatten()
        data = np.ones(len(row_ind))
        spatial_graph = csr_matrix(
            (data, (row_ind, col_ind)), shape=(n_spots, n_spots)
        )
        row_sums = np.array(spatial_graph.sum(axis=1)).flatten()
        row_sums[row_sums == 0] = 1
        spatial_weights = spatial_graph.multiply(1.0 / row_sums[:, np.newaxis])
        spatial_weights = spatial_weights.toarray()

        n_genes = apa_matrix.shape[0]
        if gene_names is None:
            gene_names = [f"Gene_{i}" for i in range(n_genes)]

        results_list = []
        for i, gene in enumerate(gene_names):
            apa_vals = apa_matrix[i, :]
            un_vals = uncertainty[i, :]
            mi_w, exp_w, var_w = self.compute_weighted_morans_i(
                apa_vals, spatial_weights, un_vals
            )
            mi, exp, var = self.compute_morans_i(apa_vals, spatial_weights)

            if np.isnan(mi_w):
                continue

            z = (mi_w - exp_w) / np.sqrt(var_w) if var_w > 0 else 0.0
            pval = 2 * (1 - stats.norm.cdf(np.abs(z)))

            results_list.append({
                'gene': gene,
                'morans_i': mi,
                'morans_i_weighted': mi_w,
                'expected_i': exp_w,
                'z_score': z,
                'pvalue': pval,
            })

        results = pd.DataFrame(results_list)
        self.morans_i_ = results
        return results

    def identify_svapa_genes_weighted(
        self,
        apa_matrix: np.ndarray,
        spatial_coords: np.ndarray,
        uncertainty: np.ndarray,
        gene_names: Optional[List[str]] = None,
        fdr_threshold: float = 0.05,
        n_neighbors: int = 6,
    ) -> List[str]:
        """
        Identify SVAPA genes with uncertainty-weighted Moran's I.

        Parameters
        ----------
        apa_matrix : np.ndarray, (n_genes, n_spots)
        spatial_coords : np.ndarray, (n_spots, 2)
        uncertainty : np.ndarray, (n_genes, n_spots)
        gene_names : list of str, optional
        fdr_threshold : float, default=0.05
        n_neighbors : int, default=6

        Returns
        -------
        svapa_genes : list of str
        """
        results = self.test_spatial_autocorrelation_weighted(
            apa_matrix, spatial_coords, uncertainty, gene_names, n_neighbors
        )
        from statsmodels.stats.multitest import multipletests
        _, padj, _, _ = multipletests(
            results['pvalue'].values, alpha=fdr_threshold, method='fdr_bh'
        )
        results['padj'] = padj
        svapa = results.loc[
            (results['padj'] < fdr_threshold) & (results['morans_i_weighted'] > 0),
            'gene'
        ].tolist()
        self.svapa_genes_ = svapa
        return svapa

    def cluster_spatial_patterns(
        self,
        apa_matrix: np.ndarray,
        spatial_coords: np.ndarray,
        n_patterns: int = 10,
        gene_names: Optional[List[str]] = None
    ) -> Tuple[np.ndarray, pd.DataFrame]:
        """
        Cluster genes by their spatial expression patterns.
        
        Parameters
        ----------
        apa_matrix : np.ndarray, shape (n_genes, n_spots)
            APA index matrix
        spatial_coords : np.ndarray, shape (n_spots, 2)
            Spatial coordinates
        n_patterns : int, default=10
            Number of spatial patterns to identify
        gene_names : list of str, optional
            Gene names
            
        Returns
        -------
        pattern_labels : np.ndarray, shape (n_genes,)
            Pattern label for each gene
        pattern_info : pd.DataFrame
            Information about each pattern
        """
        n_genes = apa_matrix.shape[0]
        
        if gene_names is None:
            gene_names = [f"Gene_{i}" for i in range(n_genes)]
            
        # Standardize APA values
        scaler = StandardScaler()
        apa_scaled = scaler.fit_transform(apa_matrix)
        
        # Hierarchical clustering
        clustering = AgglomerativeClustering(
            n_clusters=n_patterns,
            linkage='ward'
        )
        pattern_labels = clustering.fit_predict(apa_scaled)
        
        # Compute pattern statistics
        pattern_info_list = []
        
        for pattern_id in range(n_patterns):
            mask = pattern_labels == pattern_id
            pattern_genes = [gene_names[i] for i in np.where(mask)[0]]
            
            # Average spatial profile
            pattern_profile = np.mean(apa_matrix[mask, :], axis=0)
            
            # Spatial autocorrelation of pattern
            from ..spatial import build_knn_graph
            from scipy.sparse import csr_matrix
            
            distances, indices = build_knn_graph(spatial_coords, k=6)
            
            # Convert to sparse adjacency matrix
            n_spots = len(spatial_coords)
            row_ind = np.repeat(np.arange(n_spots), indices.shape[1])
            col_ind = indices.flatten()
            data = np.ones(len(row_ind))
            spatial_graph = csr_matrix(
                (data, (row_ind, col_ind)),
                shape=(n_spots, n_spots)
            )
            
            row_sums = np.array(spatial_graph.sum(axis=1)).flatten()
            row_sums[row_sums == 0] = 1
            spatial_weights = spatial_graph.multiply(1.0 / row_sums[:, np.newaxis])
            spatial_weights = spatial_weights.toarray()
            
            morans_i, _, _ = self.compute_morans_i(
                pattern_profile, spatial_weights
            )
            
            pattern_info_list.append({
                'pattern_id': pattern_id,
                'n_genes': mask.sum(),
                'morans_i': morans_i,
                'mean_apa': np.nanmean(pattern_profile),
                'std_apa': np.nanstd(pattern_profile)
            })
            
        pattern_info = pd.DataFrame(pattern_info_list)
        
        return pattern_labels, pattern_info
        
    def compute_pattern_similarity(
        self,
        pattern1: np.ndarray,
        pattern2: np.ndarray,
        method: str = 'pearson'
    ) -> float:
        """
        Compute similarity between two spatial patterns.
        
        Parameters
        ----------
        pattern1 : np.ndarray, shape (n_spots,)
            First spatial pattern
        pattern2 : np.ndarray, shape (n_spots,)
            Second spatial pattern
        method : str, default='pearson'
            Similarity metric: 'pearson', 'spearman', or 'cosine'
            
        Returns
        -------
        similarity : float
            Similarity score
        """
        # Remove NaN values
        valid_mask = ~(np.isnan(pattern1) | np.isnan(pattern2))
        if valid_mask.sum() < 3:
            return np.nan
            
        p1 = pattern1[valid_mask]
        p2 = pattern2[valid_mask]
        
        if method == 'pearson':
            similarity, _ = stats.pearsonr(p1, p2)
        elif method == 'spearman':
            similarity, _ = stats.spearmanr(p1, p2)
        elif method == 'cosine':
            similarity = np.dot(p1, p2) / (np.linalg.norm(p1) * np.linalg.norm(p2))
        else:
            raise ValueError(f"Unknown method: {method}")
            
        return similarity


def identify_svapa_genes(
    apa_matrix: np.ndarray,
    spatial_coords: np.ndarray,
    gene_names: Optional[List[str]] = None,
    fdr_threshold: float = 0.05,
    n_neighbors: int = 6,
    n_permutations: int = 1000
) -> Tuple[List[str], pd.DataFrame]:
    """
    Identify spatially variable APA genes (SVAPA).
    
    Convenience function that wraps SpatialPatternAnalyzer.
    
    Parameters
    ----------
    apa_matrix : np.ndarray, shape (n_genes, n_spots)
        APA index matrix
    spatial_coords : np.ndarray, shape (n_spots, 2)
        Spatial coordinates
    gene_names : list of str, optional
        Gene names
    fdr_threshold : float, default=0.05
        FDR threshold
    n_neighbors : int, default=6
        Number of neighbors
    n_permutations : int, default=1000
        Number of permutations for testing
        
    Returns
    -------
    svapa_genes : list of str
        List of SVAPA gene names
    results : pd.DataFrame
        Full results with Moran's I statistics
        
    Examples
    --------
    >>> import numpy as np
    >>> from spagapa.analysis import identify_svapa_genes
    >>> 
    >>> # Generate example data with spatial pattern
    >>> coords = np.random.rand(200, 2) * 100
    >>> apa_matrix = np.random.rand(50, 200)
    >>> 
    >>> # Add spatial pattern to first gene
    >>> apa_matrix[0, :] = coords[:, 0] / 100  # Gradient along x-axis
    >>> 
    >>> # Identify SVAPA genes
    >>> svapa, results = identify_svapa_genes(apa_matrix, coords)
    >>> print(f"Found {len(svapa)} SVAPA genes")
    >>> print(results.head())
    """
    analyzer = SpatialPatternAnalyzer(
        n_permutations=n_permutations,
        alpha=fdr_threshold
    )
    
    svapa_genes = analyzer.identify_svapa_genes(
        apa_matrix, spatial_coords, gene_names, fdr_threshold, n_neighbors
    )
    
    results = analyzer.morans_i_
    
    return svapa_genes, results


def cluster_spatial_patterns(
    apa_matrix: np.ndarray,
    spatial_coords: np.ndarray,
    n_patterns: int = 10,
    gene_names: Optional[List[str]] = None
) -> Tuple[np.ndarray, pd.DataFrame, Dict[int, List[str]]]:
    """
    Cluster genes by spatial expression patterns.
    
    Convenience function that wraps SpatialPatternAnalyzer.
    
    Parameters
    ----------
    apa_matrix : np.ndarray, shape (n_genes, n_spots)
        APA index matrix
    spatial_coords : np.ndarray, shape (n_spots, 2)
        Spatial coordinates
    n_patterns : int, default=10
        Number of patterns
    gene_names : list of str, optional
        Gene names
        
    Returns
    -------
    pattern_labels : np.ndarray
        Pattern label for each gene
    pattern_info : pd.DataFrame
        Pattern statistics
    pattern_genes : dict
        Dictionary mapping pattern ID to gene list
        
    Examples
    --------
    >>> import numpy as np
    >>> from spagapa.analysis import cluster_spatial_patterns
    >>> 
    >>> # Generate example data
    >>> apa_matrix = np.random.rand(100, 200)
    >>> coords = np.random.rand(200, 2) * 100
    >>> 
    >>> # Cluster patterns
    >>> labels, info, genes = cluster_spatial_patterns(
    ...     apa_matrix, coords, n_patterns=5
    ... )
    >>> print(info)
    >>> print(f"Pattern 0 has {len(genes[0])} genes")
    """
    analyzer = SpatialPatternAnalyzer()
    
    pattern_labels, pattern_info = analyzer.cluster_spatial_patterns(
        apa_matrix, spatial_coords, n_patterns, gene_names
    )
    
    # Create gene lists for each pattern
    if gene_names is None:
        gene_names = [f"Gene_{i}" for i in range(len(pattern_labels))]
        
    pattern_genes = {}
    for pattern_id in range(n_patterns):
        mask = pattern_labels == pattern_id
        pattern_genes[pattern_id] = [
            gene_names[i] for i in np.where(mask)[0]
        ]
        
    return pattern_labels, pattern_info, pattern_genes
