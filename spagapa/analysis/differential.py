"""
Differential APA analysis between spatial domains or conditions.

This module provides statistical tests to identify genes with differential
alternative polyadenylation (APA) usage between groups.
"""

import numpy as np
import pandas as pd
from typing import Optional, Union, Tuple, List, Dict
from scipy import stats
from statsmodels.stats.multitest import multipletests


class DifferentialAPAAnalyzer:
    """
    Perform differential APA analysis between groups.
    
    This class provides methods to test for differential APA usage between
    two or more groups using various statistical tests.
    
    Parameters
    ----------
    method : str, default='wilcoxon'
        Statistical test method: 'wilcoxon', 't-test', or 'permutation'
    alpha : float, default=0.05
        Significance level
    min_spots_per_group : int, default=3
        Minimum number of spots required in each group
    pseudocount : float, default=1e-6
        Pseudocount to avoid division by zero
        
    Attributes
    ----------
    results_ : pd.DataFrame
        Differential APA test results
    """
    
    def __init__(
        self,
        method: str = 'wilcoxon',
        alpha: float = 0.05,
        min_spots_per_group: int = 3,
        pseudocount: float = 1e-6
    ):
        self.method = method
        self.alpha = alpha
        self.min_spots_per_group = min_spots_per_group
        self.pseudocount = pseudocount
        
        self.results_ = None
        
    def _get_spot_weights(
        self,
        gene_idx: int,
        uncertainty: Optional[np.ndarray] = None,
    ) -> Optional[np.ndarray]:
        """Convert gene-level uncertainty to per-spot weights."""
        if uncertainty is None:
            return None
        gene_unc = uncertainty[gene_idx, :]  # (n_spots,)
        valid = ~np.isnan(gene_unc) & (gene_unc > 0)
        if valid.sum() < 3:
            return None
        w = np.ones_like(gene_unc)
        w[valid] = 1.0 / (gene_unc[valid] + 1e-8)
        return w

    def test_differential_apa(
        self,
        apa_matrix: np.ndarray,
        group1_indices: np.ndarray,
        group2_indices: np.ndarray,
        gene_names: Optional[List[str]] = None,
        uncertainty: Optional[np.ndarray] = None,
    ) -> pd.DataFrame:
        """
        Test for differential APA between two groups.
        
        Parameters
        ----------
        apa_matrix : np.ndarray, shape (n_genes, n_spots)
            APA index matrix
        group1_indices : np.ndarray
            Indices of spots in group 1
        group2_indices : np.ndarray
            Indices of spots in group 2
        gene_names : list of str, optional
            Gene names
            
        Returns
        -------
        results : pd.DataFrame
            Test results with columns:
            - gene: gene name
            - mean_group1: mean APA in group 1
            - mean_group2: mean APA in group 2
            - log2fc: log2 fold change
            - pvalue: p-value
            - statistic: test statistic
        """
        # Check group sizes
        if len(group1_indices) < self.min_spots_per_group:
            raise ValueError(
                f"Group 1 has only {len(group1_indices)} spots, "
                f"need at least {self.min_spots_per_group}"
            )
        if len(group2_indices) < self.min_spots_per_group:
            raise ValueError(
                f"Group 2 has only {len(group2_indices)} spots, "
                f"need at least {self.min_spots_per_group}"
            )
            
        n_genes = apa_matrix.shape[0]
        
        if gene_names is None:
            gene_names = [f"Gene_{i}" for i in range(n_genes)]
            
        results_list = []
        
        for i, gene in enumerate(gene_names):
            group1_values = apa_matrix[i, group1_indices]
            group2_values = apa_matrix[i, group2_indices]
            
            # Remove NaN values
            group1_values = group1_values[~np.isnan(group1_values)]
            group2_values = group2_values[~np.isnan(group2_values)]
            
            if len(group1_values) == 0 or len(group2_values) == 0:
                continue
                
            # Compute statistics
            mean1 = np.mean(group1_values)
            mean2 = np.mean(group2_values)
            
            # Log2 fold change
            log2fc = np.log2(
                (mean2 + self.pseudocount) / (mean1 + self.pseudocount)
            )
            
            # Get per-spot weights for this gene
            spot_weights = self._get_spot_weights(i, uncertainty)

            # Statistical test (weighted if uncertainty available)
            if spot_weights is not None:
                w1 = spot_weights[group1_indices]
                w2 = spot_weights[group2_indices]
                if self.method == 'wilcoxon':
                    statistic, pvalue = self._weighted_wilcoxon(
                        group1_values, group2_values, w1, w2
                    )
                else:
                    # fallback: filter by weight
                    g1 = group1_values[w1 > np.median(w1)] if len(w1) > 5 else group1_values
                    g2 = group2_values[w2 > np.median(w2)] if len(w2) > 5 else group2_values
                    statistic, pvalue = self._wilcoxon_test(g1, g2)
            elif self.method == 'wilcoxon':
                statistic, pvalue = self._wilcoxon_test(
                    group1_values, group2_values
                )
            elif self.method == 't-test':
                statistic, pvalue = self._ttest(
                    group1_values, group2_values
                )
            elif self.method == 'permutation':
                statistic, pvalue = self._permutation_test(
                    group1_values, group2_values
                )
            else:
                raise ValueError(f"Unknown method: {self.method}")
                
            results_list.append({
                'gene': gene,
                'mean_group1': mean1,
                'mean_group2': mean2,
                'log2fc': log2fc,
                'statistic': statistic,
                'pvalue': pvalue
            })
            
        results = pd.DataFrame(results_list)
        self.results_ = results
        
        return results
        
    def _wilcoxon_test(
        self,
        group1: np.ndarray,
        group2: np.ndarray
    ) -> Tuple[float, float]:
        """Perform Wilcoxon rank-sum test (Mann-Whitney U test)."""
        try:
            statistic, pvalue = stats.mannwhitneyu(
                group1, group2, alternative='two-sided'
            )
            return statistic, pvalue
        except:
            return np.nan, 1.0
            
    def _ttest(
        self,
        group1: np.ndarray,
        group2: np.ndarray
    ) -> Tuple[float, float]:
        """Perform Welch's t-test (unequal variances)."""
        try:
            statistic, pvalue = stats.ttest_ind(
                group1, group2, equal_var=False
            )
            return statistic, pvalue
        except:
            return np.nan, 1.0
            
    def _weighted_wilcoxon(
        self,
        group1: np.ndarray,
        group2: np.ndarray,
        w1: np.ndarray,
        w2: np.ndarray,
    ) -> Tuple[float, float]:
        """
        Weighted Mann-Whitney U test.

        Uses importance sampling: resamples values proportional to weights,
        then applies the standard test. This gives higher-weight (low
        uncertainty) observations more influence.
        """
        n1, n2 = len(group1), len(group2)
        # Normalize weights to sampling probabilities
        p1 = np.maximum(w1, 0) / (np.sum(np.maximum(w1, 0)) + 1e-12)
        p2 = np.maximum(w2, 0) / (np.sum(np.maximum(w2, 0)) + 1e-12)

        # Weighted resampling (bootstrap)
        n_boot = min(max(n1, n2) * 3, 500)
        rng = np.random.RandomState(42)
        g1_boot = rng.choice(group1, size=n_boot, p=p1, replace=True)
        g2_boot = rng.choice(group2, size=n_boot, p=p2, replace=True)

        try:
            stat, pval = stats.mannwhitneyu(g1_boot, g2_boot, alternative='two-sided')
            return stat, pval
        except Exception:
            return np.nan, 1.0

    def _permutation_test(
        self,
        group1: np.ndarray,
        group2: np.ndarray,
        n_permutations: int = 1000
    ) -> Tuple[float, float]:
        """Perform permutation test."""
        # Observed difference
        obs_diff = np.mean(group2) - np.mean(group1)
        
        # Combine groups
        combined = np.concatenate([group1, group2])
        n1 = len(group1)
        
        # Permutation
        null_diffs = []
        for _ in range(n_permutations):
            np.random.shuffle(combined)
            perm_group1 = combined[:n1]
            perm_group2 = combined[n1:]
            null_diff = np.mean(perm_group2) - np.mean(perm_group1)
            null_diffs.append(null_diff)
            
        null_diffs = np.array(null_diffs)
        
        # P-value (two-sided)
        pvalue = np.mean(np.abs(null_diffs) >= np.abs(obs_diff))
        
        return obs_diff, pvalue
        
    def test_all_pairwise(
        self,
        apa_matrix: np.ndarray,
        domain_labels: np.ndarray,
        gene_names: Optional[List[str]] = None,
        uncertainty: Optional[np.ndarray] = None,
    ) -> Dict[str, pd.DataFrame]:
        """
        Test all pairwise comparisons between domains.
        
        Parameters
        ----------
        apa_matrix : np.ndarray, shape (n_genes, n_spots)
            APA index matrix
        domain_labels : np.ndarray, shape (n_spots,)
            Domain labels for each spot
        gene_names : list of str, optional
            Gene names
        uncertainty : np.ndarray, optional, shape (n_genes, n_spots)
            Per-value uncertainty for weighted testing.
            
        Returns
        -------
        results : dict
            Dictionary mapping comparison names to result DataFrames
        """
        unique_domains = np.unique(domain_labels)
        results = {}
        
        for i, domain1 in enumerate(unique_domains):
            for domain2 in unique_domains[i+1:]:
                group1_indices = np.where(domain_labels == domain1)[0]
                group2_indices = np.where(domain_labels == domain2)[0]
                
                comparison_name = f"Domain{domain1}_vs_Domain{domain2}"
                
                try:
                    result = self.test_differential_apa(
                        apa_matrix,
                        group1_indices,
                        group2_indices,
                        gene_names,
                        uncertainty=uncertainty,
                    )
                    results[comparison_name] = result
                except ValueError as e:
                    print(f"Skipping {comparison_name}: {e}")
                    
        return results
        
    def adjust_pvalues(
        self,
        results: pd.DataFrame,
        method: str = 'fdr_bh'
    ) -> pd.DataFrame:
        """
        Adjust p-values for multiple testing.
        
        Parameters
        ----------
        results : pd.DataFrame
            Results from test_differential_apa
        method : str, default='fdr_bh'
            Correction method:
            - 'bonferroni': Bonferroni correction
            - 'fdr_bh': Benjamini-Hochberg FDR
            - 'fdr_by': Benjamini-Yekutieli FDR
            
        Returns
        -------
        results : pd.DataFrame
            Results with added 'padj' column
        """
        if 'pvalue' not in results.columns:
            raise ValueError("Results must contain 'pvalue' column")
            
        pvalues = results['pvalue'].values
        
        # Handle NaN p-values
        valid_mask = ~np.isnan(pvalues)
        padj = np.full(len(pvalues), np.nan)
        
        if valid_mask.sum() > 0:
            _, padj[valid_mask], _, _ = multipletests(
                pvalues[valid_mask],
                alpha=self.alpha,
                method=method
            )
            
        results = results.copy()
        results['padj'] = padj
        
        return results
        
    def rank_genes(
        self,
        results: pd.DataFrame,
        by: str = 'padj',
        ascending: bool = True,
        top_n: Optional[int] = None
    ) -> pd.DataFrame:
        """
        Rank genes by a specified metric.
        
        Parameters
        ----------
        results : pd.DataFrame
            Results from test_differential_apa
        by : str, default='padj'
            Column to rank by
        ascending : bool, default=True
            Sort order
        top_n : int, optional
            Return only top N genes
            
        Returns
        -------
        ranked : pd.DataFrame
            Ranked results
        """
        ranked = results.sort_values(by=by, ascending=ascending)
        
        if top_n is not None:
            ranked = ranked.head(top_n)
            
        return ranked
        
    def filter_results(
        self,
        results: pd.DataFrame,
        padj_threshold: float = 0.05,
        logfc_threshold: float = 0.5,
        abs_logfc: bool = True
    ) -> pd.DataFrame:
        """
        Filter results by significance and effect size.
        
        Parameters
        ----------
        results : pd.DataFrame
            Results from test_differential_apa
        padj_threshold : float, default=0.05
            Adjusted p-value threshold
        logfc_threshold : float, default=0.5
            Log2 fold change threshold
        abs_logfc : bool, default=True
            Use absolute value of log2FC
            
        Returns
        -------
        filtered : pd.DataFrame
            Filtered results
        """
        if 'padj' not in results.columns:
            results = self.adjust_pvalues(results)
            
        # Filter by adjusted p-value
        mask = results['padj'] < padj_threshold
        
        # Filter by log2 fold change
        if abs_logfc:
            mask &= np.abs(results['log2fc']) > logfc_threshold
        else:
            mask &= results['log2fc'] > logfc_threshold
            
        filtered = results[mask].copy()
        
        return filtered


def run_differential_apa_test(
    apa_matrix: np.ndarray,
    group1_indices: np.ndarray,
    group2_indices: np.ndarray,
    gene_names: Optional[List[str]] = None,
    method: str = 'wilcoxon',
    adjust_method: str = 'fdr_bh',
    padj_threshold: float = 0.05,
    logfc_threshold: float = 0.5,
    uncertainty: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """
    Test for differential APA between two groups.

    Convenience function that wraps DifferentialAPAAnalyzer.

    Parameters
    ----------
    apa_matrix : np.ndarray, shape (n_genes, n_spots)
        APA index matrix
    group1_indices : np.ndarray
        Indices of spots in group 1
    group2_indices : np.ndarray
        Indices of spots in group 2
    gene_names : list of str, optional
        Gene names
    method : str, default='wilcoxon'
        Statistical test method
    adjust_method : str, default='fdr_bh'
        Multiple testing correction method
    padj_threshold : float, default=0.05
        Adjusted p-value threshold for filtering
    logfc_threshold : float, default=0.5
        Log2 fold change threshold for filtering
        
    Returns
    -------
    results : pd.DataFrame
        Differential APA test results with adjusted p-values
        
    Examples
    --------
    >>> import numpy as np
    >>> from spagapa.analysis import run_differential_apa_test
    >>>
    >>> # Generate example data
    >>> apa_matrix = np.random.rand(100, 200)  # 100 genes, 200 spots
    >>> group1 = np.arange(100)  # First 100 spots
    >>> group2 = np.arange(100, 200)  # Last 100 spots
    >>>
    >>> # Test differential APA
    >>> results = run_differential_apa_test(
    ...     apa_matrix, group1, group2, method='wilcoxon'
    ... )
    >>> 
    >>> # Filter significant genes
    >>> sig_genes = results[results['padj'] < 0.05]
    >>> print(f"Found {len(sig_genes)} significant genes")
    """
    analyzer = DifferentialAPAAnalyzer(method=method)

    # Perform test (with uncertainty weights if available)
    results = analyzer.test_differential_apa(
        apa_matrix, group1_indices, group2_indices, gene_names,
        uncertainty=uncertainty,
    )
    
    # Adjust p-values
    results = analyzer.adjust_pvalues(results, method=adjust_method)
    
    return results


def find_domain_markers(
    apa_matrix: np.ndarray,
    domain_labels: np.ndarray,
    gene_names: Optional[List[str]] = None,
    method: str = 'wilcoxon',
    padj_threshold: float = 0.05,
    logfc_threshold: float = 0.5,
    uncertainty: Optional[np.ndarray] = None,
) -> Dict[int, pd.DataFrame]:
    """
    Find marker genes for each domain (one-vs-rest comparison).
    
    Parameters
    ----------
    apa_matrix : np.ndarray, shape (n_genes, n_spots)
        APA index matrix
    domain_labels : np.ndarray, shape (n_spots,)
        Domain labels
    gene_names : list of str, optional
        Gene names
    method : str, default='wilcoxon'
        Statistical test method
    padj_threshold : float, default=0.05
        Adjusted p-value threshold
    logfc_threshold : float, default=0.5
        Log2 fold change threshold
        
    Returns
    -------
    markers : dict
        Dictionary mapping domain ID to marker gene DataFrame
        
    Examples
    --------
    >>> import numpy as np
    >>> from spagapa.analysis import find_domain_markers
    >>> 
    >>> # Generate example data
    >>> apa_matrix = np.random.rand(100, 200)
    >>> domains = np.repeat([0, 1, 2, 3], 50)
    >>> 
    >>> # Find markers
    >>> markers = find_domain_markers(apa_matrix, domains)
    >>> 
    >>> # Print markers for domain 0
    >>> print(markers[0].head())
    """
    analyzer = DifferentialAPAAnalyzer(method=method)
    unique_domains = np.unique(domain_labels)
    
    markers = {}
    
    for domain in unique_domains:
        # One vs rest
        group1_indices = np.where(domain_labels == domain)[0]
        group2_indices = np.where(domain_labels != domain)[0]
        
        try:
            results = analyzer.test_differential_apa(
                apa_matrix, group1_indices, group2_indices, gene_names,
                uncertainty=uncertainty,
            )
            
            # Adjust p-values
            results = analyzer.adjust_pvalues(results)
            
            # Filter
            results = analyzer.filter_results(
                results, padj_threshold, logfc_threshold
            )
            
            # Rank by log2FC (descending)
            results = analyzer.rank_genes(
                results, by='log2fc', ascending=False
            )
            
            markers[domain] = results
            
        except ValueError as e:
            print(f"Skipping domain {domain}: {e}")
            
    return markers
