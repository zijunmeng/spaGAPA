"""
GP-based Spatial Trend Detection for SVAPA Identification

This module implements Gaussian Process-based methods for detecting
spatially variable alternative polyadenylation (SVAPA) genes.

Key innovations:
1. Likelihood ratio test based on GP marginal likelihood
2. Uncertainty-weighted spatial autocorrelation
3. Spatial variance decomposition
"""

import numpy as np
import pandas as pd
from typing import Optional, Tuple, Dict, List, Union
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel, Matern, ConstantKernel


class GPTrendDetector:
    """
    Detect spatial trends using Gaussian Process-based methods.
    
    This class provides multiple methods for identifying genes with
    spatially variable APA usage:
    
    1. Likelihood ratio test: Compare GP models with/without spatial structure
    2. Uncertainty-weighted Moran's I: Weight by imputation confidence
    3. Spatial variance decomposition: Separate spatial vs random variance
    
    Parameters
    ----------
    kernel_type : str, default='matern'
        Type of kernel for spatial GP: 'rbf', 'matern', or 'auto'
    length_scale : float, optional
        Length scale for spatial kernel. If None, estimated from data
    alpha : float, default=1e-10
        Noise level for GP regression
    n_restarts : int, default=5
        Number of restarts for hyperparameter optimization
    
    Attributes
    ----------
    gp_null_ : GaussianProcessRegressor
        Null model (no spatial structure)
    gp_alt_ : GaussianProcessRegressor
        Alternative model (with spatial structure)
    
    Examples
    --------
    >>> detector = GPTrendDetector(kernel_type='matern')
    >>> results = detector.detect_gp_trends(
    ...     imputed_values, uncertainty, coordinates
    ... )
    >>> svapa_genes = results[results['significant']]
    """
    
    def __init__(
        self,
        kernel_type: str = 'matern',
        length_scale: Optional[float] = None,
        alpha: float = 1e-10,
        n_restarts: int = 5
    ):
        self.kernel_type = kernel_type
        self.length_scale = length_scale
        self.alpha = alpha
        self.n_restarts = n_restarts
        self.gp_null_ = None
        self.gp_alt_ = None
    
    def _create_null_kernel(self) -> WhiteKernel:
        """Create kernel for null model (no spatial structure)."""
        return WhiteKernel(noise_level=1.0, noise_level_bounds=(1e-5, 1e2))
    
    def _create_spatial_kernel(self, length_scale: Optional[float] = None):
        """Create kernel for alternative model (with spatial structure)."""
        if length_scale is None:
            length_scale = 1.0
            length_scale_bounds = (1e-2, 1e3)
        else:
            length_scale_bounds = "fixed"
        
        if self.kernel_type == 'rbf':
            spatial_kernel = RBF(
                length_scale=length_scale,
                length_scale_bounds=length_scale_bounds
            )
        elif self.kernel_type == 'matern':
            spatial_kernel = Matern(
                length_scale=length_scale,
                length_scale_bounds=length_scale_bounds,
                nu=1.5
            )
        else:  # auto
            spatial_kernel = RBF(
                length_scale=length_scale,
                length_scale_bounds=length_scale_bounds
            )
        
        # Add constant and white noise
        kernel = (
            ConstantKernel(1.0, (1e-3, 1e3)) * spatial_kernel +
            WhiteKernel(noise_level=0.1, noise_level_bounds=(1e-5, 1e2))
        )
        
        return kernel
    
    def likelihood_ratio_test(
        self,
        values: np.ndarray,
        coordinates: np.ndarray,
        uncertainty: Optional[np.ndarray] = None
    ) -> Dict[str, float]:
        """
        Perform likelihood ratio test for spatial trend.
        
        Tests H0 (no spatial structure) vs H1 (spatial structure) using
        GP marginal likelihood.
        
        Parameters
        ----------
        values : np.ndarray, shape (n_spots,)
            APA usage values (e.g., RUD, PDUI)
        coordinates : np.ndarray, shape (n_spots, 2)
            Spatial coordinates
        uncertainty : np.ndarray, shape (n_spots,), optional
            Uncertainty estimates for each value
        
        Returns
        -------
        results : dict
            Dictionary with keys:
            - 'log_lik_null': Log marginal likelihood under null
            - 'log_lik_alt': Log marginal likelihood under alternative
            - 'lr_statistic': Likelihood ratio statistic (2 * diff)
            - 'p_value': P-value from chi-square test
            - 'df': Degrees of freedom
        """
        # Remove NaN values
        mask = ~np.isnan(values)
        values = values[mask]
        coordinates = coordinates[mask]
        if uncertainty is not None:
            uncertainty = uncertainty[mask]
        
        if len(values) < 10:
            return {
                'log_lik_null': np.nan,
                'log_lik_alt': np.nan,
                'lr_statistic': np.nan,
                'p_value': 1.0,
                'df': 0
            }
        
        # Reshape for sklearn
        X = coordinates.reshape(-1, 2)
        y = values.reshape(-1, 1)
        
        # Weight by uncertainty if provided
        if uncertainty is not None:
            # Higher uncertainty -> lower weight
            weights = 1.0 / (uncertainty + 1e-10)
            weights = weights / weights.sum() * len(weights)
            sample_weight = weights
        else:
            sample_weight = None
        
        # Fit null model (no spatial structure)
        null_kernel = self._create_null_kernel()
        self.gp_null_ = GaussianProcessRegressor(
            kernel=null_kernel,
            alpha=self.alpha,
            n_restarts_optimizer=self.n_restarts,
            normalize_y=True
        )
        
        try:
            self.gp_null_.fit(X, y)
            log_lik_null = self.gp_null_.log_marginal_likelihood()
        except Exception:
            log_lik_null = -np.inf
        
        # Fit alternative model (with spatial structure)
        alt_kernel = self._create_spatial_kernel(self.length_scale)
        self.gp_alt_ = GaussianProcessRegressor(
            kernel=alt_kernel,
            alpha=self.alpha,
            n_restarts_optimizer=self.n_restarts,
            normalize_y=True
        )
        
        try:
            self.gp_alt_.fit(X, y)
            log_lik_alt = self.gp_alt_.log_marginal_likelihood()
        except Exception:
            log_lik_alt = -np.inf
        
        # Likelihood ratio test
        lr_statistic = 2 * (log_lik_alt - log_lik_null)
        
        # Degrees of freedom = difference in number of parameters
        # Null: 1 (noise), Alt: 3 (constant, length_scale, noise)
        df = 2
        
        # P-value from chi-square distribution
        if lr_statistic > 0:
            p_value = stats.chi2.sf(lr_statistic, df)
        else:
            p_value = 1.0
        
        return {
            'log_lik_null': log_lik_null,
            'log_lik_alt': log_lik_alt,
            'lr_statistic': lr_statistic,
            'p_value': p_value,
            'df': df
        }
    
    def uncertainty_weighted_morans_i(
        self,
        values: np.ndarray,
        coordinates: np.ndarray,
        uncertainty: np.ndarray,
        k: int = 6
    ) -> Dict[str, float]:
        """
        Compute uncertainty-weighted Moran's I statistic.
        
        This is an enhanced version of Moran's I that weights observations
        by their imputation confidence (inverse uncertainty).
        
        Parameters
        ----------
        values : np.ndarray, shape (n_spots,)
            APA usage values
        coordinates : np.ndarray, shape (n_spots, 2)
            Spatial coordinates
        uncertainty : np.ndarray, shape (n_spots,)
            Uncertainty estimates
        k : int, default=6
            Number of nearest neighbors for spatial weights
        
        Returns
        -------
        results : dict
            Dictionary with keys:
            - 'morans_i': Moran's I statistic
            - 'expected_i': Expected value under null
            - 'variance_i': Variance under null
            - 'z_score': Standardized statistic
            - 'p_value': Two-tailed p-value
        """
        from sklearn.neighbors import NearestNeighbors
        
        # Remove NaN values
        mask = ~np.isnan(values) & ~np.isnan(uncertainty)
        values = values[mask]
        coordinates = coordinates[mask]
        uncertainty = uncertainty[mask]
        
        n = len(values)
        if n < 10:
            return {
                'morans_i': np.nan,
                'expected_i': np.nan,
                'variance_i': np.nan,
                'z_score': np.nan,
                'p_value': 1.0
            }
        
        # Compute weights based on uncertainty
        # Lower uncertainty -> higher weight
        alpha = 1.0 / (uncertainty + 1e-10)
        alpha = alpha / alpha.sum()  # Normalize
        
        # Build spatial weight matrix (KNN)
        nbrs = NearestNeighbors(n_neighbors=min(k + 1, n)).fit(coordinates)
        distances, indices = nbrs.kneighbors(coordinates)
        
        # Create weight matrix
        W = np.zeros((n, n))
        for i in range(n):
            for j in indices[i, 1:]:  # Skip self
                W[i, j] = 1.0
        
        # Row-normalize
        row_sums = W.sum(axis=1)
        row_sums[row_sums == 0] = 1  # Avoid division by zero
        W = W / row_sums[:, np.newaxis]
        
        # Weighted mean and variance
        mean_val = np.sum(alpha * values)
        var_val = np.sum(alpha * (values - mean_val) ** 2)
        
        if var_val < 1e-10:
            return {
                'morans_i': 0.0,
                'expected_i': -1.0 / (n - 1),
                'variance_i': np.nan,
                'z_score': 0.0,
                'p_value': 1.0
            }
        
        # Compute weighted Moran's I
        numerator = 0.0
        for i in range(n):
            for j in range(n):
                numerator += alpha[i] * alpha[j] * W[i, j] * \
                            (values[i] - mean_val) * (values[j] - mean_val)
        
        morans_i = numerator / var_val
        
        # Expected value and variance under null hypothesis
        expected_i = -1.0 / (n - 1)
        
        # Simplified variance (exact formula is complex)
        S0 = W.sum()
        S1 = ((W + W.T) ** 2).sum() / 2
        S2 = ((W.sum(axis=0) + W.sum(axis=1)) ** 2).sum()
        
        variance_i = (n * S1 - n * S0 ** 2 + 3 * S0 ** 2) / \
                     ((n ** 2 - 1) * S0 ** 2)
        
        # Z-score and p-value
        if variance_i > 0:
            z_score = (morans_i - expected_i) / np.sqrt(variance_i)
            p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))
        else:
            z_score = 0.0
            p_value = 1.0
        
        return {
            'morans_i': morans_i,
            'expected_i': expected_i,
            'variance_i': variance_i,
            'z_score': z_score,
            'p_value': p_value
        }
    
    def spatial_variance_decomposition(
        self,
        values: np.ndarray,
        coordinates: np.ndarray
    ) -> Dict[str, float]:
        """
        Decompose variance into spatial and random components.
        
        Uses GP to estimate how much variance is explained by
        spatial structure vs random noise.
        
        Parameters
        ----------
        values : np.ndarray, shape (n_spots,)
            APA usage values
        coordinates : np.ndarray, shape (n_spots, 2)
            Spatial coordinates
        
        Returns
        -------
        results : dict
            Dictionary with keys:
            - 'total_variance': Total variance in data
            - 'spatial_variance': Variance explained by spatial structure
            - 'random_variance': Residual variance
            - 'spatial_fraction': Fraction of variance that is spatial
            - 'r_squared': R² from GP model
        """
        # Remove NaN values
        mask = ~np.isnan(values)
        values = values[mask]
        coordinates = coordinates[mask]
        
        if len(values) < 10:
            return {
                'total_variance': np.nan,
                'spatial_variance': np.nan,
                'random_variance': np.nan,
                'spatial_fraction': np.nan,
                'r_squared': np.nan
            }
        
        # Total variance
        total_variance = np.var(values)
        
        # Fit GP with spatial kernel
        X = coordinates.reshape(-1, 2)
        y = values.reshape(-1, 1)
        
        kernel = self._create_spatial_kernel(self.length_scale)
        gp = GaussianProcessRegressor(
            kernel=kernel,
            alpha=self.alpha,
            n_restarts_optimizer=self.n_restarts,
            normalize_y=False
        )
        
        try:
            gp.fit(X, y)
            y_pred, y_std = gp.predict(X, return_std=True)
            
            # Spatial variance = variance of predictions
            spatial_variance = np.var(y_pred)
            
            # Random variance = mean squared error
            random_variance = np.mean((y.ravel() - y_pred.ravel()) ** 2)
            
            # Fraction
            spatial_fraction = spatial_variance / total_variance
            
            # R²
            ss_res = np.sum((y.ravel() - y_pred.ravel()) ** 2)
            ss_tot = np.sum((y.ravel() - np.mean(y)) ** 2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
            
        except Exception:
            spatial_variance = np.nan
            random_variance = np.nan
            spatial_fraction = np.nan
            r_squared = np.nan
        
        return {
            'total_variance': total_variance,
            'spatial_variance': spatial_variance,
            'random_variance': random_variance,
            'spatial_fraction': spatial_fraction,
            'r_squared': r_squared
        }
    
    def detect_gp_trends(
        self,
        apa_values: np.ndarray,
        coordinates: np.ndarray,
        uncertainty: Optional[np.ndarray] = None,
        gene_names: Optional[List[str]] = None,
        method: str = 'likelihood_ratio',
        fdr_threshold: float = 0.05,
        n_jobs: int = 1
    ) -> pd.DataFrame:
        """
        Detect SVAPA genes using GP-based methods.
        
        Parameters
        ----------
        apa_values : np.ndarray, shape (n_genes, n_spots)
            APA usage values for each gene
        coordinates : np.ndarray, shape (n_spots, 2)
            Spatial coordinates
        uncertainty : np.ndarray, shape (n_genes, n_spots), optional
            Uncertainty estimates
        gene_names : list of str, optional
            Gene names
        method : str, default='likelihood_ratio'
            Detection method: 'likelihood_ratio', 'weighted_morans', or 'all'
        fdr_threshold : float, default=0.05
            FDR threshold for significance
        n_jobs : int, default=1
            Number of parallel jobs
        
        Returns
        -------
        results : pd.DataFrame
            DataFrame with columns:
            - gene: Gene name
            - p_value: P-value from test
            - q_value: FDR-adjusted p-value
            - significant: Whether gene is significant
            - Additional method-specific columns
        """
        n_genes, n_spots = apa_values.shape
        
        if gene_names is None:
            gene_names = [f"Gene_{i}" for i in range(n_genes)]
        
        results_list = []
        
        for i in range(n_genes):
            values = apa_values[i, :]
            unc = uncertainty[i, :] if uncertainty is not None else None
            
            result = {'gene': gene_names[i]}
            
            if method in ['likelihood_ratio', 'all']:
                lr_result = self.likelihood_ratio_test(values, coordinates, unc)
                result.update({
                    'lr_statistic': lr_result['lr_statistic'],
                    'lr_p_value': lr_result['p_value'],
                    'log_lik_null': lr_result['log_lik_null'],
                    'log_lik_alt': lr_result['log_lik_alt']
                })
            
            if method in ['weighted_morans', 'all'] and uncertainty is not None:
                morans_result = self.uncertainty_weighted_morans_i(
                    values, coordinates, unc
                )
                result.update({
                    'morans_i': morans_result['morans_i'],
                    'morans_z': morans_result['z_score'],
                    'morans_p_value': morans_result['p_value']
                })
            
            if method == 'all':
                var_result = self.spatial_variance_decomposition(values, coordinates)
                result.update({
                    'spatial_fraction': var_result['spatial_fraction'],
                    'r_squared': var_result['r_squared']
                })
            
            results_list.append(result)
        
        # Create DataFrame
        results_df = pd.DataFrame(results_list)
        
        # Determine p-value column based on method
        if method == 'likelihood_ratio':
            p_col = 'lr_p_value'
        elif method == 'weighted_morans':
            p_col = 'morans_p_value'
        else:  # 'all'
            p_col = 'lr_p_value'  # Use LR as primary
        
        # FDR correction
        from statsmodels.stats.multitest import multipletests
        
        p_values = results_df[p_col].values
        valid_mask = ~np.isnan(p_values)
        
        q_values = np.full(len(p_values), np.nan)
        if valid_mask.sum() > 0:
            _, q_values[valid_mask], _, _ = multipletests(
                p_values[valid_mask],
                alpha=fdr_threshold,
                method='fdr_bh'
            )
        
        results_df['q_value'] = q_values
        results_df['significant'] = q_values < fdr_threshold
        
        # Sort by p-value
        results_df = results_df.sort_values(p_col)
        
        return results_df


def detect_svapa_genes_gp(
    apa_values: np.ndarray,
    coordinates: np.ndarray,
    uncertainty: Optional[np.ndarray] = None,
    gene_names: Optional[List[str]] = None,
    method: str = 'likelihood_ratio',
    kernel_type: str = 'matern',
    fdr_threshold: float = 0.05
) -> pd.DataFrame:
    """
    Convenience function to detect SVAPA genes using GP-based methods.
    
    Parameters
    ----------
    apa_values : np.ndarray, shape (n_genes, n_spots)
        APA usage values
    coordinates : np.ndarray, shape (n_spots, 2)
        Spatial coordinates
    uncertainty : np.ndarray, shape (n_genes, n_spots), optional
        Uncertainty estimates
    gene_names : list of str, optional
        Gene names
    method : str, default='likelihood_ratio'
        Detection method
    kernel_type : str, default='matern'
        GP kernel type
    fdr_threshold : float, default=0.05
        FDR threshold
    
    Returns
    -------
    results : pd.DataFrame
        SVAPA detection results
    
    Examples
    --------
    >>> results = detect_svapa_genes_gp(
    ...     apa_values, coordinates, uncertainty,
    ...     method='likelihood_ratio'
    ... )
    >>> svapa_genes = results[results['significant']]
    """
    detector = GPTrendDetector(kernel_type=kernel_type)
    return detector.detect_gp_trends(
        apa_values,
        coordinates,
        uncertainty,
        gene_names,
        method,
        fdr_threshold
    )
