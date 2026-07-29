"""
Gaussian Process-based imputation for spatial APA data.

This module implements Gaussian Process (GP) regression for imputing missing
or low-quality APA measurements in spatial transcriptomics data, with
uncertainty quantification.
"""

import numpy as np
from typing import Optional, Tuple, Dict, Literal
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import (
    RBF, Matern, WhiteKernel, ConstantKernel as C
)
from scipy.spatial.distance import cdist
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

try:
    from threadpoolctl import threadpool_limits
except Exception:  # pragma: no cover - optional runtime optimization
    threadpool_limits = None

logger = logging.getLogger(__name__)


class GPImputer:
    """
    Gaussian Process imputer for spatial APA data.
    
    Uses spatial coordinates to impute missing or low-quality measurements
    with uncertainty quantification. Particularly useful for sparse spatial
    transcriptomics data.
    
    Parameters
    ----------
    kernel_type : str, default='matern'
        Type of kernel: 'rbf', 'matern', or 'auto'
    length_scale : float, optional
        Length scale for the kernel. If None, will be optimized.
    nu : float, default=1.5
        Smoothness parameter for Matérn kernel (0.5, 1.5, or 2.5)
    alpha : float, default=1e-10
        Noise level (regularization parameter)
    n_restarts_optimizer : int, default=5
        Number of restarts for kernel hyperparameter optimization
    normalize_y : bool, default=True
        Whether to normalize target values
    
    Attributes
    ----------
    gp_ : GaussianProcessRegressor
        Fitted GP model
    kernel_ : Kernel
        Kernel used for GP
    
    Examples
    --------
    >>> imputer = GPImputer(kernel_type='matern', nu=1.5)
    >>> imputer.fit(spatial_coords, apa_counts)
    >>> imputed_counts, uncertainty = imputer.predict(spatial_coords)
    """
    
    def __init__(
        self,
        kernel_type: Literal['rbf', 'matern', 'auto'] = 'matern',
        length_scale: Optional[float] = None,
        nu: float = 1.5,
        alpha: float = 1e-10,
        n_restarts_optimizer: int = 5,
        normalize_y: bool = True
    ):
        self.kernel_type = kernel_type
        self.length_scale = length_scale
        self.nu = nu
        self.alpha = alpha
        self.n_restarts_optimizer = n_restarts_optimizer
        self.normalize_y = normalize_y
        
        self.gp_ = None
        self.kernel_ = None
        self._fitted_coords = None
        self._fitted_values = None
    
    def _create_kernel(self, length_scale: Optional[float] = None):
        """
        Create GP kernel based on kernel_type.
        
        Parameters
        ----------
        length_scale : float, optional
            Length scale for the kernel
        
        Returns
        -------
        Kernel
            Scikit-learn kernel object
        """
        if length_scale is None:
            length_scale = self.length_scale if self.length_scale else 1.0
        
        if self.kernel_type == 'rbf':
            # RBF (Gaussian) kernel: smooth, infinitely differentiable
            kernel = C(1.0, (1e-3, 1e3)) * RBF(
                length_scale=length_scale,
                length_scale_bounds=(1e-2, 1e2)
            )
        
        elif self.kernel_type == 'matern':
            # Matérn kernel: more flexible, controlled smoothness
            kernel = C(1.0, (1e-3, 1e3)) * Matern(
                length_scale=length_scale,
                length_scale_bounds=(1e-2, 1e2),
                nu=self.nu
            )
        
        elif self.kernel_type == 'auto':
            # Automatic kernel selection based on data
            # Use Matérn as default (good balance)
            kernel = C(1.0, (1e-3, 1e3)) * Matern(
                length_scale=length_scale,
                length_scale_bounds=(1e-2, 1e2),
                nu=1.5
            )
        
        else:
            raise ValueError(f"Unknown kernel type: {self.kernel_type}")
        
        # Add white noise kernel for numerical stability
        kernel = kernel + WhiteKernel(noise_level=self.alpha)
        
        return kernel
    
    def fit(
        self,
        coordinates: np.ndarray,
        values: np.ndarray,
        mask: Optional[np.ndarray] = None
    ):
        """
        Fit GP model to observed data.
        
        Parameters
        ----------
        coordinates : np.ndarray, shape (n_spots, 2)
            Spatial coordinates
        values : np.ndarray, shape (n_spots,) or (n_genes, n_spots)
            Observed values (can be 1D for single gene or 2D for multiple genes)
        mask : np.ndarray, optional
            Boolean mask indicating which values to use for training
            (True = use, False = ignore). If None, uses all non-zero values.
        
        Returns
        -------
        self
        """
        # Handle 2D input (multiple genes)
        if values.ndim == 2:
            # For now, fit on first gene (will be extended for batch processing)
            logger.warning("2D input detected, fitting on first gene only. Use fit_batch() for multiple genes.")
            values = values[0, :]
        
        # Create mask if not provided
        if mask is None:
            mask = values > 0
        
        # Filter to training data
        train_coords = coordinates[mask]
        train_values = values[mask]
        
        if len(train_values) == 0:
            raise ValueError("No training data available (all values are zero or masked)")
        
        # Estimate length scale from data if not provided
        if self.length_scale is None:
            # Use median pairwise distance as initial length scale
            if len(train_coords) > 1:
                dists = cdist(train_coords, train_coords)
                np.fill_diagonal(dists, np.inf)
                median_dist = np.median(dists[dists < np.inf])
                length_scale = median_dist
            else:
                length_scale = 1.0
        else:
            length_scale = self.length_scale
        
        # Create kernel
        self.kernel_ = self._create_kernel(length_scale)
        
        # Create and fit GP
        self.gp_ = GaussianProcessRegressor(
            kernel=self.kernel_,
            alpha=self.alpha,
            n_restarts_optimizer=self.n_restarts_optimizer,
            normalize_y=self.normalize_y,
            random_state=42
        )
        
        self.gp_.fit(train_coords, train_values)
        
        # Store fitted data
        self._fitted_coords = train_coords
        self._fitted_values = train_values
        
        logger.info(
            f"Fitted GP with {len(train_values)} training points, "
            f"kernel={self.kernel_type}, length_scale={self.gp_.kernel_.k1.k2.length_scale:.3f}"
        )
        
        return self
    
    def predict(
        self,
        coordinates: np.ndarray,
        return_std: bool = True,
        return_cov: bool = False
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        """
        Predict values at new coordinates.
        
        Parameters
        ----------
        coordinates : np.ndarray, shape (n_spots, 2)
            Coordinates for prediction
        return_std : bool, default=True
            Whether to return standard deviation (uncertainty)
        return_cov : bool, default=False
            Whether to return full covariance matrix
        
        Returns
        -------
        predictions : np.ndarray, shape (n_spots,)
            Predicted values
        uncertainty : np.ndarray, optional
            Standard deviation (if return_std=True) or covariance matrix (if return_cov=True)
        """
        if self.gp_ is None:
            raise ValueError("Model not fitted. Call fit() first.")
        
        if return_cov:
            predictions, cov = self.gp_.predict(coordinates, return_cov=True)
            return predictions, cov
        elif return_std:
            predictions, std = self.gp_.predict(coordinates, return_std=True)
            return predictions, std
        else:
            predictions = self.gp_.predict(coordinates, return_std=False)
            return predictions, None
    
    def impute(
        self,
        coordinates: np.ndarray,
        values: np.ndarray,
        mask: Optional[np.ndarray] = None,
        return_uncertainty: bool = True
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        """
        Impute missing values.
        
        Fits GP on observed values and predicts at all coordinates.
        
        Parameters
        ----------
        coordinates : np.ndarray, shape (n_spots, 2)
            Spatial coordinates
        values : np.ndarray, shape (n_spots,)
            Observed values (with zeros/missing)
        mask : np.ndarray, optional
            Boolean mask for training data
        return_uncertainty : bool, default=True
            Whether to return uncertainty estimates
        
        Returns
        -------
        imputed_values : np.ndarray, shape (n_spots,)
            Imputed values (including original observed values)
        uncertainty : np.ndarray, optional
            Uncertainty estimates
        """
        # Fit on observed data
        self.fit(coordinates, values, mask)
        
        # Predict at all coordinates
        imputed, uncertainty = self.predict(
            coordinates,
            return_std=return_uncertainty
        )
        
        # Replace observed values with original (don't impute where we have data)
        if mask is None:
            mask = values > 0
        imputed[mask] = values[mask]
        
        # Set uncertainty to zero for observed values
        if uncertainty is not None:
            uncertainty[mask] = 0.0
        
        return imputed, uncertainty
    
    def fit_batch(
        self,
        coordinates: np.ndarray,
        values: np.ndarray,
        mask: Optional[np.ndarray] = None,
        n_jobs: int = 1,
        verbose: bool = True
    ) -> 'GPImputerBatch':
        """
        Fit GP models for multiple genes in parallel.
        
        Parameters
        ----------
        coordinates : np.ndarray, shape (n_spots, 2)
            Spatial coordinates
        values : np.ndarray, shape (n_genes, n_spots)
            Values for multiple genes
        mask : np.ndarray, optional
            Boolean mask, shape (n_genes, n_spots)
        n_jobs : int, default=1
            Number of parallel jobs (-1 for all CPUs)
        verbose : bool, default=True
            Show progress bar
        
        Returns
        -------
        GPImputerBatch
            Batch imputer object with fitted models
        """
        return GPImputerBatch(
            base_imputer=self,
            coordinates=coordinates,
            values=values,
            mask=mask,
            n_jobs=n_jobs,
            verbose=verbose
        )
    
    def score(
        self,
        coordinates: np.ndarray,
        values: np.ndarray,
        mask: Optional[np.ndarray] = None
    ) -> float:
        """
        Compute log marginal likelihood of the model.
        
        Parameters
        ----------
        coordinates : np.ndarray
            Spatial coordinates
        values : np.ndarray
            True values
        mask : np.ndarray, optional
            Boolean mask for evaluation data
        
        Returns
        -------
        float
            Log marginal likelihood
        """
        if self.gp_ is None:
            raise ValueError("Model not fitted. Call fit() first.")
        
        if mask is None:
            mask = values > 0
        
        test_coords = coordinates[mask]
        test_values = values[mask]
        
        return self.gp_.score(test_coords, test_values)


class GPImputerBatch:
    """
    Batch GP imputer for multiple genes.
    
    Fits separate GP models for each gene in parallel.
    
    Parameters
    ----------
    base_imputer : GPImputer
        Base imputer configuration
    coordinates : np.ndarray
        Spatial coordinates
    values : np.ndarray, shape (n_genes, n_spots)
        Values for multiple genes
    mask : np.ndarray, optional
        Boolean mask
    n_jobs : int, default=1
        Number of parallel jobs
    verbose : bool, default=True
        Show progress bar
    """
    
    def __init__(
        self,
        base_imputer: GPImputer,
        coordinates: np.ndarray,
        values: np.ndarray,
        mask: Optional[np.ndarray] = None,
        n_jobs: int = 1,
        verbose: bool = True
    ):
        self.base_imputer = base_imputer
        self.coordinates = coordinates
        self.values = values
        self.mask = mask
        self.n_jobs = n_jobs
        self.verbose = verbose
        
        self.imputers_ = []
        self._fit_batch()
    
    def _fit_single_gene(self, gene_idx: int) -> Tuple[int, GPImputer]:
        """Fit GP for a single gene."""
        imputer = GPImputer(
            kernel_type=self.base_imputer.kernel_type,
            length_scale=self.base_imputer.length_scale,
            nu=self.base_imputer.nu,
            alpha=self.base_imputer.alpha,
            n_restarts_optimizer=self.base_imputer.n_restarts_optimizer,
            normalize_y=self.base_imputer.normalize_y
        )
        
        gene_values = self.values[gene_idx, :]
        gene_mask = self.mask[gene_idx, :] if self.mask is not None else None
        
        try:
            if threadpool_limits is None:
                imputer.fit(self.coordinates, gene_values, gene_mask)
            else:
                with threadpool_limits(limits=1):
                    imputer.fit(self.coordinates, gene_values, gene_mask)
            return gene_idx, imputer
        except Exception as e:
            logger.warning(f"Failed to fit gene {gene_idx}: {e}")
            return gene_idx, None
    
    def _fit_batch(self):
        """Fit all genes in parallel."""
        n_genes = self.values.shape[0]
        
        if self.n_jobs == 1:
            # Sequential processing
            iterator = range(n_genes)
            if self.verbose:
                iterator = tqdm(iterator, desc="Fitting GP models")
            
            for gene_idx in iterator:
                _, imputer = self._fit_single_gene(gene_idx)
                self.imputers_.append(imputer)
        
        else:
            # Parallel processing
            self.imputers_ = [None] * n_genes
            
            with ProcessPoolExecutor(max_workers=self.n_jobs) as executor:
                futures = {
                    executor.submit(self._fit_single_gene, i): i
                    for i in range(n_genes)
                }
                
                iterator = as_completed(futures)
                if self.verbose:
                    iterator = tqdm(iterator, total=n_genes, desc="Fitting GP models")
                
                for future in iterator:
                    gene_idx, imputer = future.result()
                    self.imputers_[gene_idx] = imputer
        
        n_success = sum(1 for imp in self.imputers_ if imp is not None)
        logger.info(f"Fitted {n_success}/{n_genes} GP models successfully")
    
    def predict(
        self,
        coordinates: Optional[np.ndarray] = None,
        return_std: bool = True
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        """
        Predict for all genes.
        
        Parameters
        ----------
        coordinates : np.ndarray, optional
            Coordinates for prediction. If None, uses training coordinates.
        return_std : bool, default=True
            Whether to return uncertainty
        
        Returns
        -------
        predictions : np.ndarray, shape (n_genes, n_spots)
            Predicted values
        uncertainty : np.ndarray, optional
            Uncertainty estimates
        """
        if coordinates is None:
            coordinates = self.coordinates
        
        n_genes = len(self.imputers_)
        n_spots = len(coordinates)
        
        predictions = np.zeros((n_genes, n_spots))
        uncertainty = np.zeros((n_genes, n_spots)) if return_std else None
        
        for i, imputer in enumerate(self.imputers_):
            if imputer is not None:
                pred, unc = imputer.predict(coordinates, return_std=return_std)
                predictions[i, :] = pred
                if uncertainty is not None:
                    uncertainty[i, :] = unc
        
        return predictions, uncertainty
    
    def impute(
        self,
        return_uncertainty: bool = True
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        """
        Impute all genes.
        
        Parameters
        ----------
        return_uncertainty : bool, default=True
            Whether to return uncertainty
        
        Returns
        -------
        imputed_values : np.ndarray, shape (n_genes, n_spots)
            Imputed values
        uncertainty : np.ndarray, optional
            Uncertainty estimates
        """
        predictions, uncertainty = self.predict(
            self.coordinates,
            return_std=return_uncertainty
        )
        
        # Replace observed values with original
        if self.mask is not None:
            predictions[self.mask] = self.values[self.mask]
            if uncertainty is not None:
                uncertainty[self.mask] = 0.0
        else:
            mask = self.values > 0
            predictions[mask] = self.values[mask]
            if uncertainty is not None:
                uncertainty[mask] = 0.0
        
        return predictions, uncertainty


def impute_spatial_apa(
    coordinates: np.ndarray,
    apa_counts: np.ndarray,
    kernel_type: str = 'matern',
    n_jobs: int = 1,
    return_uncertainty: bool = True
) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """
    Convenience function for spatial APA imputation.
    
    Parameters
    ----------
    coordinates : np.ndarray, shape (n_spots, 2)
        Spatial coordinates
    apa_counts : np.ndarray, shape (n_genes, n_spots)
        APA count matrix
    kernel_type : str, default='matern'
        Kernel type for GP
    n_jobs : int, default=1
        Number of parallel jobs
    return_uncertainty : bool, default=True
        Whether to return uncertainty
    
    Returns
    -------
    imputed_counts : np.ndarray
        Imputed APA counts
    uncertainty : np.ndarray, optional
        Uncertainty estimates
    """
    imputer = GPImputer(kernel_type=kernel_type)
    batch_imputer = imputer.fit_batch(
        coordinates,
        apa_counts,
        n_jobs=n_jobs,
        verbose=True
    )
    
    return batch_imputer.impute(return_uncertainty=return_uncertainty)
