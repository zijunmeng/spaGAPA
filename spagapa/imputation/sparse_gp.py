"""
Sparse Gaussian Process approximation for large-scale spatial data.

This module implements sparse GP methods (inducing points) to handle
large spatial transcriptomics datasets efficiently.
"""

import numpy as np
from typing import Optional, Tuple, Literal
from scipy.spatial import distance_matrix
from scipy.linalg import cho_solve, cho_factor, solve_triangular
from sklearn.cluster import KMeans
import logging

logger = logging.getLogger(__name__)


class SparseGPImputer:
    """
    Sparse Gaussian Process imputer using inducing points.
    
    Uses a subset of inducing points to approximate the full GP,
    reducing computational complexity from O(n³) to O(nm²) where
    m << n is the number of inducing points.
    
    Parameters
    ----------
    n_inducing : int, default=100
        Number of inducing points
    inducing_method : str, default='kmeans'
        Method for selecting inducing points: 'kmeans', 'random', or 'grid'
    length_scale : float, default=1.0
        Length scale for RBF kernel
    noise_level : float, default=0.1
        Observation noise level
    
    Attributes
    ----------
    inducing_points_ : np.ndarray
        Selected inducing point locations
    alpha_ : np.ndarray
        Fitted coefficients
    
    Examples
    --------
    >>> imputer = SparseGPImputer(n_inducing=100)
    >>> imputer.fit(spatial_coords, apa_counts)
    >>> imputed = imputer.predict(spatial_coords)
    """
    
    def __init__(
        self,
        n_inducing: int = 100,
        inducing_method: Literal['kmeans', 'random', 'grid'] = 'kmeans',
        length_scale: float = 1.0,
        noise_level: float = 0.1,
        kernel_type: Literal['rbf'] = 'rbf',
    ):
        if kernel_type != 'rbf':
            raise ValueError("SparseGPImputer currently supports only kernel_type='rbf'")
        self.n_inducing = n_inducing
        self.inducing_method = inducing_method
        self.length_scale = length_scale
        self.noise_level = noise_level
        self.kernel_type = kernel_type
        
        self.inducing_points_ = None
        self.alpha_ = None
        self._K_mm_inv = None
        self._train_coords = None
        self._train_values = None
    
    def _rbf_kernel(
        self,
        X1: np.ndarray,
        X2: np.ndarray,
        length_scale: Optional[float] = None
    ) -> np.ndarray:
        """
        Compute RBF (Gaussian) kernel matrix.
        
        Parameters
        ----------
        X1 : np.ndarray, shape (n1, d)
            First set of points
        X2 : np.ndarray, shape (n2, d)
            Second set of points
        length_scale : float, optional
            Length scale parameter
        
        Returns
        -------
        np.ndarray, shape (n1, n2)
            Kernel matrix
        """
        if length_scale is None:
            length_scale = self.length_scale
        
        dists = distance_matrix(X1, X2)
        return np.exp(-0.5 * (dists / length_scale) ** 2)
    
    def _select_inducing_points(
        self,
        coordinates: np.ndarray
    ) -> np.ndarray:
        """
        Select inducing points from training data.
        
        Parameters
        ----------
        coordinates : np.ndarray, shape (n, 2)
            Training coordinates
        
        Returns
        -------
        np.ndarray, shape (m, 2)
            Inducing point locations
        """
        n_points = len(coordinates)
        n_inducing = min(self.n_inducing, n_points)
        
        if self.inducing_method == 'kmeans':
            # Use k-means clustering to find representative points
            kmeans = KMeans(n_clusters=n_inducing, random_state=42, n_init=10)
            kmeans.fit(coordinates)
            inducing_points = kmeans.cluster_centers_
            logger.info(f"Selected {n_inducing} inducing points using k-means")
        
        elif self.inducing_method == 'random':
            # Random subset of training points
            indices = np.random.choice(n_points, size=n_inducing, replace=False)
            inducing_points = coordinates[indices]
            logger.info(f"Selected {n_inducing} inducing points randomly")
        
        elif self.inducing_method == 'grid':
            # Regular grid covering the data range
            x_min, x_max = coordinates[:, 0].min(), coordinates[:, 0].max()
            y_min, y_max = coordinates[:, 1].min(), coordinates[:, 1].max()
            
            # Create grid
            n_per_dim = int(np.sqrt(n_inducing))
            x_grid = np.linspace(x_min, x_max, n_per_dim)
            y_grid = np.linspace(y_min, y_max, n_per_dim)
            xx, yy = np.meshgrid(x_grid, y_grid)
            inducing_points = np.column_stack([xx.ravel(), yy.ravel()])
            logger.info(f"Selected {len(inducing_points)} inducing points on grid")
        
        else:
            raise ValueError(f"Unknown inducing method: {self.inducing_method}")
        
        return inducing_points
    
    def fit(
        self,
        coordinates: np.ndarray,
        values: np.ndarray,
        mask: Optional[np.ndarray] = None,
        inducing_points: Optional[np.ndarray] = None,
        K_mm_inv: Optional[np.ndarray] = None,
    ):
        """
        Fit sparse GP model.

        Parameters
        ----------
        coordinates : np.ndarray, shape (n_spots, 2)
            Spatial coordinates
        values : np.ndarray, shape (n_spots,)
            Observed values
        mask : np.ndarray, optional
            Boolean mask for training data
        inducing_points : np.ndarray, optional, shape (m, 2)
            Precomputed inducing points. When supplied, the (expensive)
            KMeans selection is skipped. Inducing points depend only on
            coordinates, so a single set computed across all spots can be
            reused for every gene in a batch.
        K_mm_inv : np.ndarray, optional, shape (m, m)
            Precomputed inverse of the inducing-point kernel ``K_mm``.
            Must be consistent with ``inducing_points``. Skipping the
            per-gene ``np.linalg.inv`` is the second half of the reuse
            optimization.

        Returns
        -------
        self
        """
        # Filter training data
        if mask is None:
            mask = values > 0

        train_coords = coordinates[mask]
        train_values = values[mask]

        if len(train_values) == 0:
            raise ValueError("No training data available")

        # Select inducing points (or reuse a shared precomputed set)
        if inducing_points is not None:
            self.inducing_points_ = inducing_points
        else:
            self.inducing_points_ = self._select_inducing_points(train_coords)
        m = len(self.inducing_points_)
        n = len(train_coords)

        # K_nm: kernel between training points and inducing points (n x m)
        K_nm = self._rbf_kernel(train_coords, self.inducing_points_)

        # Compute K_mm and K_mm^{-1} (or reuse a shared precomputed inverse).
        # K_mm is still needed below for Sigma = K_mm + K_mn @ Lambda^{-1} @ K_nm.
        if K_mm_inv is not None:
            self._K_mm_inv = K_mm_inv
            K_mm = self._rbf_kernel(self.inducing_points_, self.inducing_points_)
            K_mm += 1e-6 * np.eye(m)  # Jitter, consistent with the precompute path
        else:
            K_mm = self._rbf_kernel(self.inducing_points_, self.inducing_points_)
            K_mm += 1e-6 * np.eye(m)  # Jitter for numerical stability
            self._K_mm_inv = np.linalg.inv(K_mm)

        # Local alias: the math below uses K_mm_inv in two places.
        K_mm_inv = self._K_mm_inv
        
        # Compute Lambda (diagonal correction term)
        # Lambda = diag(K_nn - Q_nn) + noise
        K_nn_diag = np.ones(n)  # RBF kernel diagonal is 1
        Q_nn_diag = np.sum(K_nm @ K_mm_inv * K_nm, axis=1)
        Lambda = K_nn_diag - Q_nn_diag + self.noise_level
        
        # Compute Sigma = K_mm + K_mn @ Lambda^{-1} @ K_nm
        K_mn = K_nm.T
        Sigma = K_mm + K_mn @ np.diag(1.0 / Lambda) @ K_nm
        
        # Cholesky decomposition for efficient solving
        L = cho_factor(Sigma, lower=True)
        
        # Compute alpha = Sigma^{-1} @ K_mn @ Lambda^{-1} @ y
        temp = K_mn @ (train_values / Lambda)
        self.alpha_ = cho_solve(L, temp)
        
        # Store training data
        self._train_coords = train_coords
        self._train_values = train_values
        
        logger.info(
            f"Fitted sparse GP with {n} training points, "
            f"{m} inducing points (compression ratio: {n/m:.1f}x)"
        )
        
        return self
    
    def predict(
        self,
        coordinates: np.ndarray,
        return_std: bool = False
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        """
        Predict at new coordinates.
        
        Parameters
        ----------
        coordinates : np.ndarray, shape (n_test, 2)
            Test coordinates
        return_std : bool, default=False
            Whether to return standard deviation
        
        Returns
        -------
        predictions : np.ndarray, shape (n_test,)
            Predicted values
        std : np.ndarray, optional
            Standard deviations
        """
        if self.inducing_points_ is None or self.alpha_ is None:
            raise ValueError("Model not fitted. Call fit() first.")
        
        # K_*m: kernel between test points and inducing points
        K_star_m = self._rbf_kernel(coordinates, self.inducing_points_)
        
        # Predictions: K_*m @ alpha
        predictions = K_star_m @ self.alpha_
        
        if return_std:
            # Compute predictive variance
            # Var = K_** - K_*m @ K_mm^{-1} @ K_m*
            K_star_star_diag = np.ones(len(coordinates))  # RBF diagonal
            Q_star_star_diag = np.sum(K_star_m @ self._K_mm_inv * K_star_m, axis=1)
            variance = K_star_star_diag - Q_star_star_diag + self.noise_level
            std = np.sqrt(np.maximum(variance, 0))  # Ensure non-negative
            
            return predictions, std
        
        return predictions, None
    
    def impute(
        self,
        coordinates: np.ndarray,
        values: np.ndarray,
        mask: Optional[np.ndarray] = None,
        return_uncertainty: bool = False
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        """
        Impute missing values.
        
        Parameters
        ----------
        coordinates : np.ndarray
            Spatial coordinates
        values : np.ndarray
            Observed values
        mask : np.ndarray, optional
            Training mask
        return_uncertainty : bool, default=False
            Whether to return uncertainty
        
        Returns
        -------
        imputed_values : np.ndarray
            Imputed values
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
        
        # Keep observed values
        if mask is None:
            mask = values > 0
        imputed[mask] = values[mask]
        
        if uncertainty is not None:
            uncertainty[mask] = 0.0
        
        return imputed, uncertainty

    def fit_batch(
        self,
        coordinates: np.ndarray,
        values: np.ndarray,
        mask: Optional[np.ndarray] = None,
        n_jobs: int = 1,
        verbose: bool = True,
    ) -> 'SparseGPImputerBatch':
        """
        Fit sparse GP models for multiple genes.

        Parameters
        ----------
        coordinates : np.ndarray, shape (n_spots, 2)
            Spatial coordinates.
        values : np.ndarray, shape (n_genes, n_spots)
            APA values for multiple genes.
        mask : np.ndarray, optional, shape (n_genes, n_spots)
            Boolean training mask.
        n_jobs : int, default=1
            Reserved for API compatibility. Sparse batch currently runs
            sequentially to keep memory use predictable.
        verbose : bool, default=True
            Print progress every 100 genes.

        Returns
        -------
        SparseGPImputerBatch
            Fitted batch imputer.
        """
        return SparseGPImputerBatch(
            base_imputer=self,
            coordinates=coordinates,
            values=values,
            mask=mask,
            verbose=verbose,
        )


class SparseGPImputerBatch:
    """Batch wrapper with the same public shape contract as GPImputerBatch."""

    def __init__(
        self,
        base_imputer: SparseGPImputer,
        coordinates: np.ndarray,
        values: np.ndarray,
        mask: Optional[np.ndarray] = None,
        verbose: bool = True,
    ):
        if values.ndim != 2:
            raise ValueError("values must have shape (n_genes, n_spots)")
        if values.shape[1] != len(coordinates):
            raise ValueError(
                "values must have shape (n_genes, n_spots) with n_spots "
                "matching coordinates"
            )
        if mask is not None and mask.shape != values.shape:
            raise ValueError("mask must have the same shape as values")

        self.base_imputer = base_imputer
        self.coordinates = coordinates
        self.values = values
        self.mask = mask
        self.verbose = verbose
        self.imputers_ = []
        self._fit_batch()

    def _new_imputer(self) -> SparseGPImputer:
        return SparseGPImputer(
            n_inducing=self.base_imputer.n_inducing,
            inducing_method=self.base_imputer.inducing_method,
            length_scale=self.base_imputer.length_scale,
            noise_level=self.base_imputer.noise_level,
            kernel_type=self.base_imputer.kernel_type,
        )

    def _fit_batch(self):
        n_genes = self.values.shape[0]

        # Precompute inducing points ONCE on all coordinates (not the
        # gene-specific training mask). Inducing points depend only on the
        # spatial layout, which is identical for every gene, so running
        # KMeans(42k points, ~424 clusters, n_init=10) per gene is pure waste.
        # Reusing the same set + precomputed K_mm^{-1} turns an O(n_genes)
        # KMeans cost into O(1) and is the bulk of the speedup.
        base = self.base_imputer
        shared_inducing = base._select_inducing_points(self.coordinates)
        K_mm = base._rbf_kernel(shared_inducing, shared_inducing)
        K_mm += 1e-6 * np.eye(len(shared_inducing))
        shared_K_mm_inv = np.linalg.inv(K_mm)
        logger.info(
            f"Precomputed {len(shared_inducing)} shared inducing points "
            f"for {n_genes} genes"
        )

        for gene_idx in range(n_genes):
            if self.verbose and gene_idx > 0 and gene_idx % 100 == 0:
                logger.info(f"Fitted {gene_idx}/{n_genes} sparse GP models")

            imputer = self._new_imputer()
            gene_values = self.values[gene_idx, :]
            gene_mask = self.mask[gene_idx, :] if self.mask is not None else None
            try:
                imputer.fit(
                    self.coordinates, gene_values, gene_mask,
                    inducing_points=shared_inducing,
                    K_mm_inv=shared_K_mm_inv,
                )
                self.imputers_.append(imputer)
            except Exception as exc:
                logger.warning(f"Failed to fit sparse GP for gene {gene_idx}: {exc}")
                self.imputers_.append(None)

        n_success = sum(1 for imp in self.imputers_ if imp is not None)
        logger.info(f"Fitted {n_success}/{n_genes} sparse GP models successfully")

    def predict(
        self,
        coordinates: Optional[np.ndarray] = None,
        return_std: bool = True,
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        if coordinates is None:
            coordinates = self.coordinates

        n_genes = len(self.imputers_)
        n_spots = len(coordinates)
        predictions = np.zeros((n_genes, n_spots))
        uncertainty = np.zeros((n_genes, n_spots)) if return_std else None

        for gene_idx, imputer in enumerate(self.imputers_):
            if imputer is None:
                continue
            pred, unc = imputer.predict(coordinates, return_std=return_std)
            predictions[gene_idx, :] = pred
            if uncertainty is not None:
                uncertainty[gene_idx, :] = unc

        return predictions, uncertainty

    def impute(
        self,
        return_uncertainty: bool = True,
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        predictions, uncertainty = self.predict(
            self.coordinates,
            return_std=return_uncertainty,
        )

        if self.mask is not None:
            observed = self.mask
        else:
            observed = self.values > 0

        predictions[observed] = self.values[observed]
        if uncertainty is not None:
            uncertainty[observed] = 0.0

        return predictions, uncertainty


class BlockGPImputer:
    """
    Block-wise GP imputation for very large datasets.
    
    Divides spatial domain into blocks and fits separate GP models
    for each block, with overlap for smooth transitions.
    
    Parameters
    ----------
    block_size : int, default=500
        Maximum number of points per block
    overlap : float, default=0.1
        Overlap fraction between blocks (0-1)
    base_imputer : str, default='sparse'
        Base imputer type: 'sparse' or 'full'
    **imputer_kwargs
        Additional arguments for base imputer
    
    Examples
    --------
    >>> imputer = BlockGPImputer(block_size=500, overlap=0.1)
    >>> imputed = imputer.impute(spatial_coords, apa_counts)
    """
    
    def __init__(
        self,
        block_size: int = 500,
        overlap: float = 0.1,
        base_imputer: str = 'sparse',
        **imputer_kwargs
    ):
        self.block_size = block_size
        self.overlap = overlap
        self.base_imputer = base_imputer
        self.imputer_kwargs = imputer_kwargs
        
        self.blocks_ = []
        self.imputers_ = []
    
    def _create_blocks(
        self,
        coordinates: np.ndarray
    ) -> list:
        """
        Divide coordinates into overlapping blocks.
        
        Parameters
        ----------
        coordinates : np.ndarray
            Spatial coordinates
        
        Returns
        -------
        list
            List of block indices
        """
        n_points = len(coordinates)
        
        if n_points <= self.block_size:
            # No need for blocking
            return [np.arange(n_points)]
        
        # Compute spatial extent
        x_min, x_max = coordinates[:, 0].min(), coordinates[:, 0].max()
        y_min, y_max = coordinates[:, 1].min(), coordinates[:, 1].max()
        
        # Estimate number of blocks needed
        area = (x_max - x_min) * (y_max - y_min)
        points_per_area = n_points / area
        block_area = self.block_size / points_per_area
        block_width = np.sqrt(block_area)
        
        # Create grid of blocks with overlap
        overlap_width = block_width * self.overlap
        
        x_edges = np.arange(x_min, x_max + block_width, block_width - overlap_width)
        y_edges = np.arange(y_min, y_max + block_width, block_width - overlap_width)
        
        blocks = []
        for i in range(len(x_edges) - 1):
            for j in range(len(y_edges) - 1):
                # Find points in this block
                mask = (
                    (coordinates[:, 0] >= x_edges[i]) &
                    (coordinates[:, 0] < x_edges[i+1]) &
                    (coordinates[:, 1] >= y_edges[j]) &
                    (coordinates[:, 1] < y_edges[j+1])
                )
                indices = np.where(mask)[0]
                
                if len(indices) > 0:
                    blocks.append(indices)
        
        logger.info(f"Created {len(blocks)} blocks for {n_points} points")
        
        return blocks
    
    def impute(
        self,
        coordinates: np.ndarray,
        values: np.ndarray,
        mask: Optional[np.ndarray] = None,
        return_uncertainty: bool = False
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        """
        Impute using block-wise processing.
        
        Parameters
        ----------
        coordinates : np.ndarray
            Spatial coordinates
        values : np.ndarray
            Observed values
        mask : np.ndarray, optional
            Training mask
        return_uncertainty : bool, default=False
            Whether to return uncertainty
        
        Returns
        -------
        imputed_values : np.ndarray
            Imputed values
        uncertainty : np.ndarray, optional
            Uncertainty estimates
        """
        # Create blocks
        self.blocks_ = self._create_blocks(coordinates)
        
        # Initialize output
        imputed = np.zeros_like(values)
        uncertainty = np.zeros_like(values) if return_uncertainty else None
        weights = np.zeros_like(values)
        
        # Process each block
        for block_idx, indices in enumerate(self.blocks_):
            block_coords = coordinates[indices]
            block_values = values[indices]
            block_mask = mask[indices] if mask is not None else None
            
            # Create imputer for this block
            if self.base_imputer == 'sparse':
                imputer = SparseGPImputer(**self.imputer_kwargs)
            else:
                from .gp_imputer import GPImputer
                imputer = GPImputer(**self.imputer_kwargs)
            
            # Impute block
            try:
                block_imputed, block_unc = imputer.impute(
                    block_coords,
                    block_values,
                    block_mask,
                    return_uncertainty=return_uncertainty
                )
                
                # Accumulate results (weighted average in overlap regions)
                imputed[indices] += block_imputed
                weights[indices] += 1
                
                if uncertainty is not None and block_unc is not None:
                    uncertainty[indices] += block_unc
            
            except Exception as e:
                logger.warning(f"Failed to impute block {block_idx}: {e}")
                # Use original values for failed blocks
                imputed[indices] += block_values
                weights[indices] += 1
        
        # Average overlapping regions
        weights[weights == 0] = 1  # Avoid division by zero
        imputed /= weights
        
        if uncertainty is not None:
            uncertainty /= weights
        
        logger.info(f"Completed block-wise imputation with {len(self.blocks_)} blocks")
        
        return imputed, uncertainty
