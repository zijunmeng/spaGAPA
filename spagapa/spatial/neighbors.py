"""
Spatial neighbor finding and graph construction.

This module provides tools for finding spatial neighbors and constructing
spatial graphs from spatial coordinates.
"""

import numpy as np
from typing import Tuple, Optional, List
from scipy.spatial import distance_matrix, Delaunay
from sklearn.neighbors import NearestNeighbors, radius_neighbors_graph
import logging

logger = logging.getLogger(__name__)


class SpatialNeighbors:
    """
    Find spatial neighbors using various methods.
    
    Parameters
    ----------
    method : str, default='knn'
        Method for finding neighbors: 'knn', 'radius', or 'delaunay'
    n_neighbors : int, default=6
        Number of neighbors for KNN method
    radius : float, optional
        Radius for radius-based method
    """
    
    def __init__(
        self,
        method: str = 'knn',
        n_neighbors: int = 6,
        radius: Optional[float] = None
    ):
        self.method = method
        self.n_neighbors = n_neighbors
        self.radius = radius
        self._nbrs = None
    
    def fit(self, coordinates: np.ndarray):
        """
        Fit the neighbor finder to spatial coordinates.
        
        Parameters
        ----------
        coordinates : np.ndarray, shape (n_spots, 2)
            Spatial coordinates (x, y)
        """
        if self.method == 'knn':
            self._nbrs = NearestNeighbors(
                n_neighbors=self.n_neighbors + 1,  # +1 to exclude self
                algorithm='ball_tree'
            )
            self._nbrs.fit(coordinates)
            logger.info(f"Fitted KNN with k={self.n_neighbors}")
        
        elif self.method == 'radius':
            if self.radius is None:
                raise ValueError("radius must be specified for radius method")
            self._nbrs = NearestNeighbors(
                radius=self.radius,
                algorithm='ball_tree'
            )
            self._nbrs.fit(coordinates)
            logger.info(f"Fitted radius neighbors with r={self.radius}")
        
        elif self.method == 'delaunay':
            # Delaunay triangulation doesn't need fitting
            self._coordinates = coordinates
            logger.info("Using Delaunay triangulation")
        
        else:
            raise ValueError(f"Unknown method: {self.method}")
        
        return self
    
    def find_neighbors(
        self,
        coordinates: Optional[np.ndarray] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Find neighbors for each point.
        
        Parameters
        ----------
        coordinates : np.ndarray, optional
            Query coordinates. If None, use training coordinates.
        
        Returns
        -------
        distances : np.ndarray
            Distances to neighbors
        indices : np.ndarray
            Indices of neighbors
        """
        if self.method == 'knn':
            distances, indices = self._nbrs.kneighbors(coordinates)
            # Remove self (first neighbor)
            return distances[:, 1:], indices[:, 1:]
        
        elif self.method == 'radius':
            distances, indices = self._nbrs.radius_neighbors(coordinates)
            # Remove self from each neighbor list
            filtered_distances = []
            filtered_indices = []
            for i, (dists, inds) in enumerate(zip(distances, indices)):
                mask = inds != i
                filtered_distances.append(dists[mask])
                filtered_indices.append(inds[mask])
            return filtered_distances, filtered_indices
        
        elif self.method == 'delaunay':
            return self._delaunay_neighbors()
        
        else:
            raise ValueError(f"Unknown method: {self.method}")
    
    def _delaunay_neighbors(self) -> Tuple[List, List]:
        """
        Find neighbors using Delaunay triangulation.
        
        Returns
        -------
        distances : list of arrays
            Distances to neighbors for each point
        indices : list of arrays
            Indices of neighbors for each point
        """
        tri = Delaunay(self._coordinates)
        
        # Build neighbor lists from triangulation
        n_points = len(self._coordinates)
        neighbors = [set() for _ in range(n_points)]
        
        for simplex in tri.simplices:
            # Each simplex is a triangle, connect all vertices
            for i in range(3):
                for j in range(3):
                    if i != j:
                        neighbors[simplex[i]].add(simplex[j])
        
        # Convert to distances and indices
        distances = []
        indices = []
        for i, neighbor_set in enumerate(neighbors):
            neighbor_list = list(neighbor_set)
            if len(neighbor_list) > 0:
                dists = np.linalg.norm(
                    self._coordinates[neighbor_list] - self._coordinates[i],
                    axis=1
                )
                distances.append(dists)
                indices.append(np.array(neighbor_list))
            else:
                distances.append(np.array([]))
                indices.append(np.array([]))
        
        return distances, indices
    
    def compute_spatial_weights(
        self,
        coordinates: np.ndarray,
        weight_type: str = 'inverse_distance'
    ) -> np.ndarray:
        """
        Compute spatial weight matrix.
        
        Parameters
        ----------
        coordinates : np.ndarray
            Spatial coordinates
        weight_type : str, default='inverse_distance'
            Type of weights: 'inverse_distance', 'gaussian', or 'uniform'
        
        Returns
        -------
        np.ndarray, shape (n_spots, n_spots)
            Spatial weight matrix (sparse)
        """
        distances, indices = self.find_neighbors(coordinates)
        n_spots = len(coordinates)
        
        # Initialize weight matrix
        from scipy.sparse import lil_matrix
        W = lil_matrix((n_spots, n_spots))
        
        if self.method == 'knn':
            # Regular array format
            for i in range(n_spots):
                for j, dist in zip(indices[i], distances[i]):
                    if weight_type == 'inverse_distance':
                        W[i, j] = 1.0 / (dist + 1e-6)
                    elif weight_type == 'gaussian':
                        W[i, j] = np.exp(-dist**2 / (2 * self.radius**2)) if self.radius else np.exp(-dist**2)
                    elif weight_type == 'uniform':
                        W[i, j] = 1.0
        else:
            # List format (radius or delaunay)
            for i in range(n_spots):
                for j, dist in zip(indices[i], distances[i]):
                    if weight_type == 'inverse_distance':
                        W[i, j] = 1.0 / (dist + 1e-6)
                    elif weight_type == 'gaussian':
                        W[i, j] = np.exp(-dist**2 / (2 * self.radius**2)) if self.radius else np.exp(-dist**2)
                    elif weight_type == 'uniform':
                        W[i, j] = 1.0
        
        # Row normalize
        row_sums = np.array(W.sum(axis=1)).flatten()
        row_sums[row_sums == 0] = 1  # Avoid division by zero
        W = W.multiply(1.0 / row_sums[:, np.newaxis])
        
        logger.info(f"Computed spatial weights ({weight_type})")
        
        return W.tocsr()


def build_knn_graph(
    coordinates: np.ndarray,
    k: int = 6
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build KNN graph from spatial coordinates.
    
    Parameters
    ----------
    coordinates : np.ndarray, shape (n_spots, 2)
        Spatial coordinates
    k : int, default=6
        Number of neighbors
    
    Returns
    -------
    distances : np.ndarray
        Distance matrix
    indices : np.ndarray
        Neighbor indices
    """
    nbrs = SpatialNeighbors(method='knn', n_neighbors=k)
    nbrs.fit(coordinates)
    return nbrs.find_neighbors(coordinates)


def build_radius_graph(
    coordinates: np.ndarray,
    radius: float
) -> Tuple[List, List]:
    """
    Build radius graph from spatial coordinates.
    
    Parameters
    ----------
    coordinates : np.ndarray, shape (n_spots, 2)
        Spatial coordinates
    radius : float
        Radius for neighbor search
    
    Returns
    -------
    distances : list of arrays
        Distances to neighbors
    indices : list of arrays
        Neighbor indices
    """
    nbrs = SpatialNeighbors(method='radius', radius=radius)
    nbrs.fit(coordinates)
    return nbrs.find_neighbors(coordinates)


def build_delaunay_graph(
    coordinates: np.ndarray
) -> Tuple[List, List]:
    """
    Build Delaunay triangulation graph.
    
    Parameters
    ----------
    coordinates : np.ndarray, shape (n_spots, 2)
        Spatial coordinates
    
    Returns
    -------
    distances : list of arrays
        Distances to neighbors
    indices : list of arrays
        Neighbor indices
    """
    nbrs = SpatialNeighbors(method='delaunay')
    nbrs.fit(coordinates)
    return nbrs.find_neighbors()
