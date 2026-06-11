"""
Spatial trajectory analysis for APA dynamics.

Tasks 3.6–3.9: Manual trajectory specification, APA curve fitting,
switch point detection, and trajectory visualization.
"""

import numpy as np
import pandas as pd
from typing import Optional, Tuple, List, Dict
from scipy.spatial.distance import cdist
from scipy.signal import find_peaks
import warnings


class TrajectoryBuilder:
    """
    Build spatial trajectories between user-specified points.

    Parameters
    ----------
    n_neighbors : int, default=6
        Neighbours used for shortest-path graph
    """

    def __init__(self, n_neighbors: int = 6):
        self.n_neighbors = n_neighbors

    def build_trajectory(
        self,
        spatial_coords: np.ndarray,
        start_point: np.ndarray,
        end_point: np.ndarray
    ) -> np.ndarray:
        """
        Find the shortest path from start to end through the spatial graph.

        Parameters
        ----------
        spatial_coords : np.ndarray, shape (n_spots, 2)
        start_point : np.ndarray, shape (2,)  — x,y coordinate
        end_point   : np.ndarray, shape (2,)  — x,y coordinate

        Returns
        -------
        path_indices : np.ndarray
            Ordered spot indices along the trajectory
        """
        # Find nearest spots to start/end
        dists_start = np.linalg.norm(spatial_coords - start_point, axis=1)
        dists_end   = np.linalg.norm(spatial_coords - end_point,   axis=1)
        start_idx = int(np.argmin(dists_start))
        end_idx   = int(np.argmin(dists_end))

        if start_idx == end_idx:
            return np.array([start_idx])

        # Build KNN adjacency
        from sklearn.neighbors import NearestNeighbors
        nbrs = NearestNeighbors(n_neighbors=self.n_neighbors + 1).fit(spatial_coords)
        distances, indices = nbrs.kneighbors(spatial_coords)

        # Dijkstra shortest path
        path = self._dijkstra(start_idx, end_idx, indices[:, 1:], distances[:, 1:],
                              len(spatial_coords))
        return np.array(path)

    def _dijkstra(self, src, dst, neighbor_idx, neighbor_dist, n):
        """Simple Dijkstra on KNN graph."""
        import heapq
        dist = np.full(n, np.inf)
        prev = np.full(n, -1, dtype=int)
        dist[src] = 0.0
        heap = [(0.0, src)]

        while heap:
            d, u = heapq.heappop(heap)
            if d > dist[u]:
                continue
            if u == dst:
                break
            for v, w in zip(neighbor_idx[u], neighbor_dist[u]):
                nd = dist[u] + w
                if nd < dist[v]:
                    dist[v] = nd
                    prev[v] = u
                    heapq.heappush(heap, (nd, v))

        # Reconstruct path
        path = []
        cur = dst
        while cur != -1:
            path.append(cur)
            cur = prev[cur]
        path.reverse()
        return path if path[0] == src else []

    def order_spots_along_path(
        self,
        path_indices: np.ndarray,
        spatial_coords: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute cumulative distance (pseudospace) along path.

        Returns
        -------
        path_indices : np.ndarray
        pseudospace  : np.ndarray  — cumulative distance, normalised to [0,1]
        """
        coords_path = spatial_coords[path_indices]
        diffs = np.diff(coords_path, axis=0)
        step_dists = np.linalg.norm(diffs, axis=1)
        cumulative = np.concatenate([[0], np.cumsum(step_dists)])
        if cumulative[-1] > 0:
            cumulative = cumulative / cumulative[-1]
        return path_indices, cumulative

    def compute_trajectory_distance(self, path_indices, spatial_coords):
        """Total Euclidean length of the trajectory."""
        coords = spatial_coords[path_indices]
        return float(np.sum(np.linalg.norm(np.diff(coords, axis=0), axis=1)))

    def validate_trajectory(self, path_indices, spatial_coords, min_length=5):
        """Check trajectory is valid (non-empty, sufficient length)."""
        if len(path_indices) < min_length:
            warnings.warn(f"Trajectory has only {len(path_indices)} spots (< {min_length})")
            return False
        return True


class TrajectoryAnalyzer:
    """
    Analyse APA dynamics along a spatial trajectory.

    Parameters
    ----------
    method : str, default='spline'
        Curve fitting method: 'spline' or 'lowess'
    n_bins : int, default=20
        Number of pseudospace bins for smoothing
    """

    def __init__(self, method: str = 'spline', n_bins: int = 20):
        self.method = method
        self.n_bins = n_bins

    def fit_apa_curves(
        self,
        apa_matrix: np.ndarray,
        path_indices: np.ndarray,
        pseudospace: np.ndarray,
        gene_names: Optional[List[str]] = None
    ) -> pd.DataFrame:
        """
        Fit smooth APA curves along the trajectory.

        Parameters
        ----------
        apa_matrix  : np.ndarray, shape (n_genes, n_spots)
        path_indices: np.ndarray  — spot indices along trajectory
        pseudospace : np.ndarray  — normalised position [0,1] per spot
        gene_names  : list of str, optional

        Returns
        -------
        curves_df : pd.DataFrame
            Columns = gene names, rows = pseudospace bins
        """
        n_genes = apa_matrix.shape[0]
        if gene_names is None:
            gene_names = [f'Gene_{i}' for i in range(n_genes)]

        # Extract APA values along path
        apa_path = apa_matrix[:, path_indices]   # (n_genes, n_path)

        # Bin pseudospace
        bins = np.linspace(0, 1, self.n_bins + 1)
        bin_centers = 0.5 * (bins[:-1] + bins[1:])
        bin_idx = np.digitize(pseudospace, bins) - 1
        bin_idx = np.clip(bin_idx, 0, self.n_bins - 1)

        curves = np.full((n_genes, self.n_bins), np.nan)
        for b in range(self.n_bins):
            mask = bin_idx == b
            if mask.sum() > 0:
                curves[:, b] = np.nanmean(apa_path[:, mask], axis=1)

        # Smooth with spline or lowess
        for i in range(n_genes):
            valid = ~np.isnan(curves[i])
            if valid.sum() >= 4:
                curves[i] = self._smooth(bin_centers, curves[i], valid)

        return pd.DataFrame(
            curves.T,
            index=bin_centers,
            columns=gene_names
        )

    def _smooth(self, x, y, valid_mask):
        """Apply smoothing to a curve."""
        x_v = x[valid_mask]
        y_v = y[valid_mask]
        smoothed = y.copy()

        if self.method == 'spline':
            try:
                from scipy.interpolate import UnivariateSpline
                spl = UnivariateSpline(x_v, y_v, k=min(3, len(x_v) - 1), s=0.1)
                smoothed = spl(x)
            except Exception:
                smoothed[valid_mask] = y_v
        elif self.method == 'lowess':
            try:
                from statsmodels.nonparametric.smoothers_lowess import lowess
                result = lowess(y_v, x_v, frac=0.4, return_sorted=True)
                from scipy.interpolate import interp1d
                f = interp1d(result[:, 0], result[:, 1],
                             bounds_error=False, fill_value='extrapolate')
                smoothed = f(x)
            except Exception:
                smoothed[valid_mask] = y_v

        return np.clip(smoothed, 0, 1)

    def detect_switch_points(
        self,
        curves_df: pd.DataFrame,
        threshold: float = 0.05
    ) -> Dict[str, List[float]]:
        """
        Detect APA switch points (inflection points) along trajectory.

        Parameters
        ----------
        curves_df : pd.DataFrame
            Output of fit_apa_curves
        threshold : float
            Minimum absolute derivative to call a switch

        Returns
        -------
        dict : {gene_name: [pseudospace positions of switches]}
        """
        x = curves_df.index.values
        switches = {}

        for gene in curves_df.columns:
            y = curves_df[gene].values
            if np.all(np.isnan(y)):
                continue

            # First derivative
            dy = np.gradient(y, x)

            # Find sign changes in derivative (inflection points)
            sign_changes = np.where(np.diff(np.sign(dy)))[0]

            # Filter by magnitude
            strong = [i for i in sign_changes if abs(dy[i]) > threshold]
            switches[gene] = [float(x[i]) for i in strong]

        return switches

    def rank_dynamic_genes(
        self,
        curves_df: pd.DataFrame,
        metric: str = 'variance'
    ) -> pd.DataFrame:
        """
        Rank genes by their dynamics along the trajectory.

        Parameters
        ----------
        curves_df : pd.DataFrame
        metric : str
            'variance' | 'range' | 'n_switches'

        Returns
        -------
        pd.DataFrame with columns: gene, score, rank
        """
        scores = {}
        x = curves_df.index.values

        for gene in curves_df.columns:
            y = curves_df[gene].dropna().values
            if len(y) < 2:
                scores[gene] = 0.0
                continue

            if metric == 'variance':
                scores[gene] = float(np.var(y))
            elif metric == 'range':
                scores[gene] = float(np.ptp(y))
            elif metric == 'n_switches':
                dy = np.gradient(y)
                scores[gene] = float(len(np.where(np.diff(np.sign(dy)))[0]))
            else:
                raise ValueError(f"Unknown metric: {metric}")

        df = pd.DataFrame({'gene': list(scores.keys()),
                           'score': list(scores.values())})
        df = df.sort_values('score', ascending=False).reset_index(drop=True)
        df['rank'] = df.index + 1
        return df

    def classify_trajectory_patterns(
        self,
        curves_df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Classify each gene's trajectory pattern.

        Patterns: 'increasing', 'decreasing', 'switch', 'stable'

        Returns
        -------
        pd.DataFrame with columns: gene, pattern, score
        """
        rows = []
        for gene in curves_df.columns:
            y = curves_df[gene].dropna().values
            if len(y) < 3:
                rows.append({'gene': gene, 'pattern': 'stable', 'score': 0.0})
                continue

            # Pearson correlation with linear trend
            x = np.linspace(0, 1, len(y))
            corr = float(np.corrcoef(x, y)[0, 1])

            # Number of direction changes
            dy = np.gradient(y)
            n_changes = len(np.where(np.diff(np.sign(dy)))[0])

            if n_changes >= 2:
                pattern = 'switch'
                score = float(np.var(y))
            elif corr > 0.5:
                pattern = 'increasing'
                score = corr
            elif corr < -0.5:
                pattern = 'decreasing'
                score = -corr
            else:
                pattern = 'stable'
                score = 1.0 - abs(corr)

            rows.append({'gene': gene, 'pattern': pattern, 'score': score})

        return pd.DataFrame(rows).sort_values('score', ascending=False).reset_index(drop=True)


# ── Convenience functions ──────────────────────────────────────────────────────

def build_spatial_trajectory(
    spatial_coords: np.ndarray,
    start_point: np.ndarray,
    end_point: np.ndarray,
    n_neighbors: int = 6
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build a spatial trajectory and return (path_indices, pseudospace).

    Parameters
    ----------
    spatial_coords : np.ndarray, shape (n_spots, 2)
    start_point    : np.ndarray, shape (2,)
    end_point      : np.ndarray, shape (2,)
    n_neighbors    : int

    Returns
    -------
    path_indices : np.ndarray
    pseudospace  : np.ndarray  — normalised [0,1]

    Examples
    --------
    >>> coords = np.random.rand(200, 2) * 100
    >>> path, ps = build_spatial_trajectory(coords, coords[0], coords[-1])
    >>> print(f"Path length: {len(path)} spots")
    """
    builder = TrajectoryBuilder(n_neighbors=n_neighbors)
    path = builder.build_trajectory(spatial_coords, start_point, end_point)
    path, pseudospace = builder.order_spots_along_path(path, spatial_coords)
    return path, pseudospace


def analyse_apa_trajectory(
    apa_matrix: np.ndarray,
    spatial_coords: np.ndarray,
    start_point: np.ndarray,
    end_point: np.ndarray,
    gene_names: Optional[List[str]] = None,
    method: str = 'spline',
    n_bins: int = 20
) -> Dict:
    """
    Full trajectory analysis pipeline.

    Parameters
    ----------
    apa_matrix     : np.ndarray, shape (n_genes, n_spots)
    spatial_coords : np.ndarray, shape (n_spots, 2)
    start_point    : np.ndarray, shape (2,)
    end_point      : np.ndarray, shape (2,)
    gene_names     : list of str, optional
    method         : str  — 'spline' or 'lowess'
    n_bins         : int

    Returns
    -------
    dict with keys:
        'path_indices', 'pseudospace', 'curves', 'switches',
        'ranked_genes', 'patterns'

    Examples
    --------
    >>> import numpy as np
    >>> from spagapa.analysis.trajectory import analyse_apa_trajectory
    >>> apa = np.random.rand(50, 200)
    >>> coords = np.random.rand(200, 2) * 100
    >>> result = analyse_apa_trajectory(apa, coords, coords[0], coords[-1])
    >>> print(result['patterns'].head())
    """
    path, pseudospace = build_spatial_trajectory(
        spatial_coords, start_point, end_point
    )

    analyzer = TrajectoryAnalyzer(method=method, n_bins=n_bins)
    curves   = analyzer.fit_apa_curves(apa_matrix, path, pseudospace, gene_names)
    switches = analyzer.detect_switch_points(curves)
    ranked   = analyzer.rank_dynamic_genes(curves, metric='variance')
    patterns = analyzer.classify_trajectory_patterns(curves)

    return {
        'path_indices': path,
        'pseudospace':  pseudospace,
        'curves':       curves,
        'switches':     switches,
        'ranked_genes': ranked,
        'patterns':     patterns,
    }
