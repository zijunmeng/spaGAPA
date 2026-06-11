"""
APA quantification indices.

This module implements standard APA quantification metrics including:
- RUD (Relative Usage of Distal site)
- PDUI (Percentage of Distal Usage Index)
- WUL (Weighted 3' UTR Length)

These indices quantify alternative polyadenylation patterns from
spatial transcriptomics data.
"""

import numpy as np
import pandas as pd
from typing import Optional, Dict, Tuple, Literal
import logging

logger = logging.getLogger(__name__)


def calculate_rud(
    proximal_counts: np.ndarray,
    distal_counts: np.ndarray,
    pseudocount: float = 1.0
) -> np.ndarray:
    """
    Calculate RUD (Relative Usage of Distal site).
    
    RUD measures the relative usage of the distal polyadenylation site
    compared to the proximal site.
    
    Formula:
        RUD = distal / (proximal + distal + pseudocount)
    
    Parameters
    ----------
    proximal_counts : np.ndarray, shape (n_genes, n_spots) or (n_spots,)
        Read counts for proximal poly(A) sites
    distal_counts : np.ndarray, shape (n_genes, n_spots) or (n_spots,)
        Read counts for distal poly(A) sites
    pseudocount : float, default=1.0
        Pseudocount to avoid division by zero
    
    Returns
    -------
    np.ndarray
        RUD values (0-1), same shape as input
    
    Examples
    --------
    >>> proximal = np.array([10, 20, 5])
    >>> distal = np.array([30, 10, 15])
    >>> rud = calculate_rud(proximal, distal)
    >>> # rud ≈ [0.73, 0.32, 0.75]
    """
    if proximal_counts.shape != distal_counts.shape:
        raise ValueError("proximal_counts and distal_counts must have same shape")
    
    total = proximal_counts + distal_counts + pseudocount
    rud = distal_counts / total
    
    return rud


def calculate_pdui(
    proximal_counts: np.ndarray,
    distal_counts: np.ndarray,
    long_form_counts: Optional[np.ndarray] = None,
    pseudocount: float = 1.0
) -> np.ndarray:
    """
    Calculate PDUI (Percentage of Distal Usage Index).
    
    PDUI quantifies the percentage of transcripts using the distal
    polyadenylation site. Similar to RUD but expressed as percentage.
    
    Formula:
        PDUI = 100 * distal / (proximal + distal + pseudocount)
    
    If long_form_counts provided (reads spanning both sites):
        PDUI = 100 * (distal + long_form) / (proximal + distal + long_form + pseudocount)
    
    Parameters
    ----------
    proximal_counts : np.ndarray
        Read counts for proximal poly(A) sites
    distal_counts : np.ndarray
        Read counts for distal poly(A) sites
    long_form_counts : np.ndarray, optional
        Read counts for long-form transcripts (spanning both sites)
    pseudocount : float, default=1.0
        Pseudocount to avoid division by zero
    
    Returns
    -------
    np.ndarray
        PDUI values (0-100), same shape as input
    
    Examples
    --------
    >>> proximal = np.array([10, 20, 5])
    >>> distal = np.array([30, 10, 15])
    >>> pdui = calculate_pdui(proximal, distal)
    >>> # pdui ≈ [73.2, 32.3, 75.0]
    """
    if proximal_counts.shape != distal_counts.shape:
        raise ValueError("proximal_counts and distal_counts must have same shape")
    
    if long_form_counts is not None:
        if long_form_counts.shape != proximal_counts.shape:
            raise ValueError("long_form_counts must have same shape as other counts")
        total = proximal_counts + distal_counts + long_form_counts + pseudocount
        pdui = 100.0 * (distal_counts + long_form_counts) / total
    else:
        total = proximal_counts + distal_counts + pseudocount
        pdui = 100.0 * distal_counts / total
    
    return pdui


def calculate_wul(
    site_counts: Dict[str, np.ndarray],
    site_positions: Dict[str, int],
    normalize: bool = True,
    pseudocount: float = 1.0
) -> np.ndarray:
    """
    Calculate WUL (Weighted 3' UTR Length).
    
    WUL computes the weighted average 3' UTR length based on the usage
    of different polyadenylation sites.
    
    Formula:
        WUL = Σ(count_i * position_i) / Σ(count_i + pseudocount)
    
    If normalize=True:
        WUL_norm = (WUL - min_pos) / (max_pos - min_pos)
    
    Parameters
    ----------
    site_counts : dict
        Dictionary mapping site names to count arrays
        e.g., {'site1': array, 'site2': array}
    site_positions : dict
        Dictionary mapping site names to genomic positions
        e.g., {'site1': 1000, 'site2': 2000}
    normalize : bool, default=True
        If True, normalize WUL to [0, 1] range
    pseudocount : float, default=1.0
        Pseudocount for total counts
    
    Returns
    -------
    np.ndarray
        WUL values, shape (n_spots,)
    
    Examples
    --------
    >>> site_counts = {
    ...     'proximal': np.array([10, 20, 5]),
    ...     'distal': np.array([30, 10, 15])
    ... }
    >>> site_positions = {'proximal': 1000, 'distal': 2000}
    >>> wul = calculate_wul(site_counts, site_positions)
    """
    # Validate inputs
    if not site_counts or not site_positions:
        raise ValueError("site_counts and site_positions cannot be empty")
    
    if set(site_counts.keys()) != set(site_positions.keys()):
        raise ValueError("site_counts and site_positions must have same keys")
    
    # Get shape from first array
    first_key = list(site_counts.keys())[0]
    shape = site_counts[first_key].shape
    
    # Validate all arrays have same shape
    for key, counts in site_counts.items():
        if counts.shape != shape:
            raise ValueError(f"All count arrays must have same shape, got {counts.shape} for {key}")
    
    # Calculate weighted sum and total counts
    weighted_sum = np.zeros(shape)
    total_counts = np.zeros(shape)
    
    for site_name, counts in site_counts.items():
        position = site_positions[site_name]
        weighted_sum += counts * position
        total_counts += counts
    
    # Calculate WUL
    wul = weighted_sum / (total_counts + pseudocount)
    
    # Normalize if requested
    if normalize:
        min_pos = min(site_positions.values())
        max_pos = max(site_positions.values())
        
        if max_pos > min_pos:
            wul = (wul - min_pos) / (max_pos - min_pos)
        else:
            logger.warning("All sites at same position, cannot normalize WUL")
            wul = np.ones_like(wul) * 0.5
    
    return wul


def calculate_pai(
    proximal_counts: np.ndarray,
    distal_counts: np.ndarray,
    method: Literal['ratio', 'difference'] = 'ratio',
    pseudocount: float = 1.0
) -> np.ndarray:
    """
    Calculate PAI (Poly(A) site Index).
    
    PAI is another measure of APA site usage preference.
    
    Formula (ratio method):
        PAI = log2((distal + pseudocount) / (proximal + pseudocount))
    
    Formula (difference method):
        PAI = (distal - proximal) / (distal + proximal + pseudocount)
    
    Parameters
    ----------
    proximal_counts : np.ndarray
        Read counts for proximal poly(A) sites
    distal_counts : np.ndarray
        Read counts for distal poly(A) sites
    method : {'ratio', 'difference'}, default='ratio'
        Method for calculating PAI
    pseudocount : float, default=1.0
        Pseudocount to avoid division by zero
    
    Returns
    -------
    np.ndarray
        PAI values
    
    Examples
    --------
    >>> proximal = np.array([10, 20, 5])
    >>> distal = np.array([30, 10, 15])
    >>> pai = calculate_pai(proximal, distal, method='ratio')
    """
    if proximal_counts.shape != distal_counts.shape:
        raise ValueError("proximal_counts and distal_counts must have same shape")
    
    if method == 'ratio':
        # Log2 ratio
        pai = np.log2((distal_counts + pseudocount) / (proximal_counts + pseudocount))
    elif method == 'difference':
        # Normalized difference
        total = proximal_counts + distal_counts + pseudocount
        pai = (distal_counts - proximal_counts) / total
    else:
        raise ValueError(f"Unknown method: {method}. Use 'ratio' or 'difference'")
    
    return pai


def normalize_apa_index(
    index_values: np.ndarray,
    method: Literal['zscore', 'minmax', 'quantile'] = 'zscore',
    axis: Optional[int] = None
) -> np.ndarray:
    """
    Normalize APA index values.
    
    Parameters
    ----------
    index_values : np.ndarray
        APA index values to normalize
    method : {'zscore', 'minmax', 'quantile'}, default='zscore'
        Normalization method:
        - 'zscore': Z-score normalization (mean=0, std=1)
        - 'minmax': Min-max scaling to [0, 1]
        - 'quantile': Quantile normalization
    axis : int, optional
        Axis along which to normalize. If None, normalize entire array.
    
    Returns
    -------
    np.ndarray
        Normalized values
    
    Examples
    --------
    >>> values = np.array([10, 20, 30, 40, 50])
    >>> normalized = normalize_apa_index(values, method='zscore')
    """
    if method == 'zscore':
        mean = np.mean(index_values, axis=axis, keepdims=True)
        std = np.std(index_values, axis=axis, keepdims=True)
        normalized = (index_values - mean) / (std + 1e-10)
    
    elif method == 'minmax':
        min_val = np.min(index_values, axis=axis, keepdims=True)
        max_val = np.max(index_values, axis=axis, keepdims=True)
        normalized = (index_values - min_val) / (max_val - min_val + 1e-10)
    
    elif method == 'quantile':
        # Rank-based quantile normalization
        if axis is None:
            flat = index_values.flatten()
            ranks = np.argsort(np.argsort(flat))
            normalized = ranks.reshape(index_values.shape) / len(flat)
        else:
            # Per-axis quantile normalization
            normalized = np.apply_along_axis(
                lambda x: np.argsort(np.argsort(x)) / len(x),
                axis=axis,
                arr=index_values
            )
    
    else:
        raise ValueError(f"Unknown method: {method}")
    
    return normalized


class APAIndexCalculator:
    """
    Calculator for multiple APA indices.
    
    This class provides a unified interface for calculating various
    APA quantification indices from count data.
    
    Parameters
    ----------
    pseudocount : float, default=1.0
        Pseudocount for all calculations
    normalize : bool, default=False
        Whether to normalize indices
    normalization_method : str, default='zscore'
        Method for normalization
    
    Examples
    --------
    >>> calculator = APAIndexCalculator(pseudocount=1.0)
    >>> indices = calculator.calculate_all(proximal_counts, distal_counts)
    >>> # Returns dict with 'RUD', 'PDUI', 'PAI'
    """
    
    def __init__(
        self,
        pseudocount: float = 1.0,
        normalize: bool = False,
        normalization_method: Literal['zscore', 'minmax', 'quantile'] = 'zscore'
    ):
        self.pseudocount = pseudocount
        self.normalize = normalize
        self.normalization_method = normalization_method
    
    def calculate_all(
        self,
        proximal_counts: np.ndarray,
        distal_counts: np.ndarray,
        site_positions: Optional[Dict[str, int]] = None,
        long_form_counts: Optional[np.ndarray] = None
    ) -> Dict[str, np.ndarray]:
        """
        Calculate all APA indices.
        
        Parameters
        ----------
        proximal_counts : np.ndarray
            Proximal site counts
        distal_counts : np.ndarray
            Distal site counts
        site_positions : dict, optional
            Positions for WUL calculation
        long_form_counts : np.ndarray, optional
            Long-form transcript counts
        
        Returns
        -------
        dict
            Dictionary with keys: 'RUD', 'PDUI', 'PAI', 'WUL' (if positions provided)
        """
        indices = {}
        
        # Calculate RUD
        rud = calculate_rud(proximal_counts, distal_counts, self.pseudocount)
        indices['RUD'] = self._maybe_normalize(rud)
        
        # Calculate PDUI
        pdui = calculate_pdui(
            proximal_counts, distal_counts,
            long_form_counts, self.pseudocount
        )
        indices['PDUI'] = self._maybe_normalize(pdui)
        
        # Calculate PAI
        pai = calculate_pai(proximal_counts, distal_counts, 'ratio', self.pseudocount)
        indices['PAI'] = self._maybe_normalize(pai)
        
        # Calculate WUL if positions provided
        if site_positions is not None:
            site_counts = {
                'proximal': proximal_counts,
                'distal': distal_counts
            }
            wul = calculate_wul(site_counts, site_positions, True, self.pseudocount)
            indices['WUL'] = self._maybe_normalize(wul)
        
        logger.info(f"Calculated {len(indices)} APA indices")
        
        return indices
    
    def _maybe_normalize(self, values: np.ndarray) -> np.ndarray:
        """Apply normalization if enabled."""
        if self.normalize:
            return normalize_apa_index(values, self.normalization_method)
        return values
    
    def to_dataframe(
        self,
        indices: Dict[str, np.ndarray],
        gene_names: Optional[list] = None,
        spot_names: Optional[list] = None
    ) -> pd.DataFrame:
        """
        Convert indices to DataFrame.
        
        Parameters
        ----------
        indices : dict
            Dictionary of index arrays
        gene_names : list, optional
            Gene names for rows
        spot_names : list, optional
            Spot names for columns
        
        Returns
        -------
        pd.DataFrame
            DataFrame with indices
        """
        # Determine shape
        first_key = list(indices.keys())[0]
        shape = indices[first_key].shape
        
        if len(shape) == 1:
            # 1D array (single gene or aggregated)
            df = pd.DataFrame(indices)
            if spot_names is not None:
                df.index = spot_names
        else:
            # 2D array (genes x spots)
            # Create multi-level DataFrame
            dfs = []
            for index_name, values in indices.items():
                df_index = pd.DataFrame(
                    values,
                    index=gene_names,
                    columns=spot_names
                )
                df_index['index_type'] = index_name
                dfs.append(df_index)
            df = pd.concat(dfs)
        
        return df
