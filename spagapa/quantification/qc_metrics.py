"""
Quality control metrics for APA quantification and imputation.

This module implements cross-validation, quality metrics, and QC report
generation for evaluating APA quantification and imputation quality.
"""

import numpy as np
import pandas as pd
from typing import Optional, Dict, Tuple, List, Callable
from sklearn.model_selection import KFold
from sklearn.metrics import mean_squared_error, mean_absolute_error
from scipy.stats import pearsonr, spearmanr
import logging

logger = logging.getLogger(__name__)


def calculate_rmse(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    mask: Optional[np.ndarray] = None
) -> float:
    """
    Calculate Root Mean Squared Error.
    
    Parameters
    ----------
    y_true : np.ndarray
        True values
    y_pred : np.ndarray
        Predicted values
    mask : np.ndarray, optional
        Boolean mask for values to include
    
    Returns
    -------
    float
        RMSE value
    """
    if mask is not None:
        y_true = y_true[mask]
        y_pred = y_pred[mask]
    
    return np.sqrt(mean_squared_error(y_true, y_pred))


def calculate_mae(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    mask: Optional[np.ndarray] = None
) -> float:
    """
    Calculate Mean Absolute Error.
    
    Parameters
    ----------
    y_true : np.ndarray
        True values
    y_pred : np.ndarray
        Predicted values
    mask : np.ndarray, optional
        Boolean mask for values to include
    
    Returns
    -------
    float
        MAE value
    """
    if mask is not None:
        y_true = y_true[mask]
        y_pred = y_pred[mask]
    
    return mean_absolute_error(y_true, y_pred)


def calculate_pearson(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    mask: Optional[np.ndarray] = None
) -> Tuple[float, float]:
    """
    Calculate Pearson correlation coefficient.
    
    Parameters
    ----------
    y_true : np.ndarray
        True values
    y_pred : np.ndarray
        Predicted values
    mask : np.ndarray, optional
        Boolean mask for values to include
    
    Returns
    -------
    correlation : float
        Pearson correlation coefficient
    pvalue : float
        Two-tailed p-value
    """
    if mask is not None:
        y_true = y_true[mask]
        y_pred = y_pred[mask]
    
    if len(y_true) < 2:
        return np.nan, np.nan
    
    return pearsonr(y_true, y_pred)


def calculate_spearman(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    mask: Optional[np.ndarray] = None
) -> Tuple[float, float]:
    """
    Calculate Spearman correlation coefficient.
    
    Parameters
    ----------
    y_true : np.ndarray
        True values
    y_pred : np.ndarray
        Predicted values
    mask : np.ndarray, optional
        Boolean mask for values to include
    
    Returns
    -------
    correlation : float
        Spearman correlation coefficient
    pvalue : float
        Two-tailed p-value
    """
    if mask is not None:
        y_true = y_true[mask]
        y_pred = y_pred[mask]
    
    if len(y_true) < 2:
        return np.nan, np.nan
    
    return spearmanr(y_true, y_pred)


def calculate_r2(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    mask: Optional[np.ndarray] = None
) -> float:
    """
    Calculate R² (coefficient of determination).
    
    Parameters
    ----------
    y_true : np.ndarray
        True values
    y_pred : np.ndarray
        Predicted values
    mask : np.ndarray, optional
        Boolean mask for values to include
    
    Returns
    -------
    float
        R² value
    """
    if mask is not None:
        y_true = y_true[mask]
        y_pred = y_pred[mask]
    
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    
    if ss_tot == 0:
        return np.nan
    
    return 1 - (ss_res / ss_tot)


def cross_validate_imputation(
    imputer,
    coordinates: np.ndarray,
    values: np.ndarray,
    n_folds: int = 5,
    mask: Optional[np.ndarray] = None,
    random_state: int = 42
) -> Dict[str, float]:
    """
    Perform k-fold cross-validation for imputation.
    
    Parameters
    ----------
    imputer : object
        Imputer object with fit() and predict() methods
    coordinates : np.ndarray, shape (n_spots, 2)
        Spatial coordinates
    values : np.ndarray, shape (n_spots,)
        Observed values
    n_folds : int, default=5
        Number of cross-validation folds
    mask : np.ndarray, optional
        Boolean mask for observed values
    random_state : int, default=42
        Random seed
    
    Returns
    -------
    dict
        Dictionary with metrics: 'rmse', 'mae', 'pearson', 'r2'
    """
    if mask is None:
        mask = values > 0
    
    # Get observed indices
    observed_indices = np.where(mask)[0]
    
    if len(observed_indices) < n_folds:
        logger.warning(f"Not enough observed values ({len(observed_indices)}) for {n_folds} folds")
        n_folds = max(2, len(observed_indices) // 2)
    
    kfold = KFold(n_splits=n_folds, shuffle=True, random_state=random_state)
    
    rmse_scores = []
    mae_scores = []
    pearson_scores = []
    r2_scores = []
    
    for train_idx, test_idx in kfold.split(observed_indices):
        # Create train/test masks
        train_mask = np.zeros_like(mask)
        test_mask = np.zeros_like(mask)
        
        train_mask[observed_indices[train_idx]] = True
        test_mask[observed_indices[test_idx]] = True
        
        # Fit on training data
        try:
            imputer.fit(coordinates, values, train_mask)
            
            # Predict on test data
            predictions, _ = imputer.predict(coordinates[test_mask])
            true_values = values[test_mask]
            
            # Calculate metrics
            rmse_scores.append(calculate_rmse(true_values, predictions))
            mae_scores.append(calculate_mae(true_values, predictions))
            
            pearson_corr, _ = calculate_pearson(true_values, predictions)
            pearson_scores.append(pearson_corr)
            
            r2_scores.append(calculate_r2(true_values, predictions))
        
        except Exception as e:
            logger.warning(f"Fold failed: {e}")
            continue
    
    if not rmse_scores:
        logger.error("All folds failed")
        return {
            'rmse': np.nan,
            'mae': np.nan,
            'pearson': np.nan,
            'r2': np.nan
        }
    
    return {
        'rmse': np.mean(rmse_scores),
        'rmse_std': np.std(rmse_scores),
        'mae': np.mean(mae_scores),
        'mae_std': np.std(mae_scores),
        'pearson': np.mean(pearson_scores),
        'pearson_std': np.std(pearson_scores),
        'r2': np.mean(r2_scores),
        'r2_std': np.std(r2_scores)
    }


def evaluate_imputation_quality(
    true_values: np.ndarray,
    imputed_values: np.ndarray,
    uncertainty: Optional[np.ndarray] = None,
    mask: Optional[np.ndarray] = None
) -> Dict[str, float]:
    """
    Evaluate imputation quality with multiple metrics.
    
    Parameters
    ----------
    true_values : np.ndarray
        True values
    imputed_values : np.ndarray
        Imputed values
    uncertainty : np.ndarray, optional
        Uncertainty estimates
    mask : np.ndarray, optional
        Mask for imputed values (True = imputed, False = observed)
    
    Returns
    -------
    dict
        Dictionary with quality metrics
    """
    if mask is None:
        # Assume all non-zero true values are observed
        mask = true_values == 0
    
    # Calculate basic metrics
    rmse = calculate_rmse(true_values, imputed_values, mask)
    mae = calculate_mae(true_values, imputed_values, mask)
    pearson_corr, pearson_pval = calculate_pearson(true_values, imputed_values, mask)
    spearman_corr, spearman_pval = calculate_spearman(true_values, imputed_values, mask)
    r2 = calculate_r2(true_values, imputed_values, mask)
    
    metrics = {
        'rmse': rmse,
        'mae': mae,
        'pearson_r': pearson_corr,
        'pearson_pval': pearson_pval,
        'spearman_r': spearman_corr,
        'spearman_pval': spearman_pval,
        'r2': r2,
        'n_imputed': np.sum(mask)
    }
    
    # Add uncertainty-based metrics if available
    if uncertainty is not None:
        imputed_uncertainty = uncertainty[mask]
        metrics['mean_uncertainty'] = np.mean(imputed_uncertainty)
        metrics['median_uncertainty'] = np.median(imputed_uncertainty)
        metrics['max_uncertainty'] = np.max(imputed_uncertainty)
        
        # Calibration: correlation between uncertainty and error
        errors = np.abs(true_values[mask] - imputed_values[mask])
        if len(errors) > 1:
            calib_corr, calib_pval = pearsonr(imputed_uncertainty, errors)
            metrics['uncertainty_calibration'] = calib_corr
            metrics['calibration_pval'] = calib_pval
    
    return metrics


class QCReportGenerator:
    """
    Generate quality control reports for APA analysis.
    
    Parameters
    ----------
    dataset_name : str, optional
        Name of the dataset
    
    Examples
    --------
    >>> qc = QCReportGenerator(dataset_name="Mouse Brain")
    >>> qc.add_imputation_metrics(metrics)
    >>> qc.add_quantification_metrics(apa_indices)
    >>> report = qc.generate_report()
    """
    
    def __init__(self, dataset_name: Optional[str] = None):
        self.dataset_name = dataset_name or "Unknown"
        self.sections = {}
    
    def add_imputation_metrics(self, metrics: Dict[str, float]):
        """Add imputation quality metrics."""
        self.sections['imputation'] = metrics
    
    def add_quantification_metrics(self, indices: Dict[str, np.ndarray]):
        """Add APA quantification metrics."""
        metrics = {}
        for index_name, values in indices.items():
            metrics[f'{index_name}_mean'] = np.mean(values)
            metrics[f'{index_name}_std'] = np.std(values)
            metrics[f'{index_name}_min'] = np.min(values)
            metrics[f'{index_name}_max'] = np.max(values)
            metrics[f'{index_name}_median'] = np.median(values)
        
        self.sections['quantification'] = metrics
    
    def add_spatial_metrics(self, spatial_stats: Dict[str, float]):
        """Add spatial analysis metrics."""
        self.sections['spatial'] = spatial_stats
    
    def add_coverage_metrics(
        self,
        total_spots: int,
        observed_spots: int,
        imputed_spots: int
    ):
        """Add data coverage metrics."""
        self.sections['coverage'] = {
            'total_spots': total_spots,
            'observed_spots': observed_spots,
            'imputed_spots': imputed_spots,
            'observed_fraction': observed_spots / total_spots,
            'imputed_fraction': imputed_spots / total_spots
        }
    
    def generate_report(self, format: str = 'dict') -> Dict:
        """
        Generate QC report.
        
        Parameters
        ----------
        format : str, default='dict'
            Output format: 'dict', 'dataframe', or 'text'
        
        Returns
        -------
        dict or pd.DataFrame or str
            QC report in requested format
        """
        if format == 'dict':
            return {
                'dataset': self.dataset_name,
                **self.sections
            }
        
        elif format == 'dataframe':
            # Flatten nested dict
            flat_dict = {'dataset': self.dataset_name}
            for section, metrics in self.sections.items():
                for key, value in metrics.items():
                    flat_dict[f'{section}_{key}'] = value
            
            return pd.DataFrame([flat_dict])
        
        elif format == 'text':
            lines = [
                f"Quality Control Report",
                f"=" * 50,
                f"Dataset: {self.dataset_name}",
                ""
            ]
            
            for section, metrics in self.sections.items():
                lines.append(f"\n{section.upper()}")
                lines.append("-" * 50)
                for key, value in metrics.items():
                    if isinstance(value, float):
                        lines.append(f"  {key}: {value:.4f}")
                    else:
                        lines.append(f"  {key}: {value}")
            
            return "\n".join(lines)
        
        else:
            raise ValueError(f"Unknown format: {format}")
    
    def save_report(self, filename: str, format: str = 'text'):
        """
        Save report to file.
        
        Parameters
        ----------
        filename : str
            Output filename
        format : str, default='text'
            Output format
        """
        report = self.generate_report(format=format)
        
        if format == 'text':
            with open(filename, 'w') as f:
                f.write(report)
        elif format == 'dataframe':
            report.to_csv(filename, index=False)
        else:
            # Save dict as JSON
            import json
            with open(filename, 'w') as f:
                json.dump(report, f, indent=2)
        
        logger.info(f"Saved QC report to {filename}")


def compare_imputation_methods(
    methods: Dict[str, Callable],
    coordinates: np.ndarray,
    values: np.ndarray,
    mask: Optional[np.ndarray] = None,
    n_folds: int = 5
) -> pd.DataFrame:
    """
    Compare multiple imputation methods using cross-validation.
    
    Parameters
    ----------
    methods : dict
        Dictionary mapping method names to imputer objects
    coordinates : np.ndarray
        Spatial coordinates
    values : np.ndarray
        Observed values
    mask : np.ndarray, optional
        Boolean mask for observed values
    n_folds : int, default=5
        Number of CV folds
    
    Returns
    -------
    pd.DataFrame
        Comparison results with metrics for each method
    """
    results = []
    
    for method_name, imputer in methods.items():
        logger.info(f"Evaluating {method_name}...")
        
        try:
            metrics = cross_validate_imputation(
                imputer, coordinates, values, n_folds, mask
            )
            metrics['method'] = method_name
            results.append(metrics)
        
        except Exception as e:
            logger.error(f"Failed to evaluate {method_name}: {e}")
            continue
    
    if not results:
        logger.error("All methods failed")
        return pd.DataFrame()
    
    df = pd.DataFrame(results)
    df = df.set_index('method')
    
    return df
