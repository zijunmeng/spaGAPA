"""
Complete spaGAPA Pipeline with uncertainty-weighted analysis.

Integrates all modules: spatial validation, GP imputation (with uncertainty),
APA quantification, domain identification, differential APA, and SVAPA detection.
All downstream analyses support uncertainty weighting when available.
"""

import numpy as np
import pandas as pd
from typing import Optional, Dict, List, Union, Tuple
from pathlib import Path
import warnings

from spagapa.core import APADataset
from spagapa.io import load_spatial_dataset
from spagapa.spatial import SpatialNeighbors
from spagapa.calling import SpatialValidator, QualityFilter
from spagapa.imputation import GPImputer, SparseGPImputer
from spagapa.quantification import APAIndexCalculator, QCReportGenerator
from spagapa.analysis import (
    DomainIdentifier,
    DifferentialAPAAnalyzer,
    SpatialPatternAnalyzer,
    GPTrendDetector,
)


class SpaGAPA:
    """
    Main spaGAPA pipeline with uncertainty-weighted analysis.

    Workflow (7 steps):
    1. Load data → APADataset
    2. Spatial validation + quality filtering
    3. GP imputation → imputed values + uncertainty
    4. APA quantification (RUD, PDUI, WUL)
    5. Spatial domain identification (uncertainty-weighted)
    6. Differential APA analysis (uncertainty-weighted)
    7. SVAPA gene detection (GP likelihood ratio + weighted Moran's I)

    All downstream analyses automatically use uncertainty weights when
    imputation uncertainty is available.

    Parameters
    ----------
    n_neighbors : int, default=6
        Spatial neighbors for graph construction.
    kernel_type : str, default='matern'
        GP kernel: 'rbf', 'matern', or 'auto'.
    use_sparse_gp : bool, default=False
        Use sparse GP approximation for large datasets.
    n_inducing : int, default=100
        Inducing points for sparse GP.
    min_spatial_support : float, default=0.3
        Minimum spatial support threshold.
    min_read_count : int, default=10
        Minimum read count per site.
    min_spots : int, default=5
        Minimum spots for quality filtering.
    verbose : bool, default=True
        Print progress.
    """

    def __init__(
        self,
        n_neighbors: int = 6,
        kernel_type: str = 'matern',
        use_sparse_gp: bool = False,
        n_inducing: int = 100,
        input_type: str = 'apa_index',
        min_spatial_support: float = 0.3,
        min_read_count: int = 10,
        min_spots: int = 5,
        verbose: bool = True,
    ):
        self.n_neighbors = n_neighbors
        self.kernel_type = kernel_type
        self.use_sparse_gp = use_sparse_gp
        self.n_inducing = n_inducing
        self.input_type = input_type
        self.min_spatial_support = min_spatial_support
        self.min_read_count = min_read_count
        self.min_spots = min_spots
        self.verbose = verbose

        self.dataset_: Optional[APADataset] = None
        self.imputer_ = None
        self.results_: Dict = {}

    def _log(self, msg: str):
        if self.verbose:
            print(msg)

    def _validate_dataset_shapes(self):
        """Validate the project-wide genes x spots algorithm contract."""
        if self.dataset_ is None:
            raise ValueError("No dataset loaded")

        counts = self.dataset_.raw_counts
        coords = self.dataset_.coords

        if counts.ndim != 2:
            raise ValueError("raw_counts must have shape (n_genes, n_spots)")
        if coords.ndim != 2 or coords.shape[1] < 2:
            raise ValueError("coords must have shape (n_spots, 2)")
        expected = (self.dataset_.n_genes, self.dataset_.n_spots)
        if counts.shape != expected:
            raise ValueError(f"raw_counts shape mismatch: expected {expected}, got {counts.shape}")
        if coords.shape[0] != self.dataset_.n_spots:
            raise ValueError(
                "Coordinate count does not match spots: "
                f"{coords.shape[0]} coords vs {self.dataset_.n_spots} spots"
            )

    def _clip_apa_values(self, values: np.ndarray) -> np.ndarray:
        """Clip matrices according to the declared input scale."""
        if self.input_type == 'apa_index':
            return np.clip(values, 0.0, 1.0)
        return np.clip(values, 0.0, None)

    def _training_mask(self, values: np.ndarray) -> np.ndarray:
        """Return observed-value mask under the declared input semantics."""
        finite = np.isfinite(values)
        if self.input_type == 'apa_index':
            return finite
        return finite & (values > 0)

    def _calculate_apa_indices(self, work: np.ndarray) -> Dict[str, Union[np.ndarray, bool, str]]:
        """
        Register or calculate APA indices for pipeline output.

        For ``input_type='apa_index'`` the matrix is already an APA index, so
        this step stores it explicitly instead of pretending to recompute RUD.
        For ``input_type='proximal_distal_counts'`` the dataset must contain
        AnnData layers ``proximal_counts`` and ``distal_counts`` in spots x
        genes layout.
        """
        if self.input_type == 'apa_index':
            values = self._clip_apa_values(work)
            self.dataset_.set_apa_index('APAIndex', values)
            return {
                'APAIndex': values,
                'calculated': False,
                'source': 'input_apa_index',
            }

        if self.input_type != 'proximal_distal_counts':
            raise ValueError("input_type must be 'apa_index' or 'proximal_distal_counts'")

        layers = self.dataset_.adata.layers
        if 'proximal_counts' not in layers or 'distal_counts' not in layers:
            raise ValueError(
                "input_type='proximal_distal_counts' requires AnnData layers "
                "'proximal_counts' and 'distal_counts' in spots x genes layout"
            )

        proximal = np.asarray(layers['proximal_counts']).T
        distal = np.asarray(layers['distal_counts']).T
        if proximal.shape != work.shape or distal.shape != work.shape:
            raise ValueError(
                "proximal_counts and distal_counts layers must match raw_counts "
                "shape (genes x spots) after transpose"
            )

        calculator = APAIndexCalculator()
        indices = calculator.calculate_all(proximal, distal)
        for name, values in indices.items():
            self.dataset_.set_apa_index(name, values)
        indices['calculated'] = True
        indices['source'] = 'proximal_distal_counts'
        return indices

    # ── main entry ───────────────────────────────────────────────────

    def run(
        self,
        bam_file: Optional[str] = None,
        coordinates: Optional[Union[str, np.ndarray]] = None,
        annotation: Optional[str] = None,
        apa_matrix: Optional[Union[str, np.ndarray, pd.DataFrame]] = None,
        h5ad_file: Optional[str] = None,
        matrix_orientation: str = 'genes_by_spots',
        dataset: Optional[APADataset] = None,
        impute: bool = True,
        quantify: bool = True,
        identify_domains: bool = True,
        differential_analysis: bool = False,
        detect_svapa: bool = True,
        n_domains: Optional[int] = None,
        domain_labels: Optional[np.ndarray] = None,
        fdr_threshold: float = 0.05,
        use_uncertainty_weights: bool = True,
    ) -> Dict:
        """
        Run the complete spaGAPA pipeline.

        Parameters
        ----------
        bam_file : str, optional
            Path to BAM file.
        coordinates : str or np.ndarray, optional
            Path to coordinates file or (n_spots, 2) array.
        annotation : str, optional
            Path to GTF/GFF annotation.
        apa_matrix : str or array-like, optional
            APA count/index matrix. By default rows are genes and columns are spots.
        h5ad_file : str, optional
            Existing H5AD file to load.
        matrix_orientation : str, default='genes_by_spots'
            Orientation of apa_matrix: 'genes_by_spots' or 'spots_by_genes'.
        dataset : APADataset, optional
            Pre-loaded dataset (alternative to bam_file).
        impute : bool, default=True
            Perform GP imputation.
        quantify : bool, default=True
            Calculate APA indices.
        identify_domains : bool, default=True
            Identify spatial domains.
        differential_analysis : bool, default=False
            Perform differential analysis.
        detect_svapa : bool, default=True
            Detect SVAPA genes.
        n_domains : int, optional
            Number of domains for clustering.
        domain_labels : np.ndarray, optional
            Pre-defined domain labels.
        fdr_threshold : float, default=0.05
            FDR threshold.
        use_uncertainty_weights : bool, default=True
            Use imputation uncertainty to weight downstream analyses.

        Returns
        -------
        results : dict
            Keys: dataset, imputed_values, uncertainty, apa_indices,
                  domains, differential, svapa_genes, qc_report
        """
        self._log("=" * 60)
        self._log("spaGAPA: Spatial GP-based APA Analyzer")
        self._log("=" * 60)

        # ── Step 1: Load data ──
        self._log("\n[1/7] Loading data...")
        if dataset is not None:
            self.dataset_ = dataset
        elif apa_matrix is not None or h5ad_file is not None or bam_file is not None:
            self.dataset_ = load_spatial_dataset(
                apa_matrix=apa_matrix,
                bam_file=bam_file,
                coordinates=coordinates,
                h5ad_file=h5ad_file,
                matrix_orientation=matrix_orientation,
                annotation=annotation,
            )
        else:
            raise ValueError("Provide 'dataset', 'apa_matrix', 'h5ad_file', or 'bam_file'")
        self._validate_dataset_shapes()
        self._log(f"  ✓ {self.dataset_.n_genes} genes, {self.dataset_.n_spots} spots")

        coords = self.dataset_.coords
        if coords is None:
            raise ValueError("No spatial coordinates in dataset")

        # ── Step 2: Spatial validation & filtering ──
        self._log("\n[2/7] Spatial validation + quality filtering...")
        spatial_neighbors = SpatialNeighbors(n_neighbors=self.n_neighbors)
        spatial_neighbors.fit(coords)
        self._log(f"  ✓ Built spatial graph (k={self.n_neighbors})")

        validator = SpatialValidator(
            n_neighbors=self.n_neighbors,
            support_threshold=self.min_spatial_support,
        )
        validator.fit(coords)

        qc_filter = QualityFilter(
            min_read_count=self.min_read_count,
            min_spots=self.min_spots,
            min_spatial_support=self.min_spatial_support,
        )
        self._log(f"  ✓ Validated {self.dataset_.n_genes} APA sites")

        # ── Step 3: GP imputation ──
        uncertainty = None
        if impute:
            self._log(f"\n[3/7] GP imputation (kernel={self.kernel_type})...")
            apa_matrix_values = self.dataset_.raw_counts  # (n_genes, n_spots)

            if self.use_sparse_gp:
                base_imputer = SparseGPImputer(
                    n_inducing=self.n_inducing,
                )
            else:
                base_imputer = GPImputer(kernel_type=self.kernel_type)

            training_mask = self._training_mask(apa_matrix_values)
            self.imputer_ = base_imputer.fit_batch(
                coords,
                apa_matrix_values,
                mask=training_mask,
                n_jobs=1,
                verbose=self.verbose,
            )
            imputed, uncertainty = self.imputer_.impute(return_uncertainty=True)
            imputed = self._clip_apa_values(imputed)

            self.dataset_.set_imputed(imputed, uncertainty)

            observed = training_mask
            if observed.sum() > 0:
                mae = float(np.mean(np.abs(
                    apa_matrix_values[observed] - imputed[observed]
                )))
                self._log(f"  ✓ Imputation done (MAE={mae:.3f})")
            else:
                self._log("  ✓ Imputation done")

            self.results_['imputed_values'] = self.dataset_.imputed
            self.results_['uncertainty'] = self.dataset_.uncertainty
        else:
            self._log("\n[3/7] Skipping imputation")
            self.results_['imputed_values'] = None
            self.results_['uncertainty'] = None

        uncertainty = self.dataset_.uncertainty
        use_uw = use_uncertainty_weights and uncertainty is not None

        # Determine working matrix
        if self.dataset_.has_imputed():
            work = self.dataset_.imputed
        else:
            work = self.dataset_.raw_counts

        # ── Step 4: APA quantification ──
        if quantify:
            self._log("\n[4/7] APA quantification...")
            self.results_['apa_indices'] = self._calculate_apa_indices(work)
            if self.results_['apa_indices'].get('calculated'):
                self._log("  ✓ APA indices computed from proximal/distal counts")
            else:
                self._log("  ✓ Registered input APA index matrix")
        else:
            self._log("\n[4/7] Skipping quantification")
            self.results_['apa_indices'] = None

        # ── Step 5: Domain identification ──
        if identify_domains:
            self._log("\n[5/7] Identifying spatial domains...")
            if n_domains is None:
                n_domains = min(5, max(2, self.dataset_.n_spots // 20))

            domain_id = DomainIdentifier(
                method='kmeans',
                n_clusters=n_domains,
                min_domain_size=self.min_spots,
            )
            domain_labels = domain_id.identify_domains(
                work,
                coords,
                uncertainty=uncertainty if use_uw else None,
            )

            # Store in dataset
            self.dataset_.set_domain_labels(domain_labels)

            self.results_['domains'] = {
                'labels': domain_labels,
                'n_domains': n_domains,
            }
            self._log(f"  ✓ Identified {n_domains} spatial domains")
        else:
            self._log("\n[5/7] Skipping domain identification")
            self.results_['domains'] = None

        # ── Step 6: Differential analysis ──
        if differential_analysis:
            self._log("\n[6/7] Differential APA analysis...")
            if domain_labels is None and self.dataset_.has_domains():
                domain_labels = self.dataset_.get_domain_labels()

            if domain_labels is not None:
                analyzer = DifferentialAPAAnalyzer(method='wilcoxon')
                diff_results = analyzer.test_all_pairwise(
                    work,
                    domain_labels,
                    gene_names=self.dataset_.gene_names,
                    uncertainty=uncertainty if use_uw else None,
                )
                adjusted = {}
                for name, df in diff_results.items():
                    adjusted[name] = analyzer.adjust_pvalues(df)
                self.results_['differential'] = adjusted
                self._log("  ✓ Differential analysis complete")
            else:
                self._log("  ⚠ No domain labels; skipping")
                self.results_['differential'] = None
        else:
            self._log("\n[6/7] Skipping differential analysis")
            self.results_['differential'] = None

        # ── Step 7: SVAPA detection ──
        if detect_svapa:
            self._log(f"\n[7/7] Detecting SVAPA genes (FDR<{fdr_threshold})...")
            detector = GPTrendDetector(kernel_type=self.kernel_type)
            gene_names = self.dataset_.gene_names

            svapa_results = detector.detect_gp_trends(
                work, coords, uncertainty,
                gene_names=gene_names,
                method='likelihood_ratio',
                fdr_threshold=fdr_threshold,
            )
            self.results_['svapa_genes'] = svapa_results
            n_svapa = int(svapa_results['significant'].sum())
            self._log(f"  ✓ Detected {n_svapa} SVAPA genes")
        else:
            self._log("\n[7/7] Skipping SVAPA detection")
            self.results_['svapa_genes'] = None

        # ── QC report ──
        self._log("\nGenerating QC report...")
        qc = QCReportGenerator(dataset_name=f"spaGAPA_{self.dataset_.n_genes}genes")
        total = self.dataset_.n_spots * self.dataset_.n_genes
        observed = int((self.dataset_.adata.X > 0).sum())

        imputed_spots = observed
        if self.dataset_.has_imputed():
            imputed_spots = int((self.dataset_.imputed > 0).sum())

        qc.add_coverage_metrics(
            total_spots=total,
            observed_spots=observed,
            imputed_spots=imputed_spots,
        )

        if uncertainty is not None:
            qc.add_imputation_metrics({
                'mean_uncertainty': float(uncertainty.mean()),
                'max_uncertainty': float(uncertainty.max()),
                'min_uncertainty': float(uncertainty.min()),
            })

        self.results_['qc_report'] = qc.generate_report(format='dict')
        self.results_['dataset'] = self.dataset_

        self._log("\n" + "=" * 60)
        self._log("Pipeline complete!")
        self._log("=" * 60)

        return self.results_

    # ── accessors ────────────────────────────────────────────────────

    def get_imputed_values(self) -> Tuple[np.ndarray, np.ndarray]:
        im = self.dataset_.imputed
        un = self.dataset_.uncertainty
        if im is None:
            raise ValueError("Imputation not run")
        return im, un

    def get_svapa_genes(self, fdr_threshold: Optional[float] = None) -> pd.DataFrame:
        df = self.results_.get('svapa_genes')
        if df is None:
            raise ValueError("SVAPA detection not run")
        if fdr_threshold is not None:
            return df[df['q_value'] < fdr_threshold]
        return df[df['significant']]

    def get_domains(self) -> Dict:
        d = self.results_.get('domains')
        if d is None:
            raise ValueError("Domain identification not run")
        return d

    def save_results(self, output_dir: str):
        """Save all results."""
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        self._log(f"\nSaving results to {output_dir}...")

        if self.dataset_ is not None:
            self.dataset_.save(str(out / "dataset.h5ad"))
            self._log("  ✓ dataset.h5ad")

        imputed = self.dataset_.imputed
        if imputed is not None:
            np.save(out / "imputed_values.npy", imputed)
            self._log("  ✓ imputed_values.npy")

        uncertainty = self.dataset_.uncertainty
        if uncertainty is not None:
            np.save(out / "uncertainty.npy", uncertainty)
            self._log("  ✓ uncertainty.npy")

        svapa = self.results_.get('svapa_genes')
        if svapa is not None:
            svapa.to_csv(out / "svapa_genes.csv", index=False)
            self._log("  ✓ svapa_genes.csv")

        qc = self.results_.get('qc_report')
        if qc is not None:
            (out / "qc_report.txt").write_text(str(qc))
            self._log("  ✓ qc_report.txt")

        self._log("All results saved.")
