# Multi-caller robustness validation (Sierra × spaGAPA)

Evidence that spaGAPA is caller-agnostic: a second PAS caller (Sierra 0.99.27)
was run on the **same** GSE183456 Space Ranger BAM, its output converted to
spaGAPA's input format (usage matrix + shared coordinates, no re-tuning), and
the full spaGAPA inference re-run. Companion to the scAPAtrap baseline.

## Headline result
Split-conformal coverage on the Sierra input is exactly nominal —
**0.801 / 0.900 / 0.949** at 80/90/95% (59,457 held-out points) — matching the
scAPAtrap-input control (0.801 / 0.900 / 0.951). The framework's statistical
guarantees do not depend on the PAS caller.

Other comparisons: gene-level distal-usage concordance across callers
(median Pearson r = 0.51 over 814 shared genes with ≥10 observed spots);
PAS overlap is caller-dependent (Jaccard 0.074 at ±50 bp; 87.5% of Sierra
peaks within 500 bp of a scAPAtrap site); the "per-gene mean is a strong RMSE
baseline" pattern holds on both inputs.

## Files (runnable; absolute paths)
- `run_sierra.R` — Sierra FindPeaks/CountPeaks on the 10x BAM
- `extract_junctions.py` — pysam-based regtools-style junction BED (Sierra input)
- `convert_to_spagapa.py` — Sierra output -> spaGAPA usage-matrix format
- `compare_callers.py` — PAS overlap + gene-level usage concordance (metrics a/b)
- `gp_vs_mean.py` — GP-vs-mean RMSE on both inputs (metric d)
- `build_summary.py` — assemble summary.json
- `report.md` / `summary.json` — full numbers + a Methods-ready English paragraph

## Outputs (gitignored, on disk)
`pipeline_output/multicaller_validation/` — sierra raw+converted outputs,
conformal runs for both callers, metric JSONs.
