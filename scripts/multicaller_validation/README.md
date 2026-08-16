# Multi-caller robustness validation (Sierra × spaGAPA), 3 datasets

Evidence that spaGAPA is caller-agnostic: a second PAS caller (Sierra 0.99.27)
was run on the **same** Space Ranger BAM of each of three Visium datasets —
GSE183456 (human kidney), GSE220442/GSM6801751 (human AD brain PFC),
GSE169749/GSM5213483 (mouse colon; 2 species × 3 tissues) — its output
converted to spaGAPA's input format (usage matrix + shared coordinates, no
re-tuning), and the full spaGAPA inference re-run. Companion to the scAPAtrap
baseline of each dataset.

## Headline result (see report.md / summary.json for full numbers)
Split-conformal coverage on the Sierra input stays at nominal (80/90/95%,
global + locally adaptive) on every dataset, matching each dataset's
scAPAtrap-input control — the framework's statistical guarantees do not depend
on the PAS caller, the species, or the tissue.

## Files (runnable; absolute paths in datasets.py)
- `datasets.py` — dataset registry (BAM / reference GTF / baseline / outdir)
- `run_dataset.sh` — one-command driver for a full per-dataset chain
  (GTF decompress → whitelist → junctions → Sierra → convert → conformal ×2 →
  metrics); hostname-aware per CLAUDE.md
- `extract_junctions.py` — pysam regtools-style junction BED (Sierra input)
- `run_sierra.R` — Sierra FindPeaks/CountPeaks on the 10x BAM
- `convert_to_spagapa.py` — Sierra output -> spaGAPA usage-matrix format
  (coord = strand-aware 3' end; coord_summit = gaussian summit)
- `compare_callers.py` — PAS overlap + gene-level usage concordance (metrics a/b)
- `gp_vs_mean.py` — GP-vs-mean RMSE on both caller inputs (metric d)
- `build_summary.py` — cross-dataset summary.json + report.md
- `report.md` / `summary.json` — full numbers + a Methods-ready English paragraph

## Outputs (gitignored, on disk)
`pipeline_output/multicaller_validation/` — per dataset (`sierra/` for
GSE183456, `gse220442_gsm6801751/`, `gse169749_gsm5213483/`): junctions,
Sierra peaks/counts, converted spaGAPA tables, conformal runs for both
callers, metric JSONs; plus the cross-dataset summary.json/report.md.
