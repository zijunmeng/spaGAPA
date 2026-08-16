#!/usr/bin/env bash
# =============================================================================
# Multi-caller robustness validation: full Sierra chain for one dataset.
# Usage: bash run_dataset.sh <dataset-key>
#   keys: gse220442_gsm6801751 | gse169749_gsm5213483 | gse183456_gsm6047774
#
# Steps (identical protocol to the GSE183456 run, report.md):
#   0. decompress reference GTF + build spot whitelist (from baseline coords)
#   1. pysam splice-junction extraction (extract_junctions.py)
#   2. Sierra FindPeaks + CountPeaks UMI-deduplicated (run_sierra.R)
#   3. convert to spaGAPA format (convert_to_spagapa.py)
#   4. split-conformal calibration on the Sierra input AND on the scAPAtrap
#      baseline input (scripts/calibrate_uncertainty.py, same protocol)
#   5. caller comparison (compare_callers.py) + GP-vs-mean (gp_vs_mean.py)
#
# Hostname-aware env (CLAUDE.md S90/S91/S97/S98). Idempotent per step.
# =============================================================================
set -euo pipefail

DATASET="${1:?usage: run_dataset.sh <dataset-key>}"
HOSTNAME=$(hostname -s)
case "$HOSTNAME" in
    S90) export R_LIBS="/s1/SHARE/01_software/R_442_SeuratV5/library"; RSCRIPT="/s1/mengzijun/anaconda3/envs/r442/bin/Rscript"; PY="/s1/mengzijun/anaconda3/envs/spagapa/bin/python"; export TMPDIR="/s2/mengzijun/tmp" ;;
    S91) export R_LIBS="/s1/SHARE/01_software/R_442_SeuratV5/library"; RSCRIPT="$HOME/anaconda3/envs/r442/bin/Rscript"; PY="$HOME/anaconda3/envs/spagapa/bin/python"; export TMPDIR="/s3/mengzijun/tmp" ;;
    S97) export R_LIBS="/s1/SHARE/01_software/R_442_SeuratV5/library"; RSCRIPT="/s1/mengzijun/anaconda3/envs/r442/bin/Rscript"; PY="/s1/mengzijun/anaconda3/envs/spagapa/bin/python"; export TMPDIR="/s972/mengzijun/tmp" ;;
    S98) export R_LIBS="/s1/SHARE/01_software/R_442_SeuratV5/library"; RSCRIPT="/s1/mengzijun/anaconda3/envs/r442/bin/Rscript"; PY="/s1/mengzijun/anaconda3/envs/spagapa/bin/python"; export TMPDIR="/s982/mengzijun/tmp" ;;
    *)  echo "ERROR: unknown host $HOSTNAME" >&2; exit 1 ;;
esac
mkdir -p "$TMPDIR"
export OPENBLAS_NUM_THREADS=8

ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
SCRIPTS="$ROOT/scripts/multicaller_validation"
BAM=$("$PY" -c "import sys; sys.path.insert(0, '$SCRIPTS'); from datasets import get; print(get('$DATASET')['bam'])")
GTF_GZ=$("$PY" -c "import sys; sys.path.insert(0, '$SCRIPTS'); from datasets import get; print(get('$DATASET')['gtf_gz'])")
BASELINE=$("$PY" -c "import sys; sys.path.insert(0, '$SCRIPTS'); from datasets import get; print(get('$DATASET')['baseline'])")
OUT=$("$PY" -c "import sys; sys.path.insert(0, '$SCRIPTS'); from datasets import get; print(get('$DATASET')['outdir'])")
mkdir -p "$OUT"
cd "$ROOT"
log () { echo "[$(date '+%F %T')] $*" | tee -a "$OUT/run_dataset.log"; }
log "dataset=$DATASET bam=$BAM gtf=$GTF_GZ baseline=$BASELINE out=$OUT"

# ---- 0. GTF (decompress if needed; hardlink an existing identical decompress) ----
if [ ! -s "$OUT/genes.gtf" ]; then
    if [ -s "$ROOT/pipeline_output/multicaller_validation/sierra/genes.gtf" ] && \
       [ "$GTF_GZ" = "/s1/SHARE/00_ref_genecode/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz" ]; then
        log "hardlinking existing GRCh38-2024-A decompressed GTF"
        ln "$ROOT/pipeline_output/multicaller_validation/sierra/genes.gtf" "$OUT/genes.gtf"
    else
        log "decompressing $GTF_GZ"
        gunzip -c "$GTF_GZ" > "$OUT/genes.gtf"
    fi
fi

# ---- 0b. whitelist = baseline spot barcodes ----
if [ ! -s "$OUT/whitelist_barcodes.tsv" ]; then
    "$PY" - "$BASELINE/coordinates.csv" "$OUT/whitelist_barcodes.tsv" <<'EOF'
import sys, pandas as pd
c = pd.read_csv(sys.argv[1], index_col=0)
c.index.to_series().to_csv(sys.argv[2], index=False, header=False)
print(f"[wl] {len(c)} barcodes -> {sys.argv[2]}")
EOF
fi
log "whitelist: $(wc -l < "$OUT/whitelist_barcodes.tsv") barcodes"

# ---- 1. junctions ----
if [ ! -s "$OUT/junctions.bed" ]; then
    log "extracting splice junctions"
    "$PY" "$SCRIPTS/extract_junctions.py" --dataset "$DATASET" 2>&1 | tee "$OUT/extract_junctions.log"
else
    log "junctions.bed exists ($(wc -l < "$OUT/junctions.bed") entries), skipping"
fi

# ---- 2. Sierra FindPeaks + CountPeaks ----
if [ ! -s "$OUT/sierra_counts/matrix.mtx.gz" ]; then
    log "running Sierra"
    "$RSCRIPT" "$SCRIPTS/run_sierra.R" "$OUT" \
        "$BAM" \
        "$OUT/genes.gtf" "$OUT/junctions.bed" "$OUT/whitelist_barcodes.tsv" 16 \
        2>&1 | tee "$OUT/run_sierra.log"
else
    log "sierra_counts exists, skipping"
fi

# ---- 3. convert ----
if [ ! -s "$OUT/apa_matrix.csv" ]; then
    log "converting to spaGAPA format"
    "$PY" "$SCRIPTS/convert_to_spagapa.py" --dataset "$DATASET" 2>&1 | tee "$OUT/convert.log"
else
    log "apa_matrix.csv exists, skipping"
fi

# ---- 4. conformal: Sierra input + scAPAtrap baseline control ----
if [ ! -s "$OUT/conformal/uncertainty_calibration.json" ]; then
    log "conformal calibration (Sierra input)"
    "$PY" "$ROOT/scripts/calibrate_uncertainty.py" \
        --apa-matrix "$OUT/apa_matrix.csv" \
        --coordinates "$BASELINE/coordinates.csv" \
        --output "$OUT/conformal" > "$OUT/conformal.log" 2>&1
    tail -20 "$OUT/conformal.log"
else
    log "conformal exists, skipping"
fi
if [ ! -s "$OUT/scapatrap_conformal/uncertainty_calibration.json" ]; then
    log "conformal calibration (scAPAtrap baseline input, same protocol)"
    "$PY" "$ROOT/scripts/calibrate_uncertainty.py" \
        --apa-matrix "$BASELINE/apa_matrix.csv" \
        --coordinates "$BASELINE/coordinates.csv" \
        --output "$OUT/scapatrap_conformal" > "$OUT/scapatrap_conformal.log" 2>&1
    tail -8 "$OUT/scapatrap_conformal.log"
else
    log "scapatrap conformal exists, skipping"
fi

# ---- 5. metrics ----
log "compare_callers (metrics a+b)"
"$PY" "$SCRIPTS/compare_callers.py" --dataset "$DATASET" 2>&1 | tee "$OUT/compare_callers.log"
log "gp_vs_mean (metric d)"
"$PY" "$SCRIPTS/gp_vs_mean.py" --dataset "$DATASET" 2>&1 | tee "$OUT/gp_vs_mean.log"

log "DATASET $DATASET COMPLETE"
