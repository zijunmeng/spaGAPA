#!/usr/bin/env bash
# =============================================================================
# Acquire two NEW Stereo-seq datasets (2026-09-09 expansion):
#   GSE293464  human retinal organoids (RA± × 16/26wk, 4 samples × 8 runs)
#   GSE333693  rat thymus (1 DNBSEQ run)
#
# Both have barcodeToPos.h5 masks DEPOSITED IN GEO supplementary (no STOmics
# portal / retention-window dependency — unlike GSE269906 whose masks expired:
# github.com/STOmics/SAW/issues/268).
#
# Assets:
#   masks : GEO FTP supplementary (aria2c, same pattern as the proven
#           GSE263789 WT mask download in download_gse263789_wt_control.sh)
#   FASTQ : SRA ODP S3 (aria2c -c -x4 -s4, same as download_gse263303_sra.sh)
#
# GSM ↔ SRR ↔ chip mapping (GSE293464):
#   GSM8882884  D02266B4  RA-/BMS1uM  26wk   SRR32936139..146
#   GSM8882885  D02266C2  RA1/10uM   26wk   SRR32936131..138
#   GSM8882886  D02266D2  RA-/BMS1uM  16wk   SRR32936123..130
#   GSM8882887  D02266D4  RA1/10uM   16wk   SRR32936115..122
#   GSE333693: GSM9770943  Y01052GC  SRR38888148
#
# Hostname-aware (CLAUDE.md). aria2c with proxy DISABLED (proven method:
# ~60-70 MiB/s vs ~250 KiB/s through the system socks proxy).
# Idempotent: aria2c -c resumes partials; allow-overwrite=false skips complete.
# =============================================================================
set -uo pipefail

HOSTNAME=$(hostname -s)
case "$HOSTNAME" in
    S91) export TMPDIR="/s3/mengzijun/tmp" ;;
    *)   echo "ERROR: unexpected host $HOSTNAME (script targets S91)" >&2; exit 1 ;;
esac

# proxy DISABLED (system socks proxy throttles aria2c)
unset http_proxy https_proxy HTTP_PROXY HTTPS_PROXY all_proxy ALL_PROXY

DATA=/s1/SHARE/mengzijun/01_project/26_spaGAPA/data/raw
LOGDIR=/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/logs
DATE_TAG=20260909
mkdir -p "$DATA/gse293464/masks" "$DATA/gse293464/sra" \
         "$DATA/gse333693/masks" "$DATA/gse333693/sra" "$LOGDIR"

log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOGDIR/${DATE_TAG}_stereo_expansion_download.log"; }

# ---- Phase 1: masks (GEO FTP, ~21 GB total) ---------------------------------
MASK_URLS=(
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8882nnn/GSM8882884/suppl/GSM8882884_D02266B4.barcodeToPos.h5|$DATA/gse293464/masks|D02266B4.barcodeToPos.h5"
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8882nnn/GSM8882885/suppl/GSM8882885_D02266C2.barcodeToPos.h5|$DATA/gse293464/masks|D02266C2.barcodeToPos.h5"
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8882nnn/GSM8882886/suppl/GSM8882886_D02266D2.barcodeToPos.h5|$DATA/gse293464/masks|D02266D2.barcodeToPos.h5"
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8882nnn/GSM8882887/suppl/GSM8882887_D02266D4.barcodeToPos.h5|$DATA/gse293464/masks|D02266D4.barcodeToPos.h5"
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9770nnn/GSM9770943/suppl/GSM9770943_Y01052GC.barcodeToPos.h5|$DATA/gse333693/masks|Y01052GC.barcodeToPos.h5"
)

log "PHASE 1: masks (5 files, ~21 GB)"
for entry in "${MASK_URLS[@]}"; do
    IFS='|' read -r url dir out <<< "$entry"
    if [[ -s "$dir/$out" && ! -f "$dir/$out.aria2" ]]; then
        log "  skip (complete): $out ($(du -h "$dir/$out" | cut -f1))"
        continue
    fi
    log "  downloading $out ..."
    aria2c -c -x4 -s4 --check-certificate=false \
        --allow-overwrite=false --auto-file-renaming=false \
        --console-log-level=warn --summary-interval=60 \
        --log="$LOGDIR/${DATE_TAG}_aria2_mask_$(echo "$out" | tr '.' '_').log" \
        -d "$dir" -o "$out" "$url" 2>&1 | tail -2 | tee -a "$LOGDIR/${DATE_TAG}_stereo_expansion_download.log"
    log "  exit=$? size=$(du -h "$dir/$out" 2>/dev/null | cut -f1)"
done

# ---- Phase 2: SRA runs (ODP S3, ~775 GB total) ------------------------------
# Build input file with per-URL dir directives: SRA files land as bare SRR ids.
SRA_INPUT="$LOGDIR/${DATE_TAG}_stereo_sra_urls.txt"
: > "$SRA_INPUT"
for srr in $(seq 32936115 32936146); do
    echo "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR${srr}/SRR${srr}" >> "$SRA_INPUT"
    echo "  dir=$DATA/gse293464/sra" >> "$SRA_INPUT"
done
echo "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR38888148/SRR38888148" >> "$SRA_INPUT"
echo "  dir=$DATA/gse333693/sra" >> "$SRA_INPUT"

log "PHASE 2: SRA runs (32 + 1 files, ~775 GB; est 3-4 h at 60-70 MiB/s)"
aria2c -c -x4 -s4 -j3 --check-certificate=false \
    --allow-overwrite=false --auto-file-renaming=false \
    --max-tries=8 --retry-wait=20 \
    --console-log-level=warn --summary-interval=300 \
    --log="$LOGDIR/${DATE_TAG}_aria2_sra.log" \
    --input-file="$SRA_INPUT" >> "$LOGDIR/${DATE_TAG}_stereo_expansion_download.log" 2>&1
rc=$?
log "PHASE 2 finished rc=$rc"

# ---- Summary ----------------------------------------------------------------
log "SUMMARY:"
log "  gse293464 masks: $(ls "$DATA/gse293464/masks" 2>/dev/null | wc -l)/4, sra: $(ls "$DATA/gse293464/sra" 2>/dev/null | wc -l)/32"
log "  gse333693 masks: $(ls "$DATA/gse333693/masks" 2>/dev/null | wc -l)/1, sra: $(ls "$DATA/gse333693/sra" 2>/dev/null | wc -l)/1"
log "  disk usage: $(du -sh "$DATA/gse293464" "$DATA/gse333693" 2>/dev/null | tr '\n' ' ')"
log "DONE (rc=$rc)"
