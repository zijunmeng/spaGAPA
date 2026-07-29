#!/usr/bin/env bash
# =============================================================================
# Download 2 smallest psoriasis-lesional Visium samples from GSE206391
# (human skin) for APA analysis:
#   GSM6252925 (Patient 4 lesional psoriasis, rep 1) ~0.80 GB  - SRR19737306-9
#   GSM6252926 (Patient 4 lesional psoriasis, rep 2) ~0.58 GB  - SRR19737310-3
#
# ENA SRA read layout (verified by partial download):
#   _1.fastq.gz = 28 bp  -> Visium R1 (16 bp barcode + 12 bp UMI)
#   _2.fastq.gz = 120 bp -> Visium R2 (cDNA)
# Standard Visium layout: NO SRA-split trick needed.
#
# Uses aria2c (proxy disabled per task) with the ENA HTTPS/S3 endpoint.
# Downloads tissue image (.jpg.gz) from GEO supplementary FTP.
# Hostname-aware (CLAUDE.md). Do NOT git commit.
# =============================================================================
set -uo pipefail

# ---- hostname-aware environment (CLAUDE.md) ----
HOSTNAME=$(hostname -s)
case "$HOSTNAME" in
    S91)
        export TMPDIR="/s3/mengzijun/tmp"
        ARIA2C="/usr/bin/aria2c"
        ;;
    *)
        echo "ERROR: this script is tuned for S91 (host=$HOSTNAME)" >&2
        exit 1
        ;;
esac
mkdir -p "$TMPDIR"

ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA"
DATA_ROOT="${ROOT}/data/raw/gse206391"
FASTQ_DIR="${DATA_ROOT}/fastq_ena"
SUPP_DIR="${DATA_ROOT}/geo_supplementary"
LOG_DIR="${DATA_ROOT}/logs"
DATE_TAG="$(date +%Y%m%d)"
mkdir -p "$FASTQ_DIR" "$SUPP_DIR" "$LOG_DIR"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "${LOG_DIR}/download_${DATE_TAG}.log"; }

# ---- chosen samples ----
# GSM -> list of SRR runs
declare -A GSM_RUNS
GSM_RUNS[GSM6252925]="SRR19737306 SRR19737307 SRR19737308 SRR19737309"
GSM_RUNS[GSM6252926]="SRR19737310 SRR19737311 SRR19737312 SRR19737313"

declare -A GSM_IMAGE
GSM_IMAGE[GSM6252925]="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6252nnn/GSM6252925/suppl/GSM6252925_4-V19S23-006-V3.jpg.gz"
GSM_IMAGE[GSM6252926]="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6252nnn/GSM6252926/suppl/GSM6252926_4-V19S23-006-V4.jpg.gz"

# ---- 1. Build ENA FASTQ URL list (HTTPS/S3-backed endpoint, proxy disabled) ----
build_fastq_urls() {
    local url_file="${LOG_DIR}/ena_fastq_urls.txt"
    : > "$url_file"
    for gsm in "${!GSM_RUNS[@]}"; do
        for srr in ${GSM_RUNS[$gsm]}; do
            # ENA FASTQ path algorithm by accession length:
            #   6 chars  : vol1/fastq/ERR123/ERR123456
            #   7 chars  : vol1/fastq/ERR123/004/ERR1234567      (00 + last 1 digit)
            #   8 chars  : vol1/fastq/ERR123/0044/ERR12345678    (00 + last 2 digits)
            #   9 chars  : vol1/fastq/ERR123/00445/ERR123456789  (00 + last 3 digits)
            #   10 chars : vol1/fastq/ERR123/004-5/ERR1234567890 (0 + last 3 / last 2)
            #   11 chars : vol1/fastq/ERR123/004-56/ERR12345678901 (0 + last 4 / last 3)
            # All chosen runs here are 11 digits (SRR197373xx), middle dir = "00" + last digit.
            local prefix="${srr:0:6}"
            local len=${#srr}
            local base
            # ENA subdir rules (verified empirically):
            #   6:  ERR123/ERR123456
            #   7:  ERR123/00d       (d = last 1 digit)
            #   8:  ERR123/00dd      (dd = last 2 digits)
            #   9:  ERR123/00ddd     (ddd = last 3 digits)
            #   10: ERR123/0ddd/dd   (ddd = chars 6-8, dd = last 2)  [two subdirs]
            #   11: ERR123/0ddd/ddd  -- but empirically 11-digit SRR197373xx use
            #       a SINGLE subdir = "0" + last 2 chars, e.g. SRR19737306 -> SRR197/006,
            #       SRR19737313 -> SRR197/013. Override the generic rule accordingly.
            case "$len" in
                6)  base="https://ftp.sra.ebi.ac.uk/vol1/fastq/${prefix}/${srr}/${srr}" ;;
                7)  base="https://ftp.sra.ebi.ac.uk/vol1/fastq/${prefix}/00${srr:6:1}/${srr}/${srr}" ;;
                8)  base="https://ftp.sra.ebi.ac.uk/vol1/fastq/${prefix}/00${srr:6:2}/${srr}/${srr}" ;;
                9)  base="https://ftp.sra.ebi.ac.uk/vol1/fastq/${prefix}/00${srr:6:3}/${srr}/${srr}" ;;
                10) base="https://ftp.sra.ebi.ac.uk/vol1/fastq/${prefix}/0${srr:5:3}/${srr:8:2}/${srr}/${srr}" ;;
                11) base="https://ftp.sra.ebi.ac.uk/vol1/fastq/${prefix}/0${srr:9:2}/${srr}/${srr}" ;;
                *)  echo "ERROR: unexpected SRR length $len for $srr" >&2; exit 3 ;;
            esac
            echo "${base}_1.fastq.gz" >> "$url_file"
            echo "${base}_2.fastq.gz" >> "$url_file"
        done
    done
    sort -u "$url_file" -o "$url_file"
    echo "$url_file"
}

# ---- 2. Download FASTQ via aria2c (proxy disabled) ----
download_fastq() {
    local url_file
    url_file=$(build_fastq_urls)
    local n
    n=$(wc -l < "$url_file")
    log "FASTQ URLs ($n files): $url_file"
    cat "$url_file" | sed 's/^/  /' | head -20 | tee -a "${LOG_DIR}/download_${DATE_TAG}.log"

    # aria2c with proxy explicitly disabled via env (no http_proxy/https_proxy)
    log "Starting aria2c FASTQ download (proxy disabled)"
    env -u http_proxy -u https_proxy -u HTTP_PROXY -u HTTPS_PROXY -u all_proxy -u ALL_PROXY \
        "$ARIA2C" \
            --input-file="$url_file" \
            --dir="$FASTQ_DIR" \
            --check-certificate=false \
            --max-concurrent-downloads=8 \
            --max-connection-per-server=2 \
            --split=2 \
            --min-split-size=1M \
            --continue=true \
            --file-allocation=none \
            --auto-file-renaming=false \
            --allow-overwrite=false \
            --max-tries=5 \
            --retry-wait=15 \
            --connect-timeout=60 \
            --timeout=300 \
            --console-log-level=warn \
            --download-result=full \
            --summary-interval=30 \
            --log="${LOG_DIR}/fastq_aria2c_${DATE_TAG}.log" \
        2>&1 | tee -a "${LOG_DIR}/download_${DATE_TAG}.log"
    local rc=${PIPESTATUS[0]}
    log "aria2c FASTQ rc=${rc}"
    return $rc
}

# ---- 3. Download GEO tissue images ----
download_images() {
    local url_file="${LOG_DIR}/geo_image_urls.txt"
    : > "$url_file"
    for gsm in "${!GSM_IMAGE[@]}"; do
        echo "${GSM_IMAGE[$gsm]}" >> "$url_file"
        echo "  dir=${SUPP_DIR}/${gsm}" >> "$url_file"
        echo "  out=$(basename "${GSM_IMAGE[$gsm]}")" >> "$url_file"
    done
    log "Downloading tissue images via aria2c"
    # aria2c input with per-URL dir/out directives
    env -u http_proxy -u https_proxy -u HTTP_PROXY -u HTTPS_PROXY \
        "$ARIA2C" \
            --input-file="$url_file" \
            --dir="$SUPP_DIR" \
            --check-certificate=false \
            --max-concurrent-downloads=2 \
            --max-connection-per-server=1 \
            --split=1 \
            --continue=true \
            --file-allocation=none \
            --auto-file-renaming=false \
            --allow-overwrite=false \
            --max-tries=5 \
            --retry-wait=15 \
            --console-log-level=warn \
            --download-result=full \
            --log="${LOG_DIR}/images_aria2c_${DATE_TAG}.log" \
        2>&1 | tee -a "${LOG_DIR}/download_${DATE_TAG}.log"
    log "Image download done"
}

# ---- 4. Verify downloads (count + sizes) ----
verify() {
    log "=== FASTQ files downloaded ==="
    ls -lh "$FASTQ_DIR"/*.fastq.gz 2>/dev/null | tee -a "${LOG_DIR}/download_${DATE_TAG}.log"
    log "Expected: 16 FASTQ (2 GSM x 4 SRR x 2 reads)"
    log "=== Image files ==="
    ls -lh "$SUPP_DIR"/*.jpg.gz 2>/dev/null | tee -a "${LOG_DIR}/download_${DATE_TAG}.log"
}

MODE="${1:-all}"
case "$MODE" in
    fastq)   download_fastq ;;
    images)  download_images ;;
    verify)  verify ;;
    all)     download_fastq; download_images; verify ;;
    urls)    build_fastq_urls; cat "$(build_fastq_urls)" ;;
    *) echo "Usage: $0 {all|fastq|images|verify|urls}" >&2; exit 2 ;;
esac
