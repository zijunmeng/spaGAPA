#!/usr/bin/env bash
# Download GSE263789 WT control (GSM8199181, SS_WT18_F5) — mask + 4 SRR FASTQ.
# Same study as the AD pilot (GSM8199179) → FASTQ + mask both confirmed available.
# Hostname-aware (CLAUDE.md), aria2c with proxy DISABLED (proven GSE263789/269906 method).
set -uo pipefail
HOSTNAME=$(hostname -s)
case "$HOSTNAME" in
    S90|S91|S97|S98) ;;
    *) echo "ERROR unknown host $HOSTNAME" >&2; exit 1 ;;
esac

OUT=/s1/SHARE/mengzijun/01_project/26_spaGAPA/data/raw/gse263789/wt_control
mkdir -p "$OUT"
cd "$OUT"

# proxy DISABLED (system socks proxy throttles aria2c to ~250 KiB/s; without it ~60-70 MiB/s)
unset http_proxy https_proxy HTTP_PROXY HTTPS_PROXY all_proxy ALL_PROXY

MASK_URL="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8199nnn/GSM8199181/suppl/GSM8199181_SS200000745BL_F5.barcodeToPos.h5"
SRRS=(SRR28637878 SRR28637882 SRR28637883 SRR28637884)

echo "=== [$(date '+%F %T')] GSE263789 WT control download ==="
echo "  mask: barcodeToPos.h5 (4.2GB, chip SS200000745BL)"
echo "  FASTQ: ${SRRS[*]} (~108GB total, 1.3B reads)"

# 1. mask
echo "--- [$(date '+%T')] mask ---"
aria2c -c -x4 -s4 -d . -o barcodeToPos.h5 "$MASK_URL" 2>&1 | tail -2
echo "mask exit=$? size=$(du -h barcodeToPos.h5 2>/dev/null|cut -f1)"

# 2. 4 SRRs (SRA S3, same as pilot SRR28637909)
for SRR in "${SRRS[@]}"; do
    echo "--- [$(date '+%T')] $SRR ---"
    aria2c -c -x4 -s4 -d . "https://sra-pub-run-odp.s3.amazonaws.com/sra/${SRR}/${SRR}" 2>&1 | tail -2
    echo "$SRR exit=$? size=$(du -h ${SRR} 2>/dev/null|cut -f1)"
done

echo "=== [$(date '+%F %T')] DOWNLOAD COMPLETE ==="
ls -lh "$OUT"
echo "total: $(du -sh "$OUT" | cut -f1)"
