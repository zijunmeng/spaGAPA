#!/usr/bin/env bash
# Download benchmark datasets for spaGAPA real-data benchmark.
# Usage: bash scripts/download_benchmark_data.sh [output_dir]
#
# Downloads:
#   1. Mouse Olfactory Bulb (MOB) — 10x Visium
#   2. Human Brain Cortex        — 10x Visium
#   3. Mouse Embryo              — Stereo-seq (optional)

set -euo pipefail

OUTDIR="${1:-data/real}"
mkdir -p "$OUTDIR"

echo "=== Downloading benchmark datasets to $OUTDIR ==="

# ── 1. Mouse Olfactory Bulb (10x Visium, public) ──────────────────────────────
MOB_DIR="$OUTDIR/mob"
mkdir -p "$MOB_DIR"
echo ""
echo "[1/3] Mouse Olfactory Bulb (MOB)"

MOB_URL="https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Mouse_Olfactory_Bulb"
for f in \
  "filtered_feature_bc_matrix.h5" \
  "spatial.tar.gz"
do
  if [ ! -f "$MOB_DIR/$f" ]; then
    echo "  Downloading $f ..."
    wget -q -O "$MOB_DIR/$f" "${MOB_URL}/${f}" || \
      curl -sL -o "$MOB_DIR/$f" "${MOB_URL}/${f}"
  else
    echo "  $f already exists, skipping."
  fi
done

# Extract spatial folder
if [ -f "$MOB_DIR/spatial.tar.gz" ] && [ ! -d "$MOB_DIR/spatial" ]; then
  tar -xzf "$MOB_DIR/spatial.tar.gz" -C "$MOB_DIR"
fi

# ── 2. Human Brain Cortex (10x Visium, public) ────────────────────────────────
BRAIN_DIR="$OUTDIR/brain"
mkdir -p "$BRAIN_DIR"
echo ""
echo "[2/3] Human Brain Cortex"

BRAIN_URL="https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Brain_Section_1"
for f in \
  "filtered_feature_bc_matrix.h5" \
  "spatial.tar.gz"
do
  if [ ! -f "$BRAIN_DIR/$f" ]; then
    echo "  Downloading $f ..."
    wget -q -O "$BRAIN_DIR/$f" "${BRAIN_URL}/${f}" || \
      curl -sL -o "$BRAIN_DIR/$f" "${BRAIN_URL}/${f}"
  else
    echo "  $f already exists, skipping."
  fi
done

if [ -f "$BRAIN_DIR/spatial.tar.gz" ] && [ ! -d "$BRAIN_DIR/spatial" ]; then
  tar -xzf "$BRAIN_DIR/spatial.tar.gz" -C "$BRAIN_DIR"
fi

# ── 3. Mouse Embryo (optional, Stereo-seq) ────────────────────────────────────
echo ""
echo "[3/3] Mouse Embryo (Stereo-seq) — manual download required"
echo "  Please download from: https://db.cngb.org/stomics/mosta/"
echo "  Place files in: $OUTDIR/embryo/"

echo ""
echo "=== Download complete ==="
echo "Next step: run scripts/prepare_benchmark_data.py to extract APA matrices"
