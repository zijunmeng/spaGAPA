#!/usr/bin/env bash
set -euo pipefail

source ~/anaconda3/etc/profile.d/conda.sh
conda activate spagapa

REPO_DIR="/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
DATA_DIR="${REPO_DIR}/data/processed/gse179572_gsm5420751_scapatrap"
OUT_DIR="${REPO_DIR}/benchmark_results/real/gse179572_scapatrap_formal_v3_fast_threadlimit"

mkdir -p "${OUT_DIR}"
cd "${REPO_DIR}"

python scripts/run_stapaminer_mob_benchmark.py \
  --data-dir "${DATA_DIR}" \
  --output-dir "${OUT_DIR}" \
  --n-genes 120 \
  --min-observed-spots 120 \
  --mask-fraction 0.2 \
  --seeds 42,43,44 \
  --knn-k 10 \
  --gp-kernel matern \
  --gp-alpha 1e-10 \
  --gp-variant spatial \
  --expr-n-components 10 \
  --n-domains 8 \
  --mask-types random,spatial_block,spatial_block_large,ring_sector,low_coverage \
  --layer-column marker_weak_label \
  --include-fast-gp \
  --fast-gp-max-train-points 0 \
  --fast-gp-n-jobs 8 \
  --fast-gp-n-restarts 1 \
  --include-gp-bioml-domains \
  --include-bioml \
  --bioml-rank 8 \
  --bioml-lambda-graph 0.5 \
  --bioml-lambda-l2 1e-2 \
  --bioml-max-iter 30 \
  --bioml-n-neighbors 15 \
  --bioml-blend 0.1 \
  --bioml-domain-method spectral \
  --bioml-spatial-weight 0.1 \
  --bioml-expression-weight 0.7 \
  --bioml-apa-weight 0.2
