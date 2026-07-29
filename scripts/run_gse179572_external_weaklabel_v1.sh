#!/usr/bin/env bash
set -euo pipefail

source ~/anaconda3/etc/profile.d/conda.sh
conda activate spagapa

REPO_DIR="/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
DATA_DIR="${REPO_DIR}/data/processed/gse179572_gsm5420751_scapatrap"
OUT_DIR="${REPO_DIR}/benchmark_results/real/gse179572_scapatrap_external_weaklabel_v1"

mkdir -p "${OUT_DIR}"
cd "${REPO_DIR}"

python scripts/run_external_bioml_validation.py \
  --dataset-dir "${DATA_DIR}" \
  --output-dir "${OUT_DIR}" \
  --dataset-name gse179572_gsm5420751_scapatrap \
  --n-genes 120 \
  --min-observed-spots 120 \
  --methods raw,stapaminer_knn_expression,spagapa_gp,spagapa_bioml \
  --layer-column marker_weak_label \
  --n-domains 8 \
  --expr-n-components 10 \
  --svapa-neighbors 6 \
  --fdr 0.05 \
  --logfc-threshold 0.5 \
  --bioml-rank 8 \
  --bioml-lambda-graph 0.5 \
  --bioml-lambda-l2 1e-2 \
  --bioml-max-iter 30 \
  --bioml-n-neighbors 15 \
  --bioml-blend 0.1 \
  --bioml-domain-method spectral \
  --bioml-spatial-weight 0.1 \
  --bioml-expression-weight 0.7 \
  --bioml-apa-weight 0.2 \
  --force
