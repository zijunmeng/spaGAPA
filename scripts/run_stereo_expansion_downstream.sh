#!/usr/bin/env bash
# =============================================================================
# spaGAPA downstream for the 2026-09 stereo expansion samples.
# For each sample with binned_200 (apa_matrix/coordinates from stereo_bin200.py):
#   1. conformal uncertainty calibration (scripts/calibrate_uncertainty.py)
#   2. spaGAPA pipeline run (highres_fast preset) → domains + imputation
# Usage: bash run_stereo_expansion_downstream.sh          # all ready samples
#        bash run_stereo_expansion_downstream.sh <gsm>    # single sample
# =============================================================================
set -uo pipefail
HOSTNAME=$(hostname -s); [ "$HOSTNAME" = "S91" ] || { echo "not S91" >&2; exit 1; }
export OPENBLAS_NUM_THREADS=8
ROOT=/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA
PY=/home/mengzijun/anaconda3/envs/spagapa/bin/python
LOG=$ROOT/logs/20260912_stereo_downstream.log
ts() { date '+%F %T'; }

declare -A BIN TAG
BIN[gsm8882884]=$ROOT/pipeline_output/gse293464_retina/GSM8882884_binned
BIN[gsm8882885]=$ROOT/pipeline_output/gse293464_retina/GSM8882885_binned
BIN[gsm8882886]=$ROOT/pipeline_output/gse293464_retina/GSM8882886_binned
BIN[gsm8882887]=$ROOT/pipeline_output/gse293464_retina/GSM8882887_binned
BIN[gsm9770943]=$ROOT/pipeline_output/gse333693_thymus/GSM9770943_binned
TAG[gsm8882884]=gse293464_GSM8882884_s26_bms
TAG[gsm8882885]=gse293464_GSM8882885_s26_ra
TAG[gsm8882886]=gse293464_GSM8882886_s16_bms
TAG[gsm8882887]=gse293464_GSM8882887_s16_ra
TAG[gsm9770943]=gse333693_GSM9770943_thymus

# bin200 follower 的输出目录命名是 <sample>_scapatrap_raw/../binned_200；规范化到上表
for k in "${!BIN[@]}"; do
    d="${BIN[$k]}"
    [ -d "$d" ] || { alt="$(dirname "$d")/binned_200"; [ -d "$alt" ] && BIN[$k]="$alt"; }
done

run_sample() {
    local k=$1 d=${BIN[$1]} t=${TAG[$1]}
    [ -s "$d/apa_matrix.csv" ] && [ -s "$d/coordinates.csv" ] || { echo "$(ts) skip $k (no binned)"; return 0; }
    local out=$ROOT/pipeline_output/stereo_expansion_downstream/$t
    mkdir -p "$out/uncertainty"
    if [ ! -s "$out/uncertainty/calibration_summary.json" ]; then
        echo "$(ts) [$k] conformal calibration"
        $PY "$ROOT/scripts/calibrate_uncertainty.py" \
            --apa-matrix "$d/apa_matrix.csv" --coordinates "$d/coordinates.csv" \
            --output "$out/uncertainty" >> "$LOG" 2>&1 \
          && echo "$(ts) [$k] calibration OK" || echo "$(ts) [$k] calibration FAIL"
    fi
    if [ ! -d "$out/spagapa_run" ]; then
        echo "$(ts) [$k] spaGAPA pipeline"
        $PY - >> "$LOG" 2>&1 << PYEOF || { echo "$(ts) [$k] spaGAPA FAIL"; return 1; }
from spagapa import SpaGAPA
from spagapa.core import APADataset
ds = APADataset.from_csv(apa_matrix="$d/apa_matrix.csv", coordinates="$d/coordinates.csv")
res = SpaGAPA(accuracy="highres_fast").fit_transform(ds)
import os, json
os.makedirs("$out/spagapa_run", exist_ok=True)
with open("$out/spagapa_run/domains.csv", "w") as f:
    f.write("domain\\n")
    f.write("\\n".join(map(str, res.domains)) + "\\n")
meta = dict(n_domains=int(res.n_domains), recovered_shape=list(res.recovered.shape),
            uncertainty_shape=list(res.uncertainty.shape))
json.dump(meta, open("$out/spagapa_run/run_meta.json", "w"), indent=1)
print("spaGAPA $k done:", meta)
PYEOF
        echo "$(ts) [$k] spaGAPA OK"
    fi
}

if [ $# -ge 1 ]; then run_sample "$1"; else for k in "${!BIN[@]}"; do run_sample "$k"; done; fi
echo "$(ts) downstream pass complete"
