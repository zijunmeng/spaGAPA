#!/usr/bin/env bash
# Install stAPAminer + its missing CRAN deps into the shared R_442 lib.
# Hostname-aware (CLAUDE.md), timed, logged.
set -uo pipefail
HOSTNAME=$(hostname -s)
case "$HOSTNAME" in
    S90|S91|S97|S98) ;;
    *) echo "ERROR unknown host $HOSTNAME" >&2; exit 1 ;;
esac

RLIB=/s1/SHARE/01_software/R_442_SeuratV5/library
ZIP=/s1/SHARE/01_software/R_442_SeuratV5/library/stAPAminer-main.zip
# MUST activate the conda env so the conda-named compilers
# (x86_64-conda-linux-gnu-cc/c++/gfortran) are on PATH + build env vars set —
# running Rscript directly leaves source packages unable to compile.
source ~/anaconda3/etc/profile.d/conda.sh
conda activate r442
RSCRIPT=Rscript
echo "CC on PATH: $(command -v x86_64-conda-linux-gnu-cc || echo MISSING)"
LOGROOT=/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output/stapaminer_install_logs
mkdir -p "$LOGROOT"

echo "[$(date '+%F %T')] === STEP 1: install CRAN deps (fdm2id, clusterSim, ClusterR) ==="
T0=$(date +%s)
"$RSCRIPT" -e "
.libPaths(c('$RLIB', .libPaths()))
options(repos = c(CRAN = 'https://cloud.r-project.org'))
options(Ncpus = 8)
pkgs <- c('arules', 'arulesViz', 'fdm2id', 'clusterSim', 'ClusterR')
inst <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]
cat('to install:', paste(inst, collapse=', '), '\n')
if (length(inst) > 0) install.packages(inst, lib='$RLIB', dependencies=NA, quiet=FALSE)
for (p in pkgs) cat(sprintf('  %-12s %s\n', p, ifelse(requireNamespace(p, quietly=TRUE), paste0('OK ', as.character(packageVersion(p))), 'FAIL')))
" 2>&1 | tee "$LOGROOT/step1_deps.log"
T1=$(date +%s)
echo "[$(date '+%F %T')] STEP 1 wall=$((T1-T0))s"

echo "[$(date '+%F %T')] === STEP 2: install_local stAPAminer (dependencies=FALSE) ==="
T0=$(date +%s)
"$RSCRIPT" -e "
.libPaths(c('$RLIB', .libPaths()))
Sys.setenv(R_INSTALL_STAGED = 'TRUE')
devtools::install_local('$ZIP', lib='$RLIB', dependencies=FALSE, upgrade='never', quiet=FALSE)
" 2>&1 | tee "$LOGROOT/step2_stapaminer.log"
T1=$(date +%s)
echo "[$(date '+%F %T')] STEP 2 wall=$((T1-T0))s"

echo "[$(date '+%F %T')] === STEP 3: verify stAPAminer loads ==="
"$RSCRIPT" -e "
.libPaths(c('$RLIB', .libPaths()))
ok <- requireNamespace('stAPAminer', quietly=TRUE)
cat('stAPAminer loadable:', ok, '\n')
if (ok) {
  suppressMessages(library(stAPAminer))
  cat('version:', as.character(packageVersion('stAPAminer')), '\n')
  cat('exports (first 20):', paste(head(ls('package:stAPAminer'), 20), collapse=', '), '\n')
}
" 2>&1 | tee "$LOGROOT/step3_verify.log"
echo "[$(date '+%F %T')] === DONE ==="