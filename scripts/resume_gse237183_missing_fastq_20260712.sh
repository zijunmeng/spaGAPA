#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-/s1/SHARE/mengzijun/01_project/26_spaGAPA}"
DATA_DIR="$ROOT/data/raw/gse237183"
FASTQ_DIR="$DATA_DIR/fastq_ena"
LOG_DIR="$DATA_DIR/logs"
MANIFEST="$LOG_DIR/ena_fastq_manifest.tsv"
TARGET_URLS="$LOG_DIR/ena_fastq_resume_missing_20260712.urls.txt"
ARIA_LOG="$LOG_DIR/ena_fastq_targeted_resume_20260712.aria2c.log"
NOHUP_LOG="$ROOT/logs/20260712_gse237183_fastq_targeted_resume.nohup.log"
UNIT_NAME="spagapa_gse237183_resume_20260712"

MODE="${1:-systemd}"

usage() {
  cat <<EOF
Usage:
  $0 [systemd|nohup|foreground|check|monitor]

Modes:
  systemd     Submit aria2c through systemd-run --user. Recommended.
  nohup       Submit aria2c through nohup.
  foreground  Run aria2c in the current shell.
  check       Regenerate target list and report missing/partial files.
  monitor     Show process/log/file-status summary.
EOF
}

build_target_list() {
  mkdir -p "$FASTQ_DIR" "$LOG_DIR" "$ROOT/logs"
  python - "$MANIFEST" "$FASTQ_DIR" "$TARGET_URLS" <<'PY'
import csv
import sys
from pathlib import Path

manifest = Path(sys.argv[1])
fastq_dir = Path(sys.argv[2])
target_urls = Path(sys.argv[3])

rows = list(csv.DictReader(manifest.open(), delimiter="\t"))
targets = []
ok = 0

for row in rows:
    filename = row["filename"]
    url = row["url"]
    expected = int(row["bytes"])
    path = fastq_dir / filename
    aria2 = Path(str(path) + ".aria2")
    if not path.exists():
        targets.append((filename, url, 0, expected, False, "missing"))
    else:
        actual = path.stat().st_size
        if actual != expected or aria2.exists():
            reason = "partial_or_aria2"
            targets.append((filename, url, actual, expected, aria2.exists(), reason))
        else:
            ok += 1

target_urls.write_text(
    "\n".join(item[1] for item in targets) + ("\n" if targets else ""),
    encoding="utf-8",
)

print(f"expected_fastq={len(rows)}")
print(f"ok_fastq={ok}")
print(f"target_fastq={len(targets)}")
print(f"target_urls={target_urls}")
for filename, _url, actual, expected, has_aria2, reason in targets:
    print(f"{reason}\t{filename}\tactual={actual}\texpected={expected}\taria2={has_aria2}")
PY
}

aria2_cmd=(
  /usr/bin/aria2c
  "--input-file=$TARGET_URLS"
  "--dir=$FASTQ_DIR"
  "--check-certificate=false"
  "--max-concurrent-downloads=4"
  "--max-connection-per-server=2"
  "--split=2"
  "--min-split-size=20M"
  "--continue=true"
  "--file-allocation=none"
  "--auto-file-renaming=false"
  "--allow-overwrite=false"
  "--max-tries=0"
  "--retry-wait=30"
  "--console-log-level=warn"
  "--download-result=full"
  "--summary-interval=120"
  "--log=$ARIA_LOG"
)

monitor() {
  echo "== systemd unit =="
  systemctl --user status "$UNIT_NAME.service" --no-pager 2>/dev/null || true
  echo
  echo "== aria2c process =="
  pgrep -af "ena_fastq_resume_missing_20260712|ena_fastq_targeted_resume_20260712|aria2c" || true
  echo
  echo "== nohup log =="
  tail -n 60 "$NOHUP_LOG" 2>/dev/null || true
  echo
  echo "== aria2 log =="
  tail -n 80 "$ARIA_LOG" 2>/dev/null || true
  echo
  echo "== remaining targets =="
  build_target_list
}

case "$MODE" in
  systemd)
    build_target_list
    if [[ ! -s "$TARGET_URLS" ]]; then
      echo "No missing/partial FASTQ remains."
      exit 0
    fi
    systemd-run --user \
      --unit="$UNIT_NAME" \
      --collect \
      --working-directory="$ROOT" \
      "${aria2_cmd[@]}"
    echo "Submitted systemd user unit: $UNIT_NAME.service"
    echo "Monitor:"
    echo "  systemctl --user status $UNIT_NAME.service --no-pager"
    echo "  tail -f $ARIA_LOG"
    ;;
  nohup)
    build_target_list
    if [[ ! -s "$TARGET_URLS" ]]; then
      echo "No missing/partial FASTQ remains."
      exit 0
    fi
    nohup "${aria2_cmd[@]}" > "$NOHUP_LOG" 2>&1 &
    echo "Submitted nohup PID: $!"
    echo "Monitor:"
    echo "  tail -f $NOHUP_LOG"
    echo "  tail -f $ARIA_LOG"
    ;;
  foreground)
    build_target_list
    if [[ ! -s "$TARGET_URLS" ]]; then
      echo "No missing/partial FASTQ remains."
      exit 0
    fi
    "${aria2_cmd[@]}"
    ;;
  check)
    build_target_list
    ;;
  monitor)
    monitor
    ;;
  -h|--help|help)
    usage
    ;;
  *)
    usage >&2
    exit 2
    ;;
esac
