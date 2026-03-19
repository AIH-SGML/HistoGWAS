#!/usr/bin/env bash
set -euo pipefail

METADATA=../../data/samples_metadata/thyroid_samples_metadata.csv
OUTDIR=../../data/slides/thyroid

# Optional positional arg:
#   all (default) -> download all slides listed in metadata
#   N             -> download only the first N slides
MAX_SLIDES_ARG="${1:-all}"

CMD=(
    python stage1_download_slides.py
    --samples-metadata "$METADATA"
    --outdir "$OUTDIR"
    --rewrite
)

if [[ "$MAX_SLIDES_ARG" != "all" ]]; then
    if [[ "$MAX_SLIDES_ARG" =~ ^[0-9]+$ ]] && (( MAX_SLIDES_ARG > 0 )); then
        CMD+=(--max_slides "$MAX_SLIDES_ARG")
    else
        echo "Usage: bash stage1_download_slides_run.sh [all|N]"
        echo "  all: download every slide in metadata"
        echo "  N  : download first N slides (N > 0)"
        exit 1
    fi
fi

"${CMD[@]}"
