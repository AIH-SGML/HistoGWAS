#!/usr/bin/env bash
set -euo pipefail

# Simple direct runner for 2_gene_prediction.py.
# Defaults mirror the local data layout under ../../data.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TISSUE="${1:-Thyroid}"
MODEL="${2:-retccl}"
DATA_ROOT="${3:-$SCRIPT_DIR/../../data}"
CLUSTER_I="${4:-}"
MICRONS="192"

HFILE="$DATA_ROOT/embedding/${TISSUE}/${TISSUE}_microns_${MICRONS}/summary/${TISSUE}_img_embedding.h5ad"
EFILE="$DATA_ROOT/gene_expression/${TISSUE,,}_gene_tpm.h5ad"
OUTDIR="$DATA_ROOT/gene_prediction"
SPLIT_FILE="$DATA_ROOT/train_test_split.csv"

if [[ ! -f "$HFILE" ]]; then
  echo "Embedding file not found: $HFILE" >&2
  exit 1
fi
if [[ ! -f "$EFILE" ]]; then
  echo "Expression file not found: $EFILE" >&2
  exit 1
fi
if [[ ! -f "$SPLIT_FILE" ]]; then
  echo "Train/test split file not found: $SPLIT_FILE" >&2
  exit 1
fi

mkdir -p "$OUTDIR"

CMD=(
  python "$SCRIPT_DIR/2_gene_prediction.py"
  --tissue "$TISSUE"
  --model_type "$MODEL"
  --hfile "$HFILE"
  --efile "$EFILE"
  --outdir "$OUTDIR"
  --train_test_split "$SPLIT_FILE"
)

if [[ -n "$CLUSTER_I" ]]; then
  CMD+=(--cluster_i "$CLUSTER_I")
fi

echo "Running tissue=$TISSUE model=$MODEL"
echo "Embedding:  $HFILE"
echo "Expression: $EFILE"
echo "Output:     $OUTDIR"

"${CMD[@]}"
