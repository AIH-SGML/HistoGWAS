#!/usr/bin/env bash
set -euo pipefail

# Simple direct runner for 1_filter_gene.py.
# Defaults mirror 1_filter_gene_run.py:
#   OUTPUT_ROOT = ../../data
#   TPM files in ../../data/tpm_gene
# Example tissue: Thyroid

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TISSUE="${1:-Thyroid}"
TISSUE_LOWER="$(echo "$TISSUE" | tr '[:upper:]' '[:lower:]')"

TPM_DIR="$SCRIPT_DIR/../../data/tpm_gene"
OUTDIR="$SCRIPT_DIR/../../data/gene_expression"

mkdir -p "$OUTDIR"

GENE_FILE="$(find "$TPM_DIR" -maxdepth 1 -type f -iname "*v8_${TISSUE_LOWER}.*" | head -n 1)"

if [[ -z "${GENE_FILE:-}" ]]; then
  echo "No TPM file found for tissue '$TISSUE_LOWER' in: $TPM_DIR" >&2
  exit 1
fi

echo "Running tissue=$TISSUE_LOWER"
echo "TPM file: $GENE_FILE"
echo "Output:   $OUTDIR"

python "$SCRIPT_DIR/1_filter_gene.py" \
  --outdir "$OUTDIR" \
  --tissue "$TISSUE_LOWER" \
  --gene_file "$GENE_FILE"
