#!/usr/bin/env bash
set -euo pipefail


OUTDIR="${OUTDIR:-out}"              # or pass: OUTDIR=out ./train_pgan.sh
CONFIG="${CONFIG:-config/config.json}"
MODEL="${MODEL:-PGAN}"
EXPNAME="${EXPNAME:-HistoDecoder}"    # optional 

# ---- run ----
echo "Running $MODEL with $CONFIG -> $OUTDIR"
mkdir -p "$OUTDIR"

python train.py "$MODEL" \
  -c "$CONFIG" \
  --dir "$OUTDIR" \
  --no_vis
  ${EXPNAME:+-n "$EXPNAME"} \
  ${DIMEMB:+--dimEmb "$DIMEMB"}   # DIMEMB optional
