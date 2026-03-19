#!/bin/bash

# --- user configuration ---
OUTPUT_ROOT=../../data
TISSUE="thyroid"
MICRONS=192
MODEL_TYPE="retccl"

INPUT=$OUTPUT_ROOT/embedding/${TISSUE}/${TISSUE}_microns_${MICRONS}/runs
OUTDIR=$OUTPUT_ROOT/embedding/${TISSUE}/${TISSUE}_microns_${MICRONS}/summary

mkdir -p $OUTDIR

# --- run feature summary ---
python stage3_tile_features_summary.py \
    --input_path $INPUT \
    --outdir $OUTDIR \
    --tissue $TISSUE
