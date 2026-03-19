
OUTPUT_ROOT=../../data
TISSUE="thyroid"
MICRONS=192

SLIDES=$OUTPUT_ROOT/slides/$TISSUE/summary/summary.tsv
OUTDIR=$OUTPUT_ROOT/tiles/$TISSUE/${TISSUE}_microns_${MICRONS}

# --- run summary aggregation ---
python stage2_tiling_summary.py \
    --slides $SLIDES \
    --outdir $OUTDIR
