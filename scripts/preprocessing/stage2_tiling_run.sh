
# --- user configuration ---
OUTPUT_ROOT=../../data
TISSUE="thyroid"
MICRONS=192
FG_MIN=0.5
N_JOBS=1
TISSUE_COUNT=10

SLIDES=$OUTPUT_ROOT/slides/$TISSUE/summary/summary.tsv
# SLIDES='/lustre/groups/casale/datasets/gtex/histology/slides_/Thyroid/summary.tsv'
OUTDIR=$OUTPUT_ROOT/tiles/thyroid/${TISSUE}_microns_${MICRONS}

python stage2_tiling.py \
    --slides $SLIDES \
    --outdir $OUTDIR \
    --njobs $N_JOBS \
    --job_i 0 \
    --fract_fg_min $FG_MIN \
    --tissue_count $TISSUE_COUNT \
    --tissue_name $TISSUE \
    --attempted_microns $MICRONS \
    --export_tiles 

