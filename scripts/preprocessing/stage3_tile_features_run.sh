#!/bin/bash

# --- user configuration ---
OUTPUT_ROOT=../../data
TISSUE="thyroid"
MICRONS=192
MODEL_TYPE="retccl" # 'simclr', 'kimiaNet', 'plip', 'Autoencoder'
N_JOBS=10
DIMENSION=256
BATCH_SIZE=128
SIMCLR_CKPT=""
KIMIANET_CKPT=""
CTRANSPATH_CKPT=""
RETCCL_CKPT="../../data/pretainedModels/RetCCL/best_ckpt.pth"
AUTOENCODER_CKPT=""

TSV=$OUTPUT_ROOT/tiles/${TISSUE}/${TISSUE}_microns_${MICRONS}/tiles.tsv
OUTDIR=$OUTPUT_ROOT/embedding/${TISSUE}/${TISSUE}_microns_${MICRONS}

mkdir -p $OUTDIR

# --- run feature extraction ---
python stage3_tile_features.py \
    --tiles $TSV \
    --outdir $OUTDIR \
    --njobs $N_JOBS \
    --job_i 0 \
    --num_workers 8 \
    --dimension $DIMENSION \
    --model_type $MODEL_TYPE \
    --tissue ${TISSUE}_microns_${MICRONS} \
    --batch_size $BATCH_SIZE \
    --simclr_ckpt "$SIMCLR_CKPT" \
    --kimianet_ckpt "$KIMIANET_CKPT" \
    --ctranspath_ckpt "$CTRANSPATH_CKPT" \
    --retccl_ckpt "$RETCCL_CKPT" \
    --autoencoder_ckpt "$AUTOENCODER_CKPT"
