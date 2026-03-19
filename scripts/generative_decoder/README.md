# PGAN Decoder for HistoGWAS

This module fine-tunes a Progressive GAN (PGAN) decoder on tile embeddings derived from GTEx histology slides. The trained generator is used to visualise clusters, interpolate between embeddings, and support downstream gene-prediction analyses within the HistoGWAS release.

## Directory Guide
- `config/` – sample JSON configs providing dataset locations and training hyperparameters (e.g. `config_Thyroid.json`).
- `train.py` / `train_pgan.sh` – CLI entry points for launching PGAN, StyleGAN, or DCGAN training runs on a local GPU.
- `experiments/train_PGAN_gpu.py` – SLURM launcher that creates per-tissue configs and submits PGAN jobs to the cluster queue.
- `experiments/train_PGAN_gpu_test.py` – lightweight helper for running the same configuration interactively without queue submission.
- `datasets.py` – AnnData-backed dataset loader for histology embeddings.
- `models/`, `visualization/` – forked components from the Facebook `pytorch_GAN_zoo` implementation.
- `eval.py`, `get_trait_to_model.py`, `save_feature_extractor.py` – utilities for checkpoint inspection, trait-to-model lookup, and feature export.
- `Notebooks/` – interactive demos (e.g. `Demo_pgan.ipynb`) showing interpolation and figure generation.
- `experiments/` – default location for logs and checkpoints produced by manual runs (`train_pgan.sh`).
- `../../../../output/PGAN/` – default output root used by the SLURM launchers.

## Prerequisites
- Create and activate the conda environment defined in `HistoGWAS_2.yml` from the repository root (e.g. `conda env create -f HistoGWAS_2.yml` and `conda activate HistoGWAS_2`).
- Python environment satisfying `requirements.txt` plus the broader HistoGWAS dependencies (PyTorch ≥1.10, torchvision, h5py, AnnData/scanpy stack).
- Access to AnnData `.h5ad` files containing tile embeddings and metadata. Each file must expose:
  - `adata.X`: numeric embedding matrix (samples × embedding dimension).
  - `adata.obs`: metadata, including `path` with image tile paths for visualisation.
- Editable installations of `emb_gwas`, `mtgwas`, and `src/` from the repository root (`pip install -e ...`).

## Configure a Run
1. Copy one of the configs in `config/` and update the following keys:
   - `pathDB`: absolute path to your `.h5ad` file.
   - `dimEmb`: embedding dimensionality (must match `adata.X.shape[1]`).
   - Optional: override `config.maxIterAtScale` or other trainer parameters to match your dataset size.
2. Optional environment variables used by the launcher script:
   - `OUTDIR`: directory for checkpoints and logs (`experiments/` by default for `train_pgan.sh`).
   - `CONFIG`: path to the JSON config (defaults to `config/config.json`).
   - `MODEL`: GAN variant to train (`PGAN`, `StyleGAN`, or `DCGAN`).
   - `EXPNAME`: custom run name (propagated to checkpoints).
   - `DIMEMB`: override embedding dimension without editing the config file.

## Launch Training (Local GPU)
Use the shell wrapper when running on a workstation or interactive GPU node:
```bash
OUTDIR=experiments/thyroid \
CONFIG=config/config_Thyroid.json \
EXPNAME=ThyroidPGAN \
./train_pgan.sh
```
The wrapper script creates the output directory and calls `train.py MODEL_NAME -c CONFIG --dir OUTDIR --no_vis ...`. Use `DIMEMB=128` if you need to override the embedding dimensionality at runtime.

## Launch Training (SLURM GPU Queue)
The new SLURM utilities automate config creation, logging, and queue submission for batch runs.

### 1. Prepare the tissue list
Edit `experiments/train_PGAN_gpu.py` and populate `tissue_list` with the tissues you wish to train. Each entry triggers:
- Creation of `config/config_<TISSUE>.json` with the shared hyperparameters and `dimEmb` set to 64 by default.
- Submission of an `sbatch` job that runs `train.py PGAN` using the generated config.

Adjust `tissue_hyperparameter[...]`, `OUTPUT_ROOT`, and SLURM options (`queue`, `gpu`, `time`, `memory`, etc.) to match your cluster policy. Logs are written to `<OUTPUT_ROOT>/eval_logs/stderr_PGAN.txt` and `stdout_PGAN.txt`.

### 2. Submit the jobs
From the `experiments/` directory (or via a module load script), run:
```bash
python train_PGAN_gpu.py
```
The script creates a transient `submit.sh`, submits it with `sbatch`, and deletes the file afterwards. Ensure `condaenv` matches the environment name you created earlier.

### 3. Dry-run locally
If you want to validate the config before queueing, use:
```bash
python train_PGAN_gpu_test.py
```
This helper writes `config/config.json`, then launches `train.py` directly without SLURM.

## Monitoring & Outputs
- Checkpoints (`*.pth`) and training curves are written under the chosen `OUTDIR` (default: `experiments/` for manual runs or `../../../../output/PGAN/` for SLURM runs).
- `train.py` can resume from the latest checkpoint automatically; add `--restart` to ignore existing weights.
- Use `python eval.py --help` to list evaluation options for generating samples or computing metrics from saved checkpoints.

## Demo & Visualisation
- After training, open `Notebooks/Demo_pgan.ipynb` to reproduce interpolation figures and sample grids. Update the notebook cells with your checkpoint path and dataset configuration.

## Tips
- Ensure the AnnData file fits in memory; the supplied configs target the low-memory Thyroid embeddings used in the release. Consider sub-sampling or sharding for larger cohorts.
- If Visdom is unavailable on your cluster, keep the `--no_vis` flag or switch to `--np_vis` to log snapshots without a server.
- Record the git commit and config JSON alongside generated figures for reproducibility.
- When using the SLURM launcher, the job inherits `CUDA_VISIBLE_DEVICES` from the scheduler; adjust the `--gres` option to control GPU counts across clusters.

## Attribution
The trainers and model definitions are adapted from [facebookresearch/pytorch_GAN_zoo](https://github.com/facebookresearch/pytorch_GAN_zoo). Please cite the original project alongside the HistoGWAS manuscript when using these artifacts.
