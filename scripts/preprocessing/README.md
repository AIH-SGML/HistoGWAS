# HistoGWAS Slide Preprocessing

This directory bundles slide preprocessing pipeline that turns tissue manifests into slide tiles and tile embeddings ready for downstream association analysis. The steps are designed to run on SLURM array jobs, but can also be executed locally.

## Pipeline Overview
| Stage | Primary script(s) | Purpose | Key outputs |
| --- | --- | --- | --- |
| 1 | `stage1_download_slides.py` | Split a tissue manifest, download SVS slides, log status | `slides/<slide_id>.svs`, `summary/summary_jobXXXX.tsv` |
| 1b | `stage1_download_slides_merge_summary.py` | Merge per-job summaries | `summary/summary.tsv` |
| 2 | `stage2_tiling.py` | Generate PNG tiles and per-slide tile metadata | `tiles/<slide_id>/*.png`, `tsvs/<slide_id>.tsv`, `logs/*.txt`, optional `plots/*.png` |
| 2b | `stage2_tiling_summary.py` | Combine slide-level TSVs into a cohort table | `tiles.tsv` with per-tile metadata + file paths |
| 3 | `stage3_tile_features.py` | Embed tiles with a selected encoder | `runs/<njobs>_<job>.h5ad`, `logs/*.h5ad` |
| 3b | `stage3_tile_features_summary.py` | Concatenate embeddings and compute Scanpy QC reductions | `summary_normal.h5ad`, `summary_scanpy.h5ad` |

Companion `*_run.sh` wrappers capture the defaults (`OUTPUT_ROOT=../../data`) and common CLI flags. Edit those scripts to match your tissue and cluster configuration before launching.

## Dependencies
- Stage 1 & 2: Python ≥3.8, `pandas`, `numpy`, `openslide-python`, `opencv-python`, `scikit-image`, `tqdm`, plus OpenSlide system libraries and `wget`.
- Stage 3: Stage 1 & 2 requirements **plus** `torch`, `torchvision`, `anndata`, `scanpy`, `transformers`, `eidenalg`, `igraph` (for `plip`), and `pl-bolts` (for `simclr`). Several model options expect local checkpoints; update the hard-coded paths in `stage3_tile_features.py` if your environment differs.

Activate an environment that satisfies the dependencies listed above before running each stage.

## Stage 1 – Download Slides
**Inputs**
- A sample metadata CSV with a `Tissue Sample ID` column.
- An output directory where slides and logs can be written (e.g. `../../data/slides/thyroid`).

**Command**
```bash
bash stage1_download_slides_run.sh
```
Update `METADATA` and `OUTDIR` inside the script before launching. By default, the wrapper downloads all slides listed in the metadata file and passes `--rewrite`.

To download only a subset for smoke testing:
```bash
bash stage1_download_slides_run.sh 50
```
This downloads only the first `50` slides from the metadata file.

If you need array-job sharding, call `stage1_download_slides.py` directly with `--n_jobs` and `--job_i`:
```bash
python stage1_download_slides.py \
  --samples-metadata ../../data/samples_metadata/thyroid_samples_metadata.csv \
  --outdir ../../data/slides/thyroid \
  --n_jobs 10 \
  --job_i 2 \
  --rewrite
```

**Outputs**
- `slides/` – downloaded SVS files.
- `summary/summary_jobXXXX_of_YYYY.tsv` – job-level download status log.

### Merge Stage 1 Summaries
After all array jobs finish:
```bash
bash stage1_download_slides_merge_summary_run.sh
```
This creates `summary/summary.tsv`, which feeds Stage 2.

## Stage 2 – Tile Slides
**Inputs**
- `summary.tsv` from Stage 1 (tab-separated, with a `path` column pointing to SVS files).
- Output root (e.g. `../../data/tiles/thyroid/thyroid_microns_192`).

**Command**
```bash
bash stage2_tiling_run.sh
```
Adjust `OUTPUT_ROOT`, `TISSUE`, `MICRONS`, `FG_MIN`, `N_JOBS`, and `job_i` in the script. The wrapper already enables `--export_tiles`; add `--make_plots` directly in the script if you need diagnostics.
Key flags:
- `--fract_fg_min`: minimum foreground fraction to retain a tile (default `0.5`).
- `--attempted_microns`: tile edge size in microns (default `192`).
- `--export_tiles`: save PNG tiles under `tiles/<slide_id>/` (required for Stage 3).
- `--make_plots` (optional): emit diagnostic plots to `plots/` for the first `--tissue_count` slides.
- `--njobs`/`--job_i`: split slides across processes; run each job index once.

**Outputs**
- `logs/<slide_id>.txt`: per-slide run log.
- `tiles/<slide_id>/<tile>.png`: PNG tiles (only when `--export_tiles` is set).
- `tsvs/<slide_id>.tsv`: tile metadata including coordinates, QC metrics, and PNG paths.
- `plots/*.png`: optional overlays when `--make_plots` is enabled.

### Aggregate Stage 2 Tiles
Combine slide-level TSVs into a single table required for Stage 3:
```bash
bash stage2_tiling_summary_run.sh
```
This writes `tiles.tsv` alongside the `tiles/` directory.

## Stage 3 – Extract Tile Features
Model sources (download/checkpoint reference links):
- RetCCL: https://github.com/Xiyue-Wang/RetCCL
- SimCLR ImageNet checkpoint: https://pl-bolts-weights.s3.us-east-2.amazonaws.com/simclr/bolts_simclr_imagenet/simclr_imagenet.ckpt
- KimiaNet: https://github.com/KimiaLabMayo/KimiaNet
- CTransPath: https://github.com/Xiyue-Wang/TransPath
- PLIP: https://github.com/PathologyFoundation/plip

Set checkpoint paths in `stage3_tile_features_run.sh` using:
- `SIMCLR_CKPT`
- `KIMIANET_CKPT`
- `CTRANSPATH_CKPT`
- `RETCCL_CKPT`
- `AUTOENCODER_CKPT`


**Inputs**
- `tiles.tsv` produced by Stage 2.
- Tile PNGs in `tiles/<slide_id>/*.png`.
- Output directory (e.g. `../../data/embedding/thyroid/thyroid_microns_192`).

**Command**
```bash
bash stage3_tile_features_run.sh
```
Tune `OUTPUT_ROOT`, `TISSUE`, `MICRONS`, `MODEL_TYPE`, `N_JOBS`, `BATCH_SIZE`, and `job_i` as needed. Also set the relevant `*_CKPT` variable for your selected model.
Model choices for `--model_type`: `retccl` (default), `simclr`, `kimiaNet`, `ctranspath`, `Autoencoder`, `plip`. `plip` downloads from HuggingFace and does not require a local checkpoint path. The `--tissue` string is retained for metadata compatibility.

**Outputs**
- `runs/<njobs>_<job>.h5ad`: per-job AnnData file containing embeddings and tile metadata.
- `logs/<njobs>_<job>.h5ad`: plain-text log (legacy suffix) for the matching job index.

### Summarise Tile Embeddings
Concatenate per-job outputs and compute Scanpy reductions:
```bash
bash stage3_tile_features_summary_run.sh
```
Two files are produced:
- `summary_normal.h5ad`: raw concatenated embeddings.
- `summary_scanpy.h5ad`: PCA, neighbor graph, UMAP, and Leiden clustering annotations.

## Helper Run Scripts
`stage*_run.sh` files wrap the commands above with example defaults. Copy or modify them when setting up SLURM array jobs. Each script assumes directories relative to this folder.

## Tips
- Use the `--debug` flag in any stage to drop into `pdb` before processing begins.
- Array jobs rely on `--njobs`/`--job_i` consistency across stages; record which job indices succeeded before resuming.
