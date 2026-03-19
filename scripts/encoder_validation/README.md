# Stage 2 · Encoder Selection

Utilities in this folder drive Stage 2 of the HistoGWAS release: filtering GTEx
expression matrices and benchmarking multiple histology encoders for
tile-to-gene prediction. The jobs assume a Slurm scheduler and the
`histoGWAS_2` conda environment.

## Environment

```bash
conda env create -f HistoGWAS_2.yml
conda activate histoGWAS_2
```

If you rename the environment, update the `CONDA_ENV` constant defined near the
top of each `*_run.py` launcher so that the submitted jobs activate the correct
conda env on the compute nodes.

## Data Requirements

- **GTEx TPM matrices** (`.gct.gz`): download per tissue from the GTEx portal.
- **Tile embeddings** (`.h5ad`): expected under
  `../../data/img_embedding/<model>_embedding/<Tissue>_img_embedding.h5ad`
  for direct and launcher defaults.
- **Train/test split**: CSV at
  `../../data/train_test_split.csv`.
- Optional caches such as `*_pc.h5ad` will be generated automatically to avoid
  recomputing PCA.

Both `*_run.py` launchers are currently legacy Slurm submitters with defaults
defined in the script body. For direct one-off runs, use `run_filter_gene.sh`
and `run_gene_prediction.sh`.

## Directory Contents

| File | Purpose |
| ---- | ------- |
| `1_filter_gene.py` | Filter bulk TPM matrices and retain highly variable genes as `.h5ad` files. |
| `1_filter_gene_run.py` | Legacy Slurm submitter for Step 1 with defaults in the script body. |
| `2_gene_prediction.py` | Fits GLMM models to predict gene expression from encoder embeddings. |
| `2_gene_prediction_run.py` | Legacy Slurm submitter for Step 2 with hardcoded tissue/model/path defaults. |
| `run_filter_gene.sh` | Simple bash wrapper that directly runs `1_filter_gene.py` for one tissue. |
| `run_gene_prediction.sh` | Simple bash wrapper that directly runs `2_gene_prediction.py` for one tissue/model. |
| `3_gene_prediction_encoder_selection.ipynb` | Notebook used during encoder model selection to review gene-prediction performance across backbones. |
| `HistoGWAS_2.yml` | Conda environment specification for the Stage 2 workflow. |

## Workflow

### 1. Filter GTEx expression
Use this step to convert per-tissue GTEx TPM matrices (`.gct.gz`) into filtered
gene-expression `.h5ad` files.

Input TPM files are expected in `../../data/tpm_gene/` with names like:
`gene_tpm_2017-06-05_v8_thyroid.gct.gz`.

If a tissue file is missing, `1_filter_gene_run.py` can download it from:
`https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/tpms-by-tissue/gene_tpm_2017-06-05_v8_{tissue.lower}.gct.gz`.

- Example: `https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/tpms-by-tissue/gene_tpm_2017-06-05_v8_thyroid.gct.gz`.

- Replace `{tissue}` with one of:
  - `Skin_Sun_Exposed_Lower_leg`
  - `Esophagus_Muscularis`
  - `Stomach`
  - `Nerve_Tibial`
  - `Colon_Transverse`
  - `Esophagus_Mucosa`
  - `Artery_Tibial`
  - `Breast_Mammary_Tissue`
  - `Adipose_Subcutaneous`
  - `Muscle_Skeletal`

GTEx also provides official downloads here:
`https://www.gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression#bulk_tissue_expression-gtex_analysis_v8-rna-seq`

Simple direct run (no Slurm):

```bash
./run_filter_gene.sh Thyroid
```

By default this reads from `../../data/tpm_gene/` and writes
`../../data/gene_expression/thyroid_gene_tpm.h5ad`.

If you need Slurm batch submission for Step 1, use `1_filter_gene_run.py` and
edit tissue/path defaults in that file.

Run a single tissue without Slurm by invoking the script directly:

```bash
python 1_filter_gene.py \
  --outdir /path/to/output/gene_expression \
  --tissue thyroid \
  --gene_file /path/to/tpm_gene/gene_tpm_2017-06-05_v8_thyroid.gct.gz
```

### 2. Evaluate encoders via gene prediction

Simple bash wrapper:

```bash
./run_gene_prediction.sh [TISSUE] [MODEL] [DATA_ROOT] [CLUSTER_I]
```

Example:

```bash
./run_gene_prediction.sh Thyroid retccl ../../data
```

Defaults:
- `TISSUE=Thyroid`
- `MODEL=retccl`
- `DATA_ROOT=../../data`
- `CLUSTER_I` optional

Legacy Slurm launcher (currently uses internal defaults from
`2_gene_prediction_run.py`):

```bash
python 2_gene_prediction_run.py
```

As an alternative to local Stage 2 preprocessing outputs, you can download the
required data from Zenodo (`https://zenodo.org/records/18773562`):

- **Tissue embeddings** (already processed PCA embeddings):
  `https://zenodo.org/records/18773562/files/{tissue}_img_embedding.h5ad?download=1`
- Replace `{tissue}` with one of:
  - `Skin_Sun_Exposed_Lower_leg`
  - `Esophagus_Muscularis`
  - `Stomach`
  - `Nerve_Tibial`
  - `Colon_Transverse`
  - `Esophagus_Mucosa`
  - `Artery_Tibial`
  - `Breast_Mammary_Tissue`
  - `Adipose_Subcutaneous`
  - `Muscle_Skeletal`
- Example:
  `https://zenodo.org/records/18773562/files/Stomach_img_embedding.h5ad?download=1`
  save the downloaded embedding at location `../../data/embedding/Thyroid/Thyroid_microns_192/summary`

- **Preprocessed gene expression**:
  `https://zenodo.org/records/18773562/files/gene_expression.zip?download=1`

Each job writes `${OUTPUT_ROOT}/gene_prediction/<tissue>_<encoder>_<cluster>_glmm_nComp_<n>.csv`
containing R², Spearman correlation, and p-values per gene, plus logs in
`${OUTPUT_ROOT}/gene_prediction/eval_logs/`.

For on-demand debugging, run the driver manually:

```bash
python 2_gene_prediction.py \
  --tissue thyroid \
  --model_type retccl \
  --hfile /path/to/img_embedding/retccl_embedding/Thyroid_img_embedding.h5ad \
  --efile /path/to/gene_expression/Thyroid_gene_tpm.h5ad \
  --outdir /path/to/output/gene_prediction \
  --train_test_split /path/to/train_test_split.csv \
  --cluster_i 0
```

`--cluster_i` (and optional `--grade`) limit the run to a specific Leiden
cluster; omit them to run the full tissue cohort.

## Tips

- Resource requests (`time`, `cpus`/`memory`, `qos`) are hardcoded near the top
  of each `*_run.py` launcher—edit those values for your Slurm partition.
- Slurm launchers create temporary sbatch scripts before submission; keep a copy
  if you need to inspect resolved commands for debugging.
- The GLMM stage normalises embeddings slide-wise and caches PCA projections as
  `<embedding>_pc.h5ad`. Delete the cache if you tweak PCA settings and need to
  recompute.

Review the `stderr` logs after each run for Scanpy warnings or missing data, and
use the notebook for downstream visualisation if needed.
