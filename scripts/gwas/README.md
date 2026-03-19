# HistoGWAS Association Pipeline

This directory contains utilities for running HistoGWAS association tests on
GTEx histology embeddings. The scripts use the in-repo `histogwas` package
(including `histogwas.emb_gwas`) and assume access to Stage 3 embedding
outputs, GTEx variants, and matching covariates. Alternatively, you can use the
already preprocessed RetCCL embeddings for each tissue from the Zenodo link
listed below.

## Environment Setup

Create and activate the shared conda environment before running any step:

```bash
conda env create -f HistoGWAS_2.yml
conda activate HistoGWAS_2
```

Install the bundled `histogwas` package in editable mode so the scripts resolve
the local modules:

```bash
pip install --upgrade pip
pip install -e ../..
```

### chiscore dependency

`emb_gwas` depends on the `chiscore` and `chi2comb` wheels. On Apple-silicon
macOS you may need Homebrew's GCC toolchain:

```bash
brew install gcc
pip install chi2comb
pip install chiscore
```

On Linux clusters, install the same wheels inside the activated `HistoGWAS_2`
environment.

## Directory Contents

- `hgwas_all_cluster.py` — command-line driver that matches embeddings, genotype
  data, and covariates before invoking the `emb_gwas` association routine. It
  expects explicit inputs via `--hfile`, `--bfile`, `--pcfile`, and `--covfile`.
- `hgwas_all_clusterwise_run.py` — SLURM launcher that iterates over a
  tissue-to-cluster mapping and submits one job per cluster using the driver
  above. Update the hard-coded paths near the top of the file to reflect your
  storage layout.

## Data Requirements

- AnnData `.h5ad` embedding files exported from Stage 3. They must include the
  cluster column referenced by `--cluster_type` (default `leiden_0.5`) and a
  subject/bag identifier column (`SID` or `slide`).
- Tissue embeddings (already processed PCA embeddings) are available on Zenodo:
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
- GTEx genotype data in PLINK format (`.bed/.bim/.fam`) and the corresponding
  eigenvectors (`.pca.eigenvec`).
- GTEx phenotype covariates table (for example
  `GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt`).

## Running the CLI

Invoke the association driver directly from this directory. Example:

```bash
python hgwas_all_cluster.py \
  --outdir ../../../output/association \
  --tissue Thyroid \
  --cluster_i 0 \
  --pcfile ../../data/wgs/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.pca.eigenvec \
  --bfile ../../data/wgs/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF02 \
  --hfile ../../data/embedding/thyroid/thyroid_microns_192/summary_scanpy.h5ad \
  --covfile ../../data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt \
  --cluster_type leiden_0.5
```

The driver writes `*_hgwas.csv` summaries and (unless `--no-plot` is set)
Manhattan plots under the supplied `--outdir`.

Use the launcher to batch across clusters via Slurm:

```bash
python hgwas_all_clusterwise_run.py
```

Adjust the resource requests (`time`, `qos`, `memory`, etc.) and path templates
inside the launcher before submitting jobs.

## Input Checklist

Confirm that the embedding file:

1. Contains the cluster labels column that matches `--cluster_type`.
2. Has a subject identifier column (`SID` by default).
3. Provides the embeddings expected by the script (`adata.X` or
   `adata.obsm['X_pca']`).
4. Aligns with the individuals present in the PLINK genotype and covariate
   tables passed on the CLI.

Outputs can be organised per tissue by pointing `--outdir` to a tissue-specific
folder (the script will create it). Update the defaults in
`hgwas_all_clusterwise_run.py` to mirror your environment before submitting
batch jobs.
