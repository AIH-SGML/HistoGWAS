# HistoGWAS Downstream Characterization

Notebook-driven analyses for inspecting HistoGWAS release outputs once
preprocessing and GWAS runs have finished. Use these materials to visualise
semantic decoder effects, evaluate gene-expression prediction results, perform
pathway enrichment, and assemble inputs for colocalisation.

## Layout

- `notebooks/`
  - `README.md` – quick catalogue of the available notebooks and helper scripts.
  - `visualization_thyroid.ipynb` – semantic decoder visualisations and
    manuscript-style figures for the thyroid cohort.
  - `analysisng_gene_prediction.ipynb` – exploratory diagnostics for
    expression-prediction performance across tissues and models.
  - `pathway_enrichment_analysis.ipynb` – pathway/gene-set summary tables.
  - `colcalization_with_eqtls.ipynb` – harmonisation of GWAS and GTEx signals
    prior to colocalisation.
  - `coloc_script.R` – thin wrapper around the `coloc` R package to convert log
    Bayes factors into posterior probabilities.
  - `utils.py` – shared loaders for GTEx expression matrices and matching
    helpers consumed by the notebooks.

## Dependencies

Activate the release analysis environment defined in `histogwas.yml`, which
provides Python ≥3.8 plus `anndata`, `scanpy`, `torch`, `pandas`, `numpy`,
`matplotlib`, `seaborn`, `tqdm`, `scikit-learn`, and the internal `histogwas`
package. The colocalisation notebooks additionally rely on `limix-core`.

Install R ≥4.0 with the `coloc` package to run `notebooks/coloc_script.R`.

## Required Inputs

The notebooks assume access to the canonical release filesystem layout:

- Stage 2/Stage 3 embedding AnnData files, e.g.
  `/lustre/.../stage2/<TISSUE>/embedding/<MODEL>/summary/summary_scanpy.h5ad` or
  the low-memory exports referenced in the notebooks.
- Stage 2 tile manifests (`tiles.tsv`) with tile-level metadata and cluster IDs.
- Processed gene expression matrices (`gene_tpm_*.h5ad`) from the Emb-GWAS
  preprocessing pipeline.
- GWAS summary statistics and auxiliary tables produced by the association
  stage.

Adjust path constants in the notebooks (usually near the top of each file) if
your data reside outside the shared Casale lab hierarchy.

## Usage

Launch Jupyter from this directory after activating the `histogwas.yml`
environment:

```bash
jupyter lab notebooks/visualization_thyroid.ipynb
```

For colocalisation, first execute the notebook to export log Bayes factors, then
run the packaged R script:

```bash
Rscript notebooks/coloc_script.R input_logBF.csv output_posteriors.csv
```

Outputs are written to the directories defined inside each notebook or script;
update those paths if you prefer alternative destinations. 
