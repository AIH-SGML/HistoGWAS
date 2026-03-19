# HistoGWAS Downstream Characterization Notebooks

This directory groups interactive analyses for exploring HistoGWAS release outputs once the preprocessing, GWAS, and downstream pipelines have completed. Each notebook targets a specific slice of the post-GWAS workflow and can be run independently as long as the required summary statistics and embeddings are available.

## Contents
- `visualization_thyroid.ipynb` – original exploratory walkthrough for the thyroid cohort; visualizes genetic effects recovered by the semantic decoder, inspects tile embeddings, highlights lead associations.
- `analysisng_gene_prediction.ipynb` – evaluates gene expression prediction models built from histology-derived embeddings, including train/test splits, performance diagnostics, and QQ-style calibration plots across tissues.
- `pathway_enrichment_analysis.ipynb` – compiles pathway and gene-set enrichment results (Table 1) from significant features, producing tables suitable for manuscript export.
- `colcalization_with_eqtls.ipynb` – prepares inputs for eQTL colocalization, harmonises association statistics, and generates summary plots to compare GWAS and GTEx signals.
- `coloc_script.R` – lightweight wrapper around the `coloc` R package; consumes log Bayes factors written by the notebooks and exports posterior probabilities to CSV.
- `utils.py` – shared Python helpers for merging phenotype/variant tables and loading GTEx expression matrices used by the notebooks above.

## Environment
- Python ≥3.8 with `scanpy`, `anndata`, `pandas`, `numpy`, `torch`, `scikit-learn`, `seaborn`, `matplotlib`, `tqdm`, and the internal `histogwas` package.
- `limix-core` and related preprocessing utilities for colocalization notebooks.
- R ≥4.0 with the `coloc` package to execute `coloc_script.R`.
- The `histogwas.yml` environment file in this repository provisions the required Python dependencies; add the R `coloc` package separately if needed.

## Running the notebooks
1. Activate the HistoGWAS conda/mamba environment described by `histoGWAS_2.yml` and ensure GTEx datasets referenced in `utils.py` are accessible from your machine.
2. Launch Jupyter from this directory (or the repository root) and open the notebook of interest, e.g. `jupyter lab pathway_enrichment_analysis.ipynb`.
3. Update any configuration cells that point to Casale lab filesystem locations so they match your local workspace. Paths to embeddings, association summaries, and colocalization inputs are typically defined near the top of each notebook.
4. For colocalization workflows, run the Python notebook to assemble log Bayes factors, then execute `coloc_script.R input.csv output.csv` inside the same environment (with R) to obtain posterior probabilities.

## Data expectations and tips
- Generate GWAS summary statistics and tile-level embeddings before launching these notebooks; they assume TSV/NPY/H5AD outputs already exist.
- The helper routines in `utils.py` expect GTEx v8 expression resources in the shared `/lustre/groups/casale` hierarchy; adjust the hard-coded paths if you mirror the data elsewhere.
- Keep track of the git revision and environment details when exporting figures or tables for publication.
- When adding new analyses, drop the notebook in this directory, note its expected inputs at the top of the file, and add a short description above so collaborators can quickly navigate the available tooling.
