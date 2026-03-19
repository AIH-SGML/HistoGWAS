# Quickstart: Running the HistoGWAS Pipeline

This guide assumes you have already installed **HistoGWAS** and set up the appropriate environment.  
We start directly from preprocessing and continue through GWAS and downstream analyses.

---

## 1. Preprocessing Slides

Run preprocessing on the sample data (Stage 1–3):

```bash
cd scripts/preprocessing
python stage1_get_data_run.py --config ../../data/samplelist.csv
python stage2_tiling_run.py --input stage1_output --output stage2_output
python stage3_tile_features_run.py --input stage2_output --output stage3_output
```

Each stage produces intermediate outputs (`stage1_output/`, `stage2_output/`, etc.).

---

## 2. Encoder Selection

Filter genes and run gene-prediction evaluation:

```bash
cd ../encoder_validation
python 1_filter_gene_run.py --input ../../data/samplelist.csv
python 2_gene_prediction_run.py --input filtered_genes.csv
```

---

## 3. Genome-Wide Association

Run GWAS on the embeddings:

```bash
cd ../gwas
python hgwas_all_clusterwise_run.py --input ../../scripts/encoder_validation/results
```

---

## 4. Downstream Characterization

Example analyses are available in notebooks under [`scripts/downstream_characterization/notebooks`](../scripts/downstream_characterization/notebooks):

- `colocalization_with_eqtls.ipynb`
- `pathway_enrichment_analysis.ipynb`
- `visualization_thyroid.ipynb`

Open one in Jupyter:

```bash
jupyter notebook notebooks/colocalization_with_eqtls.ipynb
```

---

## 5. Next Steps

- Explore the [`docs/figures`](../docs/figures) for visual summaries of the pipeline.
- See per-stage `README.md` files under `scripts/` for details and CLI arguments.
- For running at scale, use the provided Slurm wrappers in each stage folder.

---

🚀 You’re now ready to run **HistoGWAS** end-to-end!
