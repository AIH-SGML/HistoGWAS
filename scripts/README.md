# HistoGWAS Scripts

This folder contains the main pipeline scripts for the HistoGWAS project.

## Modules

- **preprocessing**: Downloading slides, tiling whole-slide images, and generating embeddings stored in AnnData format.
- **encoder_validation**: Validating and selecting the best encoder for downstream analysis.
- **gwas**: Running GWAS-style association testing from slide-level embeddings.
- **downstream_characterization**: Post-association analyses and biological/clinical characterization of findings.

## Typical workflow

1. Run `preprocessing` to prepare slide embeddings.
2. Use `encoder_validation` to select the encoder.
3. Run `gwas` for association testing.
4. Use `downstream_characterization` for follow-up interpretation.
