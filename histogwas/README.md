# histogwas package

Core Python package distributed with the HistoGWAS release.

## Install

From the repository root (`HistoGWAS_release`):

```bash
pip install --upgrade pip
pip install .
```

For editable development installs:

```bash
pip install --upgrade pip
pip install -e .
```

If you are already in the `histogwas/` directory, point `pip` one level up
(`pip install ..`).

## Package layout

```
histogwas/
├── __init__.py           # exposes emb_gwas and vctest convenience imports
├── emb_gwas/             # main GLMM-based association toolkit
│   ├── emb_het_score.py  # EmbGWAS class for inference
│   ├── geno_utils.py     # PLINK readers and genotype helpers
│   ├── misc.py, utils.py
│   └── migmm/            # GLMM implementation and likelihoods
├── vctest.py             # auxiliary validation / chi-square tests
└── setup.py              # legacy setup shim (prefer root pyproject.toml)
```

## Usage

```python
from histogwas import emb_gwas, vctest

model = emb_gwas.EmbGWAS(...)
result = vctest.run_test(...)
```

Direct submodule imports continue to work, e.g.
`from histogwas.emb_gwas.migmm import GLMM`.

## Notes

- The package is shipped alongside the Stage 2–4 scripts.
- Dependencies such as `chiscore`, `chi2comb`, `limix_core`, and `pandas_plink`
  are declared in `setup.py`; ensure they are present in your environment before
  running the association drivers.
