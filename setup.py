from setuptools import setup, find_packages

setup(
    name="histogwas",
    version="0.1",
    author="Casale Lab",
    author_email="francescopaolo.casale@helmholtz-munich.de",
    description="HistoGWAS tool box",
    packages=find_packages(where="histogwas"),
    package_dir={"": "histogwas"},
    include_package_data=True,
    zip_safe=False,
    install_requires=[
        "numpy",
        "limix_core",
        "tqdm",
        "scipy",
        "statsmodels",
        "pandas",
        "anndata",
        "matplotlib",
        "pandas_plink",
        "scikit-learn",
        "chiscore",
    ],
    extras_require={
        "chi2": ["chi2comb"],
    },
)
