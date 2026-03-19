"""Convenience exports for the HistoGWAS release package."""

from importlib import metadata as _metadata

from . import emb_gwas
from . import vctest

try:
    __version__ = _metadata.version("histogwas")
except _metadata.PackageNotFoundError:  # pragma: no cover
    __version__ = "0.0.0"

__all__ = ["__version__", "emb_gwas", "vctest"]
