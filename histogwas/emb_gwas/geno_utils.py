import os
import pdb
import warnings
from os.path import basename, dirname, join

import anndata as ad
import numpy as np
import pandas as pd
import pandas_plink as pp
from anndata._core.index import _normalize_indices
from sklearn.impute import SimpleImputer


class GenoAnnData:
    def __init__(self, X, obs, var, computed):
        self.X = X
        self.obs = obs
        self.var = var
        self.computed = computed
        self.standardized = False

    @property
    def obs_names(self):
        return self.obs.index

    @property
    def var_names(self):
        return self.var.index

    @property
    def n_obs(self):
        return self.obs.shape[0]

    @property
    def n_vars(self):
        return self.var.shape[0]

    @property
    def shape(self):
        return self.X.shape

    def compute(self):
        if not self.computed:
            return GenoAnnData(self.X.compute(), self.obs, self.var, True)
        else:
            warnings.warn("Compute was already run")
            return self

    def impute(self):
        assert self.computed, "Run compute first!"
        imputer = SimpleImputer()
        self.X = imputer.fit_transform(self.X)

    def standardize(self):
        assert self.computed, "Run compute first!"
        assert not np.isnan(self.X).any(), "X contains NaN, run impute first!"
        self._mean = self.X.mean(0)
        self._std = self.X.std(0)
        self.X = (self.X - self._mean) / self._std
        self.standardized = True

    @property
    def allele_counts(self):
        if self.standardized:
            return self.X * self._std + self._mean
        return self.X

    def calc_maf(self):
        assert self.computed, "Run compute first!"
        assert not np.isnan(self.X).any(), "X contains NaN, run impute first!"
        self.var["maf"] = 0.5 * self.X.mean(0)

    def _normalize_indices(self, index):
        return _normalize_indices(index, self.obs_names, self.var_names)

    def __len__(self):
        return self.shape[0]

    def __getitem__(self, index):
        """Returns a sliced view of the object."""
        idx1, idx2 = self._normalize_indices(index)
        _X = self.X[idx1, idx2]
        _obs = self.obs.iloc[idx1].copy()
        _var = self.var.iloc[idx2].copy()
        return GenoAnnData(_X, _obs, _var, computed=self.computed)

    def __repr__(self):
        descr = f"GenoAnnData object with n_obs × n_vars = {self.n_obs} × {self.n_vars}"
        for attr in ["obs", "var"]:
            keys = getattr(self, attr).keys()
            if len(keys) > 0:
                descr += f"\n    {attr}: {str(list(keys))[1:-1]}"
        return descr

    def copy(self):
        return GenoAnnData(
            self.X.copy(), self.obs.copy(), self.var.copy(), computed=self.computed
        )


def read_plink(bfile, pcfile=None, num_pcs=None):
    # load geno
    bim, fam, bed = pp.read_plink(bfile)
    del bim["cm"]
    del bim["i"]
    bim = bim.set_index("snp")
    fam = fam.set_index("iid")[[]]

    if pcfile is not None:
        # read pcs
        usecols = None if num_pcs is None else np.arange(num_pcs + 2)
        dfpc = pd.read_csv(pcfile, sep=" ", header=None, usecols=usecols)
        dfpc.columns = ["fid", "iid"] + [f"PC{i+1}" for i in range(dfpc.shape[1] - 2)]
        del dfpc["fid"]
        dfpc = dfpc.set_index("iid")
        assert (dfpc.index == fam.index).all(), "Individuals in PC files not matching"
        fam = pd.concat([fam, dfpc], axis=1)

    return GenoAnnData(bed.T, fam, bim, computed=False)