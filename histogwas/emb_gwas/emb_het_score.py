import numpy as np
import pandas as pd
import scipy.linalg as la
from chiscore._davies import _pvalue_lambda
from limix_core.util.preprocess import regressOut
from tqdm import tqdm
import pdb

from .numpy_scatter import scatter_dot, scatter_mean, scatter_sum

class EmbGWAS:
    def __init__(self, idata, dfcovs, compute=True, idata_bag_index="bag_index", het_test_data_norm="unit_norm"):
        # assert data match
        assert (
            np.sort(np.unique(idata.obs[idata_bag_index])) == np.arange(dfcovs.shape[0])
        ).all(), "bag_index not valid"

        # set self
        self.idata = idata
        self.dfcovs = dfcovs
        self.compute = compute

        # define het gwas
        F = dfcovs.values
        X = (idata.X - idata.X.mean(0)) / idata.X.std(0)
        i = idata.obs[idata_bag_index].values
        # pdb.set_trace()

        # define Z
        self.compute_pseudobulk(X, i, norm=het_test_data_norm)

        
        # define nohet gwas
        self.gwas_test = GaussScoreTest(F, self.Z)

    # compute pseudobulk
    def compute_pseudobulk(self, X, i, norm="unit_norm"):
        # pdb.set_trace()
        Z = scatter_mean(X, i)
        if norm == "standard":
            Z = (Z - Z.mean(0)) / Z.std(0)
            Z = Z / np.sqrt(Z.shape[1])
        elif norm == "unit_norm":
            Z = Z / np.sqrt((Z**2).sum(1))[:, None]
        self.Z = Z

    def gwas(self, gdata, block_size=10000):
        # assert ids match
        assert (self.dfcovs.index == gdata.obs.index).all(), "Ids do not match"

        # perform GWAS
        i0s = np.arange(0, gdata.shape[1], block_size)
        i1s = i0s + block_size
        dfres = []
        for i0, i1 in tqdm(zip(i0s, i1s), total=len(i0s)):
            # subset data and standardize
            if not self.compute:
                _gdata = gdata[:, i0:i1]
            else:
                _gdata = gdata[:, i0:i1].compute()
            _gdata.impute()
            _gdata.calc_maf()
            _gdata.standardize()

            # run test
            _res = self.gwas_test.association(_gdata.X, verbose=False)

            # append all results to dfres
            _dfres = _gdata.var
            for key in _res.keys():
                _dfres[key] = _res[key].values
            dfres.append(_dfres)

        dfres = pd.concat(dfres, axis=0)
        return dfres


class GaussScoreTest:
    def __init__(self, F, X):
        # store
        self.W = X
        self.F = F

        # precompute eigenvals
        self.P0W = regressOut(self.W, self.F)
        Lambda0 = 0.5 * self.W.T.dot(self.P0W)
        self.lambd0 = la.eigvalsh(Lambda0)[::-1]

    def association(self, G, verbose=True):
        # Fit null
        P0y = regressOut(G, self.F)
        null_var = P0y.var(0)
        Py = P0y / null_var

        # Compute Q
        WPy = self.W.T.dot(Py)
        Q = 0.5 * np.einsum("rs,rs->s", WPy, WPy)

        # loop to compute pvals
        iterator = range(G.shape[1])
        iterator = tqdm(iterator) if verbose else iterator
        RV = []
        for s in iterator:
            lambd = self.lambd0 / null_var[s]
            _ = _pvalue_lambda(lambd, Q[[s]])
            rv = {key: _[key].ravel()[0] for key in ["p_value", "p_val_liu", "is_converge"]}
            RV.append(rv)
        return pd.DataFrame(RV)
