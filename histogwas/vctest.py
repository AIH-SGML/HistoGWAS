import numpy as np
import pandas as pd
import scipy.linalg as la
from statsmodels.stats.multitest import multipletests
import copy
from tqdm import tqdm


def vctest_solve(UtF, UtY, Sr, delta):

    # compute D
    S = delta * Sr + 1 - delta
    D = 1. / S

    # compute beta
    UtFtD = UtF.T / S[np.newaxis, :]
    Areml_inv = la.pinv(UtFtD.dot(UtF))
    beta = Areml_inv.dot(UtFtD.dot(UtY))

    # compute sigma2, lml, vg and vn
    n = UtF.shape[0]
    quad = D[:, np.newaxis] * (UtY - UtF.dot(beta))**2
    sigma2 = quad.mean(0)
    lml = - 0.5 * (n * np.log(2 * np.pi) + np.log(S).sum() + n + n * np.log(sigma2))
    vg = delta * sigma2
    vn = (1 - delta) * sigma2

    res = pd.DataFrame(
        {
            'idx': np.arange(sigma2.shape[0]),
            'beta': beta.ravel(),
            'vg': vg,
            'vn': vn,
            'delta': delta * np.ones(sigma2.shape[0]),
            'sigma2': sigma2,
            'lml': lml
        }
    )
    return res


def vctestlr_solve(X, Y, F, delta):

    # rank-by-rank  matrix from woodbury identity
    W = np.eye(X.shape[1]) + delta / (1 - delta) * np.dot(X.T, X)
    W_chol = la.cholesky(W).T
    W_inv = la.cho_solve((W_chol,True), np.eye(W_chol.shape[0]))
    W_inv = la.inv(W)
    def Vi_dot(x):
        return 1 / (1 - delta) * x - delta / (1 - delta)**2 * X.dot(W_inv.dot(X.T.dot(x)))

    # this is just for debugging purposes against efficient Vi_dot
    def Vi_dot_naive(x):
        V = delta * np.dot(X, X.T) + (1 - delta) * np.eye(X.shape[0])
        return la.inv(V).dot(x)

    # compute beta
    Areml_inv = la.pinv(F.T.dot(Vi_dot(F)))
    beta = Areml_inv.dot(F.T.dot(Vi_dot(Y)))

    # compute sigma2, lml, vg and vn
    n = X.shape[0]
    Ytilde = (Y - F.dot(beta))
    sigma2 = np.einsum('ij,ij->j', Ytilde, Vi_dot(Ytilde)) / n
    logdetV = 2 * np.log(W_chol.diagonal()).sum() + n * np.log(1 - delta)
    lml = - 0.5 * (n * np.log(2 * np.pi) + logdetV + n + n * np.log(sigma2))
    vg = delta * sigma2
    vn = (1 - delta) * sigma2
    
    res = pd.DataFrame(
        {
            'idx': np.arange(sigma2.shape[0]),
            'beta': beta.ravel(),
            'vg': vg,
            'vn': vn,
            'delta': delta * np.ones(sigma2.shape[0]),
            'sigma2': sigma2,
            'lml': lml
        }
    )
    return res


def xgower_factor_(X):
    a = np.power(X, 2).sum()
    b = X.dot(X.sum(0)).sum()
    return np.sqrt((a - b / X.shape[0]) / (X.shape[0] - 1))


class VCTEST():
    # X = #inds x 2048
    # Y = #inds x #genes
    # X -> Y
    # y_g = Y[:,g] #inds x 1
    # y_g ~ N(0, \sigma_x^2 XXt + \sigma_n^2 I_N)
    def __init__(self, ndeltas=100):
        self.ndeltas = ndeltas
        deltas = np.linspace(0, np.sqrt(1 - 1e-4), self.ndeltas)**2
        self.deltas = deltas
        self.fitted = False

    def fit(self, X, Y, normalize_X=True, compute_pvals=False, nperms=10, verbose=True):
        
        self.verbose = verbose
        self.normalize_X = normalize_X

        if self.normalize_X:
            self._xm = X.mean(0)
            self._gf = xgower_factor_(X)
            X = (X - self._xm) / self._gf

        # compute R and its eigen decomp
        R = np.dot(X, X.T)
        S_R, U_R = la.eigh(R)
        S_R = np.clip(S_R, 0, np.inf)

        # set X and F to self
        self.X = X
        self.F = np.ones([Y.shape[0], 1])
        self.U_R = U_R
        self.S_R = S_R
        self.UtF = np.dot(U_R.T, self.F)

        # set Y to self
        self._set_Y(Y)

        # to compute
        self.beta = np.zeros(Y.shape[1])
        self.vg = np.zeros(Y.shape[1])
        self.vn = np.zeros(Y.shape[1])
        self.llr = np.zeros(Y.shape[1])
        self.Ystar = np.zeros_like(Y)

        # core computations
        self.res = self._compute()

        if compute_pvals:
            self._compute_pvals(nperms=nperms)

        return self.res

    def _compute(self):

        res = []
        iter_deltas = tqdm(self.deltas) if self.verbose else self.deltas
        for _delta in iter_deltas:
            _res = vctest_solve(self.UtF, self.UtY, self.S_R, _delta)
            res.append(_res)
        res = pd.concat(res, axis=0).reset_index(drop=1)
        res0 = res.loc[res['delta']==0, ['idx', 'lml']]
        res0['lml0'] = res0.pop('lml')
        res = res.merge(res0, on='idx')
        res['llr'] = res['lml'] - res['lml0']
        res = res.iloc[res.groupby('idx')['llr'].idxmax()]
        self.fitted = True

        return res

    def _set_Y(self, Y):
        self.Y = Y
        self.UtY = np.dot(self.U_R.T, self.Y)

    def _compute_pvals(self, nperms=10, verbose=False):

        vctest0 = copy.copy(self)

        np.random.seed(0)

        llr0 = []
        iter_nperms = tqdm(range(nperms)) if verbose else range(nperms)
        for perm_i in iter_nperms:
            Y0 = np.random.randn(*self.Y.shape)
            vctest0._set_Y(Y0)
            _res0 = vctest0._compute()
            llr0.append(_res0['llr'].values)
        llr0 = np.concatenate(llr0)
        count = np.array([(llr0 >= _x).sum() for _x in self.res['llr'].values])
        self.res['pvals'] = (count + 1.) / (llr0.shape[0] + 1.)
        self.res['qvals'] = multipletests(self.res['pvals'].values, method='fdr_bh')[1]

    def predict_loo(self):
        beta = self.res['beta'].values
        vg = self.res['vg'].values
        vn = self.res['vn'].values
        Yr = self.Y - self.F.dot(beta[None, :])
        D = 1 / (vg[None, :] * self.S_R[:, None] + vn[None, :])
        HY = vg[None, :] * np.dot(self.U_R, self.S_R[:, None] * D * np.dot(self.U_R.T, Yr))
        Hdiag = vg[None, :] * np.dot(self.U_R**2, self.S_R[:, None] * D)
        return (HY - Hdiag * Yr) / (1 - Hdiag + 1e-9)

    def predict(self, Xstar):
        if self.normalize_X:
            Xstar = (Xstar - self._xm) / self._gf
        return np.dot(Xstar, self.get_beta())

    def get_beta(self):
        beta = self.res['beta'].values
        vg = self.res['vg'].values
        vn = self.res['vn'].values
        Yr = self.Y - self.F.dot(beta[None, :])
        D = 1 / (vg[None, :] * self.S_R[:, None] + vn[None, :])
        KiYr = np.dot(self.U_R, D * np.dot(self.U_R.T, Yr))
        return vg[None, :] * np.dot(self.X.T, KiYr)


class VCTESTLR(VCTEST):

    def __init__(self, ndeltas=100):
        VCTEST.__init__(self, ndeltas=ndeltas)

    def fit(self, X, Y, normalize_X=True, compute_pvals=False, nperms=10, verbose=True):

        self.verbose = verbose
        self.normalize_X = normalize_X

        if self.normalize_X:
            self._xm = X.mean(0)
            self._gf = xgower_factor_(X)
            X = (X - self._xm) / self._gf

        self.Y = Y
        self.X = X
        self.F = np.ones([Y.shape[0], 1])
        self.beta = np.zeros(Y.shape[1])
        self.vg = np.zeros(Y.shape[1])
        self.vn = np.zeros(Y.shape[1])
        self.llr = np.zeros(Y.shape[1])
        self.Ystar = np.zeros_like(Y)

        res = []
        iter_deltas = tqdm(self.deltas) if self.verbose else self.deltas 
        for _delta in iter_deltas:
            _res = vctestlr_solve(self.X, self.Y, self.F, _delta)
            res.append(_res)
        res = pd.concat(res, 0).reset_index(drop=1)
        
        res0 = res.loc[res['delta']==0, ['idx', 'lml']]
        res0['lml0'] = res0.pop('lml')
        res = res.merge(res0, on='idx')
        res['llr'] = res['lml'] - res['lml0']
        res = res.iloc[res.groupby('idx')['llr'].idxmax()]
        self.res = res
        self.fitted = True

        if compute_pvals:
            self._compute_pvals(nperms=nperms)

        return res

    def _Vi_dot(self, x, delta):
        X = self.X
        W = np.eye(X.shape[1]) + delta / (1 - delta) * np.dot(X.T, X)
        W_chol = la.cholesky(W).T
        W_inv = la.cho_solve((W_chol,True), np.eye(W_chol.shape[0]))
        return 1 / (1 - delta) * x - delta / (1 - delta)**2 * X.dot(W_inv.dot(X.T.dot(x)))
        
    def predict_loo(self):
        beta = self.res['beta'].values
        delta = self.res['delta'].values
        X = self.X
        Yr = self.Y - self.F.dot(beta[None, :])
        Vi_dot_Yr = np.concatenate([self._Vi_dot(Yr[:,[ip]], delta[ip]) for ip in range(Yr.shape[1])], 1)
        HY = delta * np.dot(X, np.einsum('ir,ip->rp', X, Vi_dot_Yr))
        Vi_dot_X = np.concatenate([self._Vi_dot(X, delta[ip])[:,:,None] for ip in range(Yr.shape[1])], 2)
        Hdiag = delta * np.einsum('ir,irp->ip', X, Vi_dot_X)
        return (HY - Hdiag * Yr) / (1 - Hdiag + 1e-9)
    
    def get_beta(self):
        beta = self.res['beta'].values
        delta = self.res['delta'].values
        Yr = self.Y - self.F.dot(beta[None, :])
        Vi_dot_Yr = np.concatenate([self._Vi_dot(Yr[:,[ip]], delta[ip]) for ip in range(Yr.shape[1])], 1)
        return delta * self.X.T.dot(Vi_dot_Yr)
 


