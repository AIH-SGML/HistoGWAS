import time

import numpy as np
import torch
import tqdm
from torch import nn

from emb_gwas.migmm.likelihoods import BinomialLik


class GLMM(nn.Module):
    def __init__(self, y, F, X, n_samples=128, seed=0, repar="xmil_std", lik=BinomialLik()):
        super().__init__()

        # reparametrize
        assert repar in ["xmil_std", "unorm_z", None], f"value '{repar}' not valid for repar"
        self.repar = repar

        # data
        self.y = y
        self._set_F(F)
        self._set_X(X)

        # for consistent typing
        self._dtype = y.dtype

        # init
        self._init_params()
        self._sample_eps(n_samples, seed)

        # set likelihood
        self.lik = lik

    def _set_F(self, F):
        self.F = F
        self._K = F.shape[1]

    def _set_X(self, X):
        if self.repar == "unorm_z":
            X = X / torch.norm(X, dim=1, keepdim=True)
        self.X = X
        self._Q = X.shape[1]

    @property
    def alpha(self):
        """
        Get the parameters of the fixed effects.

        Returns
        -------
        ndarray
            Parameters of the fixed effects.
        """
        return self._alpha

    @alpha.setter
    def alpha(self, value):
        assert value.ndim == 1, "alpha should be a 1D array."
        assert value.shape[0] == self._K, "Dimensions of alpha do not match with the number of fixed effects."
        self._alpha.data = value

    @property
    def beta_m(self):
        """
        Get the mean of the random effects.

        Returns
        -------
        ndarray
            Mean of the random effects.
        """
        return self._beta_m

    @beta_m.setter
    def beta_m(self, value):
        assert value.ndim == 1, "beta_m should be a 1D array."
        assert value.shape[0] == self._Q, "Dimensions of beta_m do not match with the number of random effects."
        self._beta_m.data = value

    @property
    def beta_s(self):
        """
        Get the standard deviation of the random effects.

        Returns
        -------
        ndarray
            Standard deviation of the random effects.
        """
        return torch.exp(self._log_beta_s)

    @beta_s.setter
    def beta_s(self, value):
        assert value.ndim == 1, "beta_s should be a 1D array."
        assert value.shape[0] == self._Q, "Dimensions of beta_s do not match with the number of random effects."
        assert (value > 0).all(), "Dimensions of beta_s do not match with the number of random effects."
        self._log_beta_s.data = torch.log(value)

    def predict(self, X=None, F=None):
        xmil = self._get_xlim(X=X)
        F = self.F if F is None else F
        logits = F.mm(self.alpha[:, None]) + xmil
        ystar = self.lik.predict(logits)
        return ystar.mean(1), ystar.std(1)

    def _init_params(self):
        self._alpha = torch.nn.Parameter(1e-3 * torch.randn(self._K, dtype=self._dtype))
        self._beta_m = torch.nn.Parameter(1e-3 * torch.randn(self._Q, dtype=self._dtype))
        self._log_beta_s = torch.nn.Parameter(1e-3 * torch.randn(self._Q, dtype=self._dtype))

    def _sample_eps(self, n_samples=128, seed=0):
        torch.manual_seed(42)
        self._eps_beta = torch.randn(self._Q, n_samples)

    def _get_xlim(self, X=None):
        # sample
        beta = self.beta_m[:, None] + self.beta_s[:, None] * self._eps_beta
        if self.repar == "xmil_std":
            beta_std = torch.sqrt((self.beta_m**2 + self.beta_s**2).mean())
            beta = beta / beta_std

        if X is None:
            X = self.X
        else:
            X = X / torch.norm(X, dim=1, keepdim=True)

        # xmil
        xmil = X.mm(beta)
        if self.repar == "xmil_std":
            xmil = beta_std * (xmil - xmil.mean(0)) / xmil.std(0)

        return xmil

    def elbo(self):
        # compute optimal beta_var and gamma_var
        beta_var = torch.mean(self.beta_m**2 + self.beta_s**2)

        # log likelihood
        logits = self.F.mm(self.alpha[:, None]) + self._get_xlim()
        loglik = self.lik.log_prob(logits, self.y)

        # kld
        kld_beta = 0.5 * torch.log(beta_var / self.beta_s**2).sum() + 0.5 * self._Q - 0.5

        # elbo
        elbo = loglik - kld_beta
        return elbo

    def optimize(self, max_iter=20, factr=1e7, verbose=False):
        """
        Optimize the model parameters.

        Returns
        -------
        float
            The minimum value of the objective function.
        ndarray
            gradient
        bool
            converged
        """
        ftol = factr * np.finfo(float).eps
        optimizer = torch.optim.LBFGS(self.parameters(), line_search_fn="strong_wolfe")

        def closure():
            optimizer.zero_grad()
            loss = -self.elbo()
            loss.backward()
            return loss

        def grad():
            optimizer.zero_grad()
            loss = -self.elbo()
            loss.backward()
            return torch.cat([p.grad.flatten() for p in self.parameters()])

        t0 = time.time()
        conv = False
        previous_loss = -self.elbo()
        iterator = (
            tqdm.tqdm(range(max_iter), desc=f"Optimize {self.__class__.__name__}") if verbose else range(max_iter)
        )
        for _iter in iterator:
            loss = optimizer.step(closure)
            if _iter > 2 and (loss.item() - previous_loss.item()) < ftol:
                conv = True
                if verbose:
                    print(f"Converged after {_iter} steps")
                break
            previous_loss = loss
        if verbose:
            print("Elapsed:", time.time() - t0)
        return loss, grad(), conv
