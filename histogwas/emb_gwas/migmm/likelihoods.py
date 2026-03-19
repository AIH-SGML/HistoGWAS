import torch
from torch import nn
from torch.nn import functional as F
import numpy as np
from scipy import stats
from scipy.optimize import minimize
from torch.autograd import Variable
from torch.distributions import Binomial, Normal, NegativeBinomial, Poisson
import time


class BinomialLik(nn.Module):
    
    def __init__(self, total_count=1):
        super().__init__()
        self.total_count = total_count
        
    def log_prob(self, logits, targets):
        likelihood = Binomial(total_count=self.total_count, logits=logits)
        return likelihood.log_prob(targets).sum(0).mean()

    def predict(self, logits):
        return self.total_count * torch.sigmoid(logits)

class PoissonLik(nn.Module):
    def __init__(self):
        super().__init__()

    def log_prob(self, logits, targets):
        rate = torch.exp(logits)
        poisson = Poisson(rate=rate)
        return poisson.log_prob(targets).sum(0).mean()

    def predict(self, logits):
        return torch.exp(logits)


class QuasiNegBinomialLik(nn.Module):
    def __init__(self):
        super().__init__()
        self._logalpha = nn.Parameter(torch.zeros(1))

    @property
    def alpha(self):
        return torch.exp(self._logalpha)

    @alpha.setter
    def alpha(self, value):
        self._logalpha = torch.log(value)

    def log_prob(self, logits, targets):
        mean = torch.exp(logits)
        var = mean + self.alpha
        std = torch.sqrt(var)
        normal = torch.distributions.Normal(mean, std)
        return normal.log_prob(targets).mean()

    def predict(self, logits):
        return torch.exp(logits)


class NormalLik(nn.Module):
    
    def __init__(self):
        super().__init__()
        self._logvar = nn.Parameter(torch.zeros(1))
        
    @property
    def std(self):
        return torch.exp(0.5 * self._logvar)
    
    @std.setter
    def std(self, value):
        self._logvar = 2 * torch.log(value)
    
    def log_prob(self, means, targets):
        likelihood = Normal(means, self.std)
        return likelihood.log_prob(targets).sum(0).mean()

    def predict(self, logits):
        return logits

