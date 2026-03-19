import numpy as np
import scipy.special as ss


def scatter_sum(x, i, keepdim=False):
    out_shape = list(x.shape)
    out_shape[0] = i.max() + 1
    out = np.zeros(out_shape)
    np.add.at(out, i, x)
    if keepdim:
        out = out[i]
    return out

def scatter_dot(x1, x2, i, keepdim=False):
    out_shape = list((x1[[0]]*x2[[0]]).shape)
    out_shape[0] = i.max() + 1
    out = np.zeros(out_shape)
    for _ in range(x1.shape[0]):
        out[i[_]] += x1[_] * x2[_]
    if keepdim:
        out = out[i]
    return out

def scatter_mean(x, i, keepdim=False):
    wt = 1 / scatter_sum(np.ones([x.shape[0], 1]), i, keepdim=True)
    out = scatter_dot(x, wt, i)
    if keepdim:
        out = out[i]
    return out

def scatter_softmax(x, i, keepdim=False):
    ex = np.exp(x)
    out = ex / scatter_sum(ex, i, keepdim=True)
    return out
