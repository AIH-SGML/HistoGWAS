import pandas as pd
import numpy as np


def df_match(dfs, keys=None):
    if keys is None:
        keys = [None] * len(dfs)
    dfidxs = []
    for _i, (df, key) in enumerate(zip(dfs, keys)):
        _id = df[key].values if key is not None else df.index.values
        _dfidx = {'ID': _id, f'id{_i}': np.arange(df.shape[0])}
        dfidxs.append(pd.DataFrame(_dfidx))
    dfidx = dfidxs[0].copy()
    for _i in range(1, len(dfidxs)):
        dfidx = dfidx.merge(dfidxs[_i], how='inner', on='ID')
    out = [dfidx[f'id{_i}'].values for _i in range(len(dfs))]
    return out
