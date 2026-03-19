import pandas as pd
from .utils import df_match
import numpy as np
import scipy.stats as st
import pdb


def match_idata_gdata_dfcov(idata, gdata, dfcov=None, idata_bag_key='bag_id'):
    """
    Inputs:
    - idata: instance level anndata
    - gdata: genoanndata
    - dfcov: dataframe with additional covariates (optional)
    - idata_bag_key: key in idata.obs denoting bag id

    Returns:
    - idata: filtered instance level anndata on common bags
    - gdata: filtered genoanndata
    - dfcov: filtered dfcov, if a dfcov is specified as input
    """

    # bag-level dataframe
    dfX = pd.DataFrame({idata_bag_key: idata.obs[idata_bag_key].unique()}).set_index(idata_bag_key)
    
    # match
    if dfcov is not None:
        idx_g, idx_c, idx_h = df_match([gdata.obs, dfcov, dfX])
        gdata = gdata[idx_g]
        dfcov = dfcov.iloc[idx_c]
        dfX = dfX.iloc[idx_h]
        assert (gdata.obs.index.values.astype(str)==dfX.index.values).all()
        assert (gdata.obs.index.values.astype(str)==dfcov.index.values).all()
        
    else:
        # match
        idx_g, idx_h = df_match([gdata.obs, dfX])
        gdata = gdata[idx_g]
        dfX = dfX.iloc[idx_h]
        assert (gdata.obs.index.values.astype(str)==dfX.index.values).all()

    # match idata
    bag_index = f'{idata_bag_key}_index'
    map_sid = dict(zip(dfX.index.values, np.arange(dfX.shape[0])))
    idata = idata[idata.obs[idata_bag_key].isin(dfX.index.values)].copy()
    idata.obs[bag_index] = idata.obs[idata_bag_key].map(map_sid).astype(int)

    # reorders instance anndata
    idata.obs['i'] = np.arange(idata.shape[0])
    idxs = idata.obs[['i', bag_index]].sort_values(by=bag_index)['i'].values
    idata = idata[idxs].copy()
    del idata.obs['i']

    if dfcov is None:
        return idata, gdata
    else:
        return idata, gdata, dfcov
