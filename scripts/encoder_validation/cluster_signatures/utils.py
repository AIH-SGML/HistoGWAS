import anndata as ad
import numpy as np
import pandas as pd
from mtgwas import GWAS
from mtgwas.utils import df_match
from limix_core.util.preprocess import gaussianize
from tqdm import tqdm
import gseapy as gp


h5ad_path = f"/lustre/groups/casale/datasets/gtex/histology/20230425_v2_tiles/stage2/%s/embedding/low_memory_scanpy.h5ad"
expr_path = f"/lustre/groups/casale/code/users/shubham.chaudhary/output/projects/gtex/pysrc_v2/emb-gwas/filtered_tpm_v1/gene_tpm_2017-06-05_v8_%s.filtered.processed.h5ad"


# This is the list of clusters for each tissue which has been used in the integration of GWAS
tissue_list = {
    'Breast_Mammary_Tissue': [0, 1, 3, 4, 8],
    'Adipose_Subcutaneous': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
    'Colon_Transverse': [0, 1, 2, 3, 5, 6, 7, 8],
    'Artery_Tibial': [1, 2, 3, 4, 5],
    'Stomach': [0, 2, 3, 4, 5, 6, 7, 8, 10],
    'Esophagus_Mucosa': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
    'Esophagus_Muscularis': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
    'Muscle_Skeletal': [1, 2, 3, 4, 5, 6, 7],
    'Skin_Sun_Exposed_Lower_leg': [0, 1, 2, 3, 4, 5, 6, 7],
    'Thyroid': [0, 1, 3, 2, 4, 5, 6],
    'Nerve_Tibial': [0, 1, 2, 3, 4, 5, 6, 8]
}

def load_anndata(tissue):
    
    print(f'.. loading {tissue} data')
    idata = ad.read_h5ad(h5ad_path % tissue)
    
    print('.. filtering based on the clusters which are considered')
    Ikeep = idata.obs['leiden_0.5'].isin(np.array(tissue_list[tissue]).astype(str))
    idata = idata[Ikeep].copy()
    
    print('.. number of slides with at least 20 tiles in each cluster')
    dfX = pd.get_dummies(idata.obs['leiden_0.5'])
    dfX['slide']  = idata.obs['slide']
    print((dfX.groupby('slide').sum()>20).sum(0))
    
    return idata

def get_expression(tissue):
    edata = ad.read_h5ad(expr_path % tissue.lower())
    edata.obs.index = edata.obs.index.str.split('-SM').str.get(0)
    return edata

def get_fractions(idata):
    dfX = pd.get_dummies(idata.obs['leiden_0.5'])
    dfX['slide'] = idata.obs['slide']
    dfX = dfX.groupby('slide').mean()
    return dfX

def get_pseudobulk(idata):
    dfW = pd.DataFrame(idata.X)
    dfW['slide'] = idata.obs['slide'].values
    dfW = dfW.groupby('slide').mean()
    return dfW

def expr_vs_clustfract(idata, edata):

    # load X
    dfX = get_fractions(idata)

    # load Y
    dfY = edata.to_df()

    # match
    idx, idy = df_match([dfX, dfY])
    dfX = dfX.iloc[idx]
    dfY = dfY.iloc[idy]
    assert (dfX.index==dfY.index).all(), 'Outch!'

    # association test
    Y = gaussianize(dfY.values)
    X = dfX.values
    gwas = GWAS(Y, F=np.ones([Y.shape[0], 1]))
    gwas.process(X)

    # get results
    BETA = gwas.getBetaSNP()
    SE = gwas.getBetaSNPste()
    ZSTAT = BETA / SE
    P = gwas.getPv()

    # compile all together
    dfout = []
    for i in range(BETA.shape[0]):
        _dfout = {}
        _dfout['cluster'] = dfX.columns[i]
        _dfout['gene_id'] = edata.var.index.values
        _dfout['gene_name'] = edata.var['Description'].values
        _dfout['beta'] = BETA[i]
        _dfout['se'] = SE[i]
        _dfout['zstat'] = ZSTAT[i]
        _dfout['pv'] = P[i]
        dfout.append(pd.DataFrame(_dfout).sort_values(by='zstat', ascending=False))
    
    return dfout


def get_enrichment(dfgene, n_top=100, top_vs_bottom=False):

    dfenr = []
    n_clusters = len(dfgene)
    for cluster_i in tqdm(range(n_clusters)):

        # define sets
        gene_list = dfgene[cluster_i].sort_values(by='zstat', ascending=False).iloc[:n_top]['gene_name'].values.astype(str).tolist()
        bottom = dfgene[cluster_i].sort_values(by='zstat', ascending=False).iloc[n_top:]['gene_name'].values.astype(str).tolist()
        if top_vs_bottom:
            background = gene_list + bottom
        else:
            background = dfgene[cluster_i]['gene_name'].values.astype(str).tolist()

        # backgound only reconigized a gene list input.
        enr_bg = gp.enrichr(gene_list=gene_list, gene_sets='MSigDB_Hallmark_2020', background=background, outdir=None)
        _dfenr = enr_bg.results
        _dfenr['cluster'] = dfgene[cluster_i]['cluster'].unique()[0]
        dfenr.append(_dfenr)
        
    return dfenr

    
    