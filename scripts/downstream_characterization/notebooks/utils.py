import numpy as np
import pandas as pd
import anndata as ad


def merge_dfs(df1, df2, key1=None, key2=None, left_index=False, right_index=False):


    if key2 is None:
        key2 = key1

    # df1 indexing
    if left_index:
        _df1 = pd.DataFrame({'x': df1.index.values})
    else:
        _df1 = pd.DataFrame({'x': df1[key1]})
    _df1['idx1'] = np.arange(len(_df1))

    # df2 indexing
    if right_index:
        _df2 = pd.DataFrame({'x': df2.index.values})
    else:
        _df2 = pd.DataFrame({'x': df2[key2]})
    _df2['idx2'] = np.arange(len(_df2))

    # match
    _dfm = _df1.merge(_df2, how='inner', on='x')
    idx1 = _dfm['idx1'].values
    idx2 = _dfm['idx2'].values

    return idx1, idx2

def get_snp(snp_id, gdata):
    _ = np.where(gdata.var.index==snp_id)
    snp_to_compute = gdata.X[:, _[0]].compute()
    snp_df = pd.DataFrame(snp_to_compute, columns=['snp'], index=gdata.obs.index.values)
    return snp_df


def get_expression(tissue):
    efile = f'/lustre/groups/casale/code/users/shubham.chaudhary/output/projects/gtex/pysrc_v2/emb-gwas/filtered_tpm_v1/gene_tpm_2017-06-05_v8_{tissue.lower()}.filtered.processed.h5ad'
    
    edata = ad.read_h5ad(efile)
    edata.obs.index = edata.obs.index.str.split('-SM').str.get(0)
    _ = edata.obs.index.str.split('-')
    SID = _.str.get(0) +'-'+_.str.get(1)
    edata.obs['SID'] = SID
    
    df_gene = pd.DataFrame(edata.X, index=edata.obs['SID'], columns=edata.var['Description'].values)
    return df_gene


def get_expression_normalized(tissue):
    expr_path = f"/lustre/groups/casale/datasets/gtex/phenotypes/GTEx_Analysis_v8_eQTL_expression_matrices/{tissue}.v8.normalized_expression.bed.gz"
    seed = 0
    dfe = pd.read_csv(expr_path, sep="\t", index_col=3)
    dfe = dfe[dfe['#chr'].str.split('chr').str[-1].isin(np.arange(23).astype(str))]
    return dfe


def get_egenes(tissue):
    # Reading Egenes
    path = f'/lustre/groups/casale/datasets/gtex/eqtls/GTEx_Analysis_v8_eQTL/{tissue}.v8.egenes.txt.gz'
    egenes = pd.read_csv(path, sep='\t')
    return egenes