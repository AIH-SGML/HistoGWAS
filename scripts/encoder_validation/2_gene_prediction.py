from math import e
from emb_gwas.migmm import GLMM
from emb_gwas.migmm.likelihoods import NormalLik
import torch
from sklearn.metrics import r2_score
import anndata as ad
import pandas as pd
import emb_gwas as eg
import numpy as np
from tqdm import tqdm
from scipy import stats
import os
from os.path import join
import argparse


def get_args():
    parser = argparse.ArgumentParser()

    
    '''split_index basically tells us from where we start taking the partition in dataset for that pertical job
        split_cout basically stores the  number of jobs that is created
    '''
    parser.add_argument("--tissue", dest="tissue", type=str)
    parser.add_argument("--model_type", dest="model_type", type=str)
    parser.add_argument("--hfile", dest="hfile", type=str, default=None)
    parser.add_argument("--efile", dest="efile", type=str, default=None)
    parser.add_argument("--cluster_i", dest="cluster_i", type=str, default=None)
    parser.add_argument("--grade", dest="grade", type=str, default=None)
    parser.add_argument("--outdir", dest="outdir", type=str, default=None)
    parser.add_argument("--pc_embed", dest="pc_embed", action="store_true", default=False)
    parser.add_argument("--train_test_split", dest="train_test_split", type=str, default="../../data/train_test_split.csv")
    
    args = parser.parse_args()
    
    return args


def main(args):
    
    # pdb.set_trace()
    
    # Reading the anndata
    tissue = args.tissue #'Thyroid'
    model_type = args.model_type
    hfile = args.hfile
    efile = args.efile
    

    n_comp = 32 if model_type == 'plip' else 64

    if os.path.exists(hfile):
        idata = ad.read_h5ad(hfile)
        if args.cluster_i is not None:
            if 'leiden_0.5' not in idata.obs.columns:
                raise ValueError("cluster_i was provided but 'leiden_0.5' is missing in embedding metadata.")
            idata = idata[idata.obs['leiden_0.5']==args.cluster_i]
        if 'SID' not in idata.obs.columns:
            if 'slide' not in idata.obs.columns:
                raise ValueError("Embedding metadata must contain either 'SID' or 'slide'.")
            _ = idata.obs['slide'].astype(str).str.split('-')
            idata.obs['SID'] = _.str.get(0) + '-' + _.str.get(1)
        idata = idata[idata.obs['SID'].isin((idata.obs['SID'].value_counts()[idata.obs['SID'].value_counts()>=10]).index)]
    else:
        # loading the embedding of tissue, which is an anndata object
        hfile = os.path.join(os.path.dirname(hfile), f'summary_scanpy.h5ad')
        adata = ad.read_h5ad(hfile)
        
        pc_path = os.path.join(os.path.dirname(hfile), f'{tissue}_img_embedding.h5ad') # this is for saving the pc_componenet

        _ = adata.obs['slide'].str.split('-')
        adata.obs['SID'] = _.str.get(0) + '-'+ _.str.get(1)
        


        # Here I am thresholding the clusters

        adata_new = adata
        if args.cluster_i is not None:
            if 'leiden_0.5' not in adata.obs.columns:
                raise ValueError("cluster_i was provided but 'leiden_0.5' is missing in embedding metadata.")
            adata_new = adata[adata.obs['leiden_0.5']==args.cluster_i]
        adata_new = adata_new[adata_new.obs['SID'].isin((adata_new.obs['SID'].value_counts()[adata_new.obs['SID'].value_counts()>=10]).index)]
        n_comp = min(adata_new.obsm['X_pca'].shape[1], n_comp)
        idata = ad.AnnData(adata_new.obsm['X_pca'][:,:n_comp], obs=adata_new.obs)
        idata.write_h5ad(pc_path)
    
    # Readin the preprocessed TPM
    # add tile id
    slide = idata.obs['slide'].astype(str)
    tnumber = idata.obs['tile_id'].astype(str).str.split('tile_').str.get(1)
    obs = idata.obs.copy()
    obs['id'] = slide + '-' + tnumber
    idata.obs = obs

    # read edata
    xmin = 0.1
   
    
    
    edata = ad.read_h5ad(efile)
            
    edata.obs.index = edata.obs.index.str.split('-SM').str.get(0)

    # match idata, gdata and dfcov
    idata, edata = eg.match_idata_gdata_dfcov(idata, edata, idata_bag_key='slide')
    idata.X = (idata.X - idata.X.mean(0)) / (idata.X.std(0) * np.sqrt(idata.shape[1]))
    
    
    # split data in train and test
    # slides_t, slides_v = train_test_split(edata.obs.index.values, train_size=0.5)
    
    df_ = pd.read_csv(args.train_test_split)
    slides_t = list(df_[df_['train_test_split'] == 'Train']['sample_id'])
    slides_v = list(df_[df_['train_test_split'] == 'Test']['sample_id'])
    
    edata_t = edata[edata.obs.index.isin(slides_t)]
    edata_v = edata[edata.obs.index.isin(slides_v)]
    
    
    # define Y, X, Z
    Yt = torch.Tensor(edata_t.X.astype(np.float32))
    Yv = torch.Tensor(edata_v.X.astype(np.float32))
    Ft = torch.ones([edata_t.X.shape[0], 1])
    Fv = torch.ones([edata_v.X.shape[0], 1])
    Xts = [torch.Tensor(idata.X[idata.obs['slide']==sid].astype(np.float32)) for sid in edata_t.obs.index.values]
    Xvs = [torch.Tensor(idata.X[idata.obs['slide']==sid].astype(np.float32)) for sid in edata_v.obs.index.values]
    Zt = torch.cat([Xt.mean(0, keepdim=True) for Xt in Xts], axis=0)
    Zv = torch.cat([Xv.mean(0, keepdim=True) for Xv in Xvs], axis=0)
    
    
    # tile ids
    tile_ids_t = np.concatenate([idata.obs.loc[idata.obs['slide']==sid, 'id'].values for sid in edata_t.obs.index.values])
    tile_ids_v = np.concatenate([idata.obs.loc[idata.obs['slide']==sid, 'id'].values for sid in edata_v.obs.index.values])
    
    
    R2 = []
    rho = []
    pv = []
    correlation = []
    p_value = []
    for i in tqdm(range(Yt.shape[1])):
        yt = Yt[:,[i]]
        yv = Yv[:,[i]]
        blmm = GLMM(yt, Ft, Zt, repar='unorm_z', lik=NormalLik())
        blmm.optimize()
        xmil0_mean, xmil0_std = blmm.predict(X=Zv, F=Fv)
        R2.append(r2_score(yv.data.numpy(), xmil0_mean.data.numpy()))
        correlation.append(stats.spearmanr(yv.data.numpy(), xmil0_mean.data.numpy())[0])
        p_value.append(stats.spearmanr(yv.data.numpy(), xmil0_mean.data.numpy())[1])
        
        
    df = pd.DataFrame(R2, columns=['R2_score'])
    
    if tissue == 'Osteoarthritis':
        df['gene_id'] = edata_t.var['gene_id'].values
        df['p_value'] = p_value
        df['correlation'] = correlation
    else:
        df['gene_id'] = edata_t.var.index
        df['gene_name'] = edata_t.var.Description.values
        df['p_value'] = p_value
        df['correlation'] = correlation
   
    os.makedirs(args.outdir, exist_ok=True)
    if args.grade is not None:
        outdir_path = join(args.outdir, f'{tissue}_{args.grade}_{model_type}_{args.cluster_i}_glmm_nComp_{n_comp}.csv')
    else:
        outdir_path = join(args.outdir, f'{tissue}_{model_type}_{args.cluster_i}_glmm_nComp_{n_comp}.csv')
    df.to_csv(outdir_path, index=True)


    
    # df.to_csv(outdir_path, index=True)


if __name__ == '__main__':
    
    args = get_args()
    main(args)
    
