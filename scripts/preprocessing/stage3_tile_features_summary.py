#Anndata processing
import anndata as ad
import glob
import re
import pandas as pd
from os.path import join
from tqdm import tqdm
import numpy as np
import anndata
import os
import argparse
import scanpy as sc
import pdb


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument("--input_path", dest="input_path", type=str)
    parser.add_argument("--outdir", dest="outdir", type=str)
    parser.add_argument("--tissue", dest="tissue", type=str)
    args = parser.parse_args()

    return args


def main(args):

    file_list = glob.glob(join(args.input_path,'*'))
    print(len(file_list))
    count=0
    X = []
    obs = []
    
    for i in file_list:
        adata = ad.read_h5ad(i)
        X.append(adata.X)
        obs.append(adata.obs)
        # count+=len(temp)
    if len(X)==1:
        X = X[0]
        obs = obs[0]
    else:
        X = np.concatenate(X, 0)
        obs = pd.concat(obs, 0)
    adata = anndata.AnnData(X, obs=obs)
    filename = join(args.outdir,'summary'+'_normal.h5ad')

    adata.write(filename, compression="gzip")



    sc.tl.pca(adata, svd_solver='arpack', n_comps=256)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
    sc.tl.umap(adata)

    for res in [0.5, 1]:
        sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}')
        filename = join(args.outdir,'summary'+'_scanpy.h5ad')
        adata.write(filename, compression="gzip")
    
if __name__ =='__main__':
    args = get_args()
    main(args)