import anndata as ad
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from os.path import join, dirname
from histogwas.emb_gwas.utils import df_match
import pylab as pl
import scipy.stats as st
import scanpy as sc
import torch
from torchvision.utils import make_grid
from PIL import Image
import argparse



# This is the list of clusters for each tissue which has been used in the integration of GWAS
tissue_list = {
        'Lung': [0, 1, 2, 3, 5, 6],
        'Artery_Aorta': [0, 1, 2, 3, 4, 5, 6, 10],
        'Pancreas': [0,2, 4],
        'Spleen': [0, 1, 2, 3, 5],
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
        'Osteoarthritis': [0, 1, 2, 3, 4, 5, 7],
    }


def make_grid_numpy(tiles, nrow):
    tiles = [torch.from_numpy(_.transpose((2, 0, 1))) for _ in tiles]
    img = make_grid(tiles, nrow=nrow)
    return img.numpy().transpose((1, 2, 0))

def make_plots_clusters(idata, tissue, cluster_i, outdir):
    # cluster_i = '1'
    idata_temp = idata[idata.obs['leiden_0.5'] == cluster_i]
    n_images = 25
    sc.pp.subsample(idata_temp, n_obs=n_images)
    tiles = [np.array(Image.open(path)) for path in idata_temp.obs['path'].values]
    img_cluster = make_grid_numpy(tiles, nrow=5)


    plotdir_concept = outdir
    os.makedirs(plotdir_concept, exist_ok=True)
    pl.figure(1, figsize=(20, 20))
    plt = pl.subplot(111)
    pl.title(f'Cluster {cluster_i} ')
    pl.imshow(img_cluster)
    plt.set_xticks([])
    plt.set_yticks([])
    pl.tight_layout()
    pl.savefig(join(plotdir_concept, f'cluster_{cluster_i}.png'), dpi=100)
    pl.close()
    
    
    
def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument("--tissue", dest="tissue", type=str)
    parser.add_argument("--outdir", dest="outdir", type=str)
    parser.add_argument("--hfile", dest="hfile", type=str, help="path to the anndata file")
    parser.add_argument("--hfile_low", dest="hfile_low", type=str, help="path to the low memory anndata file")
    args = parser.parse_args()

    return args
    
def main(args):
    
    
    
    
    
    # Saving the low memory anndata
    outdir = args.outdir
    tissue = args.tissue
    out_path = args.hfile_low
    if not os.path.exists(out_path):
        # reading the anndata

        path = args.hfile
        adata = ad.read_h5ad(path)


        dfX = pd.get_dummies(adata.obs['leiden_0.5'])
        dfX['slide']  = adata.obs['slide']
        dfX = dfX.groupby('slide').mean()


        cluster_i = ((dfX>0.05).sum(0)>20)  # using this for crtanspath
        # cluster_i = ((dfX>0.10).sum(0)>100)
        cluster_i = cluster_i[cluster_i==True].index

        adata_new = adata[adata.obs['leiden_0.5'].isin(cluster_i)]



        # creating low memory cluster
        n_comp = 64
        sample_length = len(adata_new.obs['slide'].unique())
        idata = ad.AnnData(adata_new.obsm['X_pca'][:,:n_comp], obs=adata_new.obs, obsm=adata_new.obsm, uns=adata_new.uns,  obsp=adata_new.obsp)
        _ = idata.obs['slide'].str.split('-')
        idata.obs['SID'] = _.str.get(0) + '-' + _.str.get(1)

        idata.write_h5ad(out_path)

    else:
        idata = ad.read_h5ad(out_path)
    print(f"Making plots for {tissue}")
    for cluster_i in tissue_list[tissue]:
        cluster_i = str(cluster_i)
        make_plots_clusters(idata, tissue, cluster_i, outdir)
        
if __name__ == '__main__':
    
    
    args = get_args()
    main(args)