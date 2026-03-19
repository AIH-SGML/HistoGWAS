import pandas as pd
import numpy as np
import scanpy as sc
import numpy as np
from PIL import Image
from torchvision.utils import make_grid
import pylab as plt
import matplotlib.colors as mcolors
import anndata as ad
from os.path import join
import scanpy as sc
import numpy as np
from PIL import Image
from torchvision.utils import make_grid
import pylab as plt
import matplotlib.colors as mcolors
import os
import torch
import argparse

tissue_list_cluster = {
    'Skin_Sun_Exposed_Lower_leg': [1, 7, 2, 0, 4, 3, 5],
 'Esophagus_Muscularis': [6, 7, 1, 5, 0, 8, 2, 3, 4],
 'Stomach': [0, 7, 2, 5, 10, 4],
 'Nerve_Tibial': [0, 6],
 'Colon_Transverse': [2, 8, 0, 1, 5, 6, 7, 3],
 'Esophagus_Mucosa': [7, 0, 2, 3, 4, 8, 1],
 'Artery_Tibial': [0, 2, 4, 5],
 'Breast_Mammary_Tissue': [1, 8, 0],
 'Adipose_Subcutaneous': [6, 5, 4, 1, 7, 0, 8, 2, 3],
 'Muscle_Skeletal': [5, 4, 1, 3, 2, 7],
 'Thyroid': [2, 6, 5, 3, 1, 4, 0]
                      }






h5ad_path = f"/lustre/groups/casale/datasets/gtex/histology/20230425_v2_tiles/stage2/%s/embedding/low_memory_scanpy.h5ad"
def make_grid_numpy(tiles, nrow):
    tiles = [torch.from_numpy(_.transpose((2, 0, 1))) for _ in tiles]
    img = make_grid(tiles, nrow=nrow)
    return img.numpy().transpose((1, 2, 0))

def add_contour(img, contour_size, contour_color):
    color_rgb = np.array(mcolors.to_rgb(contour_color))
    vert = np.tile(color_rgb, (contour_size, img.shape[1], 1))
    img = np.concatenate([vert, img / 255., vert], 0)
    hori = np.tile(color_rgb, (img.shape[0], contour_size, 1))
    img = np.concatenate([hori, img, hori], 1)
    return img

def load_anndata(tissue):
    
    print(f'.. loading {tissue} data')
    idata = ad.read_h5ad(h5ad_path % tissue)
    
    print('.. filtering based on the clusters which are considered')
    Ikeep = idata.obs['leiden_0.5'].isin(np.array(tissue_list_cluster[tissue]).astype(str))
    idata = idata[Ikeep].copy()
    
    print('.. number of slides with at least 20 tiles in each cluster')
    dfX = pd.get_dummies(idata.obs['leiden_0.5'])
    dfX['slide']  = idata.obs['slide']
    print((dfX.groupby('slide').sum()>20).sum(0))
    
    return idata

def get_args():
    
    parser = argparse.ArgumentParser()

    
    '''split_index basically tells us from where we start taking the partition in dataset for that pertical job
        split_cout basically stores the  number of jobs that is created
    '''
    parser.add_argument("--tissue", dest="tissue", type=str)
    parser.add_argument("--outdir", dest="outdir", type=str)
    parser.add_argument("--figdir", dest="figdir", type=str)
    args = parser.parse_args()
    
    return args


def main(args):
    
    tissue = args.tissue
    figdir = args.figdir
    outdir = args.outdir
    
    
    path = join(args.clusterFraction_dir, 'gene_prediction' , 'Tile_fraction_in_cluster_to_gene_{tissue}.csv')
    
    dfout = pd.read_csv(path)
    dfout = dfout.rename(columns={'Unnamed: 0':'gene_name'})
    print(f"********************{tissue}**********************")
    for cluster in np.sort(dfout['cluster_i'].unique()):
        Isign = dfout.loc[dfout['cluster_i']==cluster, 'p_value'] < (0.05/len(dfout))
        print(f'Custer {cluster}: {Isign.sum()} / {Isign.shape[0]} genes')
    
    
    dfgene = {}
    for key in np.sort(dfout['cluster_i'].unique()): dfgene[int(key)] = dfout.loc[dfout['cluster_i']==key]
    
    for cluster in dfgene.keys():
        Isign = dfgene[cluster]['p_value']<(0.05/len(dfgene[cluster]))
        dftop5 = dfgene[cluster].loc[Isign].sort_values('p_value').iloc[:5]
        genes = (dftop5['gene_name'].values)
        signs = (dftop5['beta'].values>0)
        print("")
        print(cluster)
        if len(dftop5)>0:
            for gene, sign in zip(genes, signs):
                _ = '+' if sign else '-'
                print(f"{gene} ({_})")
            
    
    # Reading the pathway
    enrichment_type = ['upregulated', 'downregulated']
    df_type = {}
    cluster_all = set()
    for type_ in enrichment_type:
        path = join(args.clusterFraction_dir, 'Pathway_Enrichment', f"{type_}_ClusterFraction_to_gene_pathwayEnrichment_{tissue}_{type_}.csv")
        df = pd.read_csv(path)
        df_type[type_] = df
        cluster_all.update(df['cluster_i'].unique())

    for cluster_i in list(cluster_all):
        print(f"For cluster {cluster_i}")
        print(f"***********************")
        for type_ in enrichment_type:
            df_ = df_type[type_]
            df_ = df_[df_['cluster_i']==cluster_i]
            pathway = df_[df_['Adjusted P-value'] < 0.05][:5]
            if len(pathway)> 0:
                for path_term in pathway['Term'].values:
                    _ = '+' if type_ == 'upregulated' else '-'
                    print(f"{path_term} {_}")
    
    
    # This part is for plotting the figure
    
    n_images = 6
    nrow = 3
    alpha0 = 0.01
    alpha1 = 0.1
    n_clusters = len(tissue_list_cluster[tissue])
    idata = load_anndata(tissue)
    plt.figure(1, figsize=(15, 5))
    for cl_i, cl in enumerate(tissue_list_cluster[tissue]):

        _color = 'C%d' % cl_i

        # plot umap
        ax = plt.subplot(3, n_clusters, cl_i + 1)
        plt.title('Cluster %d' % (cl ))
        Ikeep = idata.obs['leiden_0.5']==str(cl)
        idata0 = idata[~Ikeep].copy()
        idata1 = idata[Ikeep].copy()
        plt.plot(idata0.obsm['X_umap'][:,0], idata0.obsm['X_umap'][:,1], '.k', ms=0.5, alpha=alpha0)
        plt.plot(idata1.obsm['X_umap'][:,0], idata1.obsm['X_umap'][:,1], '.', ms=0.5, color=_color, alpha=alpha1)
        ax.axis('off')

        # tiles 
        ax = plt.subplot(3, n_clusters, n_clusters + cl_i + 1)
        #plt.title('C%d' % cl, color=_color, fontweight='bold')
        Ikeep = idata.obs['leiden_0.5']==str(cl)
        _idata = idata[Ikeep].copy()
        sc.pp.subsample(_idata, n_obs=n_images)
        tiles = [np.array(Image.open(path)) for path in _idata.obs['path'].values]
        img = make_grid_numpy(tiles, nrow=nrow)
        img = add_contour(img, 20, _color)
        plt.imshow(img)
        ax.axis('off')

        # qqplot 
        ax = plt.subplot(3, n_clusters, 2 * n_clusters + cl_i + 1)
        _dfout = dfout.loc[dfout['cluster_i']==cl]
        pv = _dfout['p_value'].values
        beta = _dfout['beta'].values
        idxs = np.argsort(pv)
        pv = pv[idxs]; beta = beta[idxs]
        pvo = -np.log10(np.sort(pv))
        pvo = np.clip(pvo, 0, 10)
        pve = -np.log10(np.linspace(0, 1, pvo.shape[0]+2)[1:-1])
        #pv_thr = -np.log10(0.05 / float(pvo.shape[0]))
        Isign = pv<0.05/(len(_dfout))
        plt.title(f'{Isign.sum()} genes')
        plt.plot([0, 4], [0, 4], color='Gray')
        plt.plot(pve[~Isign], pvo[~Isign], '.k')
        plt.plot(pve[Isign], pvo[Isign], '.', color='slategray', label='adj.P<0.05')
        #plt.plot([0, 4], pv_thr * np.ones(2), '--', 'r')
        plt.xlabel('Expected -log$_{10}$P')
        plt.ylabel('Observed -log$_{10}$P')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.xlim(0, 4)
        plt.ylim(0, 10.5)
        #plt.legend(loc='lower right', frameon=False, numpoints=1, prop={'size': 8})

        #if cl_i==1: break

        plt.tight_layout()
        plt.savefig(os.path.join(figdir, f'clusters_full_{tissue}.png'), dpi=300)



if __name__ == '__main__':
    args = get_args()
    main(args)