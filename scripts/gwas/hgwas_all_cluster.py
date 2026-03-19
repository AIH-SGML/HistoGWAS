from os.path import join, basename, dirname
import os
os.chdir('..')
import histogwas.emb_gwas as eg
import numpy as np
import pandas as pd
import pylab as pl
import anndata as ad
import pandas as pd
import argparse
import pdb


def load_covs(covfile):
    dfcov0 = pd.read_csv(covfile, sep='\t').dropna()
    dfcov = pd.get_dummies(dfcov0['DTHHRDY'], dtype=np.int8)
    dfcov['SEX'] = 1 * (dfcov0['SEX'].values==2)
    dfcov['AGE'] = dfcov0['AGE'].str.split('-').str.get(0).astype(int) + 5
    dfcov['SID'] = dfcov0['SUBJID']
    dfcov = dfcov.set_index('SID')
    return dfcov


def manhattan(ax, df, hgwas=None, pv_thr=None, colors=None, offset=None, callback=None):
    """
    Utility function to make manhattan plot
    Parameters
    ----------
    ax : pyplot plot
        subplot
    df : pandas.DataFrame
        pandas DataFrame with chrom, pos and pv
    colors : list
        colors to use in the manhattan plot
    offset : float
        offset between in chromosome expressed as fraction of the
        length of the longest chromosome (default is 0.2)
    callback : function
        callback function that takes as input df
    Examples
    --------
    """
    if pv_thr is not None:
        df = df[df["P"] < pv_thr]
    if colors is None:
        colors = ["k", "Gray"]
    if hgwas is not None:
        colors_h = ["b", "green"]
    if offset is None:
        offset = 0.2
    dx = offset * df["POS"].values.max()
    _x = 0
    xticks = []
    for chrom_i in np.unique(df["CHR"].values):
        _df = df[df["CHR"] == chrom_i]
        if chrom_i % 2 == 0:
            color = colors[0]
            if hgwas is not None:
                color_ = colors_h[0]
        else:
            color = colors[1]
            if hgwas is not None:
                color_ = colors_h[1]
        ax.plot(_df["POS"] + _x, -np.log10(_df["P"]), ".", color=color)
        if hgwas is not None:
            ax.plot(_df["POS"] + _x, -np.log10(_df["P_h"]), ".", color=color_)
        if callback is not None:
            callback(_df)
        xticks.append(_x + 0.5 * _df["POS"].values.max())
        _x += _df["POS"].values.max() + dx
    
    ax.set_xticks(xticks)
    ax.set_xticklabels(np.unique(df["CHR"].values))

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", dest="outdir", type=str)
    parser.add_argument("--tissue", dest="tissue", type=str)
    parser.add_argument("--cluster_i", dest="cluster_i", type=str, default=None)
    parser.add_argument("--pcfile", dest="pcfile", type=str, default=None)
    parser.add_argument("--bfile", dest="bfile", type=str, default=None)   
    parser.add_argument("--hfile", dest="hfile", type=str, default=None)
    parser.add_argument("--covfile", dest="covfile", type=str, default=None)
    parser.add_argument("--cluster_type", dest="cluster_type", type=str, default='leiden_0.5')
    args = parser.parse_args()

    return args

def main(args):

    tissue = args.tissue
    outdir = args.outdir
    hfile = args.hfile
    bfile = args.bfile
    pcfile = args.pcfile
    covfile = args.covfile
    
    # pdb.set_trace()

    pc_path = hfile.split('.')[0] + '_pc' + '.h5ad'
    if os.path.exists(pc_path):  # If the pc file already exists, then load it
        idata = ad.read_h5ad(pc_path)
        idata = idata[idata.obs[args.cluster_type]==args.cluster_i]
        idata = idata[idata.obs['SID'].isin((idata.obs['SID'].value_counts()[idata.obs['SID'].value_counts()>=10]).index)]
    else:
        # loading the embedding of tissue, which is an anndata object
        adata = ad.read_h5ad(hfile)

        _ = adata.obs['slide'].str.split('-')
        adata.obs['SID'] = _.str.get(0) + '-'+ _.str.get(1)
        


        # Here I am thresholding the clusters

        adata_new = adata[adata.obs[args.cluster_type]==args.cluster_i]
        adata_new = adata_new[adata_new.obs['SID'].isin((adata_new.obs['SID'].value_counts()[adata_new.obs['SID'].value_counts()>=10]).index)]
        n_comp = 64
        idata = ad.AnnData(adata_new.obsm['X_pca'][:,:n_comp], obs=adata_new.obs)
        idata.write_h5ad(pc_path)


    # read bfile and pcs
    num_pcs = 4
    gdata = eg.read_plink(bfile, pcfile, num_pcs=num_pcs)

    # read covs
    dfcov = load_covs(covfile)


    # match idata, gdata and dfcov
    idata, gdata, dfcov = eg.match_idata_gdata_dfcov(idata, gdata, dfcov, idata_bag_key='SID')


    # define final covs
    dfcov['AGE'] = (dfcov['AGE'] - dfcov['AGE'].mean(0)) / dfcov['AGE'].std(0)
    dfpc = (gdata.obs - gdata.obs.mean(0)) / gdata.obs.std(0)
    dfcovall = pd.concat([dfcov, dfpc], axis=1)

    

    sample_length = len(dfcovall)
    # define embgwas
    embgwas = eg.EmbGWAS(idata, dfcovall, idata_bag_index='SID_index')


    # do gwas
    dfres = embgwas.gwas(gdata, block_size=10000)



    outfile = join(outdir, f'gwas_stat_cluster_{tissue}', f'{tissue}_cluster_{args.cluster_i}_hgwas.csv')
    os.makedirs(dirname(outfile), exist_ok=True)
    dfres = dfres[~np.logical_or(dfres['maf'] < 0.05, dfres['maf']>0.95)] # Exporting the result for snp with MAF atleast 5%
    dfres.to_csv(outfile, index=True)


    # plot results
    outfile = join(outdir, 'plots_cluster', f'{tissue}_cluster_{args.cluster_i}_hgwas.png')
    os.makedirs(dirname(outfile), exist_ok=True)


    df_pv = dfres.loc[dfres['p_value']<1e-2]
    dfss1 = df_pv[['pos', 'chrom', 'p_value']]
    dfss1.columns = ['POS', 'CHR', 'P']
    dfss1['CHR']= dfss1['CHR'].astype('int32')

    pl.figure(1, figsize=(10, 5))
    plt = pl.subplot(111)
    manhattan(plt, dfss1)
    # manhattan(plt, dfss2)
    xlim = pl.xlim()
    pl.plot(xlim, -np.log10(5e-8)*np.ones(2), 'r')
    pl.plot(xlim, -np.log10(1e-6)*np.ones(2), 'Orange')
    pl.xlim(xlim)
    pl.ylim(pl.ylim()[0], 15)
    pl.title(f'for {tissue}')
    pl.xlabel('Chromosome')
    pl.ylabel('-log$_{10}$P')
    pl.tight_layout()
    pl.title(f"for {tissue} with sample {len(dfcovall)}")

    pl.savefig(outfile)



if __name__ =='__main__':
    args = get_args()
    main(args)

    