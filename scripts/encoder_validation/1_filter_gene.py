import pandas as pd
import gzip
import anndata as ad
import numpy as np
import scanpy as sc
from sklearn.preprocessing import power_transform
import os
import argparse



def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", dest="outdir", type=str)
    parser.add_argument("--tissue", dest="tissue", type=str)
    parser.add_argument("--gene_file", dest="gene_file", type=str)
    args = parser.parse_args()

    return args

def main(args):
    tissue = args.tissue





    # load expression data
    file_path = args.gene_file
    with gzip.open(file_path, 'rt') as f:
        # Read the GCT file into a pandas DataFrame
        data = pd.read_csv(f, skiprows=2, delimiter='\t')

    # make anndata
    var = data[['Name', 'Description']].set_index('Name')
    obs = pd.DataFrame({'tissueID': data.keys()[3:].values}).set_index('tissueID')
    X = data[obs.index.values].values.T.astype(np.float32)
    adata = ad.AnnData(obs=obs, var=var, X=X)

    # Only genes expressed in at least 20 samples
    is_expressed = adata.X>1
    Ikeep = (is_expressed).sum(0)>20
    adata = adata[:, Ikeep].copy()


    adata1 = adata.copy()
    Izero = adata1.X==0
    xmin = adata1.X[~Izero].min()

    adata1.X = adata1.X + xmin
    adata1.X = np.log(adata1.X)

    # highly variable genes
    sc.pp.highly_variable_genes(adata1, min_mean=0., max_mean=11, min_disp=0.5)

    adata1 = adata1[:, adata1.var.highly_variable]

    sc.pp.scale(adata1, max_value=10)
    sc.tl.pca(adata1, svd_solver='arpack')
    sc.pp.neighbors(adata1, n_neighbors=10, n_pcs=30)
    sc.tl.umap(adata1)
    sc.tl.leiden(adata1)



    # export
    file_name = f'{tissue}_gene_tpm.h5ad'
    outfile = os.path.join(args.outdir, file_name)
    adata1.write_h5ad(outfile)


if __name__ =='__main__':
    args = get_args()
    main(args)