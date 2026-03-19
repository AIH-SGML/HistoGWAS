import anndata as ad
import pandas as pd
import os
import numpy as np
from utils import get_expression, merge_dfs, get_egenes
from histogwas.emb_gwas import GWAS
from limix_core.util.preprocess import gaussianize
import matplotlib.pyplot as plt
import seaborn as sns
from os.path import join, dirname
from tqdm import tqdm
import gseapy as gp
import argparse
import pdb


# This is the list of clusters for each tissue which has been used in the integration of GWAS
tissue_list = {
    'Artery_Aorta': [0, 1, 2, 3, 4, 5, 6, 10],
    'Pancreas': [0,2, 4],
    'Spleen': [0, 1, 2, 3, 5],
    'Breast_Mammary_Tissue': [0, 1, 3, 4, 8],
    'Adipose_Subcutaneous': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
    'Colon_Transverse': [0, 1, 2, 3, 5, 6, 7, 8],
    'Artery_Tibial': [0, 1, 2, 3, 4, 5],
    'Stomach': [0, 2, 3, 4, 5, 6, 7, 8, 10],
    'Esophagus_Mucosa': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
    'Esophagus_Muscularis': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
    'Muscle_Skeletal': [1, 2, 3, 4, 5, 6, 7],
    'Skin_Sun_Exposed_Lower_leg': [0, 1, 2, 3, 4, 5, 6, 7],
    'Thyroid': [0, 1, 3, 2, 4, 5, 6],
    'Nerve_Tibial' : [0, 1, 2, 3, 4, 5, 6, 8],
}




def get_unique_gene(dfgene):
    columns_set = set()
    columns_index = []
    for idx, columns in enumerate(dfgene.columns):
        if columns not in columns_set:
            columns_set.add(columns)
            columns_index.append(idx)
    dfgene = dfgene.iloc[:,columns_index]
    return dfgene

def get_unique_egene(_egene):
    index_set = set()
    final_index = []
    for idx, index_name in enumerate(_egene['gene_name']):
        if index_name not in index_set:
            index_set.add(index_name)
            final_index.append(idx)
    _egene_2 = _egene.iloc[final_index,:]
    return _egene_2


def plot_heat_map(heat_df, heat_map_outdir, type_, tissue):
    plt.figure(figsize=(15, 15))  # Adjust the figure size as needed
    
    cbar_kws = {"label": "-log10P-value"}
    clustermap = sns.clustermap(heat_df,
    vmin=0,          
    vmax=-np.log10(1e-5), cbar_kws=cbar_kws)
    
    
    
    clustermap.fig.suptitle(f'{tissue}_{type_}', fontsize=16)

    # Add labels for the axes
    # clustermap.ax_heatmap.set_xlabel('Genes', fontsize=15)
    # clustermap.ax_heatmap.set_ylabel('Cluster', fontsize=15)
    filename = f'For_{tissue}_retccl_reach_cluster_{type_}.png'
    outfile = join(heat_map_outdir, filename)
    
    # plt.tight_layout()
    
    # Remove legend as it's not typically used in clustermaps
    plt.legend().remove()
    plt.savefig(outfile)
    
    # plt.close()

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument("--tissue", dest="tissue", type=str)
    parser.add_argument("--outdir", dest="outdir", type=str)
    parser.add_argument("--hfile", dest="hfile", type=str)
    args = parser.parse_args()

    return args


def main(args):
    
    
    tissue = args.tissue
    outdir = args.outdir
    hfile = args.hfile
    
     # If the pc file already exists, then load it
    idata = ad.read_h5ad(hfile)
    idata = idata[idata.obs[args.cluster_type]==args.cluster_i]
    idata = idata[idata.obs['SID'].isin((idata.obs['SID'].value_counts()[idata.obs['SID'].value_counts()>=10]).index)]
    
        
    # Filtering based on the clusters which are considered
    idata_ = idata[idata.obs['leiden_0.5'].isin(np.array(tissue_list[tissue]).astype(str))]
    dfX = pd.get_dummies(idata_.obs['leiden_0.5'])
    dfX['slide']  = idata_.obs['slide']
    dfX = dfX.groupby('slide').mean()
    _ = dfX.index.str.split('-').str
    new_name = _[0] + '-' +_[1]
    dfX.index = new_name
    
    # Making association between quantity of tile from each tiles in each cluster

    result_tissue = []
    # pdb.set_trace()
    for cluster_i in dfX.columns:




        snp_id = str(cluster_i)
        tissue = tissue
        snp_df = pd.DataFrame(dfX[snp_id])

        dfgene = get_expression(tissue)
        dfgene = get_unique_gene(dfgene)
        idx_l, idx_r = merge_dfs(dfgene, snp_df, left_index=True, right_index=True)
        dfgene_2 = dfgene.iloc[idx_l]
        dfgene = dfgene.iloc[idx_l]

        snp_df = snp_df.iloc[idx_r]
        df_result = pd.DataFrame(columns = ['p_value'], index=dfgene.columns)


        gwas = GWAS(np.array(snp_df))

        # Here I am gaussianizing the expression
        # U = svd(dfgene.values, full_matrices=False)[0]
        dfgene  = gaussianize(dfgene.values)

        gwas.process(dfgene)   # Outcome is expression
        df_result['p_value'] = gwas.pv[:,0]

        df_result['tissue'] = tissue
        df_result['beta'] = gwas.getBetaSNP()
        df_result['snp_id'] = snp_id
        df_result['cluster_i'] = cluster_i
        df_result = df_result.sort_values(by=['p_value'])
        temp = df_result
        # temp = df_result[df_result['p_value']<(0.05/len(dfgene))]

        # if len(temp) >0:
        _egene = get_egenes(temp['tissue'][0])
        _egene = get_unique_egene(_egene)
        _egene = _egene[_egene['gene_name'].isin(temp.index.values)][['gene_name', 'gene_chr', 'gene_start', 'gene_end']].set_index('gene_name')
        idx_l, idx_r = merge_dfs(temp, _egene, left_index=True, right_index=True)
        temp = temp.iloc[idx_l]
        _egene = _egene.iloc[idx_r]
        temp[[ 'gene_chr', 'gene_start', 'gene_end']] = _egene[[ 'gene_chr', 'gene_start', 'gene_end']]
        # temp['gene_name'] = _egene.index
        temp = temp.sort_values(by=['p_value'])
        result_tissue.append(temp)
        
    # Gene predicted from each fraction of cluster
    outdir_gene = outdir 
    outfile_gene = join(outdir_gene, 'gene_prediction', 'Tile_fraction_in_cluster_to_gene_{tissue}.csv')
    os.makedirs(dirname(outfile_gene), exist_ok=True)

    result_tissue = pd.concat(result_tissue)
    result_tissue.to_csv(outfile_gene, index=True)



# Path enrichment analysis

    # Doing pathway enrichment analysis from the gene predicted from the significant snp
    import warnings

    # Ignore all warnings
    warnings.filterwarnings("ignore")

    outdir = join(outdir,'Pathway_Enrichment')
    os.makedirs(outdir, exist_ok=True)
    enrichment_type = ['upregulated', 'downregulated']
    all_pathway = {}
    for type_ in enrichment_type:
        all_pathway[type_] = []
        for cluster_i in result_tissue['cluster_i'].unique():
            cluster_gene_df_ = result_tissue[result_tissue['cluster_i']==cluster_i]
            cluster_gene_df_['gene_name'] = cluster_gene_df_.index
            background_ = cluster_gene_df_['gene_name'].to_list()
            
            if type_ == 'upregulated':
                cluster_gene_df = cluster_gene_df_[cluster_gene_df_['beta']>0].iloc[:100]  # Predicted genes with positive effect
            else:
                cluster_gene_df = cluster_gene_df_[cluster_gene_df_['beta']<0].iloc[:100]
            gene_list = cluster_gene_df['gene_name'].tolist()


            enr_bg = gp.enrichr(gene_list=gene_list,
                        gene_sets=['MSigDB_Hallmark_2020'],
                                background=background_,
                                
                        )
            path_way = enr_bg.results
            path_way['cluster_i'] = cluster_i
            if len(path_way)>0:
                all_pathway[type_].append(path_way)
        all_pathway[type_] = pd.concat(all_pathway[type_])
        fileName = join(outdir, f'{type_}_ClusterFraction_to_gene_pathwayEnrichment_{tissue}_{type_}.csv')
        all_pathway[type_].to_csv(fileName, index=False)
        
if __name__ == '__main__':
    args = get_args()
    main(args)