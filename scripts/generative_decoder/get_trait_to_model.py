import emb_gwas as eg
import pandas as pd
import numpy as np

def get_snp(snp, gdata):
    df_gene = pd.DataFrame(gdata.X[:,np.where(gdata.var.index == snp)[0]].compute())
    df_gene.index = gdata.obs.index
    return df_gene

def get_trait(snp):
    pcfile='/lustre/groups/casale/datasets/gtex/wgs/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.pca.eigenvec'
    bfile='/lustre/groups/casale/datasets/gtex/wgs/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF02'

    # # Reading the genotype
    num_pcs = 4
    gdata = eg.read_plink(bfile, pcfile, num_pcs=num_pcs)
    df_gene = get_snp(snp, gdata)
    return df_gene, gdata