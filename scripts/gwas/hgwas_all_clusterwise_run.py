import os
from os.path import join, dirname
from pathlib import Path



# specify output where the association results will be stored

OUTPUT_ROOT = Path("../../../output/association")
os.makedirs(OUTPUT_ROOT, exist_ok=True)

# These are the cluster signature which are already defined in preprocessing stage
tissue_list_cluster = {
    # 'Skin_Sun_Exposed_Lower_leg': [0, 1, 2, 3, 4, 5, 7],
#  'Esophagus_Muscularis': [0, 1, 2, 3, 4, 5, 6, 7, 8],
#  'Stomach': [0, 2, 4, 5, 7, 10],
#  'Nerve_Tibial': [0, 6],
#  'Colon_Transverse': [0, 1, 2, 3, 5, 6, 7, 8],
#  'Esophagus_Mucosa': [0, 1, 2, 3, 4, 7, 8],
#  'Artery_Tibial': [0, 2, 4, 5],
#  'Breast_Mammary_Tissue': [0, 1, 8],
#  'Adipose_Subcutaneous': [0, 1, 2, 3, 4, 5, 6, 7, 8],
#  'Muscle_Skeletal': [1, 2, 3, 4, 5, 7],
 'Thyroid': [0, 1, 2, 3, 4, 5, 6]}


def submit_job(command, opts, cluster_i):

    os.makedirs(dirname(opts["stdout"]), exist_ok=True)
    os.makedirs(dirname(opts["stderr"]), exist_ok=True)

    f = open(f"submit_{cluster_i}.sh", "w")

    f.write(f'#!/bin/bash\n')
    f.write(f'\n')
    f.write(f'#SBATCH -J {opts["name"]}\n')
    f.write(f'#SBATCH -o {opts["stdout"]}\n')
    f.write(f'#SBATCH -e {opts["stderr"]}\n')
    f.write(f'#SBATCH -p {opts["queue"]}\n')
    f.write(f'#SBATCH -t {opts["time"]}\n')
    f.write(f'#SBATCH -c {opts["nodes"]}\n')
    f.write(f'#SBATCH --mem={opts["memory"]}G\n')
    f.write(f'#SBATCH --qos={opts["qos"]}\n')
    f.write(f'#SBATCH --nice=10000\n')
    f.write(f'\n')
    f.write(f'source $HOME/.bashrc\n')
    f.write(f'conda activate {opts["condaenv"]}\n')
    f.write(command)
    f.write(f'\n')
    f.close()

    os.system(f'sbatch submit_{cluster_i}.sh')
    os.system(f'rm submit_{cluster_i}.sh')
    
def run_jobs(tissue_hyperparameter):
    opts = {}
    opts['queue'] = 'cpu_p'
    opts['time'] = '12:00:00'
    opts['qos'] = 'cpu_normal'
    opts['nodes'] = 4
    opts['memory'] = 64
    opts['condaenv'] = 'histoGWAS_2'
    outdir = tissue_hyperparameter['outdir']
    job_name = tissue_hyperparameter['tissue']
    opts['name'] = job_name
    opts['stdout'] = join(outdir, 'eval_logs', f"{tissue_hyperparameter['tissue']}_stdout_{tissue_hyperparameter['cluster_i']}.txt")
    opts['stderr'] = join(outdir, 'eval_logs', f"{tissue_hyperparameter['tissue']}_stderr_{tissue_hyperparameter['cluster_i']}.txt")
    command = f"python hgwas_all_cluster.py --outdir {tissue_hyperparameter['outdir']} --tissue {tissue_hyperparameter['tissue']} --cluster_i {tissue_hyperparameter['cluster_i']} --pcfile {tissue_hyperparameter['pcfile']} --bfile {tissue_hyperparameter['bfile']} --hfile {tissue_hyperparameter['hfile']} --covfile {tissue_hyperparameter['covfile']}"
    # os.system(command)
    submit_job(command, opts, tissue_hyperparameter['cluster_i'])   


def main():

    tissue_hyperparameter = {}
    for tissue in tissue_list_cluster.keys():
        tissue_hyperparameter['tissue'] = tissue

        for cluster_i in tissue_list_cluster[tissue]:
            
            tissue_hyperparameter['cluster_i'] = str(cluster_i)
            tissue_hyperparameter['outdir'] = OUTPUT_ROOT
            tissue_hyperparameter['pcfile'] = '../../data/wgs/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.pca.eigenvec'
            tissue_hyperparameter['bfile'] = '../../data/wgs/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF02'
            tissue_hyperparameter['hfile'] = f"../../data/embedding/thyroid/thyroid_microns_192/summary_scanpy.h5ad"
            tissue_hyperparameter['covfile'] = '../../data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt'

            os.makedirs(tissue_hyperparameter['outdir'], exist_ok=True)
            run_jobs(tissue_hyperparameter)
    


if __name__ =='__main__':
    main()

    