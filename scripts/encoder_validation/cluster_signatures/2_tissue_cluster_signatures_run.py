import pandas as pd
import numpy as np
import glob
import os
from os.path import join, dirname
from pathlib import Path


# OUTPUT_ROOT = Path("../../output/slides")

OUTPUT_ROOT = Path('/lustre/groups/casale/code/users/shubham.chaudhary/output/histoGWAS_check/downstream_characterization')
# Conda environment that contains the dependencies for stage1_get_data.py.
CONDA_ENV = "histoGWAS_2"



def submit_job(command, opts):

    os.makedirs(dirname(opts["stdout"]), exist_ok=True)
    os.makedirs(dirname(opts["stderr"]), exist_ok=True)

    f = open(f"submit.sh", "w")

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

    os.system(f'sbatch submit.sh')
    os.system('rm submit.sh')
    
def run_jobs(tissue_hyperparameter):
    opts = {}
    opts['queue'] = 'cpu_p'
    opts['qos'] = 'cpu_normal'
    opts['time'] = '1:00:00'
    opts['nodes'] = 4
    opts['memory'] = 64
    opts['condaenv'] = CONDA_ENV
    outdir = tissue_hyperparameter['outdir']
    job_name = tissue_hyperparameter['tissue']
    opts['name'] = job_name
    opts['stdout'] = join(outdir, 'eval_logs', f"stdout_{tissue_hyperparameter['tissue']}.txt")
    opts['stderr'] = join(outdir, 'eval_logs', f"stderr_{tissue_hyperparameter['tissue']}.txt")
    command = f"python 2_tissue_cluster_signatures.py --outdir {tissue_hyperparameter['outdir']} --tissue {tissue_hyperparameter['tissue']} --figdir {tissue_hyperparameter['figdir']} --clusterFraction_dir {tissue_hyperparameter['clusterFraction_dir']} "
    submit_job(command, opts)



tissue_list = [
    # 'Breast_Mammary_Tissue',
    # 'Adipose_Subcutaneous',
    # 'Colon_Transverse',
    # 'Artery_Tibial',
    # 'Stomach',
    # 'Esophagus_Mucosa',
    # 'Esophagus_Muscularis',
    # 'Muscle_Skeletal',
    #  'Nerve_Tibial',
    # 'Skin_Sun_Exposed_Lower_leg',
    'Thyroid',
]
tissue_hyperparameter = {}



for tissue in tissue_list:
    tissue_hyperparameter['outdir'] = join(OUTPUT_ROOT, 'tissue_cluster_signatures', tissue)
    os.makedirs(tissue_hyperparameter['outdir'], exist_ok=True)
    tissue_hyperparameter['figdir'] = join(tissue_hyperparameter['outdir'], 'tissue_cluster_signatures', 'figdir')
    tissue_hyperparameter['clusterFraction_dir'] = join(OUTPUT_ROOT, 'gene_prediction', 'clusterFraction')

    os.makedirs(tissue_hyperparameter['outdir'], exist_ok=True)
    os.makedirs(tissue_hyperparameter['figdir'], exist_ok=True)
    tissue_hyperparameter['tissue'] = tissue
    run_jobs(tissue_hyperparameter)


