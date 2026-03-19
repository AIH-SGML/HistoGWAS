import pandas as pd
import numpy as np
import glob
import os
from os.path import join, dirname
from pathlib import Path


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
    os.system(f'rm submit.sh')
    
def run_jobs(tissue_hyperparameter):
    opts = {}
    opts['queue'] = 'cpu_p'
    opts['time'] = '1:00:00'
    opts['qos'] = 'cpu_normal'
    opts['nodes'] = 4
    opts['memory'] = 64
    opts['condaenv'] = 'milgan'
    outdir = tissue_hyperparameter['outdir']
    job_name = tissue_hyperparameter['tissue']
    opts['name'] = job_name
    opts['stdout'] = join(outdir, 'eval_logs', f'stdout_stage8_plot.txt')
    opts['stderr'] = join(outdir, 'eval_logs', f'stderr_stage8_plot.txt')
    
    command = f"python 1_gene_prediction_clusterFraction.py --tissue {tissue_hyperparameter['tissue']} --outdir {tissue_hyperparameter['outdir']} --hfile {tissue_hyperparameter['hfile']}"
    # os.system(command)
    submit_job(command, opts)
   

tissue_list_available = [
        # 'Artery_Aorta',
        # 'Pancreas',
        # 'Spleen',
        # 'Breast_Mammary_Tissue',
        # 'Adipose_Subcutaneous',
        # 'Colon_Transverse',
        # 'Artery_Tibial',
        # 'Stomach',
        # 'Esophagus_Mucosa',
        # 'Esophagus_Muscularis',
        # 'Muscle_Skeletal',
        # 'Skin_Sun_Exposed_Lower_leg',
        'Thyroid',
        # 'Nerve_Tibial',
    ] 


tissue_hyperparameter = {}
for tissue in tissue_list_available:
    tissue_hyperparameter['tissue'] = tissue
    tissue_hyperparameter['outdir'] = join(OUTPUT_ROOT, 'gene_prediction', 'clusterFraction')
    tissue_hyperparameter['hfile'] = f'/lustre/groups/casale/datasets/gtex/histology/20230425_v2_tiles/stage2/{tissue}/embedding/summary_scanpy_pc.h5ad'
    os.makedirs(tissue_hyperparameter['outdir'], exist_ok=True)
    run_jobs(tissue_hyperparameter)


