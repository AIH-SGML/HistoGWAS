import os
from pathlib import Path
from os.path import join, dirname

OUTPUT_ROOT = Path("../../../output/Plots/cluster_signatures")
os.makedirs(OUTPUT_ROOT, exist_ok=True)
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
    
    command = f"python make_plot_for_tissue_cluster.py --tissue {tissue_hyperparameter['tissue']} --outdir {tissue_hyperparameter['outdir']} --hfile_low {tissue_hyperparameter['hfile_low']} --hfile {tissue_hyperparameter['hfile']}"
    submit_job(command, opts)
   


tissue_list = [
    'Lung',
    'Artery_Aorta',
    'Pancreas',
    'Spleen',
    'Breast_Mammary_Tissue',
    'Adipose_Subcutaneous',
    'Colon_Transverse',
    'Artery_Tibial',
    'Stomach',
    'Esophagus_Mucosa',
    'Esophagus_Muscularis',
    'Muscle_Skeletal',
    'Skin_Sun_Exposed_Lower_leg',
    'Thyroid',
    'Nerve_Tibial',
] 


tissue_hyperparameter = {}
for tissue in tissue_list:
    tissue_hyperparameter['tissue'] = tissue
    tissue_hyperparameter['outdir'] = join(OUTPUT_ROOT, tissue)
    tissue_hyperparameter['hfile_low'] = f"/lustre/groups/casale/datasets/gtex/histology/20230425_v2_tiles/stage2/{tissue}/embedding/low_memory_scanpy.h5ad" # path to low memory anndata, with 64pc embedding
    tissue_hyperparameter['hfile'] =  f"/lustre/groups/casale/datasets/gtex/histology/20230425_v2_tiles/stage2/{tissue}/embedding/summary_scanpy.h5ad" # path to full anndata, with retccl clustering
    os.makedirs(tissue_hyperparameter['outdir'], exist_ok=True)
    run_jobs(tissue_hyperparameter)
