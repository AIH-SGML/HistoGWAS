import os
from os.path import join, dirname
from pathlib import Path



# Root directory where the downloaded SVS slides should be stored.
# OUTPUT_ROOT = Path("../../output/slides")

OUTPUT_ROOT = Path('../../data')
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
    f.write(f'#SBATCH --qos={opts["qos"]}\n')
    f.write(f'#SBATCH --mem={opts["memory"]}G\n')
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
    opts['time'] = '12:00:00'
    opts['nodes'] = 4
    opts['memory'] = 128
    opts['condaenv'] = 'milgan'
    outdir = tissue_hyperparameter['outdir']
    job_name = tissue_hyperparameter['tissue']
    opts['name'] = job_name
    opts['stdout'] = join(outdir, 'eval_logs', f'stdout_stage8.txt')
    opts['stderr'] = join(outdir, 'eval_logs', f'stderr_stage8.txt')
   # Only if different resolution is used
    
    command = f"python 2_gene_prediction.py --tissue {tissue_hyperparameter['tissue']} --model_type {tissue_hyperparameter['model_type']} --outdir {tissue_hyperparameter['outdir']} --hfile {tissue_hyperparameter['hfile']} --efile {tissue_hyperparameter['efile']}"
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
model_type = [
    # 'simclr', 
    # 'kimiaNet', 
    'retccl',
    # 'plip',
    # 'Autoencoder',
]

microns = 192
download = False
for  model_type_ in model_type:
    for tissue in tissue_list_available:    
        tissue_hyperparameter['tissue'] = tissue
        tissue_hyperparameter['outdir'] = join(OUTPUT_ROOT, f"gene_prediction")
        
        tissue_hyperparameter['hfile'] =  join(OUTPUT_ROOT, 'embedding', f'{tissue}', f'{tissue}_microns_{microns}', 'summary', f"{tissue}_img_embedding.h5ad")
        os.makedirs(dirname(tissue_hyperparameter['hfile']), exist_ok=True)
        if not (download and os.path.exists(tissue_hyperparameter['hfile'])):
            print(f"Image embedding for {tissue} is not available, downloading..")
            download_url = (
            f"https://zenodo.org/records/18773562/files/{tissue}_img_embedding.h5ad?download=1")
            os.system(f"wget -O {tissue_hyperparameter['hfile']} {download_url}")
        tissue_hyperparameter['efile'] = join(OUTPUT_ROOT,'gene_expression', f"{tissue.lower()}_gene_tpm.h5ad") 
        os.makedirs(tissue_hyperparameter['outdir'], exist_ok=True)
        tissue_hyperparameter['model_type'] = model_type_
        run_jobs(tissue_hyperparameter)