from ast import Not
import glob
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
    opts['time'] = '1:00:00'
    opts['nodes'] = 1
    opts['memory'] = 12
    opts['condaenv'] = CONDA_ENV
    outdir = tissue_hyperparameter['outdir']
    job_name = tissue_hyperparameter['tissue']
    opts['name'] = job_name
    opts['stdout'] = join(outdir, 'eval_logs', f'stdout.txt')
    opts['stderr'] = join(outdir, 'eval_logs', f'stderr.txt')
    command = f"python 1_filter_gene.py --outdir {tissue_hyperparameter['outdir']} --tissue {tissue_hyperparameter['tissue']} --gene_file {tissue_hyperparameter['gene_file']}"
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
tissue_list_available = [t.lower() for t in tissue_list_available] # since all the tpm file have annotations in lower case
tissue_hyperparameter = {}


file_path = '../../data/tpm_gene'
download_dir = Path(file_path)
download_dir.mkdir(parents=True, exist_ok=True)

for tissue in tissue_list_available:
    file = f'gene_tpm_2017-06-05_v8_{tissue}.gct.gz'
    local_file = download_dir / file
    download_url = (
        "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/tpms-by-tissue/"
        f"{file}"
    )

    if not local_file.exists():
        print(f"Gene expression does not exist for {tissue}, downloading..")
        os.system(f"wget -O {local_file} {download_url}")

    tissue_hyperparameter['outdir'] = join(OUTPUT_ROOT, 'gene_expression')
    os.makedirs(tissue_hyperparameter['outdir'], exist_ok=True)
    tissue_hyperparameter['tissue'] = tissue
    tissue_hyperparameter['gene_file'] = local_file
    run_jobs(tissue_hyperparameter)
    
