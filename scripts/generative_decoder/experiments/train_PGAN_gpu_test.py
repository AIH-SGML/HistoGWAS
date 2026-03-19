import os
from pathlib import Path
import json


OUTPUT_ROOT = Path("../../../../output/PGAN")
os.makedirs(OUTPUT_ROOT, exist_ok=True)
PROJECT_ROOT = Path(__file__).resolve().parent.parent

def create_json_file(tissue, dimEmb):

    data_tissue = {
    "pathDB": f"/lustre/groups/casale/datasets/gtex/histology/20230425_v2_tiles/stage2/{tissue}/embedding/summary_scanpy_pc.h5ad",
    "config": {
        "maxIterAtScale": [
        48000,
        96000,
        96000,
        96000,
        96000,
        96000,
        1000000
        ]
    },
    'dimEmb': dimEmb,
    }
    file_path = Path(f"{PROJECT_ROOT}/config/config.json")

    with open(file_path, "w") as json_file:
        json.dump(data_tissue, json_file, indent=4)

os.chdir('..')
tissue = 'Thyroid'
dimEmb = 64 # feature dimension
outdir = OUTPUT_ROOT
create_json_file(tissue, dimEmb)

command = f"python {PROJECT_ROOT}/train.py PGAN -c {PROJECT_ROOT}/config/config.json -n {tissue} --dir {outdir} --dimEmb {dimEmb}"
os.system(command)
