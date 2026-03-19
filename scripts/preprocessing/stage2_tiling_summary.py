import pdb
import pandas as pd
import glob
import os
import numpy as np
from tqdm import tqdm
import argparse
from os.path import join


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument("--slides", dest="slides", type=str)
    parser.add_argument("--outdir", dest="outdir", type=str)
    parser.add_argument("--debug", action="store_true", dest="debug", default=False)
    args = parser.parse_args()

    return args


def main():

    args = get_args()

    tsv_ptn = join(args.outdir, 'tsvs', '%s.tsv')
    log_ptn = join(args.outdir, 'logs', '%s.txt')
    tiles_ptn = join(args.outdir, 'tiles')


    # get indexed biopsy
    dfs = pd.read_csv(args.slides, sep='\t')
    img_paths = dfs["path"].values
    failed_biopsies = {}
    missing_tiles = {}
    dft = []
    
    for img_path in tqdm(img_paths):

        # assert if exists 
        _name = img_path.split('/')[-1].split('.svs')[0]
        _tsv = tsv_ptn % _name
        _log = log_ptn % _name
        _tiles = join(tiles_ptn, _name)
        tiles_len = len(glob.glob(_tiles))
        if not os.path.exists(_tsv) or not tiles_len>0:
            failed_biopsies[_name] = _log
            continue

        # if exists asserts if tiles are missing
        _dft = pd.read_csv(_tsv, sep='\t', index_col=0)
        _dft['tile_id'] = _dft.index
        _dft['slide'] = _name
        _dft['id'] = _dft['slide'] + '_' + _dft['tile_id']
        _dft = _dft.set_index('id', drop=True)
        dft.append(_dft)

        if not _dft['PASS'].all():
            n_missing_tiles = (1.  - _dft['PASS']).sum()
            missing_tiles[_name] = n_missing_tiles
            print(f'{_name} missing {n_missing_tiles} tiles')


    outfile = join(args.outdir, 'tiles.tsv')
    dft = pd.concat(dft, axis=0)
    dft.to_csv(outfile, sep='\t')

    # check errors in logs
    for key in failed_biopsies:
        os.system(f"cat {failed_biopsies[key]}")




if __name__=='__main__':

    main()
