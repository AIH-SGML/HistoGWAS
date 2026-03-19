from os.path import join
import os
import numpy as np
import openslide
import matplotlib.pyplot as plt
import skimage.filters as sk_filters
import pandas as pd
import torch
from torchvision.utils import make_grid
import argparse
import pdb
from skimage import morphology
from PIL import Image
import logging
import sys
from skimage.color import rgb2gray
import matplotlib.pyplot as plt
import cv2


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument("--slides", dest="slides", type=str,
                        help="TSV/CSV file containing slide metadata with 'path' column")
    parser.add_argument("--slide_idx", dest="slide_idx", type=int, default=0,
                        help="Index of a single slide to process (default: 0)")
    parser.add_argument("--outdir", dest="outdir", type=str,
                        help="Directory where outputs (logs, plots, tiles, TSVs) will be saved")
    parser.add_argument("--fract_fg_min", dest="fract_fg_min", type=float, default=0.5,
                        help="Minimum fraction of foreground required for a tile (default: 0.5)")
    parser.add_argument("--njobs", dest="njobs", type=int, default=1,
                        help="Total number of jobs to split the input file into (default: 100)")
    parser.add_argument("--job_i", dest="job_i", type=int, default=0,
                        help="Index of this job (0-based, default: 0)")
    parser.add_argument(
        "--attempted_microns", dest="attempted_microns", type=float, default=192,
        help="Target microns per tile edge (default: 192)"
    )
    parser.add_argument("--make_plots", action="store_true", dest="make_plots", default=False,
                        help="If set, generate diagnostic plots of tiling and foreground detection"
)
    parser.add_argument(
        "--export_tiles", action="store_true", dest="export_tiles", default=False,
        help="If set, export tile PNGs to disk"
    )
    parser.add_argument("--debug", action="store_true", dest="debug", default=False,
                        help="Enable debug mode (drops into pdb)")
    parser.add_argument("--tissue_count", dest="tissue_count", type=int, default=100,
                        help="Number of tissue samples to use when making plots (default: 100)")   # I am using tissue count only when the plot need to be made
    parser.add_argument("--tissue_name", dest="tissue_name", type=str,
                        help="Name of the tissue type being processed")
    args = parser.parse_args()

    return args


def init_logging(logfile):
    # initialize logging
    log_format = '%(asctime)s %(message)s'
    logging.basicConfig(stream=sys.stdout, level=logging.INFO,
        format=log_format, datefmt='%m/%d %I:%M:%S %p')
    fh = logging.FileHandler(logfile)
    fh.setFormatter(logging.Formatter(log_format))
    logging.getLogger().addHandler(fh)


def detect_foreground(slide_img):
    grayscale = 255 - np.dot(slide_img[..., :3], [0.2125, 0.7154, 0.0721])
    otsu_threshold = sk_filters.threshold_otsu(grayscale)
    binary = grayscale > otsu_threshold
    temp = binary.astype(bool)
    binary_closed = morphology.remove_small_objects(temp, 30)
    return binary_closed

def get_mask_1(slide):
    thumbnail = slide.get_thumbnail((512,512))
    gray_image = rgb2gray(np.array(thumbnail))
    img_ = cv2.resize(gray_image,(512,512))

    _,thresh = cv2.threshold(img_, np.mean(img_), 255, cv2.THRESH_BINARY_INV)
    img = thresh

    # Doing erosion followed by dialation
    # Define the kernel for erosion
    kernel = np.ones((2,2), np.uint8)

    # Perform erosion
    img = cv2.erode(img, kernel, iterations=3)
    img = cv2.dilate(img, kernel, iterations=6)
    img = cv2.erode(img, kernel, iterations=3)
    img = cv2.dilate(img, kernel, iterations=3)
    return img



def get_mask_all(slide):
    thumbnail = slide.get_thumbnail((512,512))
    gray_image = rgb2gray(np.array(thumbnail))
    img_ = cv2.resize(gray_image,(512,512))

    _,thresh = cv2.threshold(img_, np.mean(img_), 255, cv2.THRESH_BINARY_INV)
    img = thresh

    # Doing erosion followed by dialation
    # Define the kernel for erosion
    kernel = np.ones((2,2), np.uint8)
    return img


def make_tiles(slide, attempted_microns=192, level=2):

    assert slide.properties["openslide.mpp-x"] == slide.properties["openslide.mpp-y"]

    # define tile size and effective microns
    mpp = float(slide.properties["openslide.mpp-x"])
    tile_size = int(np.ceil(attempted_microns / mpp))
    microns = mpp * tile_size

    # compute the grid on level 0 and the specified level
    tile_ixs = {}
    tile_iys = {}
    tile_ixs["hr"] = np.arange(0, slide.level_dimensions[0][0], tile_size)
    tile_iys["hr"] = np.arange(0, slide.level_dimensions[0][1], tile_size)
    tile_ixs["lr"] = np.floor(tile_ixs["hr"] / slide.level_downsamples[level]).astype(
        int
    )
    tile_iys["lr"] = np.floor(tile_iys["hr"] / slide.level_downsamples[level]).astype(
        int
    )

    return tile_ixs, tile_iys, microns


def assert_tiles_enough_fg(
    tx, ty, foreground, fract_fg_min=0.9, return_masked_fg=False
):

    df = []
    if return_masked_fg:
        foreground_masked = foreground.copy()
    for i in range(len(tx["hr"]) - 1):
        for j in range(len(ty["hr"]) - 1):
            x0, x1 = tx["lr"][i], tx["lr"][i + 1]
            y0, y1 = ty["lr"][j], ty["lr"][j + 1]
            fract_fg = foreground_masked[y0:y1, x0:x1].mean()
            if fract_fg > fract_fg_min:
                _row = {}
                _row["x0"] = tx["hr"][i]
                _row["y0"] = ty["hr"][j]
                _row["w"] = tx["hr"][i + 1] - tx["hr"][i]
                _row["h"] = ty["hr"][j + 1] - ty["hr"][j]
                _row["fract_fg"] = fract_fg
                _row["index"] = f"{'%05d' % i}_{'%05d' % j}"
                df.append(_row)
            else:
                if return_masked_fg:
                    foreground_masked[y0:y1, x0:x1] = 0
                else:
                    pass
    df = pd.DataFrame(df)
    df.index = ["tile_%.4d" % _ for _ in np.arange(len(df))]

    if return_masked_fg:
        return df, foreground_masked
    else:
        return df


def tile_img(img, tile_ixs, tile_iys):
    img_tiled = img.copy()
    try:
        for tile_ix in tile_ixs:
            img_tiled[:, tile_ix] = 0
        for tile_iy in tile_iys:
            img_tiled[tile_iy, :] = 0
    except:
        print("Some error for corner case")
    return img_tiled


def make_grid_numpy(tiles, nrow=4):
    tiles = [torch.from_numpy(_.transpose((2, 0, 1))) for _ in tiles]
    img = make_grid(tiles, nrow=nrow)
    return img.numpy().transpose((1, 2, 0))


def make_debugging_plots(
    slide, slide_img, foreground_masked, tile_ixs, tile_iys, dftiles, foreground
):

    # build tiled slide img with/o foreground
    tiled_slide_img = tile_img(slide_img, tile_ixs["lr"], tile_iys["lr"])
    tiled_slide_img_fg = tiled_slide_img * foreground_masked[:, :, np.newaxis] / 255.0

    # get 16 random tiles
    idxs = np.random.permutation(len(dftiles))[:16]
    tiles = []
    for idx in idxs:
        _ = dftiles[["x0", "y0", "w", "h"]].iloc[idx]
        _tile = np.asarray(slide.read_region((_["x0"], _["y0"]), 0, (_["w"], _["h"])))[
            ..., :3
        ]
        tiles.append(_tile)
    tiles_img = make_grid_numpy(tiles)

    plt.figure(figsize=(20, 20))
    plt.subplot(221)
    plt.imshow(slide_img)
    plt.subplot(222)
    plt.imshow(tiled_slide_img)
    plt.subplot(223)
    plt.imshow(tiled_slide_img_fg)
    plt.subplot(224)
    plt.imshow(tiles_img)
    return plt


def process_biopsy(
    img_path,
    tissue,
    make_plots=False,
    export_tiles=False,
    outdir=None,
    fract_fg_min=0.9,
    attempted_microns=192,
    level=2,
):

    if args.debug:
        pdb.set_trace()

    # init logging
    name = img_path.split("/")[-1].split(".svs")[0]
    logdir = join(outdir, 'logs')
    os.makedirs(logdir, exist_ok=True)
    logfile = join(logdir, name + '.txt')
    init_logging(logfile)

    # reading low res level
    logging.info(f'reading level {level}')
    slide = openslide.open_slide(img_path)
    slide_img = np.asarray(
        slide.read_region((0, 0), level, slide.level_dimensions[level])
    )[..., :3]
    logging.info('DONE')

    # detecing foreground
    logging.info(f'detecting foreground')
    # foreground = detect_foreground(slide_img)
    if tissue == 'Adipose_Subcutaneous':
    #******************************
        mask = get_mask_1(slide)
    else:
        mask = get_mask_all(slide)
    foreground = detect_foreground(slide_img)

    # Resizing the mask to the shape of forground
    mask  = cv2.resize(mask, (foreground.shape[1], foreground.shape[0]))


    mask[mask<50] = False
    mask[mask>=50] = True
    foreground = mask

    #******************************


    logging.info('DONE')

    # make tile_grid
    logging.info(f'make tile grid')
    tile_ixs, tile_iys, microns = make_tiles(
        slide, attempted_microns=attempted_microns, level=level
    )
    dftiles, foreground_masked = assert_tiles_enough_fg(
        tile_ixs, tile_iys, foreground, fract_fg_min, return_masked_fg=True
    )
    dftiles['microns'] = microns
    logging.info('DONE')

    # make plots
    if make_plots:
        logging.info(f'making plots')
        plotdir = join(outdir, 'plots')
        os.makedirs(plotdir, exist_ok=True)
        plt = make_debugging_plots(
            slide, slide_img, foreground_masked, tile_ixs, tile_iys, dftiles, foreground
        )
        plt.tight_layout()
        plt.savefig(join(plotdir, f"{name}_fg_{fract_fg_min}.png"), dpi=300)
        plt.close()
        logging.info('DONE')

    # export tiles
    if export_tiles:
        logging.info('exporting tiles')
        tiledir = join(outdir, 'tiles', name)
        
        os.makedirs(tiledir, exist_ok=True)
        dftiles['path'] = 'NA'
        for ir in range(len(dftiles)):

            logging.info(f'\texporting tile {ir} / {len(dftiles)}')

            try:
                
                # print(dftiles.iloc[ir])
                _ = dftiles[["x0", "y0", "w", "h", "index"]].iloc[ir]
                _tile = np.asarray(slide.read_region((_["x0"], _["y0"]), 0, (_["w"], _["h"])))[
                    ..., :3
                ]
                _tile = Image.fromarray(_tile)
                _tile = _tile.resize((256,256))
                _outfile = join(tiledir, _["index"]+'.png')
                _tile.save(_outfile)
                dftiles.iat[ir, 7] = _outfile
                
            except:

                print(f"Error in saving image")
                logging.info(f'\tSomething went wrong')
                continue
        logging.info('DONE')

        # check outputs only when doing the exporting of file
        logging.info(f'checking output')
        exists = [os.path.exists(row['path']) for ir, row in dftiles.iterrows()]
        dftiles['PASS'] = exists
        logging.info('DONE')

        # export dataframe
        logging.info(f'exporting dataframe')
        tsvdir = join(outdir, 'tsvs')
        os.makedirs(tsvdir, exist_ok=True)
        outfile = join(tsvdir, f'{name}.tsv')
        dftiles.to_csv(outfile, sep='\t')
        logging.info('DONE')


if __name__ == "__main__":

    # get args
    args = get_args()
    tissue = args.tissue_name
    tissue_count = args.tissue_count

    if args.debug:
        pdb.set_trace()

    # get indexed biopsy
    if args.make_plots:
        dfs = pd.read_csv(args.slides, sep='\t')[:tissue_count]  # 
    else:
        dfs = pd.read_csv(args.slides, sep='\t')
    
    Icv = np.floor(args.njobs * np.arange(len(dfs))/ len(dfs))
    I = Icv == args.job_i
    df = dfs.loc[I]
    img_path = df["path"].values

    # process biopsy 
    
    for img_path_ in img_path:
        try: 
            process_biopsy(
                img_path_,
                tissue,
                make_plots=args.make_plots,
                outdir=args.outdir,
                export_tiles=args.export_tiles,
                attempted_microns=args.attempted_microns,
                fract_fg_min=args.fract_fg_min,
                level=2,
            )
        except:
            print(f"Problem with job {args.job_i}, for image {img_path_}")


