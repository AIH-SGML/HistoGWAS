#!/usr/bin/env python3

import glob
import pandas as pd
import os
from os.path import join
import argparse


def get_args():
    parser = argparse.ArgumentParser(description="Concatenate slide summary TSVs into a single summary.tsv")
    parser.add_argument(
        "--outdir", type=str, required=True,
        help="Directory where the merged summary.tsv will be written"
    )
    return parser.parse_args()


def main(args):
    # input pattern
    pattern = join(args.outdir, "summary_*")
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files found matching pattern: {pattern}")

    # read and concatenate
    dfs = []
    for file in files:
        print(f"[INFO] Reading {file}")
        df = pd.read_csv(file, sep="\t")
        dfs.append(df)

    df_all = pd.concat(dfs, axis=0)

    # output
    os.makedirs(args.outdir, exist_ok=True)
    outfile = join(args.outdir, "summary.tsv")
    df_all.to_csv(outfile, sep="\t", index=False)
    print(f"[INFO] Merged summary saved to {outfile}")


if __name__ == "__main__":
    args = get_args()
    main(args)
