import pandas as pd
import argparse
import os
from os.path import join
from tqdm import tqdm
import numpy as np
import pdb


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples-metadata", dest="samples_metadata", type=str, required=True,
                        help="CSV file with tissue sample metadata (e.g. skin_sun_samples_metadata.csv)")
    parser.add_argument("--outdir", dest="outdir", type=str, required=True,
                        help="Directory where outputs will be saved")
    parser.add_argument("--rewrite", action="store_true", dest="rewrite", default=False,
                        help="Rewrite existing files if found")
    parser.add_argument("--debug", action="store_true", dest="debug", default=False,
                        help="Enter debug mode")
    parser.add_argument("--n_jobs", dest="n_jobs", type=int, default=1,
                        help="Total number of jobs (split file into this many parts)")
    parser.add_argument("--job_i", dest="job_i", type=int, default=0,
                        help="Index of this job (0-based)")
    parser.add_argument("--max_slides", dest="max_slides", type=int, default=None,
                        help="Optional cap on number of slides to process from this job shard (default: all)")
    return parser.parse_args()


def main(args):
    # prepare folders
    slides_dir = join(args.outdir, "slides")
    summary_dir = join(args.outdir, "summary")
    os.makedirs(slides_dir, exist_ok=True)
    os.makedirs(summary_dir, exist_ok=True)

    # read metadata file
    df = pd.read_csv(args.samples_metadata)
    df.index = np.arange(len(df))

    # split into chunks
    if args.n_jobs > 1:
        chunks = np.array_split(df, args.n_jobs)
        if args.job_i < 0 or args.job_i >= args.n_jobs:
            raise ValueError(f"job_i must be between 0 and {args.n_jobs-1}")
        df = chunks[args.job_i].reset_index(drop=True)

    if args.max_slides is not None:
        if args.max_slides <= 0:
            raise ValueError("max_slides must be > 0")
        df = df.head(args.max_slides).reset_index(drop=True)

    # add URL + path columns
    url_base = "https://brd.nci.nih.gov/brd/imagedownload"
    df["url"] = df["Tissue Sample ID"].apply(lambda sid: join(url_base, str(sid)))
    df["path"] = df["Tissue Sample ID"].apply(lambda sid: join(slides_dir, f"{sid}.svs"))

    if args.debug:
        pdb.set_trace()

    # loop and download
    results = []
    for _, row in tqdm(df.iterrows(), total=df.shape[0], desc="Downloading slides"):
        if os.path.exists(row["path"]) and not args.rewrite:
            status = "skipped_existing"
        else:
            exit_code = os.system(f'wget -q -O "{row["path"]}" "{row["url"]}"')
            status = "success" if exit_code == 0 else "failed"
        results.append({
            "Tissue Sample ID": row["Tissue Sample ID"],
            "url": row["url"],
            "path": row["path"],
            "status": status
        })

    # save summary (per job)
    summary_name = f"summary_job{args.job_i:04d}_of_{args.n_jobs:04d}.tsv"
    df_summary = pd.DataFrame(results)
    df_summary.to_csv(join(summary_dir, summary_name), sep="\t", index=False)
    print(f"[INFO] Summary saved to {join(summary_dir, summary_name)}")


if __name__ == "__main__":
    args = get_args()
    main(args)
