#!/usr/bin/env python3

from __future__ import annotations

import argparse
import logging
import multiprocessing as mp
import sys
from pathlib import Path

from worker import process_sample
from qc import write_summary, print_summary


def parse_args():

    parser = argparse.ArgumentParser(
        description="Extract 18S rRNA regions from long-read FASTQ using Barrnap."
    )

    parser.add_argument(
        "fastq",
        nargs="+",
        help="Input FASTQ(.gz) files",
    )

    parser.add_argument(
        "-o",
        "--outdir",
        default="trimmed",
        help="Output directory",
    )

    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=8,
        help="Barrnap threads per sample",
    )

    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=1,
        help="Number of samples processed simultaneously",
    )

    parser.add_argument(
        "--kingdom",
        default="euk",
        choices=["euk", "bac", "arc", "mito"],
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
    )

    return parser.parse_args()


def setup_logging(verbose):

    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )


def main():

    args = parse_args()

    setup_logging(args.verbose)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fastqs = [Path(f) for f in args.fastq]

    missing = [f for f in fastqs if not f.exists()]

    if missing:

        for f in missing:
            logging.error("File not found: %s", f)

        sys.exit(1)

    jobs = [
        (
            fq,
            outdir,
            args.threads,
            args.kingdom,
            args.overwrite,
        )
        for fq in fastqs
    ]

    if args.processes == 1:

        results = [
            process_sample(*job)
            for job in jobs
        ]

    else:

        with mp.Pool(args.processes) as pool:

            results = pool.starmap(
                process_sample,
                jobs,
            )

    summary = outdir / "trim18s.summary.tsv"

    write_summary(
        results,
        summary,
    )

    print_summary(results)

    logging.info("Summary written to %s", summary)


if __name__ == "__main__":
    main()
