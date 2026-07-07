#!/usr/bin/env python3
"""
qc.py

Write per-sample summary statistics.
"""

from __future__ import annotations

import csv
from pathlib import Path


SUMMARY_COLUMNS = [
    # Sample
    "sample",

    # Input
    "reads",
    "bases",
    "mean_input",
    "median_input",

    # Barrnap annotation
    "reads_with_18S",
    "reads_with_5_8S",
    "reads_with_28S",
    "reads_with_all_three",
    "reads_multiple_18S",

    "18S_features",
    "5.8S_features",
    "28S_features",

    # Output
    "trimmed_reads",
    "trimmed_bases",

    # Trimmed lengths
    "mean_trimmed",
    "median_trimmed",
    "min_trimmed",
    "max_trimmed",
    "n50_trimmed",

    # Retention
    "retained_pct",

    # Runtime
    "runtime_sec",
]


def write_summary(results, outfile):
    """
    Write a TSV summary of all processed samples.

    Parameters
    ----------
    results : list[dict]
        List of dictionaries returned by worker.process_sample()

    outfile : str | Path
        Output TSV file.
    """

    outfile = Path(outfile)

    with outfile.open("w", newline="") as fh:

        writer = csv.DictWriter(
            fh,
            fieldnames=SUMMARY_COLUMNS,
            delimiter="\t",
            extrasaction="ignore",
        )

        writer.writeheader()

        for result in results:
            writer.writerow(result)


def print_summary(results):
    """
    Print a short run summary to stdout.
    """

    if not results:
        return

    total_reads = sum(r["reads"] for r in results)
    total_trimmed = sum(r["trimmed_reads"] for r in results)

    total_bases = sum(r["bases"] for r in results)
    total_trimmed_bases = sum(r["trimmed_bases"] for r in results)

    total_18S = sum(r["reads_with_18S"] for r in results)
    total_58S = sum(r["reads_with_5_8S"] for r in results)
    total_28S = sum(r["reads_with_28S"] for r in results)
    total_all = sum(r["reads_with_all_three"] for r in results)

    total_multi = sum(r["reads_multiple_18S"] for r in results)

    runtime = sum(r["runtime_sec"] for r in results)

    retained = (
        100.0 * total_trimmed_bases / total_bases
        if total_bases
        else 0.0
    )

    print("\n==============================")
    print("trim18s summary")
    print("==============================")
    print(f"Samples processed        : {len(results)}")
    print(f"Input reads             : {total_reads:,}")
    print(f"Reads with 18S          : {total_18S:,}")
    print(f"Reads with 5.8S         : {total_58S:,}")
    print(f"Reads with 28S          : {total_28S:,}")
    print(f"Reads with all three    : {total_all:,}")
    print(f"Reads with multiple 18S : {total_multi:,}")
    print(f"Trimmed reads written   : {total_trimmed:,}")
    print(f"Input bases             : {total_bases:,}")
    print(f"Trimmed bases           : {total_trimmed_bases:,}")
    print(f"Retained bases          : {retained:.2f}%")
    print(f"Aggregate runtime (s)   : {runtime:.2f}")
    print("==============================\n")
