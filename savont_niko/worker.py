#!/usr/bin/env python3
"""
worker.py

Process a single FASTQ file:

    FASTQ
        │
        ▼
    Barrnap
        │
        ▼
    Extract longest 18S
        │
        ▼
    Write *.18S.fastq.gz
        │
        ▼
    Return per-sample statistics
"""

from __future__ import annotations

import logging
import time
from pathlib import Path
from statistics import median

from barrnap import longest_feature, run_barrnap
from fastq import read_fastq, write_fastq


COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def n50(lengths):
    """Calculate N50."""

    if not lengths:
        return 0

    lengths = sorted(lengths, reverse=True)

    total = sum(lengths)
    running = 0

    for length in lengths:
        running += length
        if running >= total / 2:
            return length

    return 0


def process_sample(
    fastq,
    outdir,
    threads=8,
    kingdom="euk",
    overwrite=False,
):
    """
    Process one FASTQ.

    Returns
    -------
    dict
        Per-sample summary statistics.
    """

    fastq = Path(fastq)
    outdir = Path(outdir)

    sample = fastq.name
    stem = fastq.name

    # remove common suffixes
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break

    outfile = outdir / f"{stem}.18S.fastq.gz"

    if outfile.exists() and not overwrite:
        raise FileExistsError(
            f"{outfile} already exists. "
            "Use --overwrite to replace it."
        )

    start_time = time.time()

    annotations, feature_counts = run_barrnap(
        fastq=fastq,
        threads=threads,
        kingdom=kingdom,
    )

    stats = {
        "sample": sample,

        # input
        "reads": 0,
        "bases": 0,

        # output
        "trimmed_reads": 0,
        "trimmed_bases": 0,

        # annotation
        "reads_with_18S": 0,
        "reads_with_5_8S": 0,
        "reads_with_28S": 0,
        "reads_with_all_three": 0,
        "reads_multiple_18S": 0,

        "18S_features": feature_counts["18S_rRNA"],
        "5.8S_features": feature_counts["5.8S_rRNA"],
        "28S_features": feature_counts["28S_rRNA"],

        # lengths
        "mean_input": 0,
        "median_input": 0,
        "mean_trimmed": 0,
        "median_trimmed": 0,
        "min_trimmed": 0,
        "max_trimmed": 0,
        "n50_trimmed": 0,

        # retained
        "retained_pct": 0,

        # runtime
        "runtime_sec": 0,
    }

    input_lengths = []
    trimmed_lengths = []

    def trimmed_records():

        for name, seq, qual in read_fastq(fastq):

            readlen = len(seq)

            stats["reads"] += 1
            stats["bases"] += readlen

            input_lengths.append(readlen)

            if name not in annotations:
                continue

            feats = annotations[name]

            has18 = "18S_rRNA" in feats
            has58 = "5.8S_rRNA" in feats
            has28 = "28S_rRNA" in feats

            if has18:
                stats["reads_with_18S"] += 1

            if has58:
                stats["reads_with_5_8S"] += 1

            if has28:
                stats["reads_with_28S"] += 1

            if has18 and has58 and has28:
                stats["reads_with_all_three"] += 1

            #
            # Only trim reads containing 18S.
            #
            if not has18:
                continue

            hits = feats["18S_rRNA"]

            if len(hits) > 1:
                stats["reads_multiple_18S"] += 1

            start, end, strand = longest_feature(hits)

            #
            # Safety checks
            #
            if start < 0:
                continue

            if end > readlen:
                continue

            if start >= end:
                continue

            trimmed_seq = seq[start:end]
            trimmed_qual = qual[start:end]

            if strand == "-":
                trimmed_seq = reverse_complement(trimmed_seq)
                trimmed_qual = trimmed_qual[::-1]

            if len(trimmed_seq) != len(trimmed_qual):
                logging.warning(
                    "Quality mismatch for %s",
                    name,
                )
                continue

            stats["trimmed_reads"] += 1
            stats["trimmed_bases"] += len(trimmed_seq)

            trimmed_lengths.append(len(trimmed_seq))

            yield (
                name,
                trimmed_seq,
                trimmed_qual,
            )

    #
    # Stream output directly to compressed FASTQ.
    #
    write_fastq(outfile, trimmed_records())

    #
    # Final statistics
    #
    if input_lengths:
        stats["mean_input"] = sum(input_lengths) / len(input_lengths)
        stats["median_input"] = median(input_lengths)

    if trimmed_lengths:

        stats["mean_trimmed"] = (
            sum(trimmed_lengths) / len(trimmed_lengths)
        )

        stats["median_trimmed"] = median(trimmed_lengths)

        stats["min_trimmed"] = min(trimmed_lengths)
        stats["max_trimmed"] = max(trimmed_lengths)

        stats["n50_trimmed"] = n50(trimmed_lengths)

    if stats["bases"]:

        stats["retained_pct"] = (
            100.0
            * stats["trimmed_bases"]
            / stats["bases"]
        )

    stats["runtime_sec"] = round(
        time.time() - start_time,
        2,
    )

    logging.info(
        "%s: kept %d / %d reads",
        sample,
        stats["trimmed_reads"],
        stats["reads"],
    )

    return stats
