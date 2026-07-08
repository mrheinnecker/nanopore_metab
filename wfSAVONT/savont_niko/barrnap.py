#!/usr/bin/env python3
"""
barrnap.py

Run Barrnap and stream its GFF3 output directly into memory.

Returns
-------
annotations : dict

{
    read_id: {
        "18S_rRNA": [(start, end), ...],
        "5.8S_rRNA": [(start, end), ...],
        "28S_rRNA": [(start, end), ...],
    }
}

Coordinates follow normal Python slicing:
    start = inclusive (0-based)
    end   = exclusive
"""

from __future__ import annotations

import logging
import os
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

from fastq import read_fastq


RRNA_FEATURES = {
    "18S_rRNA",
    "5.8S_rRNA",
    "28S_rRNA",
}


def parse_gff_attributes(attributes: str) -> dict[str, str]:
    parsed = {}
    for item in attributes.split(";"):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        parsed[key] = value
    return parsed


def normalize_feature_name(feature: str, attributes: str) -> str | None:
    if feature in RRNA_FEATURES:
        return feature

    name = parse_gff_attributes(attributes).get("Name")

    if name == "5_8S_rRNA":
        name = "5.8S_rRNA"

    if name in RRNA_FEATURES:
        return name

    return None


def run_barrnap(
    fastq: str | Path,
    threads: int = 8,
    kingdom: str = "euk",
):
    """
    Run Barrnap on a FASTQ file.

    Parameters
    ----------
    fastq
        Input FASTQ(.gz)

    threads
        Threads passed to Barrnap.

    kingdom
        Barrnap kingdom.

    Returns
    -------
    annotations : dict

        Nested dictionary of feature coordinates.

    feature_counts : dict

        Total number of features reported by Barrnap.

    Raises
    ------
    RuntimeError
        If Barrnap exits with a non-zero status.
    """

    fastq = Path(fastq)

    with tempfile.NamedTemporaryFile(
        "w",
        suffix=".fasta",
        prefix=f"{fastq.stem}.barrnap.",
        delete=False,
    ) as fasta_fh:
        fasta_path = fasta_fh.name
        for name, seq, _qual in read_fastq(fastq):
            fasta_fh.write(f">{name}\n{seq}\n")

    cmd = [
        "barrnap",
        "--threads",
        str(threads),
        "--kingdom",
        kingdom,
        fasta_path,
    ]

    logging.info("Running Barrnap on %s", fastq)

    stderr_path = None

    try:
        stderr_fh = tempfile.NamedTemporaryFile(
            "w+",
            suffix=".stderr",
            prefix=f"{fastq.stem}.barrnap.",
            delete=False,
        )
        stderr_path = stderr_fh.name

        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=stderr_fh,
            text=True,
            bufsize=1,
        )

        annotations = defaultdict(lambda: defaultdict(list))

        feature_counts = {
            "18S_rRNA": 0,
            "5.8S_rRNA": 0,
            "28S_rRNA": 0,
        }

        assert proc.stdout is not None

        #
        # Stream GFF directly from stdout.
        #
        for line in proc.stdout:

            if not line or line.startswith("#"):
                continue

            cols = line.rstrip().split("\t")

            if len(cols) != 9:
                continue

            seqid = cols[0]
            feature = normalize_feature_name(cols[2], cols[8])

            if feature is None:
                continue

            #
            # Convert GFF coordinates
            # (1-based inclusive)
            # ->
            # Python slice
            # (0-based, end exclusive)
            #
            start = int(cols[3]) - 1
            end = int(cols[4])
            strand = cols[6]

            annotations[seqid][feature].append((start, end, strand))
            feature_counts[feature] += 1

        proc.wait()
        stderr_fh.close()

        with open(stderr_path) as stderr_read:
            stderr = stderr_read.read()

        if proc.returncode != 0:

            raise RuntimeError(
                f"Barrnap failed for '{fastq}'\n\n{stderr}"
            )

        if stderr.strip():
            logging.debug(stderr.strip())

        logging.info(
            "Barrnap identified %d annotated reads.",
            len(annotations),
        )

        return annotations, feature_counts
    finally:
        os.unlink(fasta_path)
        if stderr_path is not None:
            os.unlink(stderr_path)


def longest_feature(features):
    """
    Return the longest feature from a list of
    (start, end) tuples.

    Parameters
    ----------
    features : list[(start, end)]

    Returns
    -------
    (start, end)
    """

    if not features:
        return None

    return max(
        features,
        key=lambda x: x[1] - x[0],
    )
