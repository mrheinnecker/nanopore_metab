#!/usr/bin/env python3
"""
fastq.py

Simple streaming FASTQ reader/writer with transparent support for
plain text and gzip-compressed files.

Dependencies:
    pip install xopen
"""

from __future__ import annotations

from pathlib import Path
from typing import Generator, Iterable

from xopen import xopen


def read_fastq(
    filename: str | Path,
) -> Generator[tuple[str, str, str], None, None]:
    """
    Stream FASTQ records.

    Parameters
    ----------
    filename
        FASTQ or FASTQ.gz

    Yields
    ------
    (name, sequence, quality)

    Notes
    -----
    The read name is truncated at the first whitespace,
    matching Barrnap's sequence IDs.
    """

    filename = Path(filename)

    with xopen(filename, "rt") as fh:

        while True:

            header = fh.readline()

            if not header:
                break

            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()

            if not (seq and plus and qual):
                raise ValueError(
                    f"Incomplete FASTQ record in {filename}"
                )

            if not header.startswith("@"):
                raise ValueError(
                    f"Malformed FASTQ header:\n{header}"
                )

            if not plus.startswith("+"):
                raise ValueError(
                    f"Malformed FASTQ '+' line for read "
                    f"{header.strip()}"
                )

            seq = seq.rstrip("\r\n")
            qual = qual.rstrip("\r\n")

            if len(seq) != len(qual):
                raise ValueError(
                    f"Sequence/quality length mismatch "
                    f"for read {header.strip()}"
                )

            name = header[1:].split(maxsplit=1)[0]

            yield (
                name,
                seq,
                qual,
            )


def write_fastq(
    filename: str | Path,
    records: Iterable[tuple[str, str, str]],
):
    """
    Write FASTQ records.

    Parameters
    ----------
    filename
        Output FASTQ(.gz)

    records
        Iterable yielding

            (name, sequence, quality)
    """

    filename = Path(filename)

    with xopen(filename, "wt", compresslevel=6) as out:

        for name, seq, qual in records:

            if len(seq) != len(qual):
                raise ValueError(
                    f"Sequence/quality length mismatch "
                    f"for read {name}"
                )

            out.write("@")
            out.write(name)
            out.write("\n")

            out.write(seq)
            out.write("\n+\n")

            out.write(qual)
            out.write("\n")


def count_reads(filename: str | Path) -> int:
    """
    Count FASTQ records.

    Mainly useful for sanity checks.
    """

    n = 0

    for _ in read_fastq(filename):
        n += 1

    return n
