#!/usr/bin/env python3
"""Create a sequencing-run QC dashboard from Oxford Nanopore FASTQ files."""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import json
import math
import re
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, TextIO


FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
BARCODE_NAME = re.compile(r"^(barcode\d+|unclassified|unknown)$", re.IGNORECASE)
BARCODE_IN_NAME = re.compile(r"(barcode\d+|unclassified|unknown)", re.IGNORECASE)
HEADER_BARCODE_PATTERNS = (
    re.compile(r"(?:^|\s)barcode(?:_arrangement)?[=:]([^\s]+)", re.IGNORECASE),
    re.compile(r"(?:^|\s)BC:Z:([^\s]+)", re.IGNORECASE),
)
PHRED_ERROR_PROBABILITIES = tuple(10.0 ** (-phred / 10.0) for phred in range(94))


def natural_key(value: str) -> list[object]:
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", value)]


def is_fastq(path: Path) -> bool:
    return path.name.lower().endswith(FASTQ_SUFFIXES)


def fastq_stem(path: Path) -> str:
    name = path.name
    for suffix in FASTQ_SUFFIXES:
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def open_fastq(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="ascii", errors="strict")
    return path.open("rt", encoding="ascii", errors="strict")


def read_fastq(path: Path) -> Iterable[tuple[str, str, str]]:
    """Stream strict four-line FASTQ records as (header, sequence, quality)."""
    with open_fastq(path) as handle:
        record_number = 0
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            record_number += 1
            if not (sequence and plus and quality):
                raise ValueError(f"Incomplete FASTQ record {record_number} in {path}")
            if not header.startswith("@"):
                raise ValueError(f"Malformed FASTQ header at record {record_number} in {path}")
            if not plus.startswith("+"):
                raise ValueError(f"Malformed FASTQ '+' line at record {record_number} in {path}")
            sequence = sequence.rstrip("\r\n")
            quality = quality.rstrip("\r\n")
            if len(sequence) != len(quality):
                raise ValueError(
                    f"Sequence/quality length mismatch at record {record_number} in {path}: "
                    f"{len(sequence)} bases versus {len(quality)} quality values"
                )
            yield header[1:].strip(), sequence, quality


def mean_read_qscore(quality: str) -> float:
    """Calculate ONT-style mean read Q score from average error probability."""
    if not quality:
        return 0.0
    encoded = quality.encode("ascii")
    if min(encoded) < 33 or max(encoded) > 126:
        raise ValueError("FASTQ quality contains a character outside the Phred+33 range")
    error_probability = sum(PHRED_ERROR_PROBABILITIES[code - 33] for code in encoded)
    return -10.0 * math.log10(error_probability / len(quality))


def header_barcode(header: str) -> str | None:
    for pattern in HEADER_BARCODE_PATTERNS:
        match = pattern.search(header)
        if match:
            value = match.group(1).rstrip(";,|")
            if value.lower() not in {"unclassified", "unknown", "none", "-"}:
                embedded = BARCODE_IN_NAME.search(value)
                return embedded.group(1).lower() if embedded else value
            return "unclassified"
    return None


def path_barcode(path: Path, root: Path, sample_name: str | None) -> str:
    if sample_name:
        return sample_name
    if root.is_file():
        match = BARCODE_IN_NAME.search(path.name)
        return match.group(1).lower() if match else fastq_stem(path)

    relative = path.relative_to(root)
    for part in relative.parts[:-1]:
        if BARCODE_NAME.match(part):
            return part.lower()
    match = BARCODE_IN_NAME.search(path.name)
    if match:
        return match.group(1).lower()
    if len(relative.parts) > 1:
        return relative.parts[0]
    return fastq_stem(path)


def discover_fastqs(input_path: Path, sample_name: str | None = None) -> list[tuple[Path, str]]:
    input_path = input_path.resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input does not exist: {input_path}")
    if input_path.is_file():
        if not is_fastq(input_path):
            raise ValueError(f"Input file is not FASTQ/FASTQ.GZ: {input_path}")
        return [(input_path, path_barcode(input_path, input_path, sample_name))]

    files = sorted(
        (path for path in input_path.rglob("*") if path.is_file() and is_fastq(path)),
        key=lambda path: natural_key(str(path.relative_to(input_path))),
    )
    if not files:
        raise ValueError(f"No .fastq, .fq, .fastq.gz, or .fq.gz files found below {input_path}")
    return [(path, path_barcode(path, input_path, sample_name)) for path in files]


@dataclass
class BarcodeStats:
    name: str
    files: set[str] = field(default_factory=set)
    reads: int = 0
    bases: int = 0
    qscore_sum: float = 0.0
    min_length: int | None = None
    max_length: int = 0
    length_counts: Counter[int] = field(default_factory=Counter)
    qscore_tenths: Counter[int] = field(default_factory=Counter)
    qscore_pass: Counter[int] = field(default_factory=Counter)

    def add(self, length: int, qscore: float, source: Path) -> None:
        self.files.add(str(source))
        self.reads += 1
        self.bases += length
        self.qscore_sum += qscore
        self.min_length = length if self.min_length is None else min(self.min_length, length)
        self.max_length = max(self.max_length, length)
        self.length_counts[length] += 1
        self.qscore_tenths[int(round(qscore * 10.0))] += 1
        for threshold in (7, 10, 15, 20):
            if qscore >= threshold:
                self.qscore_pass[threshold] += 1


def counter_value_at_rank(counts: Counter[int], rank: int) -> int:
    cumulative = 0
    for value in sorted(counts):
        cumulative += counts[value]
        if cumulative > rank:
            return value
    raise ValueError("Cannot select a rank from an empty counter")


def counter_quantile(counts: Counter[int], fraction: float, scale: float = 1.0) -> float:
    count = sum(counts.values())
    if count == 0:
        return 0.0
    position = (count - 1) * fraction
    lower_rank = math.floor(position)
    upper_rank = math.ceil(position)
    lower = counter_value_at_rank(counts, lower_rank)
    upper = counter_value_at_rank(counts, upper_rank)
    return (lower + (upper - lower) * (position - lower_rank)) / scale


def read_n50(length_counts: Counter[int], total_bases: int) -> int:
    if total_bases == 0:
        return 0
    cumulative = 0
    threshold = total_bases / 2.0
    for length in sorted(length_counts, reverse=True):
        cumulative += length * length_counts[length]
        if cumulative >= threshold:
            return length
    return 0


def binned_distribution(
    counts: Counter[int], bin_width: float, value_scale: float = 1.0
) -> tuple[list[float], list[int]]:
    bins: Counter[int] = Counter()
    scaled_width = max(1, int(round(bin_width * value_scale)))
    for value, count in counts.items():
        bins[(value // scaled_width) * scaled_width] += count
    keys = sorted(bins)
    centers = [(key + scaled_width / 2.0) / value_scale for key in keys]
    return centers, [bins[key] for key in keys]


def metrics_for(stats: BarcodeStats, length_bin_width: int, quality_bin_width: float) -> dict[str, object]:
    length_x, length_y = binned_distribution(stats.length_counts, length_bin_width)
    quality_x, quality_y = binned_distribution(stats.qscore_tenths, quality_bin_width, 10.0)
    return {
        "barcode": stats.name,
        "files": sorted(stats.files),
        "file_count": len(stats.files),
        "reads": stats.reads,
        "bases": stats.bases,
        "mean_length": stats.bases / stats.reads if stats.reads else 0.0,
        "median_length": counter_quantile(stats.length_counts, 0.5),
        "length_q10": counter_quantile(stats.length_counts, 0.1),
        "length_q90": counter_quantile(stats.length_counts, 0.9),
        "n50": read_n50(stats.length_counts, stats.bases),
        "min_length": stats.min_length or 0,
        "max_length": stats.max_length,
        "mean_qscore": stats.qscore_sum / stats.reads if stats.reads else 0.0,
        "median_qscore": counter_quantile(stats.qscore_tenths, 0.5, 10.0),
        "q7_percent": 100.0 * stats.qscore_pass[7] / stats.reads if stats.reads else 0.0,
        "q10_percent": 100.0 * stats.qscore_pass[10] / stats.reads if stats.reads else 0.0,
        "q15_percent": 100.0 * stats.qscore_pass[15] / stats.reads if stats.reads else 0.0,
        "q20_percent": 100.0 * stats.qscore_pass[20] / stats.reads if stats.reads else 0.0,
        "length_distribution": {"x": length_x, "counts": length_y},
        "quality_distribution": {"x": quality_x, "counts": quality_y},
    }


def combine_stats(stats_by_barcode: dict[str, BarcodeStats]) -> BarcodeStats:
    combined = BarcodeStats("all barcodes")
    for stats in stats_by_barcode.values():
        combined.files.update(stats.files)
        combined.reads += stats.reads
        combined.bases += stats.bases
        combined.qscore_sum += stats.qscore_sum
        if stats.min_length is not None:
            combined.min_length = (
                stats.min_length
                if combined.min_length is None
                else min(combined.min_length, stats.min_length)
            )
        combined.max_length = max(combined.max_length, stats.max_length)
        combined.length_counts.update(stats.length_counts)
        combined.qscore_tenths.update(stats.qscore_tenths)
        combined.qscore_pass.update(stats.qscore_pass)
    return combined


def analyse_fastqs(
    inputs: list[tuple[Path, str]],
    barcode_source: str,
    max_reads_per_barcode: int,
    progress_every: int,
    per_read_handle: TextIO | None = None,
) -> dict[str, BarcodeStats]:
    stats_by_barcode: dict[str, BarcodeStats] = {}
    per_read_writer = (
        csv.writer(per_read_handle, delimiter="\t", lineterminator="\n")
        if per_read_handle is not None
        else None
    )
    if per_read_writer:
        per_read_writer.writerow(
            ("read_id", "barcode", "read_length", "mean_qscore", "source_file")
        )
    total_reads = 0
    for file_number, (path, fallback_barcode) in enumerate(inputs, start=1):
        print(f"[{file_number}/{len(inputs)}] Reading {path}", file=sys.stderr)
        for header, sequence, quality in read_fastq(path):
            detected = header_barcode(header) if barcode_source in {"auto", "header"} else None
            if barcode_source == "header" and not detected:
                detected = "unclassified"
            barcode = detected or fallback_barcode
            stats = stats_by_barcode.setdefault(barcode, BarcodeStats(barcode))
            if max_reads_per_barcode and stats.reads >= max_reads_per_barcode:
                continue
            try:
                qscore = mean_read_qscore(quality)
            except ValueError as error:
                raise ValueError(f"{error} in read {header.split(maxsplit=1)[0]} from {path}") from error
            stats.add(len(sequence), qscore, path)
            if per_read_writer:
                per_read_writer.writerow(
                    (header.split(maxsplit=1)[0], barcode, len(sequence), f"{qscore:.6f}", str(path))
                )
            total_reads += 1
            if progress_every and total_reads % progress_every == 0:
                print(f"  Processed {total_reads:,} reads", file=sys.stderr)
    if not stats_by_barcode or not any(stats.reads for stats in stats_by_barcode.values()):
        raise ValueError("No FASTQ reads were found")
    return stats_by_barcode


def write_summary_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    columns = (
        "barcode", "file_count", "reads", "bases", "mean_length", "median_length",
        "length_q10", "length_q90", "n50", "min_length", "max_length", "mean_qscore",
        "median_qscore", "q7_percent", "q10_percent", "q15_percent", "q20_percent",
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            values: list[str] = []
            for column in columns:
                value = row[column]
                values.append(f"{value:.3f}" if isinstance(value, float) else str(value))
            handle.write("\t".join(values) + "\n")


def dashboard_html(rows: list[dict[str, object]], overall: dict[str, object], input_path: Path) -> str:
    data_json = json.dumps(rows, separators=(",", ":")).replace("</", "<\\/")
    overall_json = json.dumps(overall, separators=(",", ":")).replace("</", "<\\/")
    source = html.escape(str(input_path.resolve()))
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Nanopore sequencing run QC</title><script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
:root{{--ink:#18302d;--muted:#66736c;--teal:#12685f;--gold:#c58a14;--paper:#f4f6f3;--card:#fff;--line:#dce4df}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--paper);color:var(--ink);font-family:Arial,sans-serif}}
header{{background:linear-gradient(120deg,#0c4d47,#168176);color:#fff;padding:30px 5vw}}header h1{{margin:0 0 7px;font-size:30px}}header p{{margin:0;opacity:.88;word-break:break-all}}
main{{max-width:1550px;margin:auto;padding:24px}}.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(175px,1fr));gap:14px}}
.card,.panel{{background:var(--card);border:1px solid var(--line);border-radius:12px;box-shadow:0 3px 12px #18302d12}}.card{{padding:17px}}.card span{{color:var(--muted);font-size:13px}}.card b{{display:block;color:var(--teal);font-size:27px;margin-top:6px}}
.grid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:18px}}.panel{{padding:18px;margin-top:18px}}.panel h2{{font-size:19px;margin:0 0 4px}}.hint{{color:var(--muted);font-size:13px;margin:0 0 8px}}
.table-wrap{{overflow:auto;max-height:620px}}table{{border-collapse:collapse;width:100%;font-size:12px}}th{{position:sticky;top:0;background:#e8efec;text-align:left}}th,td{{padding:8px;border-bottom:1px solid var(--line);white-space:nowrap}}tbody tr:hover{{background:#f1f7f4}}
@media(max-width:900px){{.grid{{grid-template-columns:1fr}}}}
</style></head><body>
<header><h1>Nanopore sequencing run QC</h1><p>Raw FASTQ source: {source}</p></header><main>
<section class="cards"><div class="card"><span>Barcodes</span><b id="barcodeCount"></b></div><div class="card"><span>Total reads</span><b id="totalReads"></b></div><div class="card"><span>Total yield</span><b id="totalBases"></b></div><div class="card"><span>Median read length</span><b id="medianLength"></b></div><div class="card"><span>Read N50</span><b id="readN50"></b></div><div class="card"><span>Median read quality</span><b id="medianQ"></b></div></section>
<section class="grid"><div class="panel"><h2>Read count by barcode</h2><p class="hint">Click legend entries in distribution plots to hide or isolate barcodes.</p><div id="readCounts" style="height:430px"></div></div><div class="panel"><h2>Yield by barcode</h2><p class="hint">Total sequenced bases per barcode.</p><div id="yieldChart" style="height:430px"></div></div></section>
<section class="grid"><div class="panel"><h2>Read-length summary</h2><p class="hint">Median length and length-weighted N50.</p><div id="lengthSummary" style="height:430px"></div></div><div class="panel"><h2>Read-quality summary</h2><p class="hint">Median per-read Q score and percentage of reads at or above Q15.</p><div id="qualitySummary" style="height:430px"></div></div></section>
<section class="panel"><h2>Read-length distributions</h2><p class="hint">Each barcode is normalized to percent of its reads, allowing direct comparison despite unequal yields.</p><div id="lengthDistribution" style="height:520px"></div></section>
<section class="panel"><h2>Read-quality distributions</h2><p class="hint">Per-read mean Q scores calculated from the FASTQ Phred+33 quality string.</p><div id="qualityDistribution" style="height:520px"></div></section>
<section class="panel"><h2>Barcode QC table</h2><div class="table-wrap"><table><thead><tr><th>Barcode</th><th>FASTQ files</th><th>Reads</th><th>Yield</th><th>Mean length</th><th>Median length</th><th>Length P10–P90</th><th>N50</th><th>Mean Q</th><th>Median Q</th><th>≥Q7</th><th>≥Q10</th><th>≥Q15</th><th>≥Q20</th></tr></thead><tbody id="summaryRows"></tbody></table></div></section>
</main><script>
const samples={data_json},overall={overall_json};
const fmt=n=>Number(n).toLocaleString(),bases=n=>n>=1e9?(n/1e9).toFixed(2)+' Gb':n>=1e6?(n/1e6).toFixed(2)+' Mb':n>=1e3?(n/1e3).toFixed(2)+' kb':fmt(n),pct=n=>n.toFixed(1)+'%';
document.getElementById('barcodeCount').textContent=samples.length;document.getElementById('totalReads').textContent=fmt(overall.reads);document.getElementById('totalBases').textContent=bases(overall.bases);document.getElementById('medianLength').textContent=fmt(Math.round(overall.median_length))+' bp';document.getElementById('readN50').textContent=fmt(overall.n50)+' bp';document.getElementById('medianQ').textContent='Q'+overall.median_qscore.toFixed(1);
const labels=samples.map(d=>d.barcode),config={{responsive:true,displaylogo:false}},layout={{paper_bgcolor:'white',plot_bgcolor:'white',font:{{color:'#18302d'}},margin:{{t:20,r:25,b:100,l:70}},xaxis:{{tickangle:-40}},yaxis:{{gridcolor:'#e5ebe7'}}}};
Plotly.newPlot('readCounts',[{{x:labels,y:samples.map(d=>d.reads),type:'bar',marker:{{color:'#12685f'}},hovertemplate:'%{{x}}<br>%{{y:,}} reads<extra></extra>'}}],{{...layout,yaxis:{{...layout.yaxis,title:'Reads'}}}},config);
Plotly.newPlot('yieldChart',[{{x:labels,y:samples.map(d=>d.bases/1e6),type:'bar',marker:{{color:'#c58a14'}},hovertemplate:'%{{x}}<br>%{{y:,.2f}} Mb<extra></extra>'}}],{{...layout,yaxis:{{...layout.yaxis,title:'Yield (Mb)'}}}},config);
Plotly.newPlot('lengthSummary',[{{x:labels,y:samples.map(d=>d.median_length),name:'Median',type:'bar'}},{{x:labels,y:samples.map(d=>d.n50),name:'N50',type:'bar'}}],{{...layout,barmode:'group',yaxis:{{...layout.yaxis,title:'Length (bp)'}}}},config);
Plotly.newPlot('qualitySummary',[{{x:labels,y:samples.map(d=>d.median_qscore),name:'Median Q',type:'bar',marker:{{color:'#12685f'}}}},{{x:labels,y:samples.map(d=>d.q15_percent),name:'Reads ≥ Q15 (%)',type:'scatter',mode:'lines+markers',yaxis:'y2',line:{{color:'#c58a14'}}}}],{{...layout,yaxis:{{...layout.yaxis,title:'Median Q score'}},yaxis2:{{title:'Reads ≥ Q15 (%)',overlaying:'y',side:'right',range:[0,105],showgrid:false}},legend:{{orientation:'h',y:1.12}}}},config);
const traces=(field,xTitle)=>samples.map(d=>{{const dist=d[field],total=dist.counts.reduce((a,b)=>a+b,0);return{{x:dist.x,y:dist.counts.map(n=>100*n/total),name:d.barcode,type:'scatter',mode:'lines',hovertemplate:`${{d.barcode}}<br>%{{x}}<br>%{{y:.2f}}% of reads<extra></extra>`}}}});
Plotly.newPlot('lengthDistribution',traces('length_distribution'),{{...layout,margin:{{t:25,r:25,b:65,l:70}},xaxis:{{title:'Read length (bp)'}},yaxis:{{...layout.yaxis,title:'Reads (%)'}},legend:{{orientation:'h',y:1.16}}}},config);
Plotly.newPlot('qualityDistribution',traces('quality_distribution'),{{...layout,margin:{{t:25,r:25,b:65,l:70}},xaxis:{{title:'Mean read Q score'}},yaxis:{{...layout.yaxis,title:'Reads (%)'}},legend:{{orientation:'h',y:1.16}}}},config);
document.getElementById('summaryRows').innerHTML=samples.map(d=>`<tr><td><b>${{d.barcode}}</b></td><td>${{d.file_count}}</td><td>${{fmt(d.reads)}}</td><td>${{bases(d.bases)}}</td><td>${{fmt(Math.round(d.mean_length))}} bp</td><td>${{fmt(Math.round(d.median_length))}} bp</td><td>${{fmt(Math.round(d.length_q10))}}–${{fmt(Math.round(d.length_q90))}} bp</td><td>${{fmt(d.n50)}} bp</td><td>Q${{d.mean_qscore.toFixed(1)}}</td><td>Q${{d.median_qscore.toFixed(1)}}</td><td>${{pct(d.q7_percent)}}</td><td>${{pct(d.q10_percent)}}</td><td>${{pct(d.q15_percent)}}</td><td>${{pct(d.q20_percent)}}</td></tr>`).join('');
</script></body></html>"""


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Summarize read lengths and qualities across a Nanopore FASTQ run."
    )
    parser.add_argument("input", type=Path, help="FASTQ/FASTQ.GZ file or directory searched recursively")
    parser.add_argument("-o", "--output-dir", type=Path, default=Path("run_qc_results"))
    parser.add_argument("--sample-name", help="Fallback sample name (mainly useful for one input file)")
    parser.add_argument(
        "--barcode-source", choices=("auto", "path", "header"), default="auto",
        help="Barcode assignment method (default: path, overridden by recognized header tags)",
    )
    parser.add_argument("--length-bin-width", type=int, default=100, help="Length histogram bin width in bp")
    parser.add_argument("--quality-bin-width", type=float, default=0.5, help="Q-score histogram bin width")
    parser.add_argument(
        "--max-reads-per-barcode", type=int, default=0,
        help="Optional quick-QC limit per barcode; 0 analyses all reads",
    )
    parser.add_argument("--progress-every", type=int, default=100000, help="Progress interval in reads; 0 disables")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.length_bin_width < 1:
        parser.error("--length-bin-width must be at least 1")
    if args.quality_bin_width <= 0:
        parser.error("--quality-bin-width must be greater than 0")
    if args.max_reads_per_barcode < 0:
        parser.error("--max-reads-per-barcode cannot be negative")
    if args.progress_every < 0:
        parser.error("--progress-every cannot be negative")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    per_read_path = args.output_dir / "run_qc_per_read.tsv"
    per_read_temporary_path = args.output_dir / ".run_qc_per_read.tsv.tmp"
    try:
        inputs = discover_fastqs(args.input, args.sample_name)
        print(f"Discovered {len(inputs)} FASTQ file(s)", file=sys.stderr)
        with per_read_temporary_path.open("w", encoding="utf-8", newline="") as per_read_handle:
            stats_by_barcode = analyse_fastqs(
                inputs,
                args.barcode_source,
                args.max_reads_per_barcode,
                args.progress_every,
                per_read_handle,
            )
        per_read_temporary_path.replace(per_read_path)
    except (OSError, UnicodeError, ValueError) as error:
        per_read_temporary_path.unlink(missing_ok=True)
        parser.exit(2, f"error: {error}\n")

    ordered = sorted(stats_by_barcode.values(), key=lambda stats: natural_key(stats.name))
    rows = [metrics_for(stats, args.length_bin_width, args.quality_bin_width) for stats in ordered]
    overall_stats = combine_stats(stats_by_barcode)
    overall = metrics_for(overall_stats, args.length_bin_width, args.quality_bin_width)

    summary_path = args.output_dir / "run_qc_summary.tsv"
    metrics_path = args.output_dir / "run_qc_metrics.json"
    dashboard_path = args.output_dir / "run_qc_dashboard.html"
    write_summary_tsv(summary_path, rows)
    with metrics_path.open("w", encoding="utf-8") as handle:
        json.dump({"input": str(args.input.resolve()), "overall": overall, "barcodes": rows}, handle, indent=2)
        handle.write("\n")
    dashboard_path.write_text(dashboard_html(rows, overall, args.input), encoding="utf-8")

    print(f"Analysed {overall['reads']:,} reads across {len(rows)} barcode(s)", file=sys.stderr)
    print(f"Dashboard: {dashboard_path.resolve()}")
    print(f"Summary:   {summary_path.resolve()}")
    print(f"Per read:  {per_read_path.resolve()}")
    print(f"Metrics:   {metrics_path.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
