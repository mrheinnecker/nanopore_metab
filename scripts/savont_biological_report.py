#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import html
import json
import math
from collections import defaultdict
from pathlib import Path


RANKS = [
    "Domain",
    "Supergroup",
    "Division",
    "Subdivision",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species",
]


def fmt_int(value):
    try:
        return f"{int(float(value)):,}"
    except (TypeError, ValueError):
        return "0"


def fmt_float(value, digits=2):
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return "0"


def sample_from_trim(path: Path):
    name = path.name
    name = name.replace("trim18s.summary_", "")
    return name.rsplit(".", 1)[0]


def sample_from_asv_dir(path: Path):
    return path.name.replace("savont_", "")


def sample_from_taxonomy(path: Path):
    if path.is_dir():
        return path.name.replace("pr2_", "")
    parent = path.parent.name
    return parent.replace("pr2_", "")


def read_tsv(path: Path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_trim_summaries(paths):
    records = {}
    for raw in paths:
        path = Path(raw)
        rows = read_tsv(path)
        if not rows:
            continue
        sample = sample_from_trim(path)
        records[sample] = rows[0]
    return records


def load_feature_table(path: Path):
    if not path.exists() or path.stat().st_size == 0:
        return {}
    with path.open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, [])
        if len(header) < 2:
            return {}
        counts = {}
        for row in reader:
            if len(row) < 2:
                continue
            try:
                counts[row[0]] = int(float(row[1]))
            except ValueError:
                counts[row[0]] = 0
        return counts


def count_fasta(path: Path):
    if not path.exists():
        return 0
    total = 0
    with path.open(encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                total += 1
    return total


def load_asv_dirs(paths):
    records = {}
    for raw in paths:
        path = Path(raw)
        sample = sample_from_asv_dir(path)
        counts = load_feature_table(path / "feature-table.tsv")
        records[sample] = {
            "dir": path,
            "counts": counts,
            "asv_count": count_fasta(path / "final_asvs.fasta") or len(counts),
        }
    return records


def load_taxonomy_tables(paths):
    records = {}
    for raw in paths:
        path = Path(raw)
        if path.is_dir():
            path = path / "taxonomy_table.tsv"
        sample = sample_from_taxonomy(path)
        rows = read_tsv(path) if path.exists() and path.stat().st_size else []
        records[sample] = rows
    return records


def diversity(counts):
    total = sum(counts)
    if total <= 0:
        return 0.0, 0.0
    proportions = [count / total for count in counts if count > 0]
    shannon = -sum(p * math.log(p) for p in proportions)
    simpson = 1 - sum(p * p for p in proportions)
    return shannon, simpson


def taxonomy_with_counts(tax_rows, counts):
    merged = []
    for row in tax_rows:
        otu = row.get("OTU", "")
        read_count = counts.get(otu, 0)
        merged.append({**row, "reads": read_count})
    return merged


def aggregate_rank(rows, rank):
    totals = defaultdict(int)
    for row in rows:
        taxon = row.get(rank) or "Unassigned"
        totals[taxon] += int(row.get("reads", 0))
    return sorted(totals.items(), key=lambda item: item[1], reverse=True)


def top_otus(rows, limit=12):
    return sorted(rows, key=lambda row: int(row.get("reads", 0)), reverse=True)[:limit]


def sample_metrics(sample, trim, asv, tax_rows):
    counts = list(asv["counts"].values())
    assigned_rows = taxonomy_with_counts(tax_rows, asv["counts"])
    total_asv_reads = sum(counts)
    shannon, simpson = diversity(counts)
    species = aggregate_rank(assigned_rows, "Species")
    top_species = species[0][0] if species else "Unassigned"
    top_species_reads = species[0][1] if species else 0

    reads = float(trim.get("reads", 0) or 0)
    reads18 = float(trim.get("reads_with_18S", 0) or 0)
    reads58 = float(trim.get("reads_with_5_8S", 0) or 0)
    reads28 = float(trim.get("reads_with_28S", 0) or 0)

    return {
        "sample": sample,
        "input_reads": int(reads),
        "input_bases": int(float(trim.get("bases", 0) or 0)),
        "mean_input": float(trim.get("mean_input", 0) or 0),
        "trimmed_reads": int(float(trim.get("trimmed_reads", 0) or 0)),
        "trimmed_bases": int(float(trim.get("trimmed_bases", 0) or 0)),
        "retained_pct": float(trim.get("retained_pct", 0) or 0),
        "mean_trimmed": float(trim.get("mean_trimmed", 0) or 0),
        "median_trimmed": float(trim.get("median_trimmed", 0) or 0),
        "n50_trimmed": float(trim.get("n50_trimmed", 0) or 0),
        "reads_with_18S": int(reads18),
        "reads_with_5_8S": int(reads58),
        "reads_with_28S": int(reads28),
        "rrna18_pct": reads18 / reads * 100 if reads else 0,
        "rrna58_pct": reads58 / reads * 100 if reads else 0,
        "rrna28_pct": reads28 / reads * 100 if reads else 0,
        "asv_count": int(asv["asv_count"]),
        "asv_reads": total_asv_reads,
        "shannon": shannon,
        "simpson": simpson,
        "classified_asvs": sum(1 for row in assigned_rows if row.get("Species")),
        "top_species": top_species,
        "top_species_reads": top_species_reads,
    }


def html_page(title, body, scripts=""):
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
:root {{ --bg:#f6f7f4; --surface:#ffffff; --text:#1f2522; --muted:#66736c; --line:#dce3de; --primary:#12685f; --accent:#c58a14; }}
body {{ margin:0; font-family:Inter,system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; background:var(--bg); color:var(--text); }}
.wrap {{ max-width:1180px; margin:0 auto; padding:24px; }}
header {{ display:flex; justify-content:space-between; gap:16px; align-items:flex-start; margin-bottom:18px; }}
h1 {{ margin:0; font-size:2rem; }}
h2 {{ margin:0 0 12px 0; font-size:1.15rem; }}
.sub {{ color:var(--muted); margin-top:6px; }}
.grid {{ display:grid; grid-template-columns:repeat(12,1fr); gap:16px; }}
.card {{ grid-column:span 12; background:var(--surface); border:1px solid var(--line); border-radius:8px; padding:16px; }}
.half {{ grid-column:span 6; }}
.kpis {{ display:grid; grid-template-columns:repeat(4,1fr); gap:12px; margin-bottom:16px; }}
.kpi {{ background:var(--surface); border:1px solid var(--line); border-radius:8px; padding:14px; }}
.kpi .label {{ color:var(--muted); font-size:.9rem; }}
.kpi .value {{ font-size:1.45rem; font-weight:700; margin-top:4px; }}
table {{ width:100%; border-collapse:collapse; }}
th,td {{ padding:8px 7px; border-bottom:1px solid var(--line); text-align:left; font-size:.92rem; vertical-align:top; }}
th {{ color:var(--muted); font-weight:600; }}
a {{ color:var(--primary); text-decoration:none; }}
a:hover {{ text-decoration:underline; }}
.small {{ color:var(--muted); font-size:.9rem; }}
@media (max-width:850px) {{ .kpis {{ grid-template-columns:repeat(2,1fr); }} .half {{ grid-column:span 12; }} }}
</style>
</head>
<body>
<div class="wrap">
{body}
</div>
{scripts}
</body>
</html>
"""


def render_sample(sample, metrics, rows, out_path):
    top_species = aggregate_rank(rows, "Species")[:12]
    top_genera = aggregate_rank(rows, "Genus")[:12]
    otu_rows = top_otus(rows)
    total = sum(row.get("reads", 0) for row in rows)

    kpis = f"""
<div class="kpis">
<div class="kpi"><div class="label">Input reads</div><div class="value">{fmt_int(metrics['input_reads'])}</div></div>
<div class="kpi"><div class="label">Trimmed 18S reads</div><div class="value">{fmt_int(metrics['trimmed_reads'])}</div></div>
<div class="kpi"><div class="label">ASVs</div><div class="value">{fmt_int(metrics['asv_count'])}</div></div>
<div class="kpi"><div class="label">Shannon</div><div class="value">{fmt_float(metrics['shannon'])}</div></div>
</div>
"""
    otu_html = "\n".join(
        "<tr>"
        f"<td>{html.escape(row.get('OTU', ''))}</td>"
        f"<td>{fmt_int(row.get('reads', 0))}</td>"
        f"<td>{fmt_float(row.get('Pident', 0), 1)}</td>"
        f"<td>{html.escape(row.get('Species') or 'Unassigned')}</td>"
        f"<td>{html.escape(row.get('Genus') or 'Unassigned')}</td>"
        "</tr>"
        for row in otu_rows
    )
    body = f"""
<header>
<div><h1>{html.escape(sample)}</h1><div class="sub">Savont ASV and PR2 annotation report</div></div>
<div class="small"><a href="multi_sample_dashboard.html">Multi-sample dashboard</a></div>
</header>
{kpis}
<div class="grid">
<section class="card half"><h2>Read Processing</h2><div id="readChart" style="height:320px"></div></section>
<section class="card half"><h2>Top Species</h2><div id="speciesChart" style="height:320px"></div></section>
<section class="card half"><h2>Top Genera</h2><div id="genusChart" style="height:320px"></div></section>
<section class="card half"><h2>rRNA Detection</h2><div id="rrnaChart" style="height:320px"></div></section>
<section class="card"><h2>Top ASVs</h2><table><thead><tr><th>ASV</th><th>Reads</th><th>Identity</th><th>Species</th><th>Genus</th></tr></thead><tbody>{otu_html}</tbody></table></section>
<section class="card"><h2>Metrics</h2><table><tbody>
<tr><th>Input bases</th><td>{fmt_int(metrics['input_bases'])}</td></tr>
<tr><th>Mean input length</th><td>{fmt_float(metrics['mean_input'], 1)}</td></tr>
<tr><th>Mean trimmed length</th><td>{fmt_float(metrics['mean_trimmed'], 1)}</td></tr>
<tr><th>Median trimmed length</th><td>{fmt_float(metrics['median_trimmed'], 1)}</td></tr>
<tr><th>N50 trimmed</th><td>{fmt_float(metrics['n50_trimmed'], 1)}</td></tr>
<tr><th>Retained bases</th><td>{fmt_float(metrics['retained_pct'], 2)}%</td></tr>
<tr><th>Simpson</th><td>{fmt_float(metrics['simpson'], 3)}</td></tr>
<tr><th>Top species</th><td>{html.escape(metrics['top_species'])} ({fmt_int(metrics['top_species_reads'])} ASV reads)</td></tr>
</tbody></table></section>
</div>
"""
    scripts = f"""
<script>
const topSpecies = {json.dumps(top_species)};
const topGenera = {json.dumps(top_genera)};
function bar(div, data, title) {{
  Plotly.newPlot(div, [{{x:data.map(d=>d[1]), y:data.map(d=>d[0]), type:'bar', orientation:'h', marker:{{color:'#12685f'}}}}],
    {{margin:{{l:170,r:10,t:10,b:35}}, xaxis:{{title:title}}, yaxis:{{autorange:'reversed'}}}}, {{responsive:true}});
}}
Plotly.newPlot('readChart', [{{x:['Input','18S trimmed','ASV assigned'], y:[{metrics['input_reads']},{metrics['trimmed_reads']},{total}], type:'bar', marker:{{color:['#66736c','#12685f','#c58a14']}}}}],
  {{margin:{{l:45,r:10,t:10,b:40}}, yaxis:{{title:'Reads'}}}}, {{responsive:true}});
bar('speciesChart', topSpecies, 'ASV reads');
bar('genusChart', topGenera, 'ASV reads');
Plotly.newPlot('rrnaChart', [{{labels:['18S','5.8S','28S'], values:[{metrics['reads_with_18S']},{metrics['reads_with_5_8S']},{metrics['reads_with_28S']}], type:'pie', hole:.4}}],
  {{margin:{{l:10,r:10,t:10,b:10}}}}, {{responsive:true}});
</script>
"""
    out_path.write_text(html_page(f"{sample} Savont report", body, scripts), encoding="utf-8")


def render_dashboard(metrics_rows, taxonomy_by_sample, out_path):
    rows_html = "\n".join(
        "<tr>"
        f"<td><a href=\"{html.escape(row['sample'])}_overview.html\">{html.escape(row['sample'])}</a></td>"
        f"<td>{fmt_int(row['input_reads'])}</td>"
        f"<td>{fmt_int(row['trimmed_reads'])}</td>"
        f"<td>{fmt_int(row['asv_count'])}</td>"
        f"<td>{fmt_float(row['shannon'])}</td>"
        f"<td>{fmt_float(row['retained_pct'])}%</td>"
        f"<td>{html.escape(row['top_species'])}</td>"
        "</tr>"
        for row in metrics_rows
    )
    body = f"""
<header>
<div><h1>Savont Biological Report</h1><div class="sub">ASV inference, 18S extraction, and PR2 annotation summary</div></div>
</header>
<div class="kpis">
<div class="kpi"><div class="label">Samples</div><div class="value">{len(metrics_rows)}</div></div>
<div class="kpi"><div class="label">Input reads</div><div class="value">{fmt_int(sum(r['input_reads'] for r in metrics_rows))}</div></div>
<div class="kpi"><div class="label">18S reads</div><div class="value">{fmt_int(sum(r['trimmed_reads'] for r in metrics_rows))}</div></div>
<div class="kpi"><div class="label">ASVs</div><div class="value">{fmt_int(sum(r['asv_count'] for r in metrics_rows))}</div></div>
</div>
<div class="grid">
<section class="card"><h2>Sample Summary</h2><table><thead><tr><th>Sample</th><th>Input reads</th><th>18S reads</th><th>ASVs</th><th>Shannon</th><th>Retained bases</th><th>Top species</th></tr></thead><tbody>{rows_html}</tbody></table></section>
<section class="card half"><h2>Reads By Sample</h2><div id="readsChart" style="height:340px"></div></section>
<section class="card half"><h2>ASVs And Diversity</h2><div id="divChart" style="height:340px"></div></section>
<section class="card"><h2>Top Species Across Samples</h2><div id="heatmap" style="height:460px"></div></section>
</div>
"""
    species_totals = defaultdict(int)
    for rows in taxonomy_by_sample.values():
        for species, count in aggregate_rank(rows, "Species"):
            species_totals[species] += count
    top_species = [item[0] for item in sorted(species_totals.items(), key=lambda item: item[1], reverse=True)[:20]]
    samples = [row["sample"] for row in metrics_rows]
    heatmap = []
    for species in top_species:
        heatmap.append([
            sum(int(row.get("reads", 0)) for row in taxonomy_by_sample.get(sample, []) if (row.get("Species") or "Unassigned") == species)
            for sample in samples
        ])
    scripts = f"""
<script>
const summary = {json.dumps(metrics_rows)};
const samples = {json.dumps(samples)};
const taxa = {json.dumps(top_species)};
const heatmap = {json.dumps(heatmap)};
Plotly.newPlot('readsChart', [
  {{x:samples, y:summary.map(d=>d.input_reads), name:'Input', type:'bar'}},
  {{x:samples, y:summary.map(d=>d.trimmed_reads), name:'18S trimmed', type:'bar'}}
], {{barmode:'group', margin:{{l:50,r:10,t:10,b:80}}, yaxis:{{title:'Reads'}}}}, {{responsive:true}});
Plotly.newPlot('divChart', [{{x:summary.map(d=>d.asv_count), y:summary.map(d=>d.shannon), text:samples, mode:'markers+text', textposition:'top center', marker:{{size:12,color:'#c58a14'}}}}],
  {{margin:{{l:50,r:10,t:10,b:50}}, xaxis:{{title:'ASVs'}}, yaxis:{{title:'Shannon'}}}}, {{responsive:true}});
Plotly.newPlot('heatmap', [{{z:heatmap, x:samples, y:taxa, type:'heatmap', colorscale:'Tealgrn'}}],
  {{margin:{{l:180,r:10,t:10,b:80}}}}, {{responsive:true}});
</script>
"""
    out_path.write_text(html_page("Savont Biological Report", body, scripts), encoding="utf-8")


def parse_args():
    parser = argparse.ArgumentParser(description="Generate Savont biological HTML reports")
    parser.add_argument("--trim-summaries", nargs="+", required=True)
    parser.add_argument("--asv-dirs", nargs="+", required=True)
    parser.add_argument("--taxonomy-tables", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--output-dir", default="savont_reports")
    return parser.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    trims = load_trim_summaries(args.trim_summaries)
    asvs = load_asv_dirs(args.asv_dirs)
    taxonomy = load_taxonomy_tables(args.taxonomy_tables)

    samples = sorted(set(trims) & set(asvs) & set(taxonomy))
    metrics_rows = []
    taxonomy_by_sample = {}
    for sample in samples:
        rows = taxonomy_with_counts(taxonomy[sample], asvs[sample]["counts"])
        taxonomy_by_sample[sample] = rows
        metrics = sample_metrics(sample, trims[sample], asvs[sample], taxonomy[sample])
        metrics["report_file"] = f"{sample}_overview.html"
        metrics_rows.append(metrics)
        render_sample(sample, metrics, rows, outdir / metrics["report_file"])

    dashboard = outdir / "multi_sample_dashboard.html"
    render_dashboard(metrics_rows, taxonomy_by_sample, dashboard)
    Path(args.output).write_text(dashboard.read_text(encoding="utf-8"), encoding="utf-8")


if __name__ == "__main__":
    main()
