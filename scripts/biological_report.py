#!/usr/bin/env python3

from argparse import ArgumentParser
from html import escape
from pathlib import Path
import re

import pandas as pd
import plotly.express as px
import plotly.io as pio


RANKS = [
    "domain",
    "supergroup",
    "division",
    "subdivision",
    "class",
    "order",
    "family",
    "genus",
    "species",
]


def parse_options():
    parser = ArgumentParser(description="Create a biological HTML summary report for BaNaNA outputs")
    parser.add_argument("--taxonomy", required=True, help="BaNaNA final taxonomy.tsv")
    parser.add_argument("--otu-table", required=True, help="BaNaNA final otu_table.tsv")
    parser.add_argument("--output", required=True, help="Output HTML report")
    parser.add_argument("--output-dir", default=None, help="Directory for multi-sample and per-sample reports")
    parser.add_argument("--nanoplot-dirs", nargs="*", default=[], help="NanoPlot output directories")
    parser.add_argument("--rrna-stats", nargs="*", default=[], help="rRNA extraction stats files")
    parser.add_argument("--top-n", type=int, default=15, help="Number of top taxa/OTUs to show")
    return parser.parse_args()


def clean_sample_name(name):
    name = str(name)
    name = name.replace("inputs/abundance_", "")
    name = name.replace("abundance_", "")
    return name


def load_otu_table(path):
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={df.columns[0]: "OTU"})
    sample_cols = [col for col in df.columns if col != "OTU"]
    rename = {col: clean_sample_name(col) for col in sample_cols}
    df = df.rename(columns=rename)
    sample_cols = [rename[col] for col in sample_cols]
    for col in sample_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)
    return df, sample_cols


def parse_taxonomy_field(value):
    parts = str(value).split("|")
    tax_parts = parts[4:] if len(parts) >= 13 else []
    values = {}
    for rank, taxon in zip(RANKS, tax_parts):
        values[rank] = taxon or "unclassified"
    for rank in RANKS:
        values.setdefault(rank, "unclassified")
    return values


def load_taxonomy(path):
    cols = [
        "OTU",
        "reference",
        "identity",
        "alignment_length",
        "mismatches",
        "gap_opens",
        "q_start",
        "q_end",
        "s_start",
        "s_end",
        "evalue",
        "bitscore",
    ]
    df = pd.read_csv(path, sep="\t", header=None, names=cols[:12])
    if df.empty:
        for rank in RANKS:
            df[rank] = pd.Series(dtype="object")
    else:
        parsed = df["reference"].apply(parse_taxonomy_field).apply(pd.Series)
        df = pd.concat([df, parsed], axis=1)
        for rank in RANKS:
            df[rank] = df[rank].fillna("unclassified")
    df["identity"] = pd.to_numeric(df["identity"], errors="coerce")
    return df


def sample_from_nanoplot_dir(path):
    name = Path(path).name
    return name.replace("nanoplot_", "")


def sample_from_rrna_stats(path):
    name = Path(path).name
    name = name.replace("rrna_extraction_stats_", "")
    return name.rsplit(".", 1)[0]


def load_nanostats(paths):
    records = {}
    for path in paths:
        sample = sample_from_nanoplot_dir(path)
        stats_file = Path(path) / "NanoStats.txt"
        if not stats_file.exists():
            continue
        stats = {}
        with stats_file.open() as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2 and parts[0] != "Metrics":
                    stats[parts[0]] = parts[1]
        records[sample] = stats
    return records


def parse_q_metric(value):
    if not isinstance(value, str):
        return 0.0, 0.0
    match = re.search(r"([0-9.]+)\s+\(([0-9.]+)%\)", value)
    if not match:
        return 0.0, 0.0
    return float(match.group(1)), float(match.group(2))


def load_rrna_stats(paths):
    records = {}
    for path in paths:
        sample = sample_from_rrna_stats(path)
        stats = {}
        current_total = None
        with Path(path).open() as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("All sequences:"):
                    current_total = int(line.split(":", 1)[1].strip())
                    continue
                match = re.match(r"(.+?) detected (once|multiple times|not detected):\s+(\d+)", line)
                if not match:
                    continue
                rrna, state, value = match.groups()
                rrna = rrna.replace("_rRNA", "").replace("_", ".")
                stats.setdefault(rrna, {"total": current_total or 0, "zero": 0, "one": 0, "multiple": 0})
                key = {"once": "one", "multiple times": "multiple", "not detected": "zero"}[state]
                stats[rrna][key] = int(value)
                stats[rrna]["total"] = current_total or stats[rrna]["total"]
        records[sample] = stats
    return records


def summarize_by_rank(long_df, rank, top_n):
    data = long_df.groupby(["sample", rank], as_index=False)["abundance"].sum()
    totals = data.groupby("sample", as_index=False)["abundance"].sum().rename(columns={"abundance": "total"})
    data = data.merge(totals, on="sample", how="left")
    data["relative_abundance"] = data["abundance"] / data["total"].where(data["total"] != 0, 1)

    top_taxa = (
        data.groupby(rank, as_index=False)["abundance"]
        .sum()
        .sort_values("abundance", ascending=False)
        .head(top_n)[rank]
        .tolist()
    )
    data[rank] = data[rank].where(data[rank].isin(top_taxa), "Other")
    data = data.groupby(["sample", rank], as_index=False)[["abundance", "relative_abundance"]].sum()
    return data


def fig_html(fig):
    return pio.to_html(fig, include_plotlyjs=False, full_html=False, config={"responsive": True})


def table_html(df, max_rows=25):
    return df.head(max_rows).to_html(index=False, escape=True, classes="summary-table")


def sample_report_html(sample, sample_summary, sample_otus, tax_records, nanostats, rrna_stats, top_n):
    q_labels = ["Q10", "Q15", "Q20", "Q25", "Q30"]
    q_values = [parse_q_metric(nanostats.get(f"Reads >{q}:"))[1] for q in q_labels]
    rrna_rows = []
    rrna_figs = []
    for rrna, vals in rrna_stats.items():
        total = vals.get("total", 0) or 1
        rrna_rows.append({
            "rRNA": rrna,
            "one": vals.get("one", 0),
            "multiple": vals.get("multiple", 0),
            "zero": vals.get("zero", 0),
            "detected_percent": round((vals.get("one", 0) + vals.get("multiple", 0)) / total * 100, 2),
        })
        fig = px.pie(
            names=["zero", "one", "multiple"],
            values=[vals.get("zero", 0), vals.get("one", 0), vals.get("multiple", 0)],
            title=f"{rrna} detection",
            hole=0.45,
        )
        fig.update_layout(template="plotly_white", margin=dict(l=20, r=20, t=60, b=20))
        rrna_figs.append(fig_html(fig))

    quality_fig = px.bar(
        x=q_labels,
        y=q_values,
        title="Read quality thresholds",
        labels={"x": "Threshold", "y": "Reads (%)"},
    )
    quality_fig.update_layout(template="plotly_white", margin=dict(l=40, r=20, t=60, b=50))

    species_data = sample_otus.groupby("species", as_index=False)["abundance"].sum().sort_values("abundance", ascending=False).head(top_n)
    species_fig = px.bar(
        species_data,
        x="abundance",
        y="species",
        orientation="h",
        title=f"Top {top_n} species",
        labels={"abundance": "Reads assigned to OTUs", "species": "Species"},
    )
    species_fig.update_layout(template="plotly_white", yaxis={"autorange": "reversed"}, margin=dict(l=160, r=20, t=60, b=50))

    class_data = sample_otus.groupby("class", as_index=False)["abundance"].sum().sort_values("abundance", ascending=False)
    class_fig = px.pie(class_data, names="class", values="abundance", title="Taxonomic composition by class", hole=0.35)
    class_fig.update_layout(template="plotly_white", margin=dict(l=20, r=20, t=60, b=20))

    top_otus = sample_otus.sort_values("abundance", ascending=False).head(top_n)
    top_species = species_data.iloc[0]["species"] if not species_data.empty else "none"
    mean_identity = top_otus["identity"].mean() if not top_otus.empty else 0

    rrna_table = table_html(pd.DataFrame(rrna_rows), max_rows=20) if rrna_rows else "<p>No rRNA stats available.</p>"
    top_table = table_html(top_otus[["OTU", "species", "genus", "abundance", "identity"]], max_rows=top_n)

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{escape(sample)} overview</title>
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>{REPORT_CSS}</style>
</head>
<body>
  <main class="container">
    <header><div><div class="pill">BaNaNA sample overview</div><h1>{escape(sample)}</h1><p class="muted">Sequencing quality, rRNA extraction, OTU composition, and taxonomy assignment.</p></div></header>
    <section class="grid metrics">
      <div class="metric"><div class="value">{int(float(nanostats.get("number_of_reads", 0) or 0))}</div><div>Filtered reads</div></div>
      <div class="metric"><div class="value">{nanostats.get("mean_read_length", "NA")}</div><div>Mean read length</div></div>
      <div class="metric"><div class="value">{parse_q_metric(nanostats.get("Reads >Q20:", ""))[0]:.0f}</div><div>Reads &gt; Q20</div></div>
      <div class="metric"><div class="value">{int(sample_summary.get("observed_otus", 0))}</div><div>Observed OTUs</div></div>
    </section>
    <section class="card"><h2>Sample Interpretation</h2><p>This sample contains {int(sample_summary.get("observed_otus", 0))} observed OTUs with {int(sample_summary.get("reads_assigned_to_otus", 0))} reads represented in the final OTU table. The dominant species-level assignment is <strong>{escape(str(top_species))}</strong>. Mean identity across the top OTUs shown below is {mean_identity:.2f}%.</p></section>
    <section class="grid">
      <div class="card half">{fig_html(quality_fig)}</div>
      <div class="card half">{''.join(rrna_figs) if rrna_figs else '<p>No rRNA plots available.</p>'}</div>
      <div class="card half">{fig_html(species_fig)}</div>
      <div class="card half">{fig_html(class_fig)}</div>
      <div class="card"><h2>rRNA Extraction Summary</h2>{rrna_table}</div>
      <div class="card"><h2>Top OTUs</h2>{top_table}</div>
    </section>
  </main>
</body>
</html>
"""


REPORT_CSS = """
body{font-family:Arial,sans-serif;margin:0;background:#f7f6f2;color:#22201b;line-height:1.5}
.container{max-width:1240px;margin:0 auto;padding:28px}
header{margin-bottom:20px}.pill{display:inline-block;background:#e8f1f1;color:#01696f;border-radius:999px;padding:5px 10px;font-weight:700;font-size:13px}
h1{font-size:42px;margin:8px 0 0}h2{margin:0 0 10px}.muted{color:#6f6b63}
.grid{display:grid;grid-template-columns:repeat(12,1fr);gap:16px}.metrics{margin-bottom:16px}
.metric,.card{background:#fbfbf9;border:1px solid #d8d4cc;border-radius:8px;padding:16px;box-shadow:0 6px 20px rgba(0,0,0,.04)}
.metric{grid-column:span 3}.metric .value{font-size:28px;font-weight:700}.half{grid-column:span 6}.card{grid-column:span 12}
.summary-table{border-collapse:collapse;width:100%;font-size:14px}.summary-table th,.summary-table td{border:1px solid #ddd;padding:7px 9px;text-align:left}.summary-table th{background:#f1efea}
@media(max-width:900px){.metric,.half{grid-column:span 12}}
"""


def write_dashboard(output_dir, sample_summary, long_df, nanostats_by_sample, rrna_by_sample, top_n):
    output_dir.mkdir(parents=True, exist_ok=True)
    samples = sample_summary["sample"].tolist()

    richness_fig = px.bar(sample_summary, x="sample", y="observed_otus", title="Observed OTUs per sample")
    reads_fig = px.bar(sample_summary, x="sample", y="reads_assigned_to_otus", title="Reads assigned to final OTUs")
    quality_rows = []
    for sample in samples:
        stats = nanostats_by_sample.get(sample, {})
        quality_rows.append({
            "sample": sample,
            "filtered_reads": float(stats.get("number_of_reads", 0) or 0),
            "mean_quality": float(stats.get("mean_qual", 0) or 0),
            "mean_read_length": float(stats.get("mean_read_length", 0) or 0),
            "q20_percent": parse_q_metric(stats.get("Reads >Q20:", ""))[1],
        })
    quality_df = pd.DataFrame(quality_rows)
    q20_fig = px.bar(quality_df, x="sample", y="q20_percent", title="Reads above Q20 after filtering")
    length_fig = px.bar(quality_df, x="sample", y="mean_read_length", title="Mean read length after filtering")

    genus_summary = summarize_by_rank(long_df, "genus", top_n)
    genus_fig = px.bar(genus_summary, x="sample", y="relative_abundance", color="genus", title="Relative composition by genus")
    class_summary = summarize_by_rank(long_df, "class", top_n)
    class_fig = px.bar(class_summary, x="sample", y="relative_abundance", color="class", title="Relative composition by class")

    rrna_rows = []
    for sample, rrnas in rrna_by_sample.items():
        for rrna, vals in rrnas.items():
            total = vals.get("total", 0) or 1
            rrna_rows.append({
                "sample": sample,
                "rRNA": rrna,
                "detected_percent": (vals.get("one", 0) + vals.get("multiple", 0)) / total * 100,
            })
    rrna_df = pd.DataFrame(rrna_rows)
    rrna_fig = px.bar(rrna_df, x="sample", y="detected_percent", color="rRNA", barmode="group", title="rRNA detection")

    for fig in [richness_fig, reads_fig, q20_fig, length_fig, genus_fig, class_fig, rrna_fig]:
        fig.update_layout(template="plotly_white", margin=dict(l=40, r=20, t=70, b=80))

    links = sample_summary[["sample", "reads_assigned_to_otus", "observed_otus"]].copy()
    links["report"] = links["sample"].apply(lambda sample: f'<a href="{escape(sample)}_overview.html">{escape(sample)}_overview.html</a>')
    table = links.to_html(index=False, escape=False, classes="summary-table")

    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>BaNaNA multi-sample dashboard</title>
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>{REPORT_CSS}</style>
</head>
<body>
  <main class="container">
    <header><div class="pill">BaNaNA dashboard</div><h1>Multi-Sample Dashboard</h1><p class="muted">Study-level overview of read quality, rRNA detection, OTU abundance, and taxonomy composition.</p></header>
    <section class="grid metrics">
      <div class="metric"><div class="value">{len(samples)}</div><div>Samples</div></div>
      <div class="metric"><div class="value">{int(sample_summary["reads_assigned_to_otus"].sum())}</div><div>Reads assigned to OTUs</div></div>
      <div class="metric"><div class="value">{int(sample_summary["observed_otus"].median())}</div><div>Median observed OTUs</div></div>
      <div class="metric"><div class="value">{long_df["OTU"].nunique()}</div><div>Total OTUs</div></div>
    </section>
    <section class="card"><h2>Sample Reports</h2>{table}</section>
    <section class="grid">
      <div class="card half">{fig_html(richness_fig)}</div>
      <div class="card half">{fig_html(reads_fig)}</div>
      <div class="card half">{fig_html(q20_fig)}</div>
      <div class="card half">{fig_html(length_fig)}</div>
      <div class="card half">{fig_html(rrna_fig)}</div>
      <div class="card half">{fig_html(class_fig)}</div>
      <div class="card">{fig_html(genus_fig)}</div>
    </section>
  </main>
</body>
</html>
"""
    (output_dir / "multi_sample_dashboard.html").write_text(html)


def main():
    args = parse_options()
    otu_df, sample_cols = load_otu_table(args.otu_table)
    tax_df = load_taxonomy(args.taxonomy)
    merged = otu_df.merge(tax_df, on="OTU", how="left")
    nanostats_by_sample = load_nanostats(args.nanoplot_dirs)
    rrna_by_sample = load_rrna_stats(args.rrna_stats)

    long_df = merged.melt(
        id_vars=["OTU", "identity"] + RANKS,
        value_vars=sample_cols,
        var_name="sample",
        value_name="abundance",
    )
    long_df["abundance"] = pd.to_numeric(long_df["abundance"], errors="coerce").fillna(0)
    detected = long_df[long_df["abundance"] > 0].copy()

    all_samples = pd.DataFrame({"sample": sample_cols})
    sample_totals = (
        long_df.groupby("sample", as_index=False)["abundance"]
        .sum()
        .rename(columns={"abundance": "reads_assigned_to_otus"})
    )
    sample_richness = (
        detected.groupby("sample", as_index=False)["OTU"]
        .nunique()
        .rename(columns={"OTU": "observed_otus"})
    )
    sample_summary = (
        all_samples
        .merge(sample_totals, on="sample", how="left")
        .merge(sample_richness, on="sample", how="left")
        .fillna(0)
    )
    sample_summary["observed_otus"] = sample_summary["observed_otus"].astype(int)

    otu_totals = otu_df.copy()
    otu_totals["total_abundance"] = otu_totals[sample_cols].sum(axis=1)
    otu_totals = otu_totals[["OTU", "total_abundance"]].merge(tax_df[["OTU", "identity", "genus", "species"]], on="OTU", how="left")
    otu_totals = otu_totals.sort_values("total_abundance", ascending=False)

    class_summary = summarize_by_rank(long_df, "class", args.top_n)
    genus_summary = summarize_by_rank(long_df, "genus", args.top_n)

    class_fig = px.bar(
        class_summary,
        x="sample",
        y="relative_abundance",
        color="class",
        title="Relative taxonomic composition by class",
        labels={"relative_abundance": "Relative abundance", "sample": "Sample", "class": "Class"},
    )
    genus_fig = px.bar(
        genus_summary,
        x="sample",
        y="relative_abundance",
        color="genus",
        title="Relative taxonomic composition by genus",
        labels={"relative_abundance": "Relative abundance", "sample": "Sample", "genus": "Genus"},
    )
    richness_fig = px.bar(
        sample_summary,
        x="sample",
        y="observed_otus",
        title="Observed OTUs per sample",
        labels={"observed_otus": "Observed OTUs", "sample": "Sample"},
    )
    abundance_fig = px.bar(
        sample_summary,
        x="sample",
        y="reads_assigned_to_otus",
        title="Reads assigned to final OTUs per sample",
        labels={"reads_assigned_to_otus": "Reads assigned to OTUs", "sample": "Sample"},
    )
    top_otu_fig = px.bar(
        otu_totals.head(args.top_n),
        x="OTU",
        y="total_abundance",
        color="genus",
        title=f"Top {args.top_n} OTUs by total abundance",
        labels={"total_abundance": "Total abundance", "genus": "Genus"},
    )

    for fig in [class_fig, genus_fig, richness_fig, abundance_fig, top_otu_fig]:
        fig.update_layout(template="plotly_white", legend_title_text="", margin=dict(l=40, r=20, t=70, b=80))

    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>BaNaNA Biological Report</title>
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 32px; color: #222; }}
    h1, h2 {{ margin-bottom: 0.3rem; }}
    .muted {{ color: #666; }}
    .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); gap: 16px; margin: 20px 0; }}
    .metric {{ border: 1px solid #ddd; padding: 16px; border-radius: 6px; }}
    .metric .value {{ font-size: 28px; font-weight: 700; }}
    .plot {{ margin: 26px 0; }}
    .summary-table {{ border-collapse: collapse; width: 100%; margin: 16px 0 32px; font-size: 14px; }}
    .summary-table th, .summary-table td {{ border: 1px solid #ddd; padding: 7px 9px; text-align: left; }}
    .summary-table th {{ background: #f4f4f4; }}
  </style>
</head>
<body>
  <h1>BaNaNA Biological Report</h1>
  <p class="muted">Generated from {escape(Path(args.taxonomy).name)} and {escape(Path(args.otu_table).name)}.</p>

  <div class="grid">
    <div class="metric"><div class="value">{len(sample_cols)}</div><div>Samples</div></div>
    <div class="metric"><div class="value">{len(otu_df)}</div><div>Final OTUs</div></div>
    <div class="metric"><div class="value">{int(sample_summary["reads_assigned_to_otus"].sum())}</div><div>Reads assigned to final OTUs</div></div>
    <div class="metric"><div class="value">{tax_df["identity"].mean():.1f}%</div><div>Mean taxonomy identity</div></div>
  </div>

  <h2>Sample Summary</h2>
  {table_html(sample_summary)}

  <div class="plot">{fig_html(richness_fig)}</div>
  <div class="plot">{fig_html(abundance_fig)}</div>
  <div class="plot">{fig_html(class_fig)}</div>
  <div class="plot">{fig_html(genus_fig)}</div>
  <div class="plot">{fig_html(top_otu_fig)}</div>

  <h2>Top OTUs</h2>
  {table_html(otu_totals)}
</body>
</html>
"""
    Path(args.output).write_text(html)

    if args.output_dir:
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        write_dashboard(output_dir, sample_summary, long_df, nanostats_by_sample, rrna_by_sample, args.top_n)
        for sample in sample_cols:
            sample_otus = detected[detected["sample"] == sample].copy()
            sample_stats = sample_summary[sample_summary["sample"] == sample]
            if sample_stats.empty:
                sample_stats_dict = {"observed_otus": 0, "reads_assigned_to_otus": 0}
            else:
                sample_stats_dict = sample_stats.iloc[0].to_dict()
            sample_html = sample_report_html(
                sample=sample,
                sample_summary=sample_stats_dict,
                sample_otus=sample_otus,
                tax_records=tax_df,
                nanostats=nanostats_by_sample.get(sample, {}),
                rrna_stats=rrna_by_sample.get(sample, {}),
                top_n=args.top_n,
            )
            (output_dir / f"{sample}_overview.html").write_text(sample_html)


if __name__ == "__main__":
    main()
