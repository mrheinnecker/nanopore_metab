#!/usr/bin/env python3
# run as -> python banana_sample_reporter_dashboard.py /path/to/run_root --outdir
import argparse, math, json, re, html
from pathlib import Path
import pandas as pd

HTML_TEMPLATE = '''<!doctype html>
<html lang="en" data-theme="light">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=IBM+Plex+Serif:wght@600&display=swap" rel="stylesheet">
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
:root,[data-theme="light"]{--bg:#f7f6f2;--surface:#fbfbf9;--surface2:#f1efea;--text:#22201b;--muted:#6f6b63;--primary:#01696f;--border:#d8d4cc;--accent:#d19900;--shadow:0 8px 24px rgba(0,0,0,.05)}
[data-theme="dark"]{--bg:#171614;--surface:#1d1c1a;--surface2:#252320;--text:#ece8df;--muted:#a49e93;--primary:#4f98a3;--border:#3b3832;--accent:#e8af34;--shadow:0 8px 24px rgba(0,0,0,.24)}
*{box-sizing:border-box} body{margin:0;font-family:Inter,system-ui,sans-serif;background:var(--bg);color:var(--text);line-height:1.5}
.container{max-width:1180px;margin:0 auto;padding:24px} header{display:flex;justify-content:space-between;align-items:flex-start;gap:16px;flex-wrap:wrap;margin-bottom:24px}
.brand{display:flex;gap:14px;align-items:center} .logo{width:42px;height:42px;border-radius:12px;background:linear-gradient(135deg,var(--primary),var(--accent));box-shadow:var(--shadow)}
h1{font:600 clamp(2rem,3vw,3.2rem) 'IBM Plex Serif',serif;margin:0} h2{margin:0;font-size:1.2rem} .sub{color:var(--muted);margin-top:6px}
button{border:1px solid var(--border);background:var(--surface);color:var(--text);padding:10px 14px;border-radius:10px;cursor:pointer}
.grid{display:grid;grid-template-columns:repeat(12,1fr);gap:16px} .card{grid-column:span 12;background:var(--surface);border:1px solid var(--border);border-radius:18px;padding:18px;box-shadow:var(--shadow)}
.kpi{display:grid;grid-template-columns:repeat(4,1fr);gap:14px} .kpi .mini{background:var(--surface2);border:1px solid var(--border);border-radius:14px;padding:14px} .kpi .v{font-size:1.5rem;font-weight:700} .kpi .l{color:var(--muted);font-size:.92rem}
.toolbar{display:flex;justify-content:space-between;align-items:center;gap:12px;flex-wrap:wrap;margin-bottom:10px} .seg{display:inline-flex;background:var(--surface2);border:1px solid var(--border);border-radius:12px;padding:4px;gap:4px} .seg button{padding:8px 12px;border:none;background:transparent} .seg button.active{background:var(--primary);color:white}
.helper{font-size:.9rem;color:var(--muted)} .breadcrumbs{display:flex;gap:8px;align-items:center;flex-wrap:wrap;margin:8px 0 10px 0} .crumb{padding:5px 9px;border-radius:999px;background:var(--surface2);border:1px solid var(--border);font-size:.84rem} .crumb.active{background:rgba(1,105,111,.12);color:var(--primary);border-color:rgba(1,105,111,.25)}
.pill{display:inline-block;padding:4px 10px;border-radius:999px;background:rgba(1,105,111,.1);color:var(--primary);font-size:.85rem;font-weight:600} table{width:100%;border-collapse:collapse} th,td{padding:10px 8px;border-bottom:1px solid var(--border);text-align:left;font-size:.95rem;vertical-align:top} th{color:var(--muted);font-weight:600} .footer{color:var(--muted);font-size:.9rem;margin-top:18px}
@media (max-width:900px){.kpi{grid-template-columns:repeat(2,1fr)} .half{grid-column:span 12 !important}} @media (min-width:901px){.half{grid-column:span 6}}
</style>
</head>
<body>
<div class="container">
<header>
<div class="brand"><div class="logo"></div><div><div class="pill">BaNaNA sample overview</div><h1>__SAMPLE__</h1><div class="sub">Interactive overview of sequencing quality, rRNA extraction, OTU composition, and taxonomy assignment.</div></div></div>
<div style="display:flex;gap:10px"><button onclick="toggleTheme()">Toggle theme</button><button onclick="window.print()">Print / Save PDF</button></div>
</header>
<section class="card"><h2>Key metrics</h2><div class="kpi">
<div class="mini"><div class="l">Raw reads</div><div class="v">__RAWREADS__</div></div>
<div class="mini"><div class="l">Mean read length</div><div class="v">__MEANLEN__</div></div>
<div class="mini"><div class="l">Reads > Q20</div><div class="v">__Q20__</div></div>
<div class="mini"><div class="l">Observed OTUs</div><div class="v">__OBSOTUS__</div></div>
</div></section>
<div class="grid" style="margin-top:16px">
<section class="card half"><h2>Sequencing quality</h2><div id="qualChart" style="height:340px"></div></section>
<section class="card half"><h2>rRNA extraction</h2><div id="rrnaChart" style="height:340px"></div></section>
<section class="card half"><div class="toolbar"><h2>Top species</h2><div class="seg"><button id="speciesReadsBtn" class="active" onclick="setSpeciesMode('reads')">Reads</button><button id="speciesPctBtn" onclick="setSpeciesMode('percent')">Percent</button></div></div><div class="helper">Switch between absolute read counts and relative abundance percentages.</div><div id="speciesChart" style="height:360px"></div></section>
<section class="card half"><div class="toolbar"><div><h2>Taxonomic composition</h2><div class="helper">Click a pie segment to drill down to the next taxonomic rank. Use breadcrumbs to go back.</div></div><button onclick="resetTaxonomy()">Reset</button></div><div class="breadcrumbs" id="taxBreadcrumbs"></div><div id="taxonomyChart" style="height:390px"></div></section>
<section class="card"><h2>Sample interpretation</h2><p>This sample contains __OBSOTUS__ OTUs with __TOTALREADS__ reads represented in the OTU table. Shannon diversity is __SHANNON__ and Simpson diversity is __SIMPSON__, which summarize richness and evenness across detected OTUs.</p><p>The dominant species-level assignment is <strong>__TOPSPECIES__</strong> with __TOPSPECIESREADS__ reads, and the mean identity across the five most abundant OTUs is __MEANTOP5__%.</p></section>
<section class="card half"><h2>Top OTUs</h2><table><thead><tr><th>OTU</th><th>Species</th><th>Reads</th><th>Pident</th></tr></thead><tbody>__TOPOTUS__</tbody></table></section>
<section class="card half"><h2>Overview table</h2><table><tbody>
<tr><th>Median read length</th><td>__MEDIANLEN__</td></tr>
<tr><th>N50</th><td>__N50__</td></tr>
<tr><th>Mean quality</th><td>__MEANQUAL__</td></tr>
<tr><th>18S detected</th><td>__R18__</td></tr>
<tr><th>5.8S detected</th><td>__R58__</td></tr>
<tr><th>28S detected</th><td>__R28__</td></tr>
<tr><th>Shannon</th><td>__SHANNON__</td></tr>
<tr><th>Simpson</th><td>__SIMPSON__</td></tr>
</tbody></table></section>
</div>
<div class="footer">Open this report in a browser for interactive charts, then use Print / Save PDF to export a PDF version.</div>
</div>
<script>
const speciesLabels = __SPECIES_LABELS__;
const speciesReads = __SPECIES_READS__;
const speciesPercent = __SPECIES_PCT__;
const taxonomyRecords = __TAX_RECORDS__;
const taxonomyLevels = __TAX_LEVELS__;
let speciesMode = 'reads';
let taxState = {levelIndex: 0, path: []};
function themeVars() { return {paper: 'rgba(0,0,0,0)', fontColor: getComputedStyle(document.documentElement).getPropertyValue('--text').trim(), gridColor: getComputedStyle(document.documentElement).getPropertyValue('--border').trim()}; }
function toggleTheme() { const htmlEl = document.documentElement; htmlEl.setAttribute('data-theme', htmlEl.getAttribute('data-theme') === 'dark' ? 'light' : 'dark'); drawAll(); }
function drawQuality() { const t = themeVars(); Plotly.newPlot('qualChart', [{x:['Q10','Q15','Q20','Q25','Q30'], y:[__Q10_PCT__,__Q15_PCT__,__Q20_PCT__,__Q25_PCT__,__Q30_PCT__], type:'bar', marker:{color:['#4f98a3','#4f98a3','#01696f','#d19900','#a13544']}}], {margin:{l:40,r:10,t:10,b:40}, paper_bgcolor:t.paper, plot_bgcolor:t.paper, font:{color:t.fontColor}, yaxis:{title:'% reads', gridcolor:t.gridColor}}, {responsive:true}); }
function drawRrna() { const t = themeVars(); Plotly.newPlot('rrnaChart', [{labels:['18S','5.8S','28S'], values:[__RR18__,__RR58__,__RR28__], type:'pie', hole:.45, marker:{colors:['#01696f','#d19900','#437a22']}}], {margin:{l:10,r:10,t:10,b:10}, paper_bgcolor:t.paper, font:{color:t.fontColor}}, {responsive:true}); }
function drawSpecies() { const t = themeVars(); const x = speciesMode === 'reads' ? speciesReads : speciesPercent; const xTitle = speciesMode === 'reads' ? 'Reads' : 'Percent of OTU-table reads'; const hover = speciesMode === 'reads' ? '%{y}<br>%{x} reads<extra></extra>' : '%{y}<br>%{x:.3f}%<extra></extra>'; Plotly.newPlot('speciesChart', [{x:x, y:speciesLabels, type:'bar', orientation:'h', marker:{color:'#01696f'}, hovertemplate:hover}], {margin:{l:170,r:10,t:10,b:40}, paper_bgcolor:t.paper, plot_bgcolor:t.paper, font:{color:t.fontColor}, xaxis:{title:xTitle, gridcolor:t.gridColor}, yaxis:{autorange:'reversed'}}, {responsive:true}); document.getElementById('speciesReadsBtn').classList.toggle('active', speciesMode === 'reads'); document.getElementById('speciesPctBtn').classList.toggle('active', speciesMode === 'percent'); }
function setSpeciesMode(mode) { speciesMode = mode; drawSpecies(); }
function aggregateTaxonomy(records, levelIndex, path) { const filtered = records.filter(rec => path.every((p, i) => rec[taxonomyLevels[i]] === p)); const current = taxonomyLevels[levelIndex]; const buckets = new Map(); filtered.forEach(rec => { const key = rec[current] || 'Unassigned'; buckets.set(key, (buckets.get(key) || 0) + rec.reads); }); const arr = Array.from(buckets.entries()).map(([label, value]) => ({label, value})).sort((a,b) => b.value - a.value); return {current, arr}; }
function drawBreadcrumbs() { const wrap = document.getElementById('taxBreadcrumbs'); wrap.innerHTML = taxonomyLevels.map((lv, i) => { const active = i === taxState.levelIndex ? ' active' : ''; const label = i < taxState.path.length ? `${lv}: ${taxState.path[i]}` : lv; return `<button class="crumb${active}" onclick="jumpToLevel(${i})">${label}</button>`; }).join(''); }
function drawTaxonomy() { const t = themeVars(); const agg = aggregateTaxonomy(taxonomyRecords, taxState.levelIndex, taxState.path); const labels = agg.arr.map(d => d.label); const values = agg.arr.map(d => d.value); const titleText = taxState.path.length ? `${agg.current} within ` + taxState.path.join(' › ') : agg.current; Plotly.newPlot('taxonomyChart', [{labels:labels, values:values, type:'pie', hole:.35, sort:false, textinfo:'label+percent'}], {title:{text:titleText,font:{size:16}}, margin:{l:10,r:10,t:40,b:10}, paper_bgcolor:t.paper, font:{color:t.fontColor}}, {responsive:true}); const chart = document.getElementById('taxonomyChart'); chart.on('plotly_click', function(ev) { if (taxState.levelIndex >= taxonomyLevels.length - 1) return; const clicked = ev.points[0].label; taxState.path = taxState.path.slice(0, taxState.levelIndex); taxState.path.push(clicked); taxState.levelIndex += 1; drawBreadcrumbs(); drawTaxonomy(); }); }
function jumpToLevel(i) { taxState.levelIndex = i; taxState.path = taxState.path.slice(0, i); drawBreadcrumbs(); drawTaxonomy(); }
function resetTaxonomy() { taxState = {levelIndex: 0, path: []}; drawBreadcrumbs(); drawTaxonomy(); }
function drawAll() { drawQuality(); drawRrna(); drawSpecies(); drawBreadcrumbs(); drawTaxonomy(); }
drawAll();
</script>
</body>
</html>'''

DASHBOARD_TEMPLATE = '''<!doctype html>
<html lang="en" data-theme="light">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>BaNaNA multi-sample dashboard</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=IBM+Plex+Serif:wght@600&display=swap" rel="stylesheet">
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
:root,[data-theme="light"]{--bg:#f7f6f2;--surface:#fbfbf9;--surface2:#f1efea;--text:#22201b;--muted:#6f6b63;--primary:#01696f;--border:#d8d4cc;--accent:#d19900;--shadow:0 8px 24px rgba(0,0,0,.05)}
[data-theme="dark"]{--bg:#171614;--surface:#1d1c1a;--surface2:#252320;--text:#ece8df;--muted:#a49e93;--primary:#4f98a3;--border:#3b3832;--accent:#e8af34;--shadow:0 8px 24px rgba(0,0,0,.24)}
*{box-sizing:border-box} body{margin:0;font-family:Inter,system-ui,sans-serif;background:var(--bg);color:var(--text);line-height:1.5}
.container{max-width:1280px;margin:0 auto;padding:24px}.hero{display:flex;justify-content:space-between;gap:16px;flex-wrap:wrap;align-items:flex-start;margin-bottom:18px}
.logo{width:46px;height:46px;border-radius:14px;background:linear-gradient(135deg,var(--primary),var(--accent));box-shadow:var(--shadow)} .brand{display:flex;gap:14px;align-items:center}
h1{font:600 clamp(2rem,3vw,3.2rem) 'IBM Plex Serif',serif;margin:0} h2{margin:0;font-size:1.2rem}.sub{color:var(--muted);margin-top:6px;max-width:70ch}
button{border:1px solid var(--border);background:var(--surface);color:var(--text);padding:10px 14px;border-radius:10px;cursor:pointer}
.kpis{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin:16px 0}.kpi{background:var(--surface);border:1px solid var(--border);border-radius:16px;padding:16px;box-shadow:var(--shadow)}.kpi .l{font-size:.9rem;color:var(--muted)}.kpi .v{font-size:1.6rem;font-weight:700}
.grid{display:grid;grid-template-columns:repeat(12,1fr);gap:16px}.card{grid-column:span 12;background:var(--surface);border:1px solid var(--border);border-radius:18px;padding:18px;box-shadow:var(--shadow)}.half{grid-column:span 6}
.helper{font-size:.9rem;color:var(--muted);margin-top:4px}.toolbar{display:flex;justify-content:space-between;align-items:center;gap:10px;flex-wrap:wrap;margin-bottom:10px}.seg{display:inline-flex;background:var(--surface2);border:1px solid var(--border);border-radius:12px;padding:4px;gap:4px}.seg button{padding:8px 12px;border:none;background:transparent}.seg button.active{background:var(--primary);color:#fff}
table{width:100%;border-collapse:collapse}th,td{padding:10px 8px;border-bottom:1px solid var(--border);text-align:left;font-size:.95rem}th{color:var(--muted);font-weight:600}a{color:var(--primary);text-decoration:none}a:hover{text-decoration:underline}
@media (max-width:950px){.kpis{grid-template-columns:repeat(2,1fr)}.half{grid-column:span 12}}
</style>
</head>
<body>
<div class="container">
<div class="hero">
<div class="brand"><div class="logo"></div><div><h1>Multi-sample dashboard</h1><div class="sub">Study-level overview of sequencing quality, diversity, rRNA detection, and taxonomy composition across all processed samples.</div></div></div>
<div><button onclick="toggleTheme()">Toggle theme</button></div>
</div>
<div class="kpis">
<div class="kpi"><div class="l">Samples</div><div class="v">__NSAMPLES__</div></div>
<div class="kpi"><div class="l">Total reads</div><div class="v">__TOTALREADS__</div></div>
<div class="kpi"><div class="l">Median OTUs</div><div class="v">__MEDOTUS__</div></div>
<div class="kpi"><div class="l">Median Shannon</div><div class="v">__MEDSHANNON__</div></div>
</div>
<div class="grid">
<section class="card"><div class="toolbar"><div><h2>Sample summary</h2><div class="helper">Click a sample name to open its individual HTML report.</div></div></div><div style="overflow:auto"><table><thead><tr><th>Sample</th><th>Raw reads</th><th>Mean length</th><th>Q20 reads</th><th>Observed OTUs</th><th>Shannon</th><th>Simpson</th><th>Top species</th></tr></thead><tbody>__SUMMARY_ROWS__</tbody></table></div></section>
<section class="card half"><h2>Reads per sample</h2><div id="readsChart" style="height:350px"></div></section>
<section class="card half"><h2>Diversity comparison</h2><div id="divChart" style="height:350px"></div></section>
<section class="card half"><h2>rRNA detection</h2><div id="rrnaChart" style="height:360px"></div></section>
<section class="card half"><div class="toolbar"><div><h2>Taxonomy heatmap</h2><div class="helper">Relative abundance by sample for the most abundant taxa at the selected rank.</div></div><div class="seg"><button id="rankDivision" class="active" onclick="setRank('Division')">Division</button><button id="rankClass" onclick="setRank('Class')">Class</button><button id="rankFamily" onclick="setRank('Family')">Family</button><button id="rankGenus" onclick="setRank('Genus')">Genus</button><button id="rankSpecies" onclick="setRank('Species')">Species</button></div></div><div id="taxHeatmap" style="height:420px"></div></section>
</div>
</div>
<script>
const summaryData = __SUMMARY_DATA__;
const taxonomyMatrices = __TAX_MATRICES__;
let currentRank = 'Division';
function themeVars(){return {paper:'rgba(0,0,0,0)', fontColor:getComputedStyle(document.documentElement).getPropertyValue('--text').trim(), gridColor:getComputedStyle(document.documentElement).getPropertyValue('--border').trim()};}
function toggleTheme(){const htmlEl=document.documentElement;htmlEl.setAttribute('data-theme', htmlEl.getAttribute('data-theme') === 'dark' ? 'light' : 'dark'); drawAll();}
function drawReads(){const t=themeVars(); Plotly.newPlot('readsChart',[{x:summaryData.map(d=>d.sample), y:summaryData.map(d=>d.raw_reads_num), type:'bar', marker:{color:'#01696f'}}],{margin:{l:50,r:10,t:10,b:80},paper_bgcolor:t.paper,plot_bgcolor:t.paper,font:{color:t.fontColor},xaxis:{tickangle:-40},yaxis:{title:'Raw reads',gridcolor:t.gridColor}},{responsive:true});}
function drawDiversity(){const t=themeVars(); Plotly.newPlot('divChart',[{x:summaryData.map(d=>d.observed_otus), y:summaryData.map(d=>d.shannon), text:summaryData.map(d=>d.sample), mode:'markers+text', textposition:'top center', type:'scatter', marker:{size:12,color:'#d19900',line:{color:'#01696f',width:1}}}],{margin:{l:55,r:10,t:10,b:45},paper_bgcolor:t.paper,plot_bgcolor:t.paper,font:{color:t.fontColor},xaxis:{title:'Observed OTUs',gridcolor:t.gridColor},yaxis:{title:'Shannon',gridcolor:t.gridColor}},{responsive:true});}
function drawRrna(){const t=themeVars(); Plotly.newPlot('rrnaChart',[{x:summaryData.map(d=>d.sample), y:summaryData.map(d=>d.rrna18_pct), name:'18S', type:'bar'},{x:summaryData.map(d=>d.sample), y:summaryData.map(d=>d.rrna58_pct), name:'5.8S', type:'bar'},{x:summaryData.map(d=>d.sample), y:summaryData.map(d=>d.rrna28_pct), name:'28S', type:'bar'}],{barmode:'group',margin:{l:50,r:10,t:10,b:80},paper_bgcolor:t.paper,plot_bgcolor:t.paper,font:{color:t.fontColor},xaxis:{tickangle:-40},yaxis:{title:'Detected (%)',gridcolor:t.gridColor}},{responsive:true});}
function drawHeatmap(){const t=themeVars(); const d=taxonomyMatrices[currentRank]; Plotly.newPlot('taxHeatmap',[{z:d.z,x:d.samples,y:d.taxa,type:'heatmap',colorscale:'Tealgrn',hovertemplate:'Sample: %{x}<br>Taxon: %{y}<br>Relative abundance: %{z:.3f}%<extra></extra>'}],{margin:{l:140,r:10,t:10,b:80},paper_bgcolor:t.paper,plot_bgcolor:t.paper,font:{color:t.fontColor},xaxis:{tickangle:-40}},{responsive:true}); ['Division','Class','Family','Genus','Species'].forEach(r=>document.getElementById('rank'+r).classList.toggle('active', r===currentRank));}
function setRank(rank){currentRank=rank; drawHeatmap();}
function drawAll(){drawReads(); drawDiversity(); drawRrna(); drawHeatmap();}
drawAll();
</script>
</body>
</html>'''


def parse_config_sample_names(config_path: Path):
    text = config_path.read_text(encoding='utf-8', errors='ignore')
    lines = text.splitlines()
    names = []
    in_sample_block = False
    sample_block_indent = None
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith('#'):
            continue
        if re.match(r'^sample_name\s*:\s*$', stripped):
            in_sample_block = True
            sample_block_indent = len(line) - len(line.lstrip())
            continue
        if in_sample_block:
            current_indent = len(line) - len(line.lstrip())
            if current_indent <= sample_block_indent and re.match(r'^[A-Za-z0-9_\-]+\s*:', stripped):
                break
            m = re.match(r'^["\']?([^:"\'#]+?)["\']?\s*:\s*["\']?([^"\'#]+?)["\']?\s*$', stripped)
            if m:
                names.append(m.group(1).strip())
                continue
            m2 = re.match(r'^-\s*["\']?([^"\'#]+?)["\']?\s*$', stripped)
            if m2:
                names.append(m2.group(1).strip())
    if names:
        return names
    m = re.search(r'^sample_name\s*:\s*(["\']?)([^\n#"\']+)\1\s*$', text, re.M)
    if m:
        return [m.group(2).strip()]
    return []


def clean_sample_name(name):
    name = str(name)
    name = name.replace("inputs/abundance_", "")
    name = name.replace("abundance_", "")
    return name


def sample_from_nanoplot_dir(path: Path):
    return path.name.replace("nanoplot_", "")


def sample_from_rrna_stats(path: Path):
    name = path.name.replace("rrna_extraction_stats_", "")
    return name.rsplit(".", 1)[0]


def parse_q_percent(value):
    match = re.search(r"\(([0-9.]+)%\)", str(value))
    return float(match.group(1)) if match else 0.0


def load_otu_table(path: Path):
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={df.columns[0]: "OTU"})
    rename = {col: clean_sample_name(col) for col in df.columns if col != "OTU"}
    df = df.rename(columns=rename)
    for col in rename.values():
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)
    return df, list(rename.values())


def parse_taxonomy_reference(reference):
    parts = str(reference).split("|")
    accession = parts[0] if parts else "Unassigned"
    taxa = parts[4:] if len(parts) >= 13 else []
    ranks = ["Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species"]
    values = {rank: "Unassigned" for rank in ranks}
    for rank, taxon in zip(ranks, taxa):
        values[rank] = taxon or "Unassigned"
    return accession, values


def load_taxonomy(path: Path):
    cols = [
        "OTU",
        "Reference",
        "Pident",
        "Length",
        "Mismatches",
        "Gap_opens",
        "Q_start",
        "Q_end",
        "S_start",
        "S_end",
        "Evalue",
        "Bitscore",
    ]
    if path.stat().st_size == 0:
        return pd.DataFrame(columns=["OTU", "Pident", "Accession", "Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species", "Length"])
    df = pd.read_csv(path, sep="\t", header=None, names=cols[:12])
    parsed = df["Reference"].apply(parse_taxonomy_reference)
    df["Accession"] = parsed.apply(lambda value: value[0])
    tax_df = parsed.apply(lambda value: value[1]).apply(pd.Series)
    df = pd.concat([df[["OTU", "Pident", "Accession", "Length"]], tax_df], axis=1)
    for col in ["Pident", "Length"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def parse_nanostats(path: Path):
    data = {}
    with path.open(encoding='utf-8', errors='ignore') as f:
        for line in f:
            if '\t' in line:
                k, v = line.rstrip('\n').split('\t', 1)
                data[k] = v
    return data


def parse_rrna_stats(path: Path):
    rrna = {}
    all_sequences = None
    with path.open(encoding='utf-8', errors='ignore') as f:
        for line in f:
            m = re.match(r'([^:]+):\s+(\d+)', line.strip())
            if m:
                key, val = m.group(1), int(m.group(2))
                if key == 'All sequences' and all_sequences is None:
                    all_sequences = val
                elif key != 'All sequences':
                    rrna[key] = val
    return all_sequences, rrna


def load_nanostats_by_sample(paths):
    records = {}
    for path in paths:
        path = Path(path)
        sample = sample_from_nanoplot_dir(path)
        stats_file = path / "NanoStats.txt"
        if stats_file.exists():
            records[sample] = parse_nanostats(stats_file)
    return records


def load_rrna_stats_by_sample(paths):
    records = {}
    for path in paths:
        path = Path(path)
        sample = sample_from_rrna_stats(path)
        if path.exists():
            records[sample] = parse_rrna_stats(path)
    return records


def pick_existing(paths):
    for p in paths:
        if p.exists():
            return p
    return None


def make_report(sample, df, nano, all_sequences, rrna, out_html: Path):
    levels = ['Supergroup', 'Division', 'Subdivision', 'Class', 'Order', 'Family', 'Genus', 'Species']
    total_otu_reads = int(df[sample].sum())
    observed_otus = int((df[sample] > 0).sum())
    rel = df[sample] / total_otu_reads if total_otu_reads else pd.Series(dtype=float)
    shannon = float(-(rel[rel > 0] * rel[rel > 0].map(math.log)).sum()) if total_otu_reads else 0.0
    simpson = float(1 - (rel ** 2).sum()) if total_otu_reads else 0.0
    mean_pident_top5 = float(df['Pident'].head(5).mean()) if len(df) else float('nan')
    top_species = df.groupby('Species', dropna=False)[sample].sum().sort_values(ascending=False).head(10)
    top_species_labels = [str(x) if pd.notnull(x) else 'Unassigned' for x in top_species.index]
    top_species_counts = [int(x) for x in top_species.values]
    top_species_pct = [round(v / total_otu_reads * 100, 3) if total_otu_reads else 0 for v in top_species_counts]
    taxonomy_records = []
    for _, r in df.iterrows():
        rec = {'reads': int(r[sample])}
        for lv in levels:
            rec[lv] = str(r[lv]) if pd.notnull(r[lv]) else 'Unassigned'
        taxonomy_records.append(rec)
    rrna18 = rrna.get('18S_rRNA detected once', 0) + rrna.get('18S_rRNA detected multiple times', 0)
    rrna58 = rrna.get('5_8S_rRNA detected once', 0) + rrna.get('5_8S_rRNA detected multiple times', 0)
    rrna28 = rrna.get('28S_rRNA detected once', 0) + rrna.get('28S_rRNA detected multiple times', 0)
    top_otu_html = ''.join(
        f"<tr><td>{html.escape(str(r['OTU']))}</td><td>{html.escape(str(r['Species']) if pd.notnull(r['Species']) else 'NA')}</td><td>{int(r[sample])}</td><td>{format(r['Pident'], '.1f') if pd.notnull(r['Pident']) else 'NA'}</td></tr>"
        for _, r in df.head(10).iterrows()
    )
    rep = HTML_TEMPLATE
    replacements = {
        '__TITLE__': html.escape(sample + ' report'), '__SAMPLE__': html.escape(sample), '__RAWREADS__': html.escape(str(nano.get('number_of_reads', 'NA'))),
        '__MEANLEN__': html.escape(str(nano.get('mean_read_length', 'NA'))), '__Q20__': html.escape(str(nano.get('Reads >Q20:', 'NA').split()[0] if nano.get('Reads >Q20:') else 'NA')),
        '__OBSOTUS__': str(observed_otus), '__TOTALREADS__': str(total_otu_reads), '__SHANNON__': f'{shannon:.2f}', '__SIMPSON__': f'{simpson:.3f}',
        '__TOPSPECIES__': html.escape(top_species_labels[0] if top_species_labels else 'NA'), '__TOPSPECIESREADS__': str(top_species_counts[0] if top_species_counts else 0),
        '__MEANTOP5__': f'{mean_pident_top5:.2f}', '__TOPOTUS__': top_otu_html, '__MEDIANLEN__': html.escape(str(nano.get('median_read_length', 'NA'))),
        '__N50__': html.escape(str(nano.get('n50', 'NA'))), '__MEANQUAL__': html.escape(str(nano.get('mean_qual', 'NA'))),
        '__R18__': f"{rrna18} / {all_sequences if all_sequences else 'NA'} ({round(rrna18/all_sequences*100,1) if all_sequences else 'NA'}%)",
        '__R58__': f"{rrna58} / {all_sequences if all_sequences else 'NA'} ({round(rrna58/all_sequences*100,1) if all_sequences else 'NA'}%)",
        '__R28__': f"{rrna28} / {all_sequences if all_sequences else 'NA'} ({round(rrna28/all_sequences*100,1) if all_sequences else 'NA'}%)",
        '__SPECIES_LABELS__': json.dumps(top_species_labels), '__SPECIES_READS__': json.dumps(top_species_counts), '__SPECIES_PCT__': json.dumps(top_species_pct),
        '__TAX_RECORDS__': json.dumps(taxonomy_records), '__TAX_LEVELS__': json.dumps(levels), '__RR18__': str(rrna18), '__RR58__': str(rrna58), '__RR28__': str(rrna28),
        '__Q10_PCT__': str(parse_q_percent(nano.get('Reads >Q10:', '0%'))),
        '__Q15_PCT__': str(parse_q_percent(nano.get('Reads >Q15:', '0%'))),
        '__Q20_PCT__': str(parse_q_percent(nano.get('Reads >Q20:', '0%'))),
        '__Q25_PCT__': str(parse_q_percent(nano.get('Reads >Q25:', '0%'))),
        '__Q30_PCT__': str(parse_q_percent(nano.get('Reads >Q30:', '0%'))),
    }
    for k, v in replacements.items():
        rep = rep.replace(k, v)
    out_html.write_text(rep, encoding='utf-8')
    return {
        'sample': sample,
        'raw_reads': nano.get('number_of_reads', 'NA'),
        'raw_reads_num': int(re.sub(r'[^0-9]', '', str(nano.get('number_of_reads', '0')))) if re.search(r'\d', str(nano.get('number_of_reads', '0'))) else 0,
        'mean_read_length': nano.get('mean_read_length', 'NA'),
        'reads_gt_q20': nano.get('Reads >Q20:', 'NA'),
        'observed_otus': observed_otus,
        'shannon': shannon,
        'simpson': simpson,
        'top_species': top_species_labels[0] if top_species_labels else 'NA',
        'rrna18_pct': round(rrna18/all_sequences*100, 3) if all_sequences else 0,
        'rrna58_pct': round(rrna58/all_sequences*100, 3) if all_sequences else 0,
        'rrna28_pct': round(rrna28/all_sequences*100, 3) if all_sequences else 0,
        'report_file': out_html.name,
    }


def build_dashboard(summary_rows, taxonomy_matrices, out_html: Path):
    total_reads = int(sum(r['raw_reads_num'] for r in summary_rows))
    med_otus = pd.Series([r['observed_otus'] for r in summary_rows]).median()
    med_shannon = pd.Series([r['shannon'] for r in summary_rows]).median()
    table_rows = ''.join(
        f"<tr><td><a href='{html.escape(r['report_file'])}' target='_blank' rel='noopener noreferrer'>{html.escape(r['sample'])}</a></td><td>{html.escape(str(r['raw_reads']))}</td><td>{html.escape(str(r['mean_read_length']))}</td><td>{html.escape(str(r['reads_gt_q20']))}</td><td>{r['observed_otus']}</td><td>{r['shannon']:.2f}</td><td>{r['simpson']:.3f}</td><td>{html.escape(str(r['top_species']))}</td></tr>"
        for r in summary_rows
    )
    rep = DASHBOARD_TEMPLATE
    replacements = {
        '__NSAMPLES__': str(len(summary_rows)),
        '__TOTALREADS__': str(total_reads),
        '__MEDOTUS__': str(round(float(med_otus), 1)),
        '__MEDSHANNON__': f'{float(med_shannon):.2f}',
        '__SUMMARY_ROWS__': table_rows,
        '__SUMMARY_DATA__': json.dumps(summary_rows),
        '__TAX_MATRICES__': json.dumps(taxonomy_matrices),
    }
    for k, v in replacements.items():
        rep = rep.replace(k, v)
    out_html.write_text(rep, encoding='utf-8')


def parse_options():
    ap = argparse.ArgumentParser(description='Generate BaNaNA biological HTML reports')
    ap.add_argument('--taxonomy', required=True, help='Final VSEARCH taxonomy.tsv')
    ap.add_argument('--otu-table', required=True, help='Final otu_table.tsv')
    ap.add_argument('--output', required=True, help='Main HTML report output')
    ap.add_argument('--output-dir', default='banana_reports', help='Directory for dashboard and per-sample reports')
    ap.add_argument('--nanoplot-dirs', nargs='*', default=[], help='NanoPlot output directories')
    ap.add_argument('--rrna-stats', nargs='*', default=[], help='rRNA extraction stats files')
    ap.add_argument('--top-n', type=int, default=15, help='Reserved for compatibility')
    return ap.parse_args()


def main():
    args = parse_options()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    otu, sample_names = load_otu_table(Path(args.otu_table))
    tax = load_taxonomy(Path(args.taxonomy))
    merged = otu.merge(tax, on='OTU', how='left')
    for col in ['Domain', 'Supergroup', 'Division', 'Subdivision', 'Class', 'Order', 'Family', 'Genus', 'Species']:
        if col not in merged.columns:
            merged[col] = 'Unassigned'
        merged[col] = merged[col].fillna('Unassigned')
    for col in ['Pident', 'Length']:
        if col not in merged.columns:
            merged[col] = pd.NA

    nanostats_by_sample = load_nanostats_by_sample(args.nanoplot_dirs)
    rrna_stats_by_sample = load_rrna_stats_by_sample(args.rrna_stats)
    summary = []
    tax_rank_data = {rank: {} for rank in ['Division', 'Class', 'Family', 'Genus', 'Species']}
    for sample in sample_names:
        if sample not in merged.columns:
            print(f'[WARN] Sample {sample} not found in otu_table.tsv; skipping')
            continue
        if sample not in nanostats_by_sample:
            print(f'[WARN] Missing NanoStats for {sample}; skipping')
            continue
        if sample not in rrna_stats_by_sample:
            print(f'[WARN] Missing rRNA stats for {sample}; skipping')
            continue
        all_sequences, rrna = rrna_stats_by_sample[sample]
        nano = nanostats_by_sample[sample]
        df = merged[['OTU', sample, 'Domain', 'Supergroup', 'Division', 'Subdivision', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Pident', 'Length']].copy()
        df[sample] = pd.to_numeric(df[sample], errors='coerce').fillna(0).astype(int)
        df = df[df[sample] > 0].sort_values(sample, ascending=False)
        out_html = outdir / f'{sample}_overview.html'
        row = make_report(sample, df, nano, all_sequences, rrna, out_html)
        summary.append(row)
        total_reads = int(df[sample].sum())
        for rank in tax_rank_data:
            grp = df.groupby(rank, dropna=False)[sample].sum().sort_values(ascending=False)
            tax_rank_data[rank][sample] = {str(k) if pd.notnull(k) else 'Unassigned': (float(v) / total_reads * 100 if total_reads else 0) for k, v in grp.items()}
        print(f'[OK] {sample} -> {out_html}')
    if summary:
        pd.DataFrame(summary).to_csv(outdir / 'sample_overview_summary.csv', index=False)
        print(f'[OK] Summary -> {outdir / "sample_overview_summary.csv"}')
        taxonomy_matrices = {}
        ordered_samples = [r['sample'] for r in summary]
        for rank, by_sample in tax_rank_data.items():
            totals = {}
            for mapping in by_sample.values():
                for taxon, pct in mapping.items():
                    totals[taxon] = totals.get(taxon, 0) + pct
            top_taxa = [k for k, _ in sorted(totals.items(), key=lambda kv: kv[1], reverse=True)[:20]]
            z = [[round(by_sample.get(sample, {}).get(taxon, 0), 3) for sample in ordered_samples] for taxon in top_taxa]
            taxonomy_matrices[rank] = {'samples': ordered_samples, 'taxa': top_taxa, 'z': z}
        build_dashboard(summary, taxonomy_matrices, outdir / 'multi_sample_dashboard.html')
        print(f'[OK] Dashboard -> {outdir / "multi_sample_dashboard.html"}')
        dashboard = outdir / 'multi_sample_dashboard.html'
        Path(args.output).write_text(dashboard.read_text(encoding='utf-8'), encoding='utf-8')
        print(f'[OK] Main report -> {args.output}')
    else:
        Path(args.output).write_text('<!doctype html><title>BaNaNA report</title><p>No sample reports were generated.</p>', encoding='utf-8')

if __name__ == '__main__':
    main()
