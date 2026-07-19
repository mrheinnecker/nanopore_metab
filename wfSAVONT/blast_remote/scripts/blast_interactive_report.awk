BEGIN {
    FS = "\t"
}

NR == 1 { next }

{
    rows[++row_count, "query"] = $1
    rows[row_count, "accession"] = $3
    rows[row_count, "taxid"] = $4
    rows[row_count, "identity"] = $7 + 0
    rows[row_count, "length"] = $8 + 0
    rows[row_count, "mismatches"] = $9 + 0
    rows[row_count, "gaps"] = $10 + 0
    rows[row_count, "coverage"] = $15 + 0
    rows[row_count, "evalue"] = $16
    rows[row_count, "bitscore"] = $17 + 0
    rows[row_count, "description"] = $18
}

function json(value, escaped) {
    escaped = value
    gsub(/\\/, "\\\\", escaped)
    gsub(/"/, "\\\"", escaped)
    gsub(/\r/, "", escaped)
    gsub(/\n/, "\\n", escaped)
    gsub(/</, "\\u003c", escaped)
    gsub(/>/, "\\u003e", escaped)
    return "\"" escaped "\""
}

END {
    print "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">"
    print "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">"
    print "<title>" sample " BLAST report</title>"
    print "<script src=\"https://cdn.plot.ly/plotly-2.35.2.min.js\"></script>"
    print "<style>"
    print ":root{--ink:#18302d;--muted:#66736c;--teal:#12685f;--gold:#c58a14;--paper:#f4f6f3;--card:#fff;--line:#dce4df}"
    print "*{box-sizing:border-box}body{margin:0;background:var(--paper);color:var(--ink);font-family:Arial,sans-serif}"
    print "header{background:linear-gradient(120deg,#0c4d47,#168176);color:#fff;padding:30px 5vw}header h1{margin:0 0 6px;font-size:30px}header p{margin:0;opacity:.86}"
    print "main{max-width:1500px;margin:0 auto;padding:24px}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:14px;margin-bottom:18px}"
    print ".card,.panel{background:var(--card);border:1px solid var(--line);border-radius:12px;box-shadow:0 3px 12px #18302d12}.card{padding:17px}.card b{display:block;font-size:27px;color:var(--teal);margin-top:6px}.card span{color:var(--muted);font-size:13px}"
    print ".panel{padding:18px;margin:18px 0}.panel h2{margin:0 0 4px;font-size:20px}.hint{color:var(--muted);font-size:13px;margin:0 0 10px}"
    print "#filter{width:min(520px,100%);padding:10px 12px;border:1px solid #b8c7c0;border-radius:8px;margin:8px 0 12px}"
    print ".table-wrap{overflow:auto;max-height:620px}table{border-collapse:collapse;width:100%;font-size:12px}th{position:sticky;top:0;background:#e8efec;color:var(--ink);text-align:left}th,td{padding:8px;border-bottom:1px solid var(--line);white-space:nowrap}td.desc{white-space:normal;min-width:420px}"
    print "</style></head><body>"
    print "<header><h1>" sample " · remote BLAST overview</h1><p>Interactive review of retained nucleotide BLAST hits</p></header>"
    print "<main><section class=\"cards\">"
    print "<div class=\"card\"><span>ASVs with hits</span><b id=\"asvCount\">–</b></div>"
    print "<div class=\"card\"><span>Retained hits</span><b id=\"hitCount\">–</b></div>"
    print "<div class=\"card\"><span>Best identity</span><b id=\"bestIdentity\">–</b></div>"
    print "<div class=\"card\"><span>Median coverage</span><b id=\"medianCoverage\">–</b></div>"
    print "</section>"
    print "<section class=\"panel\"><h2>Hit quality overview</h2><p class=\"hint\">High-quality hits tend toward the upper-right. Point size reflects alignment length; color reflects mismatches.</p><div id=\"scatter\" style=\"height:520px\"></div></section>"
    print "<section class=\"panel\"><h2>Identity by ASV and hit rank</h2><p class=\"hint\">Hover over a cell for complete alignment metrics and the subject description.</p><div id=\"heatmap\" style=\"height:480px\"></div></section>"
    print "<section class=\"panel\"><h2>BLAST hits</h2><input id=\"filter\" placeholder=\"Filter ASV, accession, taxon ID, or description…\"><div class=\"table-wrap\"><table><thead><tr>"
    print "<th>Query</th><th>Rank</th><th>Accession</th><th>Taxon ID</th><th>Identity %</th><th>Coverage %</th><th>Alignment</th><th>Mismatches</th><th>Gaps</th><th>E-value</th><th>Bit score</th><th>Description</th>"
    print "</tr></thead><tbody id=\"hitRows\"></tbody></table></div></section></main>"
    print "<script>const hits=["
    for (i = 1; i <= row_count; i++) {
        if (i > 1) printf ","
        printf("{q:%s,a:%s,t:%s,id:%g,len:%d,mm:%d,gaps:%d,cov:%g,e:%s,bit:%g,d:%s}", \
            json(rows[i, "query"]), json(rows[i, "accession"]), json(rows[i, "taxid"]), \
            rows[i, "identity"], rows[i, "length"], rows[i, "mismatches"], rows[i, "gaps"], \
            rows[i, "coverage"], json(rows[i, "evalue"]), rows[i, "bitscore"], json(rows[i, "description"]))
    }
    print "];"
    print "const ranks={};hits.forEach(h=>h.rank=(ranks[h.q]=(ranks[h.q]||0)+1));"
    print "const queries=[...new Set(hits.map(h=>h.q))],maxRank=Math.max(1,...hits.map(h=>h.rank));"
    print "const median=a=>{if(!a.length)return null;const b=[...a].sort((x,y)=>x-y),m=Math.floor(b.length/2);return b.length%2?b[m]:(b[m-1]+b[m])/2};"
    print "document.getElementById('asvCount').textContent=queries.length;document.getElementById('hitCount').textContent=hits.length;"
    print "document.getElementById('bestIdentity').textContent=hits.length?Math.max(...hits.map(h=>h.id)).toFixed(2)+'%':'–';const mc=median(hits.map(h=>h.cov));document.getElementById('medianCoverage').textContent=mc===null?'–':mc.toFixed(1)+'%';"
    print "const hover=h=>`<b>${h.q}</b><br>Hit ${h.rank}: ${h.a}<br>${h.d}<br>Identity: ${h.id}%<br>Query coverage: ${h.cov}%<br>Alignment length: ${h.len}<br>Mismatches: ${h.mm}<br>Gaps: ${h.gaps}<br>E-value: ${h.e}<br>Bit score: ${h.bit}<extra></extra>`;"
    print "Plotly.newPlot('scatter',[{x:hits.map(h=>h.cov),y:hits.map(h=>h.id),customdata:hits,mode:'markers',type:'scatter',marker:{size:hits.map(h=>Math.max(9,Math.min(25,h.len/90))),color:hits.map(h=>h.mm),colorscale:'YlOrRd',reversescale:true,showscale:true,colorbar:{title:'Mismatches'},line:{color:'#18302d',width:.5}},hovertemplate:hits.map(hover)}],{margin:{t:20},xaxis:{title:'Query coverage (%)',range:[0,102]},yaxis:{title:'Percent identity',range:[Math.max(0,Math.min(95,...hits.map(h=>h.id))-2),100.5]},hoverlabel:{align:'left'}},{responsive:true});"
    print "const z=queries.map(q=>Array.from({length:maxRank},(_,i)=>{const h=hits.find(x=>x.q===q&&x.rank===i+1);return h?h.id:null}));"
    print "const custom=queries.map(q=>Array.from({length:maxRank},(_,i)=>hits.find(x=>x.q===q&&x.rank===i+1)||null));"
    print "Plotly.newPlot('heatmap',[{z,x:Array.from({length:maxRank},(_,i)=>i+1),y:queries,customdata:custom,type:'heatmap',zmin:70,zmax:100,colorscale:[[0,'#d94f46'],[.5,'#d6a21f'],[1,'#12685f']],hoverongaps:false,hovertemplate:custom.map(row=>row.map(h=>h?hover(h):'<extra></extra>'))}],{margin:{t:15,l:Math.min(320,Math.max(100,...queries.map(q=>q.length*7))),b:55},xaxis:{title:'Hit rank',dtick:1},yaxis:{automargin:true},hoverlabel:{align:'left'}},{responsive:true});"
    print "const tbody=document.getElementById('hitRows');function render(term=''){term=term.toLowerCase();tbody.innerHTML=hits.filter(h=>!term||`${h.q} ${h.a} ${h.t} ${h.d}`.toLowerCase().includes(term)).map(h=>`<tr><td>${h.q}</td><td>${h.rank}</td><td>${h.a}</td><td>${h.t}</td><td>${h.id.toFixed(3)}</td><td>${h.cov}</td><td>${h.len}</td><td>${h.mm}</td><td>${h.gaps}</td><td>${h.e}</td><td>${h.bit}</td><td class=\"desc\">${h.d}</td></tr>`).join('')||'<tr><td colspan=\"12\">No matching hits</td></tr>'}render();document.getElementById('filter').addEventListener('input',e=>render(e.target.value));"
    print "</script></body></html>"
}
