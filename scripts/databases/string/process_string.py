"""
process_string.py
- Carga 9606.protein.links.detailed.v12.0.txt.gz (streaming, no carga en RAM)
- Filtra combined_score >= 400
- Mapea ENSP -> ENSG via aliases file
- Genera TSV procesado + dashboard HTML
"""
import gzip, csv, json, collections, heapq

IN_LINKS   = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/raw/string/9606.protein.links.detailed.v12.0.txt.gz"
IN_ALIASES = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/raw/string/9606.protein.aliases.v12.0.txt.gz"
OUT_TSV    = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed/string_ppi.tsv"
OUT_HTML   = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed/dashboard_string.html"

SCORE_THRESHOLD = 400

# ── Cargar mapeo ENSP -> ENSG desde aliases ──────────────────────────────────
print("Cargando aliases ENSP -> ENSG...")
ensp2ensg = {}
with gzip.open(IN_ALIASES, "rt") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        ensp   = parts[0].replace("9606.", "")
        alias  = parts[1]
        source = parts[2]
        # Buscar aliases que sean IDs de Ensembl gene (ENSG)
        if "Ensembl_gene" in source or alias.startswith("ENSG"):
            if ensp not in ensp2ensg:
                ensp2ensg[ensp] = alias
print(f"  {len(ensp2ensg):,} ENSPs mapeados a ENSG")

# ── Stream links file ─────────────────────────────────────────────────────────
print(f"Procesando STRING (combined_score >= {SCORE_THRESHOLD})...")

n_total    = 0
n_filtered = 0
n_nomatch  = 0

score_bins = collections.Counter()
channel_sums = {"experimental": 0, "coexpression": 0, "database": 0, "textmining": 0}
top10_heap = []  # (score, row)

fieldnames = ["ensg_a","ensg_b","ensp_a","ensp_b",
              "neighborhood","fusion","cooccurence","coexpression",
              "experimental","database","textmining","combined_score"]

with open(OUT_TSV, "w", newline="") as out_f:
    writer = csv.DictWriter(out_f, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()

    with gzip.open(IN_LINKS, "rt") as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.rstrip("\n").split(" ")
            if len(parts) < 10:
                continue

            n_total += 1
            combined = int(parts[9])

            if combined < SCORE_THRESHOLD:
                continue

            ensp_a_full = parts[0]
            ensp_b_full = parts[1]
            ensp_a = ensp_a_full.replace("9606.", "")
            ensp_b = ensp_b_full.replace("9606.", "")

            ensg_a = ensp2ensg.get(ensp_a, "")
            ensg_b = ensp2ensg.get(ensp_b, "")

            if not ensg_a or not ensg_b:
                n_nomatch += 1

            row = {
                "ensg_a":        ensg_a,
                "ensg_b":        ensg_b,
                "ensp_a":        ensp_a,
                "ensp_b":        ensp_b,
                "neighborhood":  parts[2],
                "fusion":        parts[3],
                "cooccurence":   parts[4],
                "coexpression":  parts[5],
                "experimental":  parts[6],
                "database":      parts[7],
                "textmining":    parts[8],
                "combined_score":parts[9],
            }
            writer.writerow(row)
            n_filtered += 1

            # Score bins
            if combined < 500:   score_bins["400-499"] += 1
            elif combined < 700: score_bins["500-699"] += 1
            elif combined < 900: score_bins["700-899"] += 1
            else:                score_bins["900-1000"] += 1

            # Channel sums (para gráfica de fuentes de evidencia)
            for ch in channel_sums:
                channel_sums[ch] += int(parts[{"experimental":6,"coexpression":5,"database":7,"textmining":8}[ch]])

            # Top10 heap
            if len(top10_heap) < 10:
                heapq.heappush(top10_heap, (combined, n_filtered, row.copy()))
            elif combined > top10_heap[0][0]:
                heapq.heapreplace(top10_heap, (combined, n_filtered, row.copy()))

print(f"  Total pares en fichero: {n_total:,}")
print(f"  Pares filtrados (>= {SCORE_THRESHOLD}): {n_filtered:,}")
print(f"  Sin mapeo ENSG: {n_nomatch:,}")
print(f"TSV guardado: {OUT_TSV}")

# Dashboard
top10_sorted = sorted(top10_heap, key=lambda x: x[0], reverse=True)
top10_html = "".join(
    f"<tr><td class='mono'>{r['ensg_a'] or r['ensp_a']}</td>"
    f"<td class='mono'>{r['ensg_b'] or r['ensp_b']}</td>"
    f"<td>{r['experimental']}</td><td>{r['coexpression']}</td>"
    f"<td>{r['database']}</td><td>{r['textmining']}</td>"
    f"<td><strong>{r['combined_score']}</strong></td></tr>"
    for _, _, r in top10_sorted
)

bin_order = ["400-499","500-699","700-899","900-1000"]
bin_colors = ["#3b82f6","#8b5cf6","#06b6d4","#10b981"]

html = f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="UTF-8">
<title>STRING v12.0 — Dashboard</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
<link href="https://fonts.googleapis.com/css2?family=Space+Grotesk:wght@300;400;500;700&family=JetBrains+Mono:wght@400;600&display=swap" rel="stylesheet">
<style>
  :root{{--bg:#03050f;--surface:#080d1f;--surface2:#0f1629;--blue:#3b82f6;--purple:#8b5cf6;--cyan:#06b6d4;--green:#10b981;--text:#e2e8f0;--muted:#475569;--border:#1e293b}}
  *{{box-sizing:border-box;margin:0;padding:0}}
  body{{font-family:'Space Grotesk',sans-serif;background:var(--bg);color:var(--text);padding:36px}}
  .header h1{{font-size:2rem;font-weight:700;background:linear-gradient(90deg,var(--blue),var(--cyan));-webkit-background-clip:text;-webkit-text-fill-color:transparent}}
  .header .meta{{color:var(--muted);font-size:.82rem;margin-top:8px;font-family:'JetBrains Mono',monospace}}
  .filter-note{{background:rgba(59,130,246,.06);border-left:3px solid var(--blue);padding:10px 16px;font-size:.82rem;color:#93c5fd;margin:20px 0}}
  .kpi-grid{{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin-bottom:26px}}
  .kpi{{background:var(--surface);border:1px solid var(--border);border-radius:6px;padding:20px}}
  .kpi .val{{font-family:'JetBrains Mono',monospace;font-size:1.9rem;font-weight:600;background:linear-gradient(135deg,var(--blue),var(--cyan));-webkit-background-clip:text;-webkit-text-fill-color:transparent}}
  .kpi .lbl{{font-size:.73rem;color:var(--muted);margin-top:5px;text-transform:uppercase;letter-spacing:.05em}}
  .charts-row{{display:grid;grid-template-columns:1fr 1fr;gap:18px;margin-bottom:22px}}
  .chart-wrap,.section{{background:var(--surface);border:1px solid var(--border);border-radius:6px;padding:22px;margin-bottom:22px}}
  h2{{font-size:.85rem;font-weight:600;color:var(--blue);margin-bottom:14px;text-transform:uppercase;letter-spacing:.06em}}
  table{{width:100%;border-collapse:collapse;font-size:.79rem}}
  th{{background:var(--surface2);color:var(--cyan);padding:9px 12px;text-align:left;font-size:.71rem;text-transform:uppercase;border-bottom:1px solid var(--border)}}
  td{{padding:8px 12px;border-bottom:1px solid var(--border)}}
  tr:hover td{{background:var(--surface2)}}
  .mono{{font-family:'JetBrains Mono',monospace;font-size:.76rem}}
  canvas{{max-height:250px}}
</style>
</head>
<body>
<div class="header">
  <h1>STRING v12.0 — PPI Network</h1>
  <div class="meta">Homo sapiens (taxID 9606) &nbsp;·&nbsp; Libre (CC-BY)</div>
</div>
<div class="filter-note">
  <strong>Filtro:</strong> combined_score ≥ {SCORE_THRESHOLD} (escala 0–1000).
  Scores: 400–499 = confianza media; 500–699 = confianza alta; 700–899 = muy alta; 900–1000 = extrema.
  El combined_score integra experimental, coexpresión, co-ocurrencia, bases de datos curadas y text mining.
</div>
<div class="kpi-grid">
  <div class="kpi"><div class="val">{n_filtered:,}</div><div class="lbl">Pares PPI filtrados</div></div>
  <div class="kpi"><div class="val">{score_bins.get("900-1000",0):,}</div><div class="lbl">Score 900–1000 (extrema)</div></div>
  <div class="kpi"><div class="val">{score_bins.get("700-899",0):,}</div><div class="lbl">Score 700–899 (muy alta)</div></div>
  <div class="kpi"><div class="val">{n_nomatch:,}</div><div class="lbl">Sin mapeo ENSG (ENSP disponible)</div></div>
</div>
<div class="charts-row">
  <div class="chart-wrap">
    <h2>Distribución combined_score</h2>
    <canvas id="scoreChart"></canvas>
  </div>
  <div class="chart-wrap">
    <h2>Suma de scores por canal de evidencia</h2>
    <canvas id="channelChart"></canvas>
  </div>
</div>
<div class="section">
  <h2>Top 10 interacciones — mayor combined_score</h2>
  <table>
    <tr><th>ENSG A</th><th>ENSG B</th><th>Experimental</th><th>Coexpr.</th><th>Database</th><th>Textmining</th><th>Combined</th></tr>
    {top10_html}
  </table>
</div>
<script>
new Chart(document.getElementById('scoreChart'),{{
  type:'bar',
  data:{{
    labels:{json.dumps(bin_order)},
    datasets:[{{data:{json.dumps([score_bins.get(b,0) for b in bin_order])},backgroundColor:{json.dumps(bin_colors)},borderRadius:4}}]
  }},
  options:{{plugins:{{legend:{{display:false}}}},scales:{{y:{{beginAtZero:true,ticks:{{color:'#475569'}},grid:{{color:'#1e293b'}}}},x:{{ticks:{{color:'#94a3b8'}},grid:{{display:false}}}}}}}}
}});
new Chart(document.getElementById('channelChart'),{{
  type:'bar',
  data:{{
    labels:{json.dumps(list(channel_sums.keys()))},
    datasets:[{{data:{json.dumps(list(channel_sums.values()))},backgroundColor:['#3b82f6','#8b5cf6','#06b6d4','#10b981'],borderRadius:4}}]
  }},
  options:{{plugins:{{legend:{{display:false}}}},scales:{{y:{{beginAtZero:true,ticks:{{color:'#475569'}},grid:{{color:'#1e293b'}}}},x:{{ticks:{{color:'#94a3b8'}},grid:{{display:false}}}}}}}}
}});
</script>
</body></html>"""

with open(OUT_HTML, "w") as f:
    f.write(html)
print(f"Dashboard guardado: {OUT_HTML}")
