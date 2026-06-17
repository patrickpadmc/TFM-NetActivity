"""
process_jensen.py
- Carga human_tissue_knowledge_full.tsv
- Columnas: ENSP, gene_symbol, BTO_id, tissue_name, source, evidence_type, score
- Filtra score >= 4 (CURATED de alta confianza)
- Genera TSV procesado + dashboard HTML
"""
import csv, json, collections

IN_TSV   = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/raw/jensen_tissues/human_tissue_knowledge_full.tsv"
MAP_FILE = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed/entrez2ensg.tsv"
OUT_TSV  = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed/jensen_gene_tissue.tsv"
OUT_HTML = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed/dashboard_jensen.html"

SCORE_THRESHOLD = 4

# Jensen usa ENSP — necesitamos mapear a ENSG via ens2uniprot de Bioteque
# Pero ENSP -> ENSG: STRING aliases file o directamente desde Ensembl
# Para este script usamos ENSP como ID primario (mappeable externamente)
# y añadimos gene symbol que ya está en el fichero

print("Procesando Jensen TISSUES...")
rows = []
skipped_score = 0
all_evidence = collections.Counter()
score_values = []

with open(IN_TSV, newline="") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 7:
            continue
        ensp        = parts[0].strip()
        gene_symbol = parts[1].strip()
        bto_id      = parts[2].strip()
        tissue_name = parts[3].strip()
        source      = parts[4].strip()
        evidence    = parts[5].strip()
        score_str   = parts[6].strip()

        all_evidence[evidence] += 1

        try:
            score = int(score_str)
        except ValueError:
            continue

        score_values.append(score)

        if score < SCORE_THRESHOLD:
            skipped_score += 1
            continue

        rows.append({
            "ensp_id":     ensp,
            "gene_symbol": gene_symbol,
            "bto_id":      bto_id,
            "tissue_name": tissue_name,
            "source":      source,
            "evidence":    evidence,
            "score":       score_str,
        })

fieldnames = ["ensp_id","gene_symbol","bto_id","tissue_name","source","evidence","score"]
with open(OUT_TSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(rows)

print(f"  Total asociaciones filtradas (score >= {SCORE_THRESHOLD}): {len(rows):,}")
print(f"  Descartadas (score < {SCORE_THRESHOLD}): {skipped_score:,}")
print(f"TSV guardado: {OUT_TSV}")

# Stats
ev_counter   = collections.Counter(r["evidence"] for r in rows)
score_counter= collections.Counter(r["score"] for r in rows)
unique_genes = {r["gene_symbol"] for r in rows}
unique_tissues= {r["tissue_name"] for r in rows if r["tissue_name"] != r["bto_id"]}

top10_tissues = collections.Counter(r["tissue_name"] for r in rows if r["tissue_name"] != r["bto_id"]).most_common(10)
top10_html = "".join(
    f"<tr><td>{t}</td><td>{c:,}</td></tr>"
    for t, c in top10_tissues
)

score_labels = ["4","5"]
score_vals   = [score_counter.get(s,0) for s in score_labels]

html = f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="UTF-8">
<title>Jensen TISSUES — Dashboard</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
<link href="https://fonts.googleapis.com/css2?family=Raleway:wght@300;400;700;900&family=JetBrains+Mono:wght@400;600&display=swap" rel="stylesheet">
<style>
  :root{{--bg:#0a0f1e;--surface:#111827;--surface2:#1f2937;--accent:#6366f1;--green:#10b981;--yellow:#f59e0b;--text:#f1f5f9;--muted:#64748b;--border:#1f2937}}
  *{{box-sizing:border-box;margin:0;padding:0}}
  body{{font-family:'Raleway',sans-serif;background:var(--bg);color:var(--text);padding:36px}}
  .header h1{{font-size:2rem;font-weight:900;color:var(--accent);letter-spacing:-.02em}}
  .header .meta{{color:var(--muted);font-size:.82rem;margin-top:8px;font-family:'JetBrains Mono',monospace}}
  .filter-note{{background:rgba(99,102,241,.08);border-left:3px solid var(--accent);padding:10px 16px;font-size:.82rem;color:#a5b4fc;margin:20px 0}}
  .kpi-grid{{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin-bottom:26px}}
  .kpi{{background:var(--surface);border:1px solid var(--border);border-radius:6px;padding:20px;border-top:2px solid var(--accent)}}
  .kpi.green{{border-top-color:var(--green)}}
  .kpi.yellow{{border-top-color:var(--yellow)}}
  .kpi .val{{font-family:'JetBrains Mono',monospace;font-size:1.9rem;font-weight:600;color:var(--accent)}}
  .kpi.green .val{{color:var(--green)}}
  .kpi.yellow .val{{color:var(--yellow)}}
  .kpi .lbl{{font-size:.73rem;color:var(--muted);margin-top:5px;text-transform:uppercase;letter-spacing:.05em}}
  .charts-row{{display:grid;grid-template-columns:1fr 1fr;gap:18px;margin-bottom:22px}}
  .chart-wrap,.section{{background:var(--surface);border:1px solid var(--border);border-radius:6px;padding:22px;margin-bottom:22px}}
  h2{{font-size:.85rem;font-weight:700;color:var(--accent);margin-bottom:14px;text-transform:uppercase;letter-spacing:.06em}}
  table{{width:100%;border-collapse:collapse;font-size:.82rem}}
  th{{background:var(--surface2);color:var(--accent);padding:9px 12px;text-align:left;font-size:.72rem;text-transform:uppercase;border-bottom:1px solid var(--border)}}
  td{{padding:8px 12px;border-bottom:1px solid var(--border)}}
  tr:hover td{{background:var(--surface2)}}
  canvas{{max-height:250px}}
</style>
</head>
<body>
<div class="header">
  <h1>Jensen TISSUES — Curated</h1>
  <div class="meta">Actualización semanal (JensenLab / SIB Zurich) &nbsp;·&nbsp; Libre</div>
</div>
<div class="filter-note">
  <strong>Filtro:</strong> score ≥ {SCORE_THRESHOLD} (escala 1–5; 4=curado literatura, 5=curado + experimental).
  Solo asociaciones gen–tejido con respaldo bibliográfico directo.
</div>
<div class="kpi-grid">
  <div class="kpi"><div class="val">{len(rows):,}</div><div class="lbl">Asociaciones filtradas</div></div>
  <div class="kpi green"><div class="val">{len(unique_genes):,}</div><div class="lbl">Genes únicos (ENSP/Symbol)</div></div>
  <div class="kpi yellow"><div class="val">{len(unique_tissues):,}</div><div class="lbl">Tejidos únicos (BTO)</div></div>
  <div class="kpi"><div class="val">{score_counter.get("5",0):,}</div><div class="lbl">Score 5 (máx. confianza)</div></div>
</div>
<div class="charts-row">
  <div class="chart-wrap">
    <h2>Por tipo de evidencia</h2>
    <canvas id="evChart"></canvas>
  </div>
  <div class="chart-wrap">
    <h2>Por score (4 vs 5)</h2>
    <canvas id="scoreChart"></canvas>
  </div>
</div>
<div class="section">
  <h2>Top 10 tejidos — mayor número de genes asociados</h2>
  <table><tr><th>Tejido</th><th>Nº genes</th></tr>{top10_html}</table>
</div>
<script>
new Chart(document.getElementById('evChart'),{{
  type:'doughnut',
  data:{{
    labels:{json.dumps(list(ev_counter.keys()))},
    datasets:[{{data:{json.dumps(list(ev_counter.values()))},backgroundColor:['#6366f1','#10b981','#f59e0b','#e11d48'],borderWidth:0}}]
  }},
  options:{{plugins:{{legend:{{position:'bottom',labels:{{color:'#f1f5f9',font:{{size:11}}}}}}}}}}
}});
new Chart(document.getElementById('scoreChart'),{{
  type:'bar',
  data:{{
    labels:{json.dumps(score_labels)},
    datasets:[{{data:{json.dumps(score_vals)},backgroundColor:['#6366f1','#10b981'],borderRadius:4}}]
  }},
  options:{{plugins:{{legend:{{display:false}}}},scales:{{y:{{beginAtZero:true,ticks:{{color:'#64748b'}},grid:{{color:'#1f2937'}}}},x:{{ticks:{{color:'#64748b'}},grid:{{display:false}}}}}}}}
}});
</script>
</body></html>"""

with open(OUT_HTML, "w") as f:
    f.write(html)
print(f"Dashboard guardado: {OUT_HTML}")
