"""
process_ctd.py  — versión memory-safe
Procesa CTD_genes_diseases.tsv.gz línea a línea sin cargar en RAM.
Escribe el TSV filtrado directamente a disco y acumula solo estadísticas.
"""
import gzip, csv, json, collections, heapq

# ── Paths ──────────────────────────────────────────────────────────────────
CTD_GZ   = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/raw/ctd/CTD_genes_diseases.tsv.gz"
MAP_FILE = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed/entrez2ensg.tsv"
OUT_TSV  = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed/ctd_gene_disease.tsv"
OUT_HTML = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed/dashboard_ctd.html"

INFER_THRESHOLD = 10.0

# ── Load mapping (small, ~19K rows) ─────────────────────────────────────────
print("Cargando mapeo Entrez->ENSG...")
entrez2ensg = {}
with open(MAP_FILE) as f:
    next(f)
    for line in f:
        p = line.strip().split("\t")
        if len(p) == 2:
            entrez2ensg[p[0]] = p[1]
print(f"  {len(entrez2ensg):,} entradas")

# ── Stream CTD, write TSV on the fly, accumulate only counters ───────────────
print("Procesando CTD linea a linea (memory-safe)...")

fieldnames = ["ensembl_gene_id","gene_symbol","entrez_gene_id","disease_name",
              "disease_id","evidence_type","direct_evidence","inference_score",
              "omim_ids","pubmed_ids"]

n_curated        = 0
n_inferred       = 0
n_skipped        = 0
genes_curated    = set()
genes_inferred   = set()
diseases_curated = set()
diseases_inferred= set()
de_counter       = collections.Counter()
score_bins       = collections.Counter()
top10_heap       = []

def pubmed_count(pubmed_str):
    return len(pubmed_str.split("|")) if pubmed_str.strip() else 0

with open(OUT_TSV, "w", newline="") as out_f:
    writer = csv.DictWriter(out_f, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()

    with gzip.open(CTD_GZ, "rt") as gz_f:
        for raw_line in gz_f:
            if raw_line.startswith("#"):
                continue
            parts = raw_line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            gene_symbol  = parts[0]
            gene_id      = parts[1]
            disease_name = parts[2]
            disease_id   = parts[3]
            direct_ev    = parts[4].strip()
            infer_score  = parts[6].strip()
            omim_ids     = parts[7].strip()
            pubmed_ids   = parts[8].strip()

            ensg = entrez2ensg.get(gene_id)
            if not ensg:
                n_skipped += 1
                continue

            if direct_ev:
                row = {
                    "ensembl_gene_id": ensg,
                    "gene_symbol":     gene_symbol,
                    "entrez_gene_id":  gene_id,
                    "disease_name":    disease_name,
                    "disease_id":      disease_id,
                    "evidence_type":   "curated",
                    "direct_evidence": direct_ev,
                    "inference_score": "",
                    "omim_ids":        omim_ids,
                    "pubmed_ids":      pubmed_ids,
                }
                writer.writerow(row)
                n_curated += 1
                genes_curated.add(ensg)
                diseases_curated.add(disease_id)
                de_counter[direct_ev] += 1

                pc = pubmed_count(pubmed_ids)
                if len(top10_heap) < 10:
                    heapq.heappush(top10_heap, (pc, n_curated, row.copy()))
                elif pc > top10_heap[0][0]:
                    heapq.heapreplace(top10_heap, (pc, n_curated, row.copy()))

            else:
                try:
                    score = float(infer_score)
                except ValueError:
                    continue
                if score < INFER_THRESHOLD:
                    continue

                row = {
                    "ensembl_gene_id": ensg,
                    "gene_symbol":     gene_symbol,
                    "entrez_gene_id":  gene_id,
                    "disease_name":    disease_name,
                    "disease_id":      disease_id,
                    "evidence_type":   "inferred",
                    "direct_evidence": "",
                    "inference_score": infer_score,
                    "omim_ids":        omim_ids,
                    "pubmed_ids":      pubmed_ids,
                }
                writer.writerow(row)
                n_inferred += 1
                genes_inferred.add(ensg)
                diseases_inferred.add(disease_id)

                if score < 20:   score_bins["10-20"] += 1
                elif score < 30: score_bins["20-30"] += 1
                elif score < 40: score_bins["30-40"] += 1
                elif score < 50: score_bins["40-50"] += 1
                else:            score_bins["50+"]   += 1

print(f"  Curadas:  {n_curated:,}")
print(f"  Inferidas (>={INFER_THRESHOLD}): {n_inferred:,}")
print(f"  Sin ENSG (descartadas): {n_skipped:,}")
print(f"TSV guardado: {OUT_TSV}")

# ── Build dashboard ──────────────────────────────────────────────────────────
print("Generando dashboard HTML...")

top10_sorted = sorted(top10_heap, key=lambda x: x[0], reverse=True)

top10_rows_html = ""
for pc, _, r in top10_sorted:
    badge_cls = r['direct_evidence'].replace('/','-').replace(' ','-').lower()
    top10_rows_html += (
        f"<tr><td>{r['ensembl_gene_id']}</td><td>{r['gene_symbol']}</td>"
        f"<td>{r['disease_name']}</td><td>{r['disease_id']}</td>"
        f"<td><span class='badge badge-{badge_cls}'>{r['direct_evidence']}</span></td>"
        f"<td>{pc}</td></tr>\n"
    )

bin_order = ["10-20","20-30","30-40","40-50","50+"]
de_labels = json.dumps(list(de_counter.keys()))
de_values = json.dumps(list(de_counter.values()))
sb_labels = json.dumps(bin_order)
sb_values = json.dumps([score_bins[b] for b in bin_order])

html = f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="UTF-8">
<title>CTD — Dashboard</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
<link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans:wght@300;400;600&display=swap" rel="stylesheet">
<style>
  :root{{
    --bg:#0f1117;--surface:#1a1d27;--surface2:#22263a;
    --accent:#00e5ff;--accent2:#ff6b35;--accent3:#7fff6b;
    --text:#e8eaf0;--muted:#6b7280;--border:#2d3148;
  }}
  *{{box-sizing:border-box;margin:0;padding:0}}
  body{{font-family:'IBM Plex Sans',sans-serif;background:var(--bg);color:var(--text);padding:32px;min-height:100vh}}
  .header{{margin-bottom:40px;border-left:4px solid var(--accent);padding-left:20px}}
  .header h1{{font-family:'IBM Plex Mono',monospace;font-size:1.8rem;color:var(--accent);letter-spacing:-.02em}}
  .header .meta{{color:var(--muted);font-size:.85rem;margin-top:6px;font-family:'IBM Plex Mono',monospace}}
  .kpi-grid{{display:grid;grid-template-columns:repeat(4,1fr);gap:16px;margin-bottom:32px}}
  .kpi{{background:var(--surface);border:1px solid var(--border);border-radius:6px;padding:20px;position:relative;overflow:hidden}}
  .kpi::before{{content:'';position:absolute;top:0;left:0;right:0;height:2px;background:var(--accent)}}
  .kpi.orange::before{{background:var(--accent2)}}
  .kpi.green::before{{background:var(--accent3)}}
  .kpi .val{{font-family:'IBM Plex Mono',monospace;font-size:2rem;font-weight:600;color:var(--accent)}}
  .kpi.orange .val{{color:var(--accent2)}}
  .kpi.green .val{{color:var(--accent3)}}
  .kpi .lbl{{font-size:.78rem;color:var(--muted);margin-top:6px;text-transform:uppercase;letter-spacing:.06em}}
  .section{{background:var(--surface);border:1px solid var(--border);border-radius:6px;padding:24px;margin-bottom:24px}}
  .section h2{{font-family:'IBM Plex Mono',monospace;font-size:1rem;color:var(--accent);margin-bottom:18px;text-transform:uppercase;letter-spacing:.08em}}
  .charts-row{{display:grid;grid-template-columns:1fr 1fr;gap:20px;margin-bottom:24px}}
  .chart-wrap{{background:var(--surface);border:1px solid var(--border);border-radius:6px;padding:20px}}
  .chart-wrap h2{{font-family:'IBM Plex Mono',monospace;font-size:.9rem;color:var(--accent);margin-bottom:14px;text-transform:uppercase;letter-spacing:.08em}}
  table{{width:100%;border-collapse:collapse;font-size:.83rem}}
  th{{background:var(--surface2);color:var(--accent);padding:10px 12px;text-align:left;font-family:'IBM Plex Mono',monospace;font-size:.75rem;text-transform:uppercase;letter-spacing:.06em;border-bottom:1px solid var(--border)}}
  td{{padding:9px 12px;border-bottom:1px solid var(--border);color:var(--text)}}
  tr:hover td{{background:var(--surface2)}}
  .badge{{display:inline-block;padding:2px 8px;border-radius:3px;font-size:.72rem;font-family:'IBM Plex Mono',monospace;font-weight:600}}
  .badge-marker-mechanism{{background:rgba(0,229,255,.12);color:var(--accent);border:1px solid var(--accent)}}
  .badge-therapeutic{{background:rgba(127,255,107,.12);color:var(--accent3);border:1px solid var(--accent3)}}
  .threshold-note{{background:var(--surface2);border-left:3px solid var(--accent2);padding:12px 16px;border-radius:0 4px 4px 0;font-size:.82rem;color:var(--muted);margin-bottom:24px}}
  .threshold-note strong{{color:var(--accent2)}}
  canvas{{max-height:280px}}
</style>
</head>
<body>
<div class="header">
  <h1>CTD — Comparative Toxicogenomics Database</h1>
  <div class="meta">Reporte: 31 Mar 2026 &nbsp;·&nbsp; Homo sapiens &nbsp;·&nbsp; Filtro: DirectEvidence != vacio | InferenceScore >= {INFER_THRESHOLD}</div>
</div>
<div class="threshold-note">
  <strong>Criterio de evidencia solida:</strong>
  Asociaciones <strong>curadas</strong> (DirectEvidence = marker/mechanism o therapeutic) +
  asociaciones <strong>inferidas</strong> con InferenceScore >= {INFER_THRESHOLD}
  (numero de publicaciones que conectan el quimico intermediario con el gen y la enfermedad).
</div>
<div class="kpi-grid">
  <div class="kpi"><div class="val">{n_curated + n_inferred:,}</div><div class="lbl">Total asociaciones filtradas</div></div>
  <div class="kpi orange"><div class="val">{n_curated:,}</div><div class="lbl">Curadas (DirectEvidence)</div></div>
  <div class="kpi green"><div class="val">{n_inferred:,}</div><div class="lbl">Inferidas (score >= {INFER_THRESHOLD})</div></div>
  <div class="kpi"><div class="val">{len(genes_curated | genes_inferred):,}</div><div class="lbl">Genes unicos (ENSG)</div></div>
</div>
<div class="kpi-grid">
  <div class="kpi"><div class="val">{len(diseases_curated | diseases_inferred):,}</div><div class="lbl">Enfermedades unicas (MeSH)</div></div>
  <div class="kpi orange"><div class="val">{len(genes_curated):,}</div><div class="lbl">Genes en curadas</div></div>
  <div class="kpi green"><div class="val">{len(diseases_curated):,}</div><div class="lbl">Enfermedades en curadas</div></div>
  <div class="kpi"><div class="val">{len(diseases_inferred):,}</div><div class="lbl">Enfermedades en inferidas</div></div>
</div>
<div class="charts-row">
  <div class="chart-wrap">
    <h2>Tipo de evidencia directa (curadas)</h2>
    <canvas id="deChart"></canvas>
  </div>
  <div class="chart-wrap">
    <h2>Distribucion InferenceScore (inferidas, score >= {INFER_THRESHOLD})</h2>
    <canvas id="scoreChart"></canvas>
  </div>
</div>
<div class="section">
  <h2>Top 10 asociaciones curadas — mayor respaldo bibliografico (no PubMed IDs)</h2>
  <table>
    <tr><th>ENSG</th><th>Gen</th><th>Enfermedad</th><th>MeSH ID</th><th>DirectEvidence</th><th>No PubMed</th></tr>
    {top10_rows_html}
  </table>
</div>
<script>
new Chart(document.getElementById('deChart'),{{
  type:'doughnut',
  data:{{labels:{de_labels},datasets:[{{data:{de_values},backgroundColor:['#00e5ff','#7fff6b','#ff6b35','#a78bfa'],borderColor:'#1a1d27',borderWidth:2}}]}},
  options:{{plugins:{{legend:{{position:'bottom',labels:{{color:'#e8eaf0',font:{{family:'IBM Plex Mono',size:11}}}}}}}}}}
}});
new Chart(document.getElementById('scoreChart'),{{
  type:'bar',
  data:{{labels:{sb_labels},datasets:[{{label:'Asociaciones inferidas',data:{sb_values},backgroundColor:'#ff6b35',borderRadius:3}}]}},
  options:{{plugins:{{legend:{{display:false}}}},scales:{{y:{{beginAtZero:true,ticks:{{color:'#6b7280'}},grid:{{color:'#2d3148'}}}},x:{{ticks:{{color:'#6b7280'}},grid:{{display:false}}}}}}}}
}});
</script>
</body></html>"""

with open(OUT_HTML, "w") as f:
    f.write(html)
print(f"Dashboard guardado: {OUT_HTML}")
