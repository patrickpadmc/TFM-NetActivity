"""
process_hpa.py
- Carga proteinatlas.tsv (formato ancho, 107 columnas)
- Filtra por RNA tissue specificity y Tau score >= 0.5
- Extrae: ENSG, gene_symbol, rna_tissue_specificity, rna_tissue_distribution, tau_score
- Genera TSV procesado + dashboard HTML
"""
import csv, json, collections

IN_TSV   = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/raw/hpa/proteinatlas.tsv"
OUT_TSV  = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed/hpa_gene_tissue.tsv"
OUT_HTML = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed/dashboard_hpa.html"

HPA_VERSION = "v25.0 (2025-11-11)"
TAU_THRESHOLD = 0.5
KEEP_SPECIFICITY = {"Tissue enriched", "Group enriched", "Tissue enhanced"}

print("Procesando HPA...")
rows = []
skipped_spec = 0
skipped_tau  = 0
skipped_ensg = 0
all_specificity = collections.Counter()

with open(IN_TSV, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        ensg  = row.get("Ensembl", "").strip()
        gene  = row.get("Gene", "").strip()
        spec  = row.get("RNA tissue specificity", "").strip().strip('"')
        dist  = row.get("RNA tissue distribution", "").strip().strip('"')
        tau_s = row.get("RNA tissue specificity score", "").strip().strip('"')
        ntpm  = row.get("RNA tissue specific nTPM", "").strip().strip('"')

        all_specificity[spec] += 1

        if not ensg:
            skipped_ensg += 1
            continue
        if spec not in KEEP_SPECIFICITY:
            skipped_spec += 1
            continue
        try:
            tau = float(tau_s)
        except ValueError:
            skipped_tau += 1
            continue
        if tau < TAU_THRESHOLD:
            skipped_tau += 1
            continue

        rows.append({
            "ensembl_gene_id":         ensg,
            "gene_symbol":             gene,
            "rna_tissue_specificity":  spec,
            "rna_tissue_distribution": dist,
            "tau_score":               tau_s,
            "rna_tissue_specific_nTPM":ntpm,
        })

fieldnames = ["ensembl_gene_id","gene_symbol","rna_tissue_specificity",
              "rna_tissue_distribution","tau_score","rna_tissue_specific_nTPM"]
with open(OUT_TSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(rows)

print(f"  Asociaciones filtradas: {len(rows):,}")
print(f"  Descartadas (sin ENSG): {skipped_ensg:,}")
print(f"  Descartadas (especificidad): {skipped_spec:,}")
print(f"  Descartadas (Tau < {TAU_THRESHOLD}): {skipped_tau:,}")
print(f"TSV guardado: {OUT_TSV}")

# Stats
spec_counter = collections.Counter(r["rna_tissue_specificity"] for r in rows)
unique_genes = {r["ensembl_gene_id"] for r in rows}

tau_bins = collections.Counter()
for r in rows:
    t = float(r["tau_score"])
    if t < 0.6:   tau_bins["0.5-0.6"] += 1
    elif t < 0.7: tau_bins["0.6-0.7"] += 1
    elif t < 0.8: tau_bins["0.7-0.8"] += 1
    elif t < 0.9: tau_bins["0.8-0.9"] += 1
    else:         tau_bins["0.9-1.0"] += 1

top10 = sorted(rows, key=lambda r: float(r["tau_score"]), reverse=True)[:10]
top10_html = "".join(
    f"<tr><td>{r['ensembl_gene_id']}</td><td>{r['gene_symbol']}</td>"
    f"<td><span class='badge badge-{r['rna_tissue_specificity'].lower().replace(' ','-')}'>"
    f"{r['rna_tissue_specificity']}</span></td>"
    f"<td>{r['tau_score']}</td><td>{r['rna_tissue_specific_nTPM']}</td></tr>"
    for r in top10
)

bin_order = ["0.5-0.6","0.6-0.7","0.7-0.8","0.8-0.9","0.9-1.0"]

html = f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="UTF-8">
<title>HPA — Dashboard</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
<link href="https://fonts.googleapis.com/css2?family=Plus+Jakarta+Sans:wght@300;400;600;800&family=JetBrains+Mono:wght@400;600&display=swap" rel="stylesheet">
<style>
  :root{{--bg:#f7f7f2;--surface:#fff;--surface2:#f0f0e8;--navy:#1a1f3c;--teal:#00897b;--amber:#f59e0b;--rose:#e11d48;--text:#1a1f3c;--muted:#6b7280;--border:#e5e5dc}}
  *{{box-sizing:border-box;margin:0;padding:0}}
  body{{font-family:'Plus Jakarta Sans',sans-serif;background:var(--bg);color:var(--text);padding:36px}}
  .header{{margin-bottom:36px}}
  .header h1{{font-size:2rem;font-weight:800;color:var(--navy)}}
  .header h1 span{{color:var(--teal)}}
  .header .meta{{color:var(--muted);font-size:.82rem;margin-top:8px}}
  .filter-note{{background:#ecfdf5;border-left:3px solid var(--teal);padding:10px 16px;border-radius:0 6px 6px 0;font-size:.82rem;color:#065f46;margin-bottom:28px}}
  .kpi-grid{{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin-bottom:26px}}
  .kpi{{background:var(--surface);border:1px solid var(--border);border-radius:8px;padding:20px;border-top:3px solid var(--teal)}}
  .kpi.amber{{border-top-color:var(--amber)}}
  .kpi.rose{{border-top-color:var(--rose)}}
  .kpi .val{{font-family:'JetBrains Mono',monospace;font-size:1.9rem;font-weight:600;color:var(--navy)}}
  .kpi .lbl{{font-size:.73rem;color:var(--muted);margin-top:5px;text-transform:uppercase;letter-spacing:.05em}}
  .charts-row{{display:grid;grid-template-columns:1fr 1fr;gap:18px;margin-bottom:22px}}
  .chart-wrap{{background:var(--surface);border:1px solid var(--border);border-radius:8px;padding:22px}}
  .chart-wrap h2{{font-size:.85rem;font-weight:700;color:var(--navy);margin-bottom:14px;text-transform:uppercase;letter-spacing:.05em}}
  .section{{background:var(--surface);border:1px solid var(--border);border-radius:8px;padding:22px;margin-bottom:22px}}
  .section h2{{font-size:.85rem;font-weight:700;color:var(--navy);margin-bottom:14px;text-transform:uppercase;letter-spacing:.05em}}
  table{{width:100%;border-collapse:collapse;font-size:.82rem}}
  th{{background:var(--surface2);color:var(--navy);padding:9px 12px;text-align:left;font-size:.72rem;text-transform:uppercase;letter-spacing:.05em;border-bottom:2px solid var(--border)}}
  td{{padding:8px 12px;border-bottom:1px solid var(--border)}}
  tr:hover td{{background:var(--surface2)}}
  .badge{{display:inline-block;padding:2px 8px;border-radius:20px;font-size:.7rem;font-weight:600}}
  .badge-tissue-enriched{{background:#d1fae5;color:#065f46}}
  .badge-group-enriched{{background:#dbeafe;color:#1e40af}}
  .badge-tissue-enhanced{{background:#fef3c7;color:#92400e}}
  canvas{{max-height:250px}}
</style>
</head>
<body>
<div class="header">
  <h1>Human Protein Atlas <span>HPA</span></h1>
  <div class="meta">{HPA_VERSION} &nbsp;·&nbsp; Fichero: proteinatlas.tsv &nbsp;·&nbsp; Libre (CC-BY)</div>
</div>
<div class="filter-note">
  <strong>Filtro aplicado:</strong> RNA tissue specificity ∈ {{Tissue enriched, Group enriched, Tissue enhanced}} + Tau score ≥ {TAU_THRESHOLD}.
  El Tau score (0–1) mide la especificidad tejido-única del gen: 1 = expresión exclusiva en un tejido, 0 = expresión ubicua.
</div>
<div class="kpi-grid">
  <div class="kpi"><div class="val">{len(rows):,}</div><div class="lbl">Genes filtrados (Tau ≥ {TAU_THRESHOLD})</div></div>
  <div class="kpi amber"><div class="val">{len(unique_genes):,}</div><div class="lbl">Genes únicos (ENSG)</div></div>
  <div class="kpi rose"><div class="val">{spec_counter.get("Tissue enriched",0):,}</div><div class="lbl">Tissue enriched</div></div>
  <div class="kpi"><div class="val">{spec_counter.get("Group enriched",0):,}</div><div class="lbl">Group enriched</div></div>
</div>
<div class="charts-row">
  <div class="chart-wrap">
    <h2>Por categoría de especificidad</h2>
    <canvas id="specChart"></canvas>
  </div>
  <div class="chart-wrap">
    <h2>Distribución Tau score (genes filtrados)</h2>
    <canvas id="tauChart"></canvas>
  </div>
</div>
<div class="section">
  <h2>Top 10 genes — mayor Tau score (mayor especificidad tisular)</h2>
  <table>
    <tr><th>ENSG</th><th>Gen</th><th>Especificidad</th><th>Tau score</th><th>nTPM tejido específico</th></tr>
    {top10_html}
  </table>
</div>
<script>
new Chart(document.getElementById('specChart'),{{
  type:'doughnut',
  data:{{
    labels:{json.dumps(list(spec_counter.keys()))},
    datasets:[{{data:{json.dumps(list(spec_counter.values()))},backgroundColor:['#00897b','#1e40af','#f59e0b'],borderWidth:0}}]
  }},
  options:{{plugins:{{legend:{{position:'bottom',labels:{{color:'#1a1f3c',font:{{size:11}}}}}}}}}}
}});
new Chart(document.getElementById('tauChart'),{{
  type:'bar',
  data:{{
    labels:{json.dumps(bin_order)},
    datasets:[{{data:{json.dumps([tau_bins[b] for b in bin_order])},backgroundColor:'#00897b',borderRadius:4}}]
  }},
  options:{{plugins:{{legend:{{display:false}}}},scales:{{y:{{beginAtZero:true,ticks:{{color:'#6b7280'}},grid:{{color:'#f0f0e8'}}}},x:{{ticks:{{color:'#6b7280'}},grid:{{display:false}}}}}}}}
}});
</script>
</body></html>"""

with open(OUT_HTML, "w") as f:
    f.write(html)
print(f"Dashboard guardado: {OUT_HTML}")
