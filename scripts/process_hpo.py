"""
process_hpo.py
- Loads genes_to_disease.txt  (gene->disease associations)
- Loads phenotype.hpoa         (disease->HPO phenotype annotations, for evidence codes)
- Maps Entrez Gene ID -> ENSG
- Filters by association_type (MENDELIAN > POLYGENIC) and evidence code
- Computes stats
- Generates dashboard HTML
"""
import csv, json, collections

# ── Paths ──────────────────────────────────────────────────────────────────
G2D_FILE  = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/raw/hpo/genes_to_disease.txt"
HPOA_FILE = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/raw/hpo/phenotype.hpoa"
MAP_FILE  = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed/entrez2ensg.tsv"
OUT_TSV   = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed/hpo_gene_disease.tsv"
OUT_HTML  = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed/dashboard_hpo.html"

HPO_VERSION = "2026-02-16"

# Evidence code hierarchy (higher = more reliable)
EVIDENCE_RANK = {"PCS": 3, "TAS": 2, "IEA": 1}
EVIDENCE_LABEL = {
    "PCS": "PCS — Publicación con casos clínicos (alta confianza)",
    "TAS": "TAS — Statement de experto (confianza media)",
    "IEA": "IEA — Inferido electrónicamente (menor confianza)",
}

# ── Load mapping ────────────────────────────────────────────────────────────
print("Cargando mapeo Entrez->ENSG...")
entrez2ensg = {}
with open(MAP_FILE) as f:
    next(f)
    for line in f:
        p = line.strip().split("\t")
        if len(p) == 2:
            entrez2ensg[p[0]] = p[1]
print(f"  {len(entrez2ensg):,} entradas")

# ── Load HPOA: disease -> best evidence code ─────────────────────────────────
print("Cargando phenotype.hpoa para evidencia...")
disease_best_evidence = {}   # disease_id -> best evidence code seen
disease_hpo_count     = {}   # disease_id -> number of HPO annotations

with open(HPOA_FILE) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 6:
            continue
        if parts[0] == "database_id":
            continue
        disease_id = parts[0]
        evidence   = parts[5].strip()
        rank       = EVIDENCE_RANK.get(evidence, 0)
        prev_rank  = EVIDENCE_RANK.get(disease_best_evidence.get(disease_id, ""), 0)
        if rank > prev_rank:
            disease_best_evidence[disease_id] = evidence
        disease_hpo_count[disease_id] = disease_hpo_count.get(disease_id, 0) + 1

print(f"  Enfermedades con evidencia: {len(disease_best_evidence):,}")

# ── Parse genes_to_disease.txt ───────────────────────────────────────────────
print("Procesando genes_to_disease.txt...")
rows = []
skipped = 0

with open(G2D_FILE) as f:
    for line in f:
        if line.startswith("#") or line.startswith("ncbi_gene_id"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 5:
            continue
        entrez_raw   = parts[0].replace("NCBIGene:", "").strip()
        gene_symbol  = parts[1].strip()
        assoc_type   = parts[2].strip()
        disease_id   = parts[3].strip()
        source       = parts[4].strip()

        ensg = entrez2ensg.get(entrez_raw)
        if not ensg:
            skipped += 1
            continue

        best_ev     = disease_best_evidence.get(disease_id, "IEA")
        hpo_count   = disease_hpo_count.get(disease_id, 0)

        rows.append({
            "ensembl_gene_id": ensg,
            "gene_symbol":     gene_symbol,
            "entrez_gene_id":  entrez_raw,
            "disease_id":      disease_id,
            "association_type":assoc_type,
            "best_evidence":   best_ev,
            "hpo_term_count":  hpo_count,
            "source":          source,
        })

print(f"  Total asociaciones mapeadas: {len(rows):,}")
print(f"  Sin mapeo ENSG (descartadas): {skipped:,}")

# ── Write TSV ────────────────────────────────────────────────────────────────
fieldnames = ["ensembl_gene_id","gene_symbol","entrez_gene_id","disease_id",
              "association_type","best_evidence","hpo_term_count","source"]
with open(OUT_TSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(rows)
print(f"TSV guardado: {OUT_TSV}")

# ── Stats ─────────────────────────────────────────────────────────────────────
type_counter = collections.Counter(r["association_type"] for r in rows)
ev_counter   = collections.Counter(r["best_evidence"]    for r in rows)

mendelian    = [r for r in rows if r["association_type"] == "MENDELIAN"]
polygenic    = [r for r in rows if r["association_type"] == "POLYGENIC"]
pcs_rows     = [r for r in rows if r["best_evidence"] == "PCS"]

unique_genes_total    = {r["ensembl_gene_id"] for r in rows}
unique_diseases_total = {r["disease_id"]      for r in rows}
unique_genes_mend     = {r["ensembl_gene_id"] for r in mendelian}
unique_genes_pcs      = {r["ensembl_gene_id"] for r in pcs_rows}

# HPO term count distribution for mendelian PCS
hpo_bins = collections.Counter()
for r in mendelian:
    c = r["hpo_term_count"]
    if c == 0:       hpo_bins["0"] += 1
    elif c <= 5:     hpo_bins["1-5"] += 1
    elif c <= 15:    hpo_bins["6-15"] += 1
    elif c <= 30:    hpo_bins["16-30"] += 1
    else:            hpo_bins["30+"] += 1
hpo_bin_order = ["0","1-5","6-15","16-30","30+"]

# Top 10 MENDELIAN + PCS by HPO term count
top10 = sorted(
    [r for r in rows if r["association_type"]=="MENDELIAN" and r["best_evidence"]=="PCS"],
    key=lambda r: r["hpo_term_count"], reverse=True
)[:10]

top10_html = ""
for r in top10:
    top10_html += (
        f"<tr><td>{r['ensembl_gene_id']}</td><td>{r['gene_symbol']}</td>"
        f"<td>{r['disease_id']}</td>"
        f"<td><span class='badge badge-{r['association_type'].lower()}'>{r['association_type']}</span></td>"
        f"<td><span class='badge badge-{r['best_evidence'].lower()}'>{r['best_evidence']}</span></td>"
        f"<td>{r['hpo_term_count']}</td></tr>\n"
    )

# ── HTML ──────────────────────────────────────────────────────────────────────
type_labels  = json.dumps(list(type_counter.keys()))
type_values  = json.dumps(list(type_counter.values()))
ev_labels    = json.dumps([k for k in ["PCS","TAS","IEA"] if k in ev_counter])
ev_values    = json.dumps([ev_counter.get(k,0) for k in ["PCS","TAS","IEA"] if k in ev_counter])
hpo_labels   = json.dumps(hpo_bin_order)
hpo_values   = json.dumps([hpo_bins[b] for b in hpo_bin_order])

html = f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="UTF-8">
<title>HPO — Dashboard</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
<link href="https://fonts.googleapis.com/css2?family=Syne:wght@400;700;800&family=DM+Sans:wght@300;400;500&display=swap" rel="stylesheet">
<style>
  :root{{
    --bg:#fafaf8;--surface:#ffffff;--surface2:#f3f4f0;
    --accent:#1a1a2e;--cyan:#0a9396;--orange:#ee9b00;--red:#ae2012;
    --text:#1a1a2e;--muted:#6b7280;--border:#e5e7eb;
  }}
  *{{box-sizing:border-box;margin:0;padding:0}}
  body{{font-family:'DM Sans',sans-serif;background:var(--bg);color:var(--text);padding:36px}}
  .header{{margin-bottom:36px}}
  .header h1{{font-family:'Syne',sans-serif;font-size:2rem;font-weight:800;color:var(--accent);line-height:1.1}}
  .header h1 span{{color:var(--cyan)}}
  .header .meta{{color:var(--muted);font-size:.82rem;margin-top:8px}}
  .kpi-grid{{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin-bottom:28px}}
  .kpi{{background:var(--surface);border:1px solid var(--border);border-radius:8px;padding:20px;border-top:3px solid var(--cyan)}}
  .kpi.orange{{border-top-color:var(--orange)}}
  .kpi.red{{border-top-color:var(--red)}}
  .kpi.dark{{border-top-color:var(--accent)}}
  .kpi .val{{font-family:'Syne',sans-serif;font-size:1.9rem;font-weight:700;color:var(--accent)}}
  .kpi .lbl{{font-size:.75rem;color:var(--muted);margin-top:5px;text-transform:uppercase;letter-spacing:.05em}}
  .section{{background:var(--surface);border:1px solid var(--border);border-radius:8px;padding:24px;margin-bottom:22px}}
  .section h2{{font-family:'Syne',sans-serif;font-size:.95rem;font-weight:700;color:var(--accent);margin-bottom:16px;text-transform:uppercase;letter-spacing:.06em}}
  .charts-row{{display:grid;grid-template-columns:1fr 1fr 1fr;gap:18px;margin-bottom:22px}}
  .chart-wrap{{background:var(--surface);border:1px solid var(--border);border-radius:8px;padding:20px}}
  .chart-wrap h2{{font-family:'Syne',sans-serif;font-size:.85rem;font-weight:700;color:var(--accent);margin-bottom:12px;text-transform:uppercase;letter-spacing:.05em}}
  table{{width:100%;border-collapse:collapse;font-size:.82rem}}
  th{{background:var(--surface2);color:var(--accent);padding:9px 12px;text-align:left;font-family:'Syne',sans-serif;font-size:.72rem;text-transform:uppercase;letter-spacing:.05em;border-bottom:2px solid var(--border)}}
  td{{padding:8px 12px;border-bottom:1px solid var(--border)}}
  tr:hover td{{background:var(--surface2)}}
  .badge{{display:inline-block;padding:2px 8px;border-radius:20px;font-size:.7rem;font-weight:500}}
  .badge-mendelian{{background:#dcfce7;color:#166534}}
  .badge-polygenic{{background:#fef3c7;color:#92400e}}
  .badge-pcs{{background:#dbeafe;color:#1e40af}}
  .badge-tas{{background:#f3e8ff;color:#6b21a8}}
  .badge-iea{{background:#f1f5f9;color:#475569}}
  .ev-legend{{display:grid;grid-template-columns:1fr;gap:8px;margin-top:14px}}
  .ev-item{{padding:10px 14px;border-radius:6px;font-size:.8rem;border-left:3px solid}}
  .ev-pcs{{background:#dbeafe;border-color:#1e40af;color:#1e40af}}
  .ev-tas{{background:#f3e8ff;border-color:#6b21a8;color:#6b21a8}}
  .ev-iea{{background:#f1f5f9;border-color:#475569;color:#475569}}
  canvas{{max-height:240px}}
</style>
</head>
<body>

<div class="header">
  <h1>Human Phenotype Ontology <span>HPO</span></h1>
  <div class="meta">Versión: {HPO_VERSION} &nbsp;·&nbsp; Fuente: genes_to_disease.txt + phenotype.hpoa &nbsp;·&nbsp; Filtro: MENDELIAN + PCS como evidencia sólida</div>
</div>

<div class="kpi-grid">
  <div class="kpi"><div class="val">{len(rows):,}</div><div class="lbl">Asociaciones totales</div></div>
  <div class="kpi dark"><div class="val">{len(mendelian):,}</div><div class="lbl">MENDELIAN</div></div>
  <div class="kpi orange"><div class="val">{len(polygenic):,}</div><div class="lbl">POLYGENIC</div></div>
  <div class="kpi red"><div class="val">{len(pcs_rows):,}</div><div class="lbl">Con evidencia PCS</div></div>
</div>
<div class="kpi-grid">
  <div class="kpi"><div class="val">{len(unique_genes_total):,}</div><div class="lbl">Genes únicos (ENSG)</div></div>
  <div class="kpi dark"><div class="val">{len(unique_genes_mend):,}</div><div class="lbl">Genes en MENDELIAN</div></div>
  <div class="kpi orange"><div class="val">{len(unique_diseases_total):,}</div><div class="lbl">Enfermedades únicas</div></div>
  <div class="kpi red"><div class="val">{len(unique_genes_pcs):,}</div><div class="lbl">Genes con evidencia PCS</div></div>
</div>

<div class="charts-row">
  <div class="chart-wrap">
    <h2>Por tipo de asociación</h2>
    <canvas id="typeChart"></canvas>
  </div>
  <div class="chart-wrap">
    <h2>Por código de evidencia</h2>
    <canvas id="evChart"></canvas>
    <div class="ev-legend">
      <div class="ev-item ev-pcs">PCS — Publicación con casos clínicos (alta confianza)</div>
      <div class="ev-item ev-tas">TAS — Statement de experto (confianza media)</div>
      <div class="ev-item ev-iea">IEA — Inferido electrónicamente (menor confianza)</div>
    </div>
  </div>
  <div class="chart-wrap">
    <h2>Riqueza fenotípica (nº términos HPO por enfermedad, MENDELIAN)</h2>
    <canvas id="hpoChart"></canvas>
  </div>
</div>

<div class="section">
  <h2>Top 10 — MENDELIAN + PCS · Mayor riqueza fenotípica (nº términos HPO)</h2>
  <table>
    <tr><th>ENSG</th><th>Gen</th><th>Enfermedad (ID)</th><th>Tipo</th><th>Evidencia</th><th>Términos HPO</th></tr>
    {top10_html}
  </table>
</div>

<script>
new Chart(document.getElementById('typeChart'),{{
  type:'doughnut',
  data:{{
    labels:{type_labels},
    datasets:[{{data:{type_values},backgroundColor:['#0a9396','#ee9b00'],borderWidth:0}}]
  }},
  options:{{plugins:{{legend:{{position:'bottom',labels:{{color:'#1a1a2e',font:{{family:'DM Sans',size:11}}}}}}}}}}
}});
new Chart(document.getElementById('evChart'),{{
  type:'bar',
  data:{{
    labels:{ev_labels},
    datasets:[{{data:{ev_values},backgroundColor:['#1e40af','#6b21a8','#475569'],borderRadius:4}}]
  }},
  options:{{
    plugins:{{legend:{{display:false}}}},
    scales:{{y:{{beginAtZero:true,ticks:{{color:'#6b7280'}},grid:{{color:'#f3f4f0'}}}},x:{{ticks:{{color:'#6b7280'}},grid:{{display:false}}}}}}
  }}
}});
new Chart(document.getElementById('hpoChart'),{{
  type:'bar',
  data:{{
    labels:{hpo_labels},
    datasets:[{{label:'Enfermedades',data:{hpo_values},backgroundColor:'#0a9396',borderRadius:4}}]
  }},
  options:{{
    plugins:{{legend:{{display:false}}}},
    scales:{{y:{{beginAtZero:true,ticks:{{color:'#6b7280'}},grid:{{color:'#f3f4f0'}}}},x:{{ticks:{{color:'#6b7280'}},grid:{{display:false}}}}}}
  }}
}});
</script>
</body></html>"""

with open(OUT_HTML, "w") as f:
    f.write(html)
print(f"Dashboard guardado: {OUT_HTML}")
