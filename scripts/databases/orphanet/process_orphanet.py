"""
process_orphanet.py
- Loads orphanet_gene_disease.tsv (already parsed, ENSG native)
- Filters by association_status = Assessed
- Computes stats by association_type
- Generates dashboard HTML
"""
import csv, json, collections

# ── Paths ──────────────────────────────────────────────────────────────────
IN_TSV   = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/raw/orphanet/orphanet_gene_disease.tsv"
OUT_TSV  = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed/orphanet_gene_disease.tsv"
OUT_HTML = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed/dashboard_orphanet.html"

ORPHA_VERSION = "v1.3.42 / 4.1.8 [2025-03-03]"
ORPHA_DATE    = "2025-12-09"

# Association types and their reliability tier
TYPE_TIER = {
    "Disease-causing germline mutation(s) in": 1,
    "Disease-causing somatic mutation(s) in":  1,
    "Modifying germline mutation in":           2,
    "Major susceptibility factor in":           2,
    "Role in the phenotype of":                 3,
    "Part of a fusion gene in":                 3,
    "Candidate gene tested in":                 4,
}

def get_tier(assoc_type):
    for key, tier in TYPE_TIER.items():
        if key in assoc_type:
            return tier
    return 5

# ── Load and filter ─────────────────────────────────────────────────────────
print("Cargando orphanet_gene_disease.tsv...")
all_rows      = []
assessed_rows = []
skipped_not_assessed = 0
skipped_no_ensg      = 0

with open(IN_TSV) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        all_rows.append(row)
        if not row.get("ensembl_id"):
            skipped_no_ensg += 1
            continue
        if row.get("association_status", "").strip() != "Assessed":
            skipped_not_assessed += 1
            continue
        assessed_rows.append(row)

print(f"  Total filas en TSV:       {len(all_rows):,}")
print(f"  Assessed:                 {len(assessed_rows):,}")
print(f"  Not yet validated (omit): {skipped_not_assessed:,}")
print(f"  Sin ENSG (omit):          {skipped_no_ensg:,}")

# ── Write filtered TSV ───────────────────────────────────────────────────────
with open(OUT_TSV, "w", newline="") as f:
    if assessed_rows:
        writer = csv.DictWriter(f, fieldnames=assessed_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(assessed_rows)
print(f"TSV filtrado guardado: {OUT_TSV}")

# ── Stats ─────────────────────────────────────────────────────────────────────
type_counter = collections.Counter(r["association_type"] for r in assessed_rows)
unique_genes    = {r["ensembl_id"]  for r in assessed_rows}
unique_diseases = {r["orpha_code"]  for r in assessed_rows}

# Tier breakdown
tier_counter = collections.Counter(get_tier(r["association_type"]) for r in assessed_rows)

# Disease burden: diseases with most genes
disease_gene_count = collections.Counter(r["orpha_code"] for r in assessed_rows)
top10_diseases = disease_gene_count.most_common(10)

# Disease names lookup
disease_names = {}
for r in assessed_rows:
    disease_names[r["orpha_code"]] = r["disease_name"]

# Top 10 genes by number of diseases
gene_disease_count = collections.Counter(r["ensembl_id"] for r in assessed_rows)
top10_genes = gene_disease_count.most_common(10)
gene_symbols = {}
for r in assessed_rows:
    gene_symbols[r["ensembl_id"]] = r["gene_symbol"]

# ── Build HTML ────────────────────────────────────────────────────────────────
# Shorten long association type labels for chart
def short_type(t):
    if "Disease-causing germline" in t: return "Causativa germinal"
    if "Disease-causing somatic"  in t: return "Causativa somática"
    if "Modifying germline"       in t: return "Modificadora"
    if "Major susceptibility"     in t: return "Factor susceptibilidad"
    if "Role in the phenotype"    in t: return "Rol en fenotipo"
    if "Part of a fusion"         in t: return "Fusión génica"
    if "Candidate gene"           in t: return "Gen candidato"
    return t[:30]

type_labels_short = [short_type(k) for k in type_counter.keys()]
type_values_list  = list(type_counter.values())

tier_labels = ["Tier 1 — Causativa", "Tier 2 — Moduladora/Susceptibilidad",
               "Tier 3 — Rol/Fusión", "Tier 4 — Candidato", "Tier 5 — Otro"]
tier_values = [tier_counter.get(i, 0) for i in range(1, 6)]
# trim trailing zeros
while tier_values and tier_values[-1] == 0:
    tier_labels.pop()
    tier_values.pop()

top10_disease_html = ""
for code, count in top10_diseases:
    name = disease_names.get(code, "")
    top10_disease_html += f"<tr><td>ORPHA:{code}</td><td>{name}</td><td>{count}</td></tr>\n"

top10_gene_html = ""
for ensg, count in top10_genes:
    sym = gene_symbols.get(ensg, "")
    top10_gene_html += f"<tr><td>{ensg}</td><td>{sym}</td><td>{count}</td></tr>\n"

html = f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="UTF-8">
<title>Orphanet — Dashboard</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
<link href="https://fonts.googleapis.com/css2?family=Fraunces:wght@400;700;900&family=Inter:wght@300;400;500&display=swap" rel="stylesheet">
<style>
  :root{{
    --bg:#0d0d0d;--surface:#161616;--surface2:#1f1f1f;
    --gold:#c9a84c;--teal:#2dd4bf;--rose:#fb7185;
    --text:#f0ede8;--muted:#737373;--border:#2a2a2a;
  }}
  *{{box-sizing:border-box;margin:0;padding:0}}
  body{{font-family:'Inter',sans-serif;background:var(--bg);color:var(--text);padding:36px}}
  .header{{margin-bottom:38px;display:flex;align-items:flex-start;gap:24px}}
  .header-text h1{{font-family:'Fraunces',serif;font-size:2.2rem;font-weight:900;color:var(--gold);line-height:1}}
  .header-text .sub{{font-size:1rem;color:var(--teal);font-family:'Fraunces',serif;margin-top:4px}}
  .header-text .meta{{color:var(--muted);font-size:.8rem;margin-top:10px}}
  .badge-assessed{{display:inline-block;padding:3px 10px;background:rgba(45,212,191,.1);color:var(--teal);border:1px solid var(--teal);border-radius:20px;font-size:.72rem;margin-top:8px}}
  .kpi-grid{{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin-bottom:26px}}
  .kpi{{background:var(--surface);border:1px solid var(--border);border-radius:6px;padding:20px;position:relative}}
  .kpi::after{{content:'';position:absolute;bottom:0;left:0;right:0;height:2px;background:var(--gold)}}
  .kpi.teal::after{{background:var(--teal)}}
  .kpi.rose::after{{background:var(--rose)}}
  .kpi .val{{font-family:'Fraunces',serif;font-size:2rem;font-weight:700;color:var(--gold)}}
  .kpi.teal .val{{color:var(--teal)}}
  .kpi.rose .val{{color:var(--rose)}}
  .kpi .lbl{{font-size:.73rem;color:var(--muted);margin-top:6px;text-transform:uppercase;letter-spacing:.05em}}
  .charts-row{{display:grid;grid-template-columns:1fr 1fr;gap:18px;margin-bottom:22px}}
  .chart-wrap{{background:var(--surface);border:1px solid var(--border);border-radius:6px;padding:20px}}
  .chart-wrap h2{{font-family:'Fraunces',serif;font-size:.9rem;color:var(--gold);margin-bottom:14px}}
  .tables-row{{display:grid;grid-template-columns:1fr 1fr;gap:18px}}
  .section{{background:var(--surface);border:1px solid var(--border);border-radius:6px;padding:22px}}
  .section h2{{font-family:'Fraunces',serif;font-size:.9rem;color:var(--gold);margin-bottom:14px}}
  table{{width:100%;border-collapse:collapse;font-size:.81rem}}
  th{{background:var(--surface2);color:var(--teal);padding:9px 12px;text-align:left;font-size:.72rem;text-transform:uppercase;letter-spacing:.05em;border-bottom:1px solid var(--border)}}
  td{{padding:8px 12px;border-bottom:1px solid var(--border);color:var(--text)}}
  tr:hover td{{background:var(--surface2)}}
  .tier-legend{{margin-top:14px;display:grid;gap:6px;font-size:.78rem}}
  .tier-item{{padding:8px 12px;border-radius:4px;border-left:3px solid}}
  .t1{{background:rgba(201,168,76,.08);border-color:var(--gold);color:var(--gold)}}
  .t2{{background:rgba(45,212,191,.08);border-color:var(--teal);color:var(--teal)}}
  .t3{{background:rgba(251,113,133,.08);border-color:var(--rose);color:var(--rose)}}
  canvas{{max-height:260px}}
</style>
</head>
<body>

<div class="header">
  <div class="header-text">
    <h1>Orphanet</h1>
    <div class="sub">Rare Disease Gene Associations</div>
    <div class="meta">{ORPHA_VERSION} &nbsp;·&nbsp; Datos: {ORPHA_DATE} &nbsp;·&nbsp; CC-BY 4.0</div>
    <span class="badge-assessed">✓ Filtro: association_status = Assessed</span>
  </div>
</div>

<div class="kpi-grid">
  <div class="kpi"><div class="val">{len(assessed_rows):,}</div><div class="lbl">Asociaciones Assessed</div></div>
  <div class="kpi teal"><div class="val">{len(unique_genes):,}</div><div class="lbl">Genes únicos (ENSG)</div></div>
  <div class="kpi rose"><div class="val">{len(unique_diseases):,}</div><div class="lbl">Enfermedades únicas (ORPHAcode)</div></div>
  <div class="kpi"><div class="val">{skipped_not_assessed:,}</div><div class="lbl">Not yet validated (omitidas)</div></div>
</div>

<div class="charts-row">
  <div class="chart-wrap">
    <h2>Por tipo de asociación gen–enfermedad</h2>
    <canvas id="typeChart"></canvas>
  </div>
  <div class="chart-wrap">
    <h2>Por tier de confianza</h2>
    <canvas id="tierChart"></canvas>
    <div class="tier-legend">
      <div class="tier-item t1">Tier 1 — Causativa (germinal / somática) · máxima confianza</div>
      <div class="tier-item t2">Tier 2 — Modificadora / Factor de susceptibilidad</div>
      <div class="tier-item t3">Tier 3–4 — Rol en fenotipo / Fusión / Candidato</div>
    </div>
  </div>
</div>

<div class="tables-row">
  <div class="section">
    <h2>Top 10 enfermedades — mayor número de genes asociados</h2>
    <table>
      <tr><th>ORPHAcode</th><th>Enfermedad</th><th>Genes</th></tr>
      {top10_disease_html}
    </table>
  </div>
  <div class="section">
    <h2>Top 10 genes — mayor número de enfermedades asociadas</h2>
    <table>
      <tr><th>ENSG</th><th>Símbolo</th><th>Enfermedades</th></tr>
      {top10_gene_html}
    </table>
  </div>
</div>

<script>
new Chart(document.getElementById('typeChart'),{{
  type:'bar',
  data:{{
    labels:{json.dumps(type_labels_short)},
    datasets:[{{
      data:{json.dumps(type_values_list)},
      backgroundColor:['#c9a84c','#2dd4bf','#fb7185','#a78bfa','#34d399','#60a5fa','#f97316'],
      borderRadius:4
    }}]
  }},
  options:{{
    indexAxis:'y',
    plugins:{{legend:{{display:false}}}},
    scales:{{
      x:{{beginAtZero:true,ticks:{{color:'#737373'}},grid:{{color:'#2a2a2a'}}}},
      y:{{ticks:{{color:'#f0ede8',font:{{size:11}}}},grid:{{display:false}}}}
    }}
  }}
}});
new Chart(document.getElementById('tierChart'),{{
  type:'doughnut',
  data:{{
    labels:{json.dumps(tier_labels)},
    datasets:[{{
      data:{json.dumps(tier_values)},
      backgroundColor:['#c9a84c','#2dd4bf','#fb7185','#a78bfa','#737373'],
      borderColor:'#161616',borderWidth:3
    }}]
  }},
  options:{{plugins:{{legend:{{position:'bottom',labels:{{color:'#f0ede8',font:{{size:11}}}}}}}}}}
}});
</script>
</body></html>"""

with open(OUT_HTML, "w") as f:
    f.write(html)
print(f"Dashboard guardado: {OUT_HTML}")
