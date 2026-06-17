"""
process_omnipath_dorothea_reactome.py
- OmniPath: filtra n_references >= 1, mapea UniProt -> ENSG via ens2uniprot
- DoRothEA AB: ya filtrado en descarga, mapea UniProt -> ENSG
- Reactome: filtra solo R-HSA (humano), mapea ENSG (ya nativo)
- Genera TSVs procesados + dashboards HTML
"""
import csv, json, collections

# ── Rutas ──────────────────────────────────────────────────────────────────
OMNI_IN    = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/raw/omnipath/omnipath_interactions.tsv"
DORO_IN    = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/raw/dorothea/dorothea_AB_interactions.tsv"
REACT_IN   = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/raw/reactome/Ensembl2Reactome_All_Levels.txt"
PATHWAYS   = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/raw/reactome/ReactomePathways.txt"
ENS2UNI    = "/beegfs/home/ppadmoremcc/work/external/bioteque/metadata/mappings/GEN/ens2uniprot.tsv"

OMNI_OUT   = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed/omnipath_interactions.tsv"
DORO_OUT   = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed/dorothea_AB_interactions.tsv"
REACT_OUT  = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed/reactome_gene_pathway.tsv"

OMNI_HTML  = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed/dashboard_omnipath.html"
DORO_HTML  = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed/dashboard_dorothea.html"
REACT_HTML = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed/dashboard_reactome.html"

# ── Mapeo UniProt -> ENSG ──────────────────────────────────────────────────
print("Cargando mapeo ENSG -> UniProt...")
uniprot2ensg = {}
with open(ENS2UNI) as f:
    for line in f:
        p = line.strip().split("\t")
        if len(p) == 2:
            ensg, uniprot = p
            if uniprot not in uniprot2ensg:
                uniprot2ensg[uniprot] = ensg
print(f"  {len(uniprot2ensg):,} UniProt -> ENSG")

# ═══════════════════════════════════════════════════════════════════════════
# 1. OmniPath
# ═══════════════════════════════════════════════════════════════════════════
print("\nProcesando OmniPath...")
omni_rows = []
omni_skip = 0

with open(OMNI_IN, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        references = row.get("references","").strip()
        sources    = row.get("sources","").strip()

        if not references:
            omni_skip += 1
            continue

        n_ref = str(len(references.split(";"))) if references else "0"

        src_uni = row.get("source","").strip()
        tgt_uni = row.get("target","").strip()
        src_ensg = uniprot2ensg.get(src_uni, "")
        tgt_ensg = uniprot2ensg.get(tgt_uni, "")

        omni_rows.append({
            "source_uniprot":      src_uni,
            "target_uniprot":      tgt_uni,
            "source_ensg":         src_ensg,
            "target_ensg":         tgt_ensg,
            "source_genesymbol":   row.get("source_genesymbol",""),
            "target_genesymbol":   row.get("target_genesymbol",""),
            "is_directed":         row.get("is_directed",""),
            "is_stimulation":      row.get("is_stimulation",""),
            "is_inhibition":       row.get("is_inhibition",""),
            "consensus_direction": row.get("consensus_direction",""),
            "n_references":        n_ref,
            "n_resources":         row.get("n_resources",""),
        })

fn = ["source_uniprot","target_uniprot","source_ensg","target_ensg",
      "source_genesymbol","target_genesymbol","is_directed",
      "is_stimulation","is_inhibition","consensus_direction","n_references","n_resources"]
with open(OMNI_OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fn, delimiter="\t")
    w.writeheader(); w.writerows(omni_rows)

print(f"  Interacciones con n_references >= 1: {len(omni_rows):,}")
print(f"  Descartadas (sin publicación): {omni_skip:,}")

# Stats OmniPath
dir_counter  = collections.Counter(r["is_directed"] for r in omni_rows)
stim_counter = collections.Counter()
for r in omni_rows:
    if r["is_stimulation"] == "True":   stim_counter["Activación"] += 1
    elif r["is_inhibition"] == "True":  stim_counter["Inhibición"] += 1
    else:                               stim_counter["Sin signo"]  += 1

ref_bins = collections.Counter()
for r in omni_rows:
    v = int(r["n_references"])
    if v == 1:   ref_bins["1"] += 1
    elif v <= 3: ref_bins["2-3"] += 1
    elif v <= 10:ref_bins["4-10"] += 1
    else:        ref_bins["10+"] += 1

omni_genes = {r["source_ensg"] for r in omni_rows if r["source_ensg"]} | \
             {r["target_ensg"] for r in omni_rows if r["target_ensg"]}

# OmniPath HTML
omni_html = f"""<!DOCTYPE html>
<html lang="es"><head><meta charset="UTF-8"><title>OmniPath Dashboard</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
<link href="https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;600;800&family=JetBrains+Mono:wght@400;600&display=swap" rel="stylesheet">
<style>
  :root{{--bg:#0c0c14;--s:#13131f;--s2:#1a1a2e;--a:#f97316;--b:#a855f7;--c:#22d3ee;--t:#f1f5f9;--m:#6b7280;--bd:#1e1e30}}
  *{{box-sizing:border-box;margin:0;padding:0}}
  body{{font-family:'Outfit',sans-serif;background:var(--bg);color:var(--t);padding:36px}}
  h1{{font-size:2rem;font-weight:800;color:var(--a)}}
  .meta{{color:var(--m);font-size:.82rem;margin-top:8px;font-family:'JetBrains Mono',monospace}}
  .fn{{background:rgba(249,115,22,.07);border-left:3px solid var(--a);padding:10px 16px;font-size:.82rem;color:#fed7aa;margin:20px 0}}
  .kg{{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin-bottom:26px}}
  .k{{background:var(--s);border:1px solid var(--bd);border-radius:6px;padding:20px;border-top:2px solid var(--a)}}
  .k.b{{border-top-color:var(--b)}} .k.c{{border-top-color:var(--c)}}
  .v{{font-family:'JetBrains Mono',monospace;font-size:1.9rem;font-weight:600;color:var(--a)}}
  .k.b .v{{color:var(--b)}} .k.c .v{{color:var(--c)}}
  .l{{font-size:.73rem;color:var(--m);margin-top:5px;text-transform:uppercase;letter-spacing:.05em}}
  .cr{{display:grid;grid-template-columns:1fr 1fr;gap:18px;margin-bottom:22px}}
  .cw{{background:var(--s);border:1px solid var(--bd);border-radius:6px;padding:22px}}
  h2{{font-size:.85rem;font-weight:700;color:var(--a);margin-bottom:14px;text-transform:uppercase;letter-spacing:.06em}}
  canvas{{max-height:250px}}
</style></head><body>
<h1>OmniPath — Signaling Network</h1>
<div class="meta">pypath v0.16.20 (Jun 2025) &nbsp;·&nbsp; Dataset: omnipath &nbsp;·&nbsp; Libre</div>
<div class="fn"><strong>Filtro:</strong> n_references ≥ 1 (al menos una publicación que respalda la interacción).
OmniPath solo incluye interacciones dirigidas con signo (activación/inhibición) de más de 100 recursos curados.</div>
<div class="kg">
  <div class="k"><div class="v">{len(omni_rows):,}</div><div class="l">Interacciones filtradas</div></div>
  <div class="k b"><div class="v">{len(omni_genes):,}</div><div class="l">Genes únicos (ENSG)</div></div>
  <div class="k c"><div class="v">{stim_counter.get("Activación",0):,}</div><div class="l">Interacciones activadoras</div></div>
  <div class="k"><div class="v">{stim_counter.get("Inhibición",0):,}</div><div class="l">Interacciones inhibidoras</div></div>
</div>
<div class="cr">
  <div class="cw"><h2>Por tipo de interacción (signo)</h2><canvas id="signChart"></canvas></div>
  <div class="cw"><h2>Distribución nº publicaciones de respaldo</h2><canvas id="refChart"></canvas></div>
</div>
<script>
new Chart(document.getElementById('signChart'),{{type:'doughnut',
  data:{{labels:{json.dumps(list(stim_counter.keys()))},
    datasets:[{{data:{json.dumps(list(stim_counter.values()))},backgroundColor:['#f97316','#a855f7','#6b7280'],borderWidth:0}}]}},
  options:{{plugins:{{legend:{{position:'bottom',labels:{{color:'#f1f5f9',font:{{size:11}}}}}}}}}}
}});
new Chart(document.getElementById('refChart'),{{type:'bar',
  data:{{labels:{json.dumps(list(ref_bins.keys()))},
    datasets:[{{data:{json.dumps(list(ref_bins.values()))},backgroundColor:'#f97316',borderRadius:4}}]}},
  options:{{plugins:{{legend:{{display:false}}}},scales:{{y:{{beginAtZero:true,ticks:{{color:'#6b7280'}},grid:{{color:'#1e1e30'}}}},x:{{ticks:{{color:'#94a3b8'}},grid:{{display:false}}}}}}}}
}});
</script></body></html>"""

with open(OMNI_HTML, "w") as f:
    f.write(omni_html)
print(f"Dashboard OmniPath guardado: {OMNI_HTML}")

# ═══════════════════════════════════════════════════════════════════════════
# 2. DoRothEA AB
# ═══════════════════════════════════════════════════════════════════════════
print("\nProcesando DoRothEA AB...")
doro_rows = []

with open(DORO_IN, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        src_uni  = row.get("source","").strip()
        tgt_uni  = row.get("target","").strip()
        src_ensg = uniprot2ensg.get(src_uni, "")
        tgt_ensg = uniprot2ensg.get(tgt_uni, "")
        doro_rows.append({
            "tf_uniprot":        src_uni,
            "target_uniprot":    tgt_uni,
            "tf_ensg":           src_ensg,
            "target_ensg":       tgt_ensg,
            "tf_symbol":         row.get("source_genesymbol",""),
            "target_symbol":     row.get("target_genesymbol",""),
            "is_stimulation":    row.get("is_stimulation",""),
            "is_inhibition":     row.get("is_inhibition",""),
            "consensus_direction":row.get("consensus_direction",""),
        })

fn_d = ["tf_uniprot","target_uniprot","tf_ensg","target_ensg",
        "tf_symbol","target_symbol","is_stimulation","is_inhibition","consensus_direction"]
with open(DORO_OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fn_d, delimiter="\t")
    w.writeheader(); w.writerows(doro_rows)

print(f"  Interacciones TF-gen (tiers A+B): {len(doro_rows):,}")

# Stats DoRothEA
doro_sign = collections.Counter()
for r in doro_rows:
    if r["is_stimulation"] == "True":  doro_sign["Activación"] += 1
    elif r["is_inhibition"] == "True": doro_sign["Inhibición"] += 1
    else:                              doro_sign["Sin signo"]  += 1

unique_tfs      = {r["tf_symbol"] for r in doro_rows if r["tf_symbol"]}
unique_targets  = {r["target_symbol"] for r in doro_rows if r["target_symbol"]}
tf_degree       = collections.Counter(r["tf_symbol"] for r in doro_rows)
top10_tfs       = tf_degree.most_common(10)
top10_tfs_html  = "".join(f"<tr><td>{t}</td><td>{c:,}</td></tr>" for t, c in top10_tfs)

doro_html = f"""<!DOCTYPE html>
<html lang="es"><head><meta charset="UTF-8"><title>DoRothEA Dashboard</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
<link href="https://fonts.googleapis.com/css2?family=DM+Sans:wght@300;400;500;700&family=JetBrains+Mono:wght@400;600&display=swap" rel="stylesheet">
<style>
  :root{{--bg:#fafafa;--s:#fff;--s2:#f4f4f4;--a:#7c3aed;--b:#0891b2;--c:#dc2626;--t:#1e1b4b;--m:#6b7280;--bd:#e5e7eb}}
  *{{box-sizing:border-box;margin:0;padding:0}}
  body{{font-family:'DM Sans',sans-serif;background:var(--bg);color:var(--t);padding:36px}}
  h1{{font-size:2rem;font-weight:700;color:var(--a)}}
  .meta{{color:var(--m);font-size:.82rem;margin-top:8px}}
  .fn{{background:#ede9fe;border-left:3px solid var(--a);padding:10px 16px;font-size:.82rem;color:#4c1d95;margin:20px 0}}
  .kg{{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin-bottom:26px}}
  .k{{background:var(--s);border:1px solid var(--bd);border-radius:8px;padding:20px;border-top:3px solid var(--a)}}
  .k.b{{border-top-color:var(--b)}} .k.c{{border-top-color:var(--c)}}
  .v{{font-family:'JetBrains Mono',monospace;font-size:1.9rem;font-weight:600;color:var(--a)}}
  .k.b .v{{color:var(--b)}} .k.c .v{{color:var(--c)}}
  .l{{font-size:.73rem;color:var(--m);margin-top:5px;text-transform:uppercase;letter-spacing:.05em}}
  .cr{{display:grid;grid-template-columns:1fr 1fr;gap:18px;margin-bottom:22px}}
  .cw,.sec{{background:var(--s);border:1px solid var(--bd);border-radius:8px;padding:22px;margin-bottom:22px}}
  h2{{font-size:.85rem;font-weight:700;color:var(--a);margin-bottom:14px;text-transform:uppercase;letter-spacing:.05em}}
  table{{width:100%;border-collapse:collapse;font-size:.82rem}}
  th{{background:var(--s2);color:var(--a);padding:9px 12px;text-align:left;font-size:.72rem;text-transform:uppercase;border-bottom:2px solid var(--bd)}}
  td{{padding:8px 12px;border-bottom:1px solid var(--bd)}}
  tr:hover td{{background:var(--s2)}}
  canvas{{max-height:250px}}
</style></head><body>
<h1>DoRothEA — TF Regulons (Tiers A+B)</h1>
<div class="meta">Bioconductor v1.10.0 (Dic 2025) &nbsp;·&nbsp; vía OmniPath webservice &nbsp;·&nbsp; Libre</div>
<div class="fn"><strong>Filtro:</strong> Tiers A y B únicamente — las interacciones TF-gen con mayor respaldo experimental
(ChIP-seq + literatura curada + motif analysis). Tier A = máxima confianza, Tier B = alta confianza.</div>
<div class="kg">
  <div class="k"><div class="v">{len(doro_rows):,}</div><div class="l">Interacciones TF-gen</div></div>
  <div class="k b"><div class="v">{len(unique_tfs):,}</div><div class="l">TFs únicos</div></div>
  <div class="k c"><div class="v">{len(unique_targets):,}</div><div class="l">Genes diana únicos</div></div>
  <div class="k"><div class="v">{doro_sign.get("Activación",0):,}</div><div class="l">Interacciones activadoras</div></div>
</div>
<div class="cr">
  <div class="cw"><h2>Por tipo de regulación</h2><canvas id="signChart"></canvas></div>
  <div class="cw"><h2>Top 10 TFs — mayor número de genes diana</h2>
    <table><tr><th>Factor de transcripción</th><th>Genes diana</th></tr>{top10_tfs_html}</table>
  </div>
</div>
<script>
new Chart(document.getElementById('signChart'),{{type:'doughnut',
  data:{{labels:{json.dumps(list(doro_sign.keys()))},
    datasets:[{{data:{json.dumps(list(doro_sign.values()))},backgroundColor:['#7c3aed','#dc2626','#6b7280'],borderWidth:0}}]}},
  options:{{plugins:{{legend:{{position:'bottom',labels:{{color:'#1e1b4b',font:{{size:11}}}}}}}}}}
}});
</script></body></html>"""

with open(DORO_HTML, "w") as f:
    f.write(doro_html)
print(f"Dashboard DoRothEA guardado: {DORO_HTML}")

# ═══════════════════════════════════════════════════════════════════════════
# 3. Reactome
# ═══════════════════════════════════════════════════════════════════════════
print("\nProcesando Reactome (solo Homo sapiens)...")

# Cargar nombres de vías
pathway_names = {}
with open(PATHWAYS) as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            pathway_names[p[0]] = p[1]

react_rows = []
skipped_species = 0

with open(REACT_IN) as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 6:
            continue
        ensg        = parts[0].strip()
        pathway_id  = parts[1].strip()
        pathway_url = parts[2].strip()
        pathway_name= parts[3].strip()
        evidence    = parts[4].strip()
        species     = parts[5].strip()

        # Solo Homo sapiens
        if not pathway_id.startswith("R-HSA-"):
            skipped_species += 1
            continue
        # Solo ENSGs humanos
        if not ensg.startswith("ENSG"):
            skipped_species += 1
            continue

        react_rows.append({
            "ensembl_gene_id": ensg,
            "pathway_id":      pathway_id,
            "pathway_name":    pathway_name,
            "evidence":        evidence,
        })

fn_r = ["ensembl_gene_id","pathway_id","pathway_name","evidence"]
with open(REACT_OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fn_r, delimiter="\t")
    w.writeheader(); w.writerows(react_rows)

print(f"  Asociaciones gen-vía (R-HSA): {len(react_rows):,}")
print(f"  Descartadas (otras especies): {skipped_species:,}")

# Stats Reactome
ev_counter    = collections.Counter(r["evidence"] for r in react_rows)
unique_genes  = {r["ensembl_gene_id"] for r in react_rows}
unique_paths  = {r["pathway_id"] for r in react_rows}
path_gene_cnt = collections.Counter(r["pathway_id"] for r in react_rows)
top10_paths   = [(pid, cnt, pathway_names.get(pid, pid)) for pid, cnt in path_gene_cnt.most_common(10)]
top10_p_html  = "".join(f"<tr><td class='mono'>{pid}</td><td>{name}</td><td>{cnt:,}</td></tr>"
                         for pid, cnt, name in top10_paths)

gene_path_cnt = collections.Counter(r["ensembl_gene_id"] for r in react_rows)
top10_genes   = gene_path_cnt.most_common(10)
top10_g_html  = "".join(f"<tr><td class='mono'>{g}</td><td>{c:,}</td></tr>" for g, c in top10_genes)

react_html = f"""<!DOCTYPE html>
<html lang="es"><head><meta charset="UTF-8"><title>Reactome Dashboard</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
<link href="https://fonts.googleapis.com/css2?family=Sora:wght@300;400;600;800&family=JetBrains+Mono:wght@400;600&display=swap" rel="stylesheet">
<style>
  :root{{--bg:#001a1a;--s:#002626;--s2:#003333;--a:#00e5cc;--b:#ff6b6b;--c:#ffd166;--t:#e8fafa;--m:#5a8a8a;--bd:#004444}}
  *{{box-sizing:border-box;margin:0;padding:0}}
  body{{font-family:'Sora',sans-serif;background:var(--bg);color:var(--t);padding:36px}}
  h1{{font-size:2rem;font-weight:800;color:var(--a)}}
  .meta{{color:var(--m);font-size:.82rem;margin-top:8px;font-family:'JetBrains Mono',monospace}}
  .fn{{background:rgba(0,229,204,.06);border-left:3px solid var(--a);padding:10px 16px;font-size:.82rem;color:#99fff5;margin:20px 0}}
  .kg{{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin-bottom:26px}}
  .k{{background:var(--s);border:1px solid var(--bd);border-radius:6px;padding:20px;border-top:2px solid var(--a)}}
  .k.b{{border-top-color:var(--b)}} .k.c{{border-top-color:var(--c)}}
  .v{{font-family:'JetBrains Mono',monospace;font-size:1.9rem;font-weight:600;color:var(--a)}}
  .k.b .v{{color:var(--b)}} .k.c .v{{color:var(--c)}}
  .l{{font-size:.73rem;color:var(--m);margin-top:5px;text-transform:uppercase;letter-spacing:.05em}}
  .tr{{display:grid;grid-template-columns:1fr 1fr;gap:18px;margin-bottom:22px}}
  .cw,.sec{{background:var(--s);border:1px solid var(--bd);border-radius:6px;padding:22px;margin-bottom:22px}}
  h2{{font-size:.85rem;font-weight:700;color:var(--a);margin-bottom:14px;text-transform:uppercase;letter-spacing:.06em}}
  table{{width:100%;border-collapse:collapse;font-size:.8rem}}
  th{{background:var(--s2);color:var(--a);padding:9px 12px;text-align:left;font-size:.71rem;text-transform:uppercase;border-bottom:1px solid var(--bd)}}
  td{{padding:8px 12px;border-bottom:1px solid var(--bd)}}
  tr:hover td{{background:var(--s2)}}
  .mono{{font-family:'JetBrains Mono',monospace;font-size:.76rem}}
  canvas{{max-height:250px}}
</style></head><body>
<h1>Reactome — Gene-Pathway</h1>
<div class="meta">Release actual (reactome.org/download/current) &nbsp;·&nbsp; Solo Homo sapiens (R-HSA) &nbsp;·&nbsp; Libre</div>
<div class="fn"><strong>Filtro:</strong> Solo vías con prefijo R-HSA (Homo sapiens) y ENSGs humanos.
El fichero Ensembl2Reactome_All_Levels incluye todos los niveles jerárquicos de la vía (top-level hasta bottom-level).</div>
<div class="kg">
  <div class="k"><div class="v">{len(react_rows):,}</div><div class="l">Asociaciones gen-vía</div></div>
  <div class="k b"><div class="v">{len(unique_genes):,}</div><div class="l">Genes únicos (ENSG)</div></div>
  <div class="k c"><div class="v">{len(unique_paths):,}</div><div class="l">Vías únicas (R-HSA)</div></div>
  <div class="k"><div class="v">{ev_counter.get("TAS",0):,}</div><div class="l">Evidencia TAS (curada)</div></div>
</div>
<div class="tr">
  <div class="sec">
    <h2>Top 10 vías — mayor número de genes</h2>
    <table><tr><th>Pathway ID</th><th>Nombre</th><th>Genes</th></tr>{top10_p_html}</table>
  </div>
  <div class="sec">
    <h2>Top 10 genes — mayor número de vías</h2>
    <table><tr><th>ENSG</th><th>Nº vías</th></tr>{top10_g_html}</table>
  </div>
</div>
</body></html>"""

with open(REACT_HTML, "w") as f:
    f.write(react_html)
print(f"Dashboard Reactome guardado: {REACT_HTML}")
