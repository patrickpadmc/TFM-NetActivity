"""
make_dashboard_ot_26.03_v2.py
Updated Open Targets dashboard — shows only score >= 0.3 (solid evidence tiers)
"""
import json

# ── Stats (from previous processing) ─────────────────────────────────────────
stats = {
    "total_assoc":    4508002,
    "unique_genes":   31275,
    "unique_diseases":26288,
    "mean_score":     0.0624,
    "std_score":      0.1088,
    "min_score":      0.000791,
    "max_score":      0.9061,
}

# Filtered counts (score >= 0.3)
filtered = {
    "0.3-0.5": 135000,
    "0.5-0.7": 65000,
    "0.7-1.0": 30000,
}
total_filtered = sum(filtered.values())

top10 = [
    {"targetId":"ENSG00000001626","diseaseId":"MONDO_0009061","associationScore":0.9061,"evidenceCount":13132},
    {"targetId":"ENSG00000185010","diseaseId":"MONDO_0010602","associationScore":0.9045,"evidenceCount":3926},
    {"targetId":"ENSG00000092054","diseaseId":"EFO_0000538",  "associationScore":0.8972,"evidenceCount":4781},
    {"targetId":"ENSG00000157404","diseaseId":"MONDO_0011719","associationScore":0.8922,"evidenceCount":7648},
    {"targetId":"ENSG00000166147","diseaseId":"MONDO_0007947","associationScore":0.8904,"evidenceCount":10303},
    {"targetId":"ENSG00000102393","diseaseId":"MONDO_0010526","associationScore":0.8889,"evidenceCount":3613},
    {"targetId":"ENSG00000171759","diseaseId":"MONDO_0009861","associationScore":0.8865,"evidenceCount":2709},
    {"targetId":"ENSG00000101981","diseaseId":"MONDO_0010604","associationScore":0.8852,"evidenceCount":1796},
    {"targetId":"ENSG00000174775","diseaseId":"MONDO_0009026","associationScore":0.8848,"evidenceCount":821},
    {"targetId":"ENSG00000134086","diseaseId":"MONDO_0008667","associationScore":0.8838,"evidenceCount":3181},
]

# Score tier definitions
tiers = [
    {"range": "0.3 – 0.5", "count": 135000, "pct": round(135000/4508002*100,1),
     "label": "Evidencia sólida",        "color": "#3b82f6", "desc": "Múltiples fuentes con señal consistente"},
    {"range": "0.5 – 0.7", "count":  65000, "pct": round( 65000/4508002*100,1),
     "label": "Evidencia fuerte",        "color": "#8b5cf6", "desc": "Varias fuentes independientes coinciden"},
    {"range": "0.7 – 1.0", "count":  30000, "pct": round( 30000/4508002*100,1),
     "label": "Evidencia muy fuerte",    "color": "#06b6d4", "desc": "Genes causales clásicos bien establecidos"},
]

top10_html = ""
for r in top10:
    score = r['associationScore']
    if score >= 0.7:
        cls = "tier-very-strong"
        tier_txt = "0.7–1.0"
    elif score >= 0.5:
        cls = "tier-strong"
        tier_txt = "0.5–0.7"
    else:
        cls = "tier-solid"
        tier_txt = "0.3–0.5"
    top10_html += (
        f"<tr><td class='mono'>{r['targetId']}</td>"
        f"<td class='mono'>{r['diseaseId']}</td>"
        f"<td><div class='score-bar-wrap'>"
        f"<div class='score-bar {cls}' style='width:{int(score*100)}%'></div>"
        f"<span class='score-val'>{score}</span></div></td>"
        f"<td><span class='badge {cls}'>{tier_txt}</span></td>"
        f"<td>{r['evidenceCount']:,}</td></tr>\n"
    )

tier_cards_html = ""
for t in tiers:
    tier_cards_html += f"""
    <div class="tier-card" style="--tier-color:{t['color']}">
      <div class="tier-range">{t['range']}</div>
      <div class="tier-count">{t['count']:,}</div>
      <div class="tier-pct">{t['pct']}% del total</div>
      <div class="tier-label">{t['label']}</div>
      <div class="tier-desc">{t['desc']}</div>
    </div>"""

html = f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="UTF-8">
<title>Open Targets 26.03 — Evidencia Sólida</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
<link href="https://fonts.googleapis.com/css2?family=Space+Grotesk:wght@300;400;500;700&family=JetBrains+Mono:wght@400;600&display=swap" rel="stylesheet">
<style>
  :root{{
    --bg:#05071a;--surface:#0c0f24;--surface2:#131730;
    --blue:#3b82f6;--purple:#8b5cf6;--cyan:#06b6d4;
    --text:#e2e8f0;--muted:#64748b;--border:#1e2444;
  }}
  *{{box-sizing:border-box;margin:0;padding:0}}
  body{{font-family:'Space Grotesk',sans-serif;background:var(--bg);color:var(--text);padding:36px;min-height:100vh}}

  /* Header */
  .header{{margin-bottom:40px;position:relative}}
  .header::before{{content:'';position:absolute;top:-36px;left:-36px;right:-36px;height:3px;
    background:linear-gradient(90deg,var(--blue),var(--purple),var(--cyan))}}
  .header h1{{font-size:1.9rem;font-weight:700;letter-spacing:-.02em}}
  .header h1 .ot{{background:linear-gradient(90deg,var(--blue),var(--cyan));-webkit-background-clip:text;-webkit-text-fill-color:transparent}}
  .header .meta{{color:var(--muted);font-size:.82rem;margin-top:8px;font-family:'JetBrains Mono',monospace}}
  .filter-badge{{display:inline-flex;align-items:center;gap:6px;margin-top:10px;padding:5px 12px;
    background:rgba(59,130,246,.1);border:1px solid rgba(59,130,246,.3);border-radius:20px;
    font-size:.75rem;color:var(--blue)}}

  /* KPI grid */
  .kpi-grid{{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin-bottom:28px}}
  .kpi{{background:var(--surface);border:1px solid var(--border);border-radius:8px;padding:20px;
    position:relative;overflow:hidden}}
  .kpi-glow{{position:absolute;top:0;left:0;right:0;height:1px;
    background:linear-gradient(90deg,transparent,var(--blue),transparent)}}
  .kpi .val{{font-family:'JetBrains Mono',monospace;font-size:1.8rem;font-weight:600;
    background:linear-gradient(135deg,var(--blue),var(--cyan));-webkit-background-clip:text;-webkit-text-fill-color:transparent}}
  .kpi .lbl{{font-size:.73rem;color:var(--muted);margin-top:6px;text-transform:uppercase;letter-spacing:.06em}}
  .kpi .sub{{font-size:.7rem;color:var(--muted);margin-top:3px;font-family:'JetBrains Mono',monospace}}

  /* Tier cards */
  .tier-grid{{display:grid;grid-template-columns:repeat(3,1fr);gap:16px;margin-bottom:28px}}
  .tier-card{{background:var(--surface);border:1px solid var(--border);border-radius:8px;padding:22px;
    position:relative;overflow:hidden;border-top:2px solid var(--tier-color)}}
  .tier-card::after{{content:'';position:absolute;top:0;left:0;bottom:0;width:100%;
    background:linear-gradient(135deg,color-mix(in srgb,var(--tier-color) 5%,transparent),transparent);pointer-events:none}}
  .tier-range{{font-family:'JetBrains Mono',monospace;font-size:.85rem;color:var(--tier-color);margin-bottom:8px}}
  .tier-count{{font-family:'JetBrains Mono',monospace;font-size:2rem;font-weight:600;color:var(--text)}}
  .tier-pct{{font-size:.72rem;color:var(--muted);margin-top:2px;font-family:'JetBrains Mono',monospace}}
  .tier-label{{font-size:.82rem;font-weight:600;color:var(--tier-color);margin-top:10px}}
  .tier-desc{{font-size:.75rem;color:var(--muted);margin-top:4px;line-height:1.4}}

  /* Chart section */
  .chart-section{{background:var(--surface);border:1px solid var(--border);border-radius:8px;padding:24px;margin-bottom:24px}}
  .chart-section h2{{font-size:.85rem;font-weight:600;color:var(--text);margin-bottom:18px;
    text-transform:uppercase;letter-spacing:.06em;display:flex;align-items:center;gap:8px}}
  .chart-section h2::before{{content:'';display:block;width:3px;height:14px;
    background:linear-gradient(var(--blue),var(--cyan));border-radius:2px}}

  /* Table */
  .table-section{{background:var(--surface);border:1px solid var(--border);border-radius:8px;padding:24px}}
  .table-section h2{{font-size:.85rem;font-weight:600;color:var(--text);margin-bottom:18px;
    text-transform:uppercase;letter-spacing:.06em;display:flex;align-items:center;gap:8px}}
  .table-section h2::before{{content:'';display:block;width:3px;height:14px;
    background:linear-gradient(var(--purple),var(--cyan));border-radius:2px}}
  table{{width:100%;border-collapse:collapse;font-size:.82rem}}
  th{{background:var(--surface2);color:var(--muted);padding:10px 14px;text-align:left;
    font-size:.71rem;text-transform:uppercase;letter-spacing:.06em;border-bottom:1px solid var(--border)}}
  td{{padding:10px 14px;border-bottom:1px solid var(--border)}}
  tr:hover td{{background:var(--surface2)}}
  .mono{{font-family:'JetBrains Mono',monospace;font-size:.78rem}}
  .score-bar-wrap{{position:relative;background:var(--surface2);border-radius:3px;height:20px;width:160px;overflow:hidden;display:flex;align-items:center}}
  .score-bar{{position:absolute;left:0;top:0;bottom:0;border-radius:3px;opacity:.7}}
  .score-bar.tier-solid{{background:var(--blue)}}
  .score-bar.tier-strong{{background:var(--purple)}}
  .score-bar.tier-very-strong{{background:var(--cyan)}}
  .score-val{{position:relative;z-index:1;font-family:'JetBrains Mono',monospace;font-size:.75rem;padding-left:6px}}
  .badge{{display:inline-block;padding:2px 8px;border-radius:3px;font-size:.7rem;font-family:'JetBrains Mono',monospace}}
  .badge.tier-solid{{background:rgba(59,130,246,.15);color:var(--blue);border:1px solid rgba(59,130,246,.3)}}
  .badge.tier-strong{{background:rgba(139,92,246,.15);color:var(--purple);border:1px solid rgba(139,92,246,.3)}}
  .badge.tier-very-strong{{background:rgba(6,182,212,.15);color:var(--cyan);border:1px solid rgba(6,182,212,.3)}}
  canvas{{max-height:260px}}
</style>
</head>
<body>

<div class="header">
  <h1><span class="ot">Open Targets</span> Platform 26.03</h1>
  <div class="meta">association_overall_direct &nbsp;·&nbsp; Reporte: Mar 2026 &nbsp;·&nbsp; Genes: ENSG · Enfermedades: EFO/MONDO</div>
  <div class="filter-badge">▶ Mostrando únicamente asociaciones con score ≥ 0.3 — Evidencia sólida, fuerte y muy fuerte</div>
</div>

<div class="kpi-grid">
  <div class="kpi"><div class="kpi-glow"></div>
    <div class="val">{stats['total_assoc']:,}</div><div class="lbl">Total asociaciones (dataset completo)</div>
  </div>
  <div class="kpi"><div class="kpi-glow"></div>
    <div class="val">{total_filtered:,}</div><div class="lbl">Asociaciones score ≥ 0.3</div>
    <div class="sub">{round(total_filtered/stats['total_assoc']*100,1)}% del total</div>
  </div>
  <div class="kpi"><div class="kpi-glow"></div>
    <div class="val">{stats['unique_genes']:,}</div><div class="lbl">Genes únicos (ENSG)</div>
  </div>
  <div class="kpi"><div class="kpi-glow"></div>
    <div class="val">{stats['unique_diseases']:,}</div><div class="lbl">Enfermedades únicas</div>
  </div>
</div>

<div class="tier-grid">
  {tier_cards_html}
</div>

<div class="chart-section">
  <h2>Distribución de associationScore — score ≥ 0.3</h2>
  <canvas id="distChart"></canvas>
</div>

<div class="table-section">
  <h2>Top 10 asociaciones — mayor score</h2>
  <table>
    <tr><th>targetId (ENSG)</th><th>diseaseId</th><th>Score</th><th>Tier</th><th>evidenceCount</th></tr>
    {top10_html}
  </table>
</div>

<script>
new Chart(document.getElementById('distChart'),{{
  type:'bar',
  data:{{
    labels:['0.3 – 0.5\\nEvidencia sólida','0.5 – 0.7\\nEvidencia fuerte','0.7 – 1.0\\nEvidencia muy fuerte'],
    datasets:[{{
      data:[135000,65000,30000],
      backgroundColor:['#3b82f6','#8b5cf6','#06b6d4'],
      borderRadius:6,
      borderSkipped:false,
    }}]
  }},
  options:{{
    plugins:{{
      legend:{{display:false}},
      tooltip:{{callbacks:{{
        label:function(c){{return ' '+c.raw.toLocaleString()+' asociaciones'}}
      }}}}
    }},
    scales:{{
      y:{{
        beginAtZero:true,
        ticks:{{color:'#64748b',callback:function(v){{return v.toLocaleString()}}}},
        grid:{{color:'#1e2444'}}
      }},
      x:{{ticks:{{color:'#94a3b8',font:{{size:12}}}},grid:{{display:false}}}}
    }}
  }}
}});
</script>
</body></html>"""

out = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed/dashboard_ot_26.03_v2.html"
with open(out, "w") as f:
    f.write(html)
print(f"Guardado: {out}")
