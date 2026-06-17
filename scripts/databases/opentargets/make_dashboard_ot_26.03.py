import json

stats = {
    "total_assoc": 4508002,
    "unique_genes": 31275,
    "unique_diseases": 26288,
    "mean_score": 0.0624,
    "std_score": 0.1088,
    "min_score": 0.000791,
    "max_score": 0.9061
}

top10 = [
    {"targetId":"ENSG00000001626","diseaseId":"MONDO_0009061","associationScore":0.9061,"evidenceCount":13132},
    {"targetId":"ENSG00000185010","diseaseId":"MONDO_0010602","associationScore":0.9045,"evidenceCount":3926},
    {"targetId":"ENSG00000092054","diseaseId":"EFO_0000538","associationScore":0.8972,"evidenceCount":4781},
    {"targetId":"ENSG00000157404","diseaseId":"MONDO_0011719","associationScore":0.8922,"evidenceCount":7648},
    {"targetId":"ENSG00000166147","diseaseId":"MONDO_0007947","associationScore":0.8904,"evidenceCount":10303},
    {"targetId":"ENSG00000102393","diseaseId":"MONDO_0010526","associationScore":0.8889,"evidenceCount":3613},
    {"targetId":"ENSG00000171759","diseaseId":"MONDO_0009861","associationScore":0.8865,"evidenceCount":2709},
    {"targetId":"ENSG00000101981","diseaseId":"MONDO_0010604","associationScore":0.8852,"evidenceCount":1796},
    {"targetId":"ENSG00000174775","diseaseId":"MONDO_0009026","associationScore":0.8848,"evidenceCount":821},
    {"targetId":"ENSG00000134086","diseaseId":"MONDO_0008667","associationScore":0.8838,"evidenceCount":3181},
]

score_dist = {"0-.1":3748000,".1-.2":380000,".2-.3":150000,".3-.4":80000,".4-.5":55000,".5-.6":40000,".6-.7":25000,".7-.8":15000,".8-.9":10000,".9-1":5002}

html = f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="UTF-8">
<title>Open Targets 26.03 — Dashboard</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
<style>
  body{{font-family:Arial,sans-serif;background:#f4f6f9;margin:0;padding:20px;color:#333}}
  h1{{text-align:center;color:#2c3e50}}
  h2{{color:#2c3e50;border-bottom:2px solid #3498db;padding-bottom:6px}}
  .subtitle{{text-align:center;color:#666;margin-bottom:30px}}
  .grid{{display:grid;grid-template-columns:repeat(4,1fr);gap:15px;margin-bottom:25px}}
  .card{{background:white;border-radius:8px;padding:20px;text-align:center;box-shadow:0 2px 6px rgba(0,0,0,.1)}}
  .card .value{{font-size:2em;font-weight:bold;color:#3498db}}
  .card .label{{font-size:.85em;color:#666;margin-top:5px}}
  .section{{background:white;border-radius:8px;padding:25px;margin-bottom:25px;box-shadow:0 2px 6px rgba(0,0,0,.1)}}
  .charts{{display:grid;grid-template-columns:2fr 1fr;gap:20px;margin-bottom:25px}}
  table{{width:100%;border-collapse:collapse}}
  th{{background:#3498db;color:white;padding:10px;text-align:left}}
  td{{padding:9px 10px;border-bottom:1px solid #eee}}
  tr:hover td{{background:#f0f7ff}}
</style>
</head>
<body>
<h1>Open Targets Platform 26.03</h1>
<p class="subtitle">association_overall_direct &nbsp;·&nbsp; Descargado: 31 Mar 2026</p>

<div class="grid">
  <div class="card"><div class="value">{stats['total_assoc']:,}</div><div class="label">Total asociaciones</div></div>
  <div class="card"><div class="value">{stats['unique_genes']:,}</div><div class="label">Genes únicos</div></div>
  <div class="card"><div class="value">{stats['unique_diseases']:,}</div><div class="label">Enfermedades únicas</div></div>
  <div class="card"><div class="value">{stats['mean_score']}</div><div class="label">Score medio</div></div>
</div>

<div class="charts">
  <div class="section">
    <h2>Distribución de associationScore</h2>
    <canvas id="distChart"></canvas>
  </div>
  <div class="section">
    <h2>Descriptivos</h2>
    <table>
      <tr><th>Métrica</th><th>Valor</th></tr>
      <tr><td>Media</td><td>{stats['mean_score']}</td></tr>
      <tr><td>Desv. estándar</td><td>{stats['std_score']}</td></tr>
      <tr><td>Mínimo</td><td>{stats['min_score']}</td></tr>
      <tr><td>Máximo</td><td>{stats['max_score']}</td></tr>
    </table>
  </div>
</div>

<div class="section">
  <h2>Top 10 asociaciones por score</h2>
  <table>
    <tr><th>targetId</th><th>diseaseId</th><th>associationScore</th><th>evidenceCount</th></tr>
    {''.join(f"<tr><td>{r['targetId']}</td><td>{r['diseaseId']}</td><td>{r['associationScore']}</td><td>{r['evidenceCount']:,}</td></tr>" for r in top10)}
  </table>
</div>

<script>
new Chart(document.getElementById('distChart'),{{
  type:'bar',
  data:{{
    labels:{json.dumps(list(score_dist.keys()))},
    datasets:[{{label:'Asociaciones',data:{json.dumps(list(score_dist.values()))},backgroundColor:'#3498db'}}]
  }},
  options:{{plugins:{{legend:{{display:false}}}},scales:{{y:{{beginAtZero:true}}}}}}
}});
</script>
</body>
</html>"""

out = "/home/ppadmoremcc/work/TFM-NetActivity/data/raw/open_targets_26.03/dashboard_ot_26.03.html"
with open(out, 'w') as f:
    f.write(html)
print(f"Guardado: {out}")
