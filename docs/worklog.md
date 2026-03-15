## 2026-03-15
### Objective
Inspect and process a Bioteque-compatible dataset for TFM-NetActivity.

### Dataset selected
Open Targets (`opentargets`), release 22.04.

### Why this dataset
DisGeNET curated download script from Bioteque was outdated and no longer downloaded the real compressed TSV file. Open Targets FTP access worked correctly.

### Commands executed
- Cloned Bioteque in `/beegfs/home/ppadmoremcc/work/external/bioteque`
- Inspected `README.txt`, `get_data.sh`, and `script.py`
- Downloaded:
  - `diseases/`
  - `associationByOverallIndirect/`
- Created custom script:
  - `scripts/opentargets_gene_summary.py`

### Results
- Total associations: 7,544,477
- Unique genes: 29,162
- Gene ID type: Ensembl Gene ID (`targetId`)
- Per-gene summary saved to:
  - `data/processed/opentargets_gene_summary.tsv`

### Notes
- Open Targets was processed from Bioteque external clone, not from inside the thesis repo.
- TFM outputs were saved inside `TFM-NetActivity/data/`.
