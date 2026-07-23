DATABASE:     OmniPath — omnipath
SOURCE:       https://omnipathdb.org / https://archive.omnipathdb.org
ARCHIVOS:
  omnipath_webservice_interactions__latest.tsv.gz   -> PPI, fosforilación, DoRothEA
  omnipath_webservice_enz_sub__latest.tsv.gz        -> enzyme-substrate (defosforilación)
FORMAT:       TSV comprimido (gzip)
DATE:         JUN 25 2026
PUBLICATION:  https://doi.org/10.1093/nar/gkaf1126

RELATIONS:
  GEN-ppi-GEN   Gene -> ppi            -> Gene
  GEN-pho-GEN   Gene -> phosphorylates -> Gene
  GEN-dph-GEN   Gene -> dephosphorylates -> Gene

GENE ID:      Gene symbol (requiere mapeo a ENSG)

NOTAS:
  - Este directorio contiene los archivos FUENTE compartidos por omnipath,
    dorothea_AB y dorothea_CD. Evita duplicar datos en disco.
  - 'latest' puede cambiar; verificar versión descargada en los logs.
