DATABASE:     COXPRESdb — coexpressdb
VERSION:      v22-05 (Hsa-r.c6-0)
SOURCE:       https://coxpresdb.jp
ARCHIVOS:
  Hsa-r.v22-05.G16651-S235187.combat_pca.subagging.z.d.zip  → coexpresión
  supportability.ver8-0.tar.bz2                              → supportability scores

RELATIONS:
  GEN-cex-GEN   Gene → coexpresses → Gene

GENE ID:      Entrez Gene ID (requiere mapeo a ENSG)
ASSOC ID:     Entrez Gene ID del gen coexpresado

⚠️  VERIFICAR:
  El URL del archivo de coexpresión fue inferido del checksum .md5 encontrado
  en el directorio. Confirmar que el archivo ZIP es descargable antes de
  lanzar el job de SLURM.
  URL a verificar: https://coxpresdb.jp/download/Hsa-r.c6-0/coex_md5/

NOTAS:
  - Bioteque usó Hsa-r.c4-0 (2014). Esta es la versión más reciente (2022).
  - El archivo de supportability indica el número de datasets que apoyan
    cada par de coexpresión (score de confianza).
