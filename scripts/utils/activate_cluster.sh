#!/bin/bash
set -e

# Cargar Python del cluster
module load Python/3.11.3-GCCcore-12.3.0

# Ir al directorio de trabajo
cd "$(dirname "$0")"

# Activar ambiente virtual
source .venv/bin/activate

# Verificar
echo "✓ Python activado:"
python --version
echo "✓ Ambiente: $VIRTUAL_ENV"

