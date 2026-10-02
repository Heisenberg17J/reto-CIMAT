#!/usr/bin/env bash
# Crea el entorno conda "radiomica" para el reto CIMAT.
# Uso:  bash crear_env.sh
set -euo pipefail

ENV_NAME="${1:-radiomica}"

echo ">> Creando entorno conda '$ENV_NAME' con Python 3.11 + numpy 1.26"
conda create -y -n "$ENV_NAME" -c conda-forge python=3.11 "numpy=1.26" pip setuptools wheel

PY="$(conda info --base)/envs/$ENV_NAME/bin/python"

echo ">> Instalando dependencias (numpy queda fijado en 1.x)"
"$PY" -m pip install -r requirements.txt

echo ">> Instalando pyradiomics desde codigo fuente"
# versioneer es dependencia de build; --no-build-isolation obliga a compilar
# contra el numpy 1.26 ya instalado en el entorno.
"$PY" -m pip install versioneer==0.29
"$PY" -m pip install --no-build-isolation --no-cache-dir pyradiomics

echo ">> Registrando kernel de Jupyter"
"$PY" -m ipykernel install --user --name "$ENV_NAME" --display-name "Python ($ENV_NAME)"

echo ">> Verificando"
"$PY" -c "import radiomics, SimpleITK, nibabel, numpy; print('pyradiomics', radiomics.__version__, '| numpy', numpy.__version__)"

echo
echo "Listo. Activa con:  conda activate $ENV_NAME"
