"""
VARIANTE DE SEGMENTACION (bloques 4-7)
Reto CIMAT / BraTS

El mismo pipeline de radiomica corre sobre dos segmentaciones distintas:
    manual  -> la segmentacion de BraTS (por defecto; lo de siempre)
    pred    -> las mascaras fuera de fold de nnU-Net (segmentaciones_pred/, D24)

Se elige con la variable de entorno SEGMENTACION, que vale para toda la cadena:
    SEGMENTACION=pred python scripts/regiones.py && python scripts/extraccion.py ...

Las imagenes normalizadas son las mismas en las dos variantes; solo cambian la
mascara, el manifest de regiones, la carpeta de resultados y el nombre de los logs,
para que una variante nunca sobrescriba a la otra.
"""

import os
from pathlib import Path

SEGMENTACION = os.environ.get("SEGMENTACION", "manual")
if SEGMENTACION not in ("manual", "pred"):
    raise SystemExit(f"SEGMENTACION={SEGMENTACION!r} no es valida: usa 'manual' o 'pred'")

SEG_PRED = Path("segmentaciones_pred")

if SEGMENTACION == "manual":
    MANIFEST_REGIONES = Path("data_normalizada/manifest_regiones.csv")
    RESULTADOS = Path("resultados")
    PREFIJO_LOG = ""
    DESCRIPCION = "manual (BraTS)"
else:
    MANIFEST_REGIONES = SEG_PRED / "manifest_regiones.csv"
    RESULTADOS = Path("resultados/pred")
    PREFIJO_LOG = "pred_"
    DESCRIPCION = "predicha fuera de fold (nnU-Net, D24)"


def ruta_seg(pac):
    """Segmentacion de la que salen las regiones WT/TC/ET."""
    if SEGMENTACION == "manual":
        return pac.seg
    return SEG_PRED / f"{pac.paciente_id}_seg.nii.gz"


def ruta_mascara(pac, region):
    if SEGMENTACION == "manual":
        return pac.seg.parent / f"{pac.paciente_id}_mask_{region}.nii.gz"
    return SEG_PRED / "mascaras" / f"{pac.paciente_id}_mask_{region}.nii.gz"
