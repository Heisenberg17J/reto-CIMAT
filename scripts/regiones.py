"""
BLOQUE 4 - REGIONES DE INTERES
Reto CIMAT / BraTS

Genera, a partir de seg, las tres regiones compuestas estandar de BraTS
como mascaras binarias (valor 1, igual que `label: 1` en el params):

    WT (tumor completo)   = etiquetas 1 + 2 + 4
    TC (nucleo tumoral)   = etiquetas 1 + 4
    ET (tumor con realce) = etiqueta 4

Ademas define el plan de extraccion que usara el Bloque 5, evitando duplicar
las caracteristicas de forma: como las 4 modalidades estan co-registradas, la
forma de una region es la misma en las 4. Por eso:
    shape                       -> 1 vez por region            (3 por paciente)
    firstorder + texturas       -> 1 vez por modalidad x region (12 por paciente)

Entrada  (no se modifica):
    data_normalizada/manifest.csv
    params_brats2018_v0.yaml      <- de aqui salen minimumROISize/Dimensions

Salida:
    data_normalizada/<paciente>/<paciente>_mask_<WT|TC|ET>.nii.gz
    data_normalizada/manifest_regiones.csv   <- una fila por paciente x region
    logs/regiones_<fecha>.csv

Una region vacia o demasiado pequena no se descarta en silencio: queda en el
manifest con estado "vacia" / "pequena" y el Bloque 5 la salta (sus
caracteristicas quedaran como NaN). Que hacer con ellas es DECISION ABIERTA.

Uso (desde la raiz del repo):
    python scripts/regiones.py
    python scripts/regiones.py --sobrescribir
"""

import argparse
from datetime import datetime
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
import yaml

import paciente as p

MANIFEST_ENTRADA = Path("data_normalizada/manifest.csv")
MANIFEST_SALIDA = Path("data_normalizada/manifest_regiones.csv")
PARAMS = Path("params_brats2018_v0.yaml")
LOGS = Path("logs")

MODALIDADES = ["t1", "t1ce", "t2", "flair"]

# Regiones compuestas de BraTS -> etiquetas de seg que incluyen
REGIONES = {
    "WT": (1, 2, 4),
    "TC": (1, 4),
    "ET": (4,),
}

# Imagen que se le pasa a PyRadiomics para shape. La forma solo depende de la
# mascara, pero la API exige una imagen; cualquiera co-registrada sirve.
MODALIDAD_FORMA = "t1"


def leer_umbrales(params=PARAMS):
    """Lee del params los mismos umbrales que aplicara PyRadiomics."""
    with open(params) as f:
        setting = yaml.safe_load(f).get("setting", {})
    return (int(setting.get("minimumROISize", 1)),
            int(setting.get("minimumROIDimensions", 2)))


def ruta_mascara(pac, region):
    return pac.seg.parent / f"{pac.paciente_id}_mask_{region}.nii.gz"


def estado_region(mascara, min_voxeles, min_dims):
    """Clasifica la mascara igual que lo haria PyRadiomics, pero sin excepcion."""
    n = int(mascara.sum())
    if n == 0:
        return n, 0, "vacia"
    idx = np.argwhere(mascara)
    tam_bbox = idx.max(axis=0) - idx.min(axis=0) + 1
    dims = int((tam_bbox > 1).sum())
    if n < min_voxeles or dims < min_dims:
        return n, dims, "pequena"
    return n, dims, "ok"


def procesar_paciente(pac, umbrales, sobrescribir):
    """Crea las 3 mascaras de un paciente. Devuelve una fila por region."""
    seg_img = nib.load(pac.seg)
    seg = np.asanyarray(seg_img.dataobj).astype(np.int16)
    vol_voxel = float(np.prod(seg_img.header.get_zooms()[:3]))
    filas = []

    for region, etiquetas in REGIONES.items():
        salida = ruta_mascara(pac, region)
        fila = {"paciente_id": pac.paciente_id, "region": region,
                "etiquetas": "+".join(map(str, etiquetas)),
                "mascara": salida.as_posix()}
        try:
            mascara = np.isin(seg, etiquetas)
            n, dims, estado = estado_region(mascara, *umbrales)

            if salida.exists() and not sobrescribir:
                mensaje = "ya existe"
            else:
                img = nib.Nifti1Image(mascara.astype(np.uint8), seg_img.affine, seg_img.header)
                img.set_data_dtype(np.uint8)
                nib.save(img, salida)
                mensaje = ""

            fila.update(n_voxeles=n, volumen_cm3=n * vol_voxel / 1000,
                        dims_roi=dims, estado=estado, mensaje=mensaje)
        except Exception as e:
            fila.update(estado="error", mensaje=str(e))

        filas.append(fila)
        print(f"  {region}  {fila.get('n_voxeles', 0):>8,} vox  {fila['estado']}"
              f"{'  (' + fila['mensaje'] + ')' if fila['mensaje'] else ''}")

    # Coherencia de la jerarquia: ET dentro de TC dentro de WT
    n = {f["region"]: f.get("n_voxeles", 0) for f in filas}
    if not (n["ET"] <= n["TC"] <= n["WT"]):
        raise ValueError(f"jerarquia ET <= TC <= WT rota: {n}")

    return filas


# ---------------------------------------------------------------------------
# Plan de extraccion (lo consume el Bloque 5)
# ---------------------------------------------------------------------------

def crear_extractores(params=PARAMS):
    """Devuelve (extractor_forma, extractor_intensidad) a partir del mismo params.

    - forma:       solo la clase shape.
    - intensidad:  todas las clases del params excepto shape.
    Asi ningun ajuste (binWidth, label, etc.) puede divergir entre los dos.
    """
    from radiomics import featureextractor

    forma = featureextractor.RadiomicsFeatureExtractor(str(params))
    forma.disableAllFeatures()
    forma.enableFeatureClassByName("shape")

    intensidad = featureextractor.RadiomicsFeatureExtractor(str(params))
    intensidad.enabledFeatures.pop("shape", None)

    return forma, intensidad


def plan_extraccion(manifest_regiones=MANIFEST_SALIDA, manifest_imagenes=MANIFEST_ENTRADA):
    """Lista de tareas: una fila por llamada a PyRadiomics.

    Columnas: paciente_id, region, tipo ('forma' | 'intensidad'), modalidad,
    imagen, mascara. Solo incluye regiones con estado 'ok'.
    """
    regiones = pd.read_csv(manifest_regiones)
    imagenes = pd.read_csv(manifest_imagenes).set_index("paciente_id")

    tareas = []
    for _, r in regiones[regiones["estado"] == "ok"].iterrows():
        rutas = imagenes.loc[r["paciente_id"]]
        base = {"paciente_id": r["paciente_id"], "region": r["region"], "mascara": r["mascara"]}
        tareas.append({**base, "tipo": "forma", "modalidad": None,
                       "imagen": rutas[MODALIDAD_FORMA]})
        for mod in MODALIDADES:
            tareas.append({**base, "tipo": "intensidad", "modalidad": mod,
                           "imagen": rutas[mod]})
    return pd.DataFrame(tareas, columns=["paciente_id", "region", "tipo",
                                         "modalidad", "imagen", "mascara"])


def main():
    parser = argparse.ArgumentParser(description="Mascaras WT/TC/ET a partir de seg")
    parser.add_argument("--sobrescribir", action="store_true",
                        help="vuelve a generar mascaras que ya existen")
    args = parser.parse_args()

    umbrales = leer_umbrales()
    print(f"Umbrales (de {PARAMS}): minimumROISize={umbrales[0]}, "
          f"minimumROIDimensions={umbrales[1]}")

    pacientes = p.cargar_pacientes(MANIFEST_ENTRADA)
    print(f"Pacientes en {MANIFEST_ENTRADA}: {len(pacientes)}")

    filas = []
    for pac in pacientes:
        print(f"\n{pac.paciente_id}")
        filas.extend(procesar_paciente(pac, umbrales, args.sobrescribir))

    df = pd.DataFrame(filas)
    df.to_csv(MANIFEST_SALIDA, index=False)

    LOGS.mkdir(exist_ok=True)
    ruta_log = LOGS / f"regiones_{datetime.now():%Y%m%d_%H%M%S}.csv"
    df_log = df.copy()
    df_log.insert(0, "fecha", datetime.now().isoformat(timespec="seconds"))
    df_log.to_csv(ruta_log, index=False)

    plan = plan_extraccion()
    n_forma = int((plan["tipo"] == "forma").sum())
    n_int = int((plan["tipo"] == "intensidad").sum())

    print(f"\n{'=' * 62}\nRESUMEN\n{'=' * 62}")
    print(pd.crosstab(df["region"], df["estado"]).reindex(list(REGIONES)).to_string())
    print(f"\nPlan de extraccion: {n_forma} llamadas de forma + "
          f"{n_int} de intensidad = {len(plan)}")
    print(f"  (sin deduplicar serian {n_forma * len(MODALIDADES)} extracciones de forma)")
    print(f"\nManifest: {MANIFEST_SALIDA}")
    print(f"Log:      {ruta_log}")


if __name__ == "__main__":
    main()
