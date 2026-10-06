"""
OBJETIVO 1 - SEGMENTACION: metricas oficiales de BraTS 2018 sobre las mascaras finales
Reto CIMAT / BraTS 2018

Compara las 285 mascaras finales fuera de fold (D24) con la segmentacion manual,
para WT, TC y ET, con las metricas del reto (Bakas et al. 2018, arXiv:1811.02629):
    Dice, Hausdorff 95 (mm), sensibilidad y especificidad.

Hausdorff 95: el maximo de los dos percentiles 95 dirigidos entre las superficies
(prediccion -> real y real -> prediccion), en mm (voxel de 1 mm). Convencion para
regiones vacias, como la herramienta de evaluacion de BraTS: si no hay region ni
real ni predicha, Dice = 1 y HD95 = 0; si falta en una sola, Dice = 0 y
HD95 = 373.13 mm (la diagonal del volumen 240 x 240 x 155).

Entrada:
    datos/segmentaciones_pred/<paciente>_seg.nii.gz   (etiquetas BraTS 0/1/2/4)
    datos/brats2018/manifest.csv                      (segmentacion manual)
    particiones/folds.csv

Salida:
    resultados/segmentacion/metricas_oof.csv   <- una fila por paciente x region
    resultados/segmentacion/metricas_resumen.csv

Uso (desde la raiz del repo):
    python segmentacion/metricas_segmentacion.py
"""

import os
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy import ndimage

PRED = Path("datos/segmentaciones_pred")
MANIFEST = Path("datos/brats2018/manifest.csv")
FOLDS = Path("particiones/folds.csv")
SALIDA = Path("resultados/segmentacion")

REGIONES = {"WT": (1, 2, 4), "TC": (1, 4), "ET": (4,)}
HD_VACIO = 373.13


def superficie(m):
    return m & ~ndimage.binary_erosion(m)


def hd95(a, b):
    """Hausdorff 95 simetrico entre dos mascaras no vacias (voxel de 1 mm)."""
    # Recorte a la caja que contiene a las dos, con margen, para acelerar
    idx = np.argwhere(a | b)
    lo, hi = np.maximum(idx.min(0) - 2, 0), idx.max(0) + 3
    a, b = a[tuple(slice(l, h) for l, h in zip(lo, hi))], b[tuple(slice(l, h) for l, h in zip(lo, hi))]
    sa, sb = superficie(a), superficie(b)
    d_a = ndimage.distance_transform_edt(~sb)[sa]     # superficie de a -> superficie de b
    d_b = ndimage.distance_transform_edt(~sa)[sb]
    return float(max(np.percentile(d_a, 95), np.percentile(d_b, 95)))


def metricas(pid, ruta_gt):
    p = np.asanyarray(nib.load(PRED / f"{pid}_seg.nii.gz").dataobj)
    g = np.asanyarray(nib.load(ruta_gt).dataobj)
    filas = []
    for region, etiquetas in REGIONES.items():
        mp, mg = np.isin(p, etiquetas), np.isin(g, etiquetas)
        tp, n_p, n_g = int((mp & mg).sum()), int(mp.sum()), int(mg.sum())
        if n_p == 0 and n_g == 0:
            dice, hd = 1.0, 0.0
        elif n_p == 0 or n_g == 0:
            dice, hd = 0.0, HD_VACIO
        else:
            dice, hd = 2 * tp / (n_p + n_g), hd95(mp, mg)
        filas.append({"paciente_id": pid, "region": region, "voxeles_real": n_g, "voxeles_pred": n_p,
                      "dice": dice, "hd95_mm": hd,
                      "sensibilidad": tp / n_g if n_g else np.nan,
                      "especificidad": float(((~mp) & (~mg)).sum() / (~mg).sum())})
    return filas


def main():
    m = pd.read_csv(MANIFEST).set_index("paciente_id")
    grupo = pd.read_csv(FOLDS).set_index("paciente_id")["grupo"]
    res = Parallel(n_jobs=max(1, (os.cpu_count() or 2) - 2))(
        delayed(metricas)(pid, m.loc[pid, "seg"]) for pid in m.index)
    t = pd.DataFrame([f for r in res for f in r])
    t["grupo"] = t["paciente_id"].map(grupo)
    SALIDA.mkdir(parents=True, exist_ok=True)
    t.to_csv(SALIDA / "metricas_oof.csv", index=False)

    cols = ["dice", "hd95_mm", "sensibilidad", "especificidad"]
    filas = []
    for nombre, sub in (("todos (285)", t), ("HGG con supervivencia (163)", t[t.grupo == "HGG_superv"])):
        for region in REGIONES:
            x = sub[sub.region == region]
            filas.append({"conjunto": nombre, "region": region,
                          **{f"{c}_media": x[c].mean() for c in cols},
                          "dice_mediana": x.dice.median(), "hd95_mediana": x.hd95_mm.median()})
    resumen = pd.DataFrame(filas)
    resumen.to_csv(SALIDA / "metricas_resumen.csv", index=False)
    print(resumen.round(3).to_string(index=False))
    print(f"\nHD95 = {HD_VACIO} mm (region ausente en una sola de las dos) en "
          f"{int((t.hd95_mm == HD_VACIO).sum())} casos")
    print(f"\ntablas: {SALIDA}/metricas_oof.csv, metricas_resumen.csv")


if __name__ == "__main__":
    os.chdir(Path(__file__).resolve().parent.parent)
    main()
