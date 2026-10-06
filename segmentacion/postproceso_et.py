"""
OBJETIVO 1 - SEGMENTACION: posprocesado de ET y mascaras finales fuera de fold
Reto CIMAT / BraTS 2018

nnU-Net a veces predice un poco de realce (ET) en pacientes que no lo tienen
(LGG sin ET). En BraTS eso vale Dice 0 y hunde el promedio. La correccion
estandar: si el ET predicho tiene menos de T voxeles, se reetiqueta como
necrosis/no realzado. TC y WT no cambian; solo ET.

El umbral T se elige con VALIDACION ANIDADA, para no ajustarlo a los mismos
pacientes con los que se mide:
    para cada fold k: T_k = el mejor umbral en los otros 4 folds,
                      y se aplica (y se evalua) en el fold k.
Asi cada mascara final se corrige con un umbral que no se eligio mirandola.
El umbral elegido con los 5 folds juntos queda como el de casos nuevos.

Convencion de Dice (la oficial de BraTS): si no hay ET real ni predicho, Dice = 1.
Tambien se reporta la de nnU-Net (ese caso no cuenta).

Entrada:
    datos/predicciones_oof/fold_k/validation/<paciente>.nii.gz   (etiquetas nnU-Net 0/1/2/3)
    datos/nnunet_raw/Dataset501_BraTS2018/labelsTr/               (segmentacion real)
    particiones/folds.csv

Salida:
    datos/segmentaciones_pred/<paciente>_seg.nii.gz   <- mascara final, etiquetas BraTS 0/1/2/4
    resultados/segmentacion/segmentacion_oof.csv             <- Dice por paciente, antes y despues

Uso (desde la raiz del repo):
    python segmentacion/postproceso_et.py
"""

import glob
import json
import os
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

from convertir_nnunet import NNUNET_A_BRATS, reasignar

PRED = Path("datos/predicciones_oof")
GT = Path("datos/nnunet_raw/Dataset501_BraTS2018/labelsTr")
FOLDS = Path("particiones/folds.csv")
SALIDA = Path("datos/segmentaciones_pred")
TABLA = Path("resultados/segmentacion/segmentacion_oof.csv")
DECISION = Path("resultados/segmentacion/postproceso_et.json")

UMBRALES = [0, 25, 50, 100, 150, 200, 300, 400, 500, 750, 1000]


def dice(a, b, vacio):
    """Dice; 'vacio' es el valor cuando la region no existe ni real ni predicha."""
    s = a.sum() + b.sum()
    return vacio if s == 0 else 2 * np.logical_and(a, b).sum() / s


def medir():
    """Una fila por paciente: Dice de WT/TC y de ET conservando o descartando el ET predicho."""
    folds = pd.read_csv(FOLDS).set_index("paciente_id")["fold"]
    rutas = {Path(p).name[:-7]: p for p in glob.glob(str(PRED / "fold_*/validation/*.nii.gz"))}
    if set(rutas) != set(folds.index):
        raise SystemExit(f"predicciones: {len(rutas)}, pacientes: {len(folds)}; no coinciden")

    filas = []
    for i, (pid, ruta) in enumerate(sorted(rutas.items()), 1):
        p = np.asanyarray(nib.load(ruta).dataobj)
        g = np.asanyarray(nib.load(GT / f"{pid}.nii.gz").dataobj)
        et_g, et_p = g == 3, p == 3
        vacio_et = np.zeros_like(et_p)
        filas.append({
            "paciente_id": pid, "fold": int(folds[pid]),
            "et_real": int(et_g.sum()), "et_pred": int(et_p.sum()),
            "dice_WT": dice(g >= 1, p >= 1, 1.0),
            "dice_TC": dice(g >= 2, p >= 2, 1.0),
            # BraTS: ET vacio real y predicho = 1;  nnU-Net: ese caso es NaN
            "et_conservar": dice(et_g, et_p, 1.0),
            "et_descartar": dice(et_g, vacio_et, 1.0),
            "et_conservar_nnunet": dice(et_g, et_p, np.nan),
            "et_descartar_nnunet": dice(et_g, vacio_et, np.nan),
        })
        if i % 50 == 0:
            print(f"  {i}/{len(rutas)}")
    return pd.DataFrame(filas).set_index("paciente_id")


def dice_et(t, umbral, nnunet=False):
    """Dice ET por paciente si se descarta el ET predicho con menos de 'umbral' voxeles."""
    suf = "_nnunet" if nnunet else ""
    descartar = t["et_pred"] < umbral
    return t[f"et_descartar{suf}"].where(descartar, t[f"et_conservar{suf}"])


def mejor_umbral(t):
    medias = {u: dice_et(t, u).mean() for u in UMBRALES}
    # Empate: el umbral mas chico (descarta menos ET reales)
    return max(medias, key=lambda u: (round(medias[u], 6), -u)), medias


def main():
    t = medir()

    print(f"\n{'=' * 70}\nPACIENTES SIN ET REAL\n{'=' * 70}")
    sin = t[t.et_real == 0]
    fp = sin[sin.et_pred > 0]
    print(f"{len(sin)} sin ET real; el modelo predice ET en {len(fp)}:")
    print(fp[["fold", "et_pred"]].sort_values("et_pred").to_string())
    chicos = t[(t.et_real > 0) & (t.et_pred < 1000)]
    print(f"\npacientes CON ET real y ET predicho < 1000 voxeles (riesgo de perderlos): {len(chicos)}")
    print(chicos[["fold", "et_real", "et_pred", "et_conservar"]].round(3).to_string())

    print(f"\n{'=' * 70}\nVALIDACION ANIDADA DEL UMBRAL\n{'=' * 70}")
    t["umbral"] = 0
    filas = []
    for k in sorted(t.fold.unique()):
        u, _ = mejor_umbral(t[t.fold != k])
        t.loc[t.fold == k, "umbral"] = u
        f = t[t.fold == k]
        filas.append({"fold": k, "umbral elegido en los otros 4": u,
                      "ET antes": dice_et(f, 0).mean(), "ET despues": dice_et(f, u).mean(),
                      "ET antes (nnU-Net)": dice_et(f, 0, True).mean(),
                      "ET despues (nnU-Net)": dice_et(f, u, True).mean()})
    anidada = pd.DataFrame(filas).set_index("fold")
    print(anidada.round(3).to_string())

    t["dice_ET"] = t["et_descartar"].where(t.et_pred < t.umbral, t["et_conservar"])
    t["dice_ET_nnunet"] = t["et_descartar_nnunet"].where(t.et_pred < t.umbral, t["et_conservar_nnunet"])
    final_u, medias = mejor_umbral(t)
    print("\nDice ET medio (BraTS) de los 285 segun el umbral:",
          {u: round(m, 3) for u, m in medias.items()})
    print(f"umbral para casos nuevos (elegido con los 5 folds): {final_u} voxeles")

    print(f"\n{'=' * 70}\nRESULTADO FUERA DE FOLD (285 pacientes, umbral anidado)\n{'=' * 70}")
    res = pd.DataFrame({
        "antes": [t.dice_WT.mean(), t.dice_TC.mean(), dice_et(t, 0).mean(), dice_et(t, 0, True).mean()],
        "despues": [t.dice_WT.mean(), t.dice_TC.mean(), t.dice_ET.mean(), t.dice_ET_nnunet.mean()],
    }, index=["WT", "TC", "ET (BraTS)", "ET (nnU-Net)"])
    print(res.round(3).to_string())
    perdidos = t[(t.et_real > 0) & (t.et_pred > 0) & (t.et_pred < t.umbral)]
    corregidos = t[(t.et_real == 0) & (t.et_pred > 0) & (t.et_pred < t.umbral)]
    print(f"\nET falsos eliminados: {len(corregidos)} | ET reales perdidos: {len(perdidos)}")

    # Mascaras finales en etiquetas BraTS
    SALIDA.mkdir(exist_ok=True)
    for pid, fila in t.iterrows():
        img = nib.load(PRED / f"fold_{int(fila.fold)}" / "validation" / f"{pid}.nii.gz")
        seg = np.asanyarray(img.dataobj).astype(np.int16)
        if fila.et_pred < fila.umbral:
            seg[seg == 3] = 2          # ET -> necrosis/no realzado: TC y WT no cambian
        nueva = nib.Nifti1Image(reasignar(seg, NNUNET_A_BRATS), img.affine, img.header)
        nueva.set_data_dtype(np.uint8)
        nib.save(nueva, SALIDA / f"{pid}_seg.nii.gz")

    TABLA.parent.mkdir(exist_ok=True)
    t[["fold", "et_real", "et_pred", "umbral", "dice_WT", "dice_TC", "dice_ET", "dice_ET_nnunet"]] \
        .to_csv(TABLA)
    DECISION.write_text(json.dumps({
        "umbral_casos_nuevos": int(final_u),
        "umbral_por_fold": {int(k): int(v) for k, v in anidada["umbral elegido en los otros 4"].items()},
        "dice_oof": {k: round(float(v), 4) for k, v in res["despues"].items()},
    }, indent=2))
    print(f"\nmascaras: {SALIDA}/ ({len(t)})\ntabla:    {TABLA}\numbral:   {DECISION}")


if __name__ == "__main__":
    os.chdir(Path(__file__).resolve().parent.parent)
    main()
