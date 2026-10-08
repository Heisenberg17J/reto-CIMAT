"""
OBJETIVO 1 - SEGMENTACION: figura de mascara manual frente a mascara predicha
Reto CIMAT / BraTS 2018

Por paciente, un corte axial con tres paneles: FLAIR | FLAIR + segmentacion manual
(radiologo) | FLAIR + mascara final fuera de fold (D24, ya con el posprocesado de ET).
Son las mismas mascaras con las que se calcularon las metricas de
resultados/segmentacion/metricas_oof.csv, y el titulo lleva su Dice.

Sin --pacientes elige tres casos representativos por el Dice medio (WT, TC, ET):
el mejor, el de la mediana y el peor del grupo pedido.

Entrada:
    datos/segmentaciones_pred/<paciente>_seg.nii.gz   (etiquetas BraTS 0/1/2/4)
    datos/brats2018/manifest.csv                      (FLAIR y segmentacion manual)
    resultados/segmentacion/metricas_oof.csv          (Dice y grupo)

Salida:
    resultados/segmentacion/figuras/<paciente>.png   <- una por paciente
    resultados/segmentacion/figuras/comparacion.png  <- todos juntos, una fila por paciente

Uso (desde la raiz del repo):
    python segmentacion/figura_comparacion.py                       # mejor, mediana, peor (285)
    python segmentacion/figura_comparacion.py --grupo HGG_superv     # solo HGG con supervivencia
    python segmentacion/figura_comparacion.py --pacientes Brats18_2013_10_1 Brats18_CBICA_AAP_1
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

PRED = Path("datos/segmentaciones_pred")
MANIFEST = Path("datos/brats2018/manifest.csv")
METRICAS = Path("resultados/segmentacion/metricas_oof.csv")
SALIDA = Path("resultados/segmentacion/figuras")

GRUPOS = ("HGG_superv", "HGG_sin_superv", "LGG")
# Etiquetas BraTS -> indice de color: 1 necrosis/no realzado, 2 edema, 4 realce (ET)
A_COLOR = np.zeros(5, dtype=np.uint8)
A_COLOR[[1, 2, 4]] = [1, 2, 3]
COLORES = ListedColormap([(0, 0, 0, 0), (0.9, 0.1, 0.1, 1), (0.1, 0.8, 0.2, 1), (1.0, 0.85, 0.0, 1)])
LEYENDA = [Patch(color=COLORES(2), label="edema"),
           Patch(color=COLORES(1), label="necrosis / no realzado"),
           Patch(color=COLORES(3), label="realce (ET)")]


def cargar(ruta):
    return np.asanyarray(nib.load(ruta).dataobj)


def elegir(metricas, grupo):
    """Mejor, mediana y peor paciente por el Dice medio de WT, TC y ET."""
    d = metricas if grupo is None else metricas[metricas.grupo == grupo]
    media = d.pivot(index="paciente_id", columns="region", values="dice").mean(axis=1).sort_values()
    return [media.index[-1], media.index[len(media) // 2], media.index[0]]


def preparar(pid, manifest):
    flair = cargar(manifest.loc[pid, "flair"]).astype(np.float32)
    gt, pred = cargar(manifest.loc[pid, "seg"]), cargar(PRED / f"{pid}_seg.nii.gz")
    ref = gt if gt.any() else pred
    z = int(np.argmax((ref > 0).sum(axis=(0, 1)))) if ref.any() else flair.shape[2] // 2
    # Recorte a la caja del cerebro en ese corte, con margen
    idx = np.argwhere(flair[:, :, z] > 0)
    (x0, y0), (x1, y1) = np.maximum(idx.min(0) - 5, 0), idx.max(0) + 6
    corte = lambda a: np.rot90(a[x0:x1, y0:y1, z])
    f = corte(flair)
    lo, hi = np.percentile(f[f > 0], [1, 99.5]) if (f > 0).any() else (0, 1)
    return z, np.clip((f - lo) / (hi - lo), 0, 1), A_COLOR[corte(gt)], A_COLOR[corte(pred)]


def dibujar_fila(ejes, pid, datos, dice, grupo, alfa):
    z, f, gt, pred = datos
    for eje, (nombre, seg) in zip(ejes, [("FLAIR", None), ("Manual (radiólogo)", gt),
                                         ("Predicha (nnU-Net, fuera de fold)", pred)]):
        eje.imshow(f, cmap="gray", vmin=0, vmax=1)
        if seg is not None:
            eje.imshow(np.ma.masked_equal(seg, 0), cmap=COLORES, vmin=0, vmax=3, alpha=alfa,
                       interpolation="nearest")
        eje.set_title(nombre, fontsize=11)
        eje.axis("off")
    ejes[0].text(-0.04, 0.5, f"{pid}\n{grupo} · corte {z}\nDice WT {dice['WT']:.2f} · "
                 f"TC {dice['TC']:.2f} · ET {dice['ET']:.2f}",
                 transform=ejes[0].transAxes, ha="right", va="center", fontsize=10)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pacientes", nargs="+", help="ids de paciente (por defecto: mejor, mediana y peor)")
    ap.add_argument("--grupo", choices=GRUPOS, help="grupo del que elegir los casos (por defecto: los 285)")
    ap.add_argument("--alfa", type=float, default=0.5, help="opacidad de la segmentacion")
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()

    metricas = pd.read_csv(METRICAS)
    manifest = pd.read_csv(MANIFEST).set_index("paciente_id")
    pacientes = args.pacientes or elegir(metricas, args.grupo)
    faltan = [p for p in pacientes if p not in manifest.index or not (PRED / f"{p}_seg.nii.gz").exists()]
    if faltan:
        raise SystemExit(f"sin segmentacion manual o predicha: {faltan}")
    SALIDA.mkdir(parents=True, exist_ok=True)

    info = {}
    for pid in pacientes:
        m = metricas[metricas.paciente_id == pid]
        info[pid] = (preparar(pid, manifest), m.set_index("region").dice, m.grupo.iloc[0])
        fig, ejes = plt.subplots(1, 3, figsize=(13, 4.6))
        dibujar_fila(ejes, pid, *info[pid], args.alfa)
        fig.legend(handles=LEYENDA, loc="lower center", ncol=3, frameon=False)
        fig.tight_layout(rect=(0, 0.07, 1, 1))
        fig.savefig(SALIDA / f"{pid}.png", dpi=args.dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"{SALIDA / pid}.png")

    fig, ejes = plt.subplots(len(pacientes), 3, figsize=(13, 4.3 * len(pacientes)), squeeze=False)
    for fila, pid in zip(ejes, pacientes):
        dibujar_fila(fila, pid, *info[pid], args.alfa)
    fig.legend(handles=LEYENDA, loc="lower center", ncol=3, frameon=False)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(SALIDA / "comparacion.png", dpi=args.dpi, bbox_inches="tight")
    print(SALIDA / "comparacion.png")


if __name__ == "__main__":
    main()
