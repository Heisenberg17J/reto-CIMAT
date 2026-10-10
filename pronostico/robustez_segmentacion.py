"""
OBJETIVO 2 - PRONOSTICO, fase B: robustez de la radiomica frente a la segmentacion
Reto CIMAT / BraTS 2018

Compara cada caracteristica calculada con la mascara manual (BraTS) y con la
predicha por nnU-Net (D24, D25), en los 163 HGG con supervivencia: los pacientes
del pronostico. Una caracteristica que cambia mucho segun quien dibujo la mascara
no es confiable para el modelo.

Dos medidas por caracteristica (corte habitual en radiomica: 0.85):
    spearman  rho de Spearman: se conserva el orden de los pacientes.
              ES EL CRITERIO DEL FILTRO (D26): entran al modelo las caracteristicas
              con spearman >= 0.85 (935 de 1146). Lo aplica datos.cargar_cohorte.
    ccc       concordancia de Lin: baja tambien si la mascara predicha desplaza o
              cambia la escala de los valores. Es DESCRIPTIVA: mide el sesgo de la
              mascara predicha y prohibe mezclar variantes, pero no filtra. La
              columna 'robusta' del CSV es ccc >= 0.85 y no la usa el modelo.
Por que Spearman: el modelo se entrena y evalua siempre con la variante predicha,
asi que un desplazamiento sistematico no le hace dano; si se lo hace que el
contorno reordene a los pacientes.

No usa la supervivencia: es un filtro sin etiqueta, como el de columnas
constantes (D17), y no filtra informacion del objetivo hacia la validacion.

Entrada:
    resultados/radiomica/manual/features.csv        (mascara manual)
    resultados/radiomica/pred/features.csv   (mascara predicha)
    particiones/folds.csv          (grupo HGG_superv)

Salida:
    resultados/pronostico/robustez_segmentacion.csv   <- una fila por caracteristica

Uso (desde la raiz del repo):
    python pronostico/robustez_segmentacion.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

MANUAL = Path("resultados/radiomica/manual/features.csv")
PRED = Path("resultados/radiomica/pred/features.csv")
FOLDS = Path("particiones/folds.csv")
SALIDA = Path("resultados/pronostico/robustez_segmentacion.csv")
CORTE = 0.85


def ccc(x, y):
    """Coeficiente de concordancia de Lin."""
    mx, my = x.mean(), y.mean()
    vx, vy = x.var(), y.var()            # poblacional (ddof=0), como en la definicion
    cov = ((x - mx) * (y - my)).mean()
    denom = vx + vy + (mx - my) ** 2
    return np.nan if denom == 0 else 2 * cov / denom


def columnas_robustas(manual, pred, corte=CORTE):
    """Columnas con Spearman (manual vs predicha) >= corte en estos pacientes.

    Mismo calculo que la columna 'spearman' de main(), pero sobre las filas que
    se le pasen: evaluar.py --filtro-por-fold la llama solo con el entrenamiento
    externo de cada fold. Un Spearman indefinido (columna constante) no pasa.
    """
    if list(manual.columns) != list(pred.columns) or list(manual.index) != list(pred.index):
        raise ValueError("las tablas manual y predicha no estan alineadas")
    return [col for col in manual.columns
            if spearmanr(manual[col].to_numpy(float), pred[col].to_numpy(float))[0] >= corte]


def main():
    grupo = pd.read_csv(FOLDS).set_index("paciente_id")["grupo"]
    ids = grupo.index[grupo == "HGG_superv"]
    m = pd.read_csv(MANUAL).set_index("paciente_id").loc[ids]
    p = pd.read_csv(PRED).set_index("paciente_id").loc[ids]
    if list(m.columns) != list(p.columns):
        raise SystemExit("las tablas manual y predicha no tienen las mismas columnas")
    if m.isna().any().any() or p.isna().any().any():
        raise SystemExit("hay NaN en los 163 HGG con supervivencia; revisar D25")

    filas = []
    for col in m.columns:
        modalidad, region, clase, nombre = col.split("_")
        x, y = m[col].to_numpy(float), p[col].to_numpy(float)
        filas.append({"columna": col, "modalidad": modalidad, "region": region, "clase": clase,
                      "nombre": nombre, "ccc": ccc(x, y), "spearman": spearmanr(x, y)[0],
                      "dif_media_rel": (y.mean() - x.mean()) / abs(x.mean()) if x.mean() else np.nan})
    r = pd.DataFrame(filas)
    r["robusta"] = r["ccc"] >= CORTE
    SALIDA.parent.mkdir(exist_ok=True)
    r.to_csv(SALIDA, index=False)

    print(f"{len(ids)} HGG con supervivencia | {len(r)} caracteristicas | corte CCC >= {CORTE}")
    print(f"robustas: {int(r.robusta.sum())} ({r.robusta.mean():.0%})  |  CCC mediano: {r.ccc.median():.3f}")
    print(f"filtro del modelo (D26), Spearman >= {CORTE}: {int((r.spearman >= CORTE).sum())} caracteristicas")

    print("\n% de caracteristicas robustas por region y clase:")
    tabla = r.pivot_table(index="clase", columns="region", values="robusta", aggfunc="mean")
    tabla = tabla.reindex(columns=["WT", "TC", "ET"])
    print((tabla * 100).round(0).astype(int).to_string())

    print("\nCCC mediano por region y modalidad (forma = mask):")
    print(r.pivot_table(index="modalidad", columns="region", values="ccc", aggfunc="median")
          .reindex(index=["mask", "t1", "t1ce", "t2", "flair"], columns=["WT", "TC", "ET"])
          .round(3).to_string())

    print("\nforma (no depende de la modalidad):")
    forma = r[r.clase == "shape"].pivot_table(index="nombre", columns="region", values="ccc")
    print(forma.reindex(columns=["WT", "TC", "ET"]).round(3).to_string())

    print("\nlas 10 menos robustas:")
    print(r.nsmallest(10, "ccc")[["columna", "ccc", "spearman", "dif_media_rel"]].round(3)
          .to_string(index=False))
    print(f"\ntabla: {SALIDA}")


if __name__ == "__main__":
    main()
