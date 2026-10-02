"""
EXPERIMENTO - VARIANTE B DE NORMALIZACION
Reto CIMAT / BraTS  (rama prueba/variante-b; no forma parte del pipeline)

Compara la tabla actual (A) contra la variante B:
    A: z-score sobre voxeles del cerebro,          binWidth 0.1, voxelArrayShift 0
    B: el mismo z-score multiplicado por 100,      binWidth 25,  voxelArrayShift 300

Reutiliza sin modificar: zscore() del Bloque 3, las mascaras y el plan del
Bloque 4 y extraer_paciente() del Bloque 5. Solo cambia la imagen y el params.

Entrada:
    data/manifest.csv                     <- imagenes originales
    resultados/caracteristicas_qc.csv     <- tabla A
    experimentos/params_b.yaml

Salida (solo variante B; no toca nada de A):
    resultados/variante_b/imagenes/<paciente>/<paciente>_<mod>.nii.gz
    resultados/variante_b/caracteristicas_b.csv

Uso (desde la raiz del repo):
    python experimentos/variante_b.py
"""

import logging
import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
import radiomics

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import extraccion as e   # noqa: E402
import normalizar as n   # noqa: E402
import paciente as p     # noqa: E402
import regiones as r     # noqa: E402

ESCALA = 100
PARAMS_B = Path("experimentos/params_b.yaml")
SALIDA = Path("resultados/variante_b")
TABLA_A = Path("resultados/caracteristicas_qc.csv")
BIN_A, BIN_B = 0.1, 25
UMBRAL_RHO = 0.99


def normalizar_b():
    """Z-score del Bloque 3 x100. Devuelve {(paciente, modalidad): ruta}."""
    rutas = {}
    for pac in p.cargar_pacientes(n.MANIFEST_ENTRADA):
        destino = SALIDA / "imagenes" / pac.paciente_id
        destino.mkdir(parents=True, exist_ok=True)
        for mod in r.MODALIDADES:
            origen = pac.rutas()[mod]
            img = nib.load(origen)
            z, _ = n.zscore(img.get_fdata())
            nueva = nib.Nifti1Image((z * ESCALA).astype(np.float32), img.affine, img.header)
            nueva.set_data_dtype(np.float32)
            salida = destino / origen.name
            nib.save(nueva, salida)
            rutas[(pac.paciente_id, mod)] = salida.as_posix()
    return rutas


def extraer_b(rutas):
    forma, intensidad = r.crear_extractores(PARAMS_B)
    extractores = {"forma": forma, "intensidad": intensidad}
    plan = r.plan_extraccion()
    # Mismas tareas y mascaras que A; solo cambia la imagen de intensidad.
    # (forma usa t1 normalizada de A, pero shape solo depende de la mascara)
    es_int = plan["tipo"] == "intensidad"
    plan.loc[es_int, "imagen"] = [rutas[(pid, mod)] for pid, mod
                                  in plan.loc[es_int, ["paciente_id", "modalidad"]].values]
    filas = []
    for pid, tareas in plan.groupby("paciente_id", sort=False):
        caracteristicas, _ = e.extraer_paciente(tareas, extractores)
        filas.append({"paciente_id": pid, **caracteristicas})
        print(f"  {pid}: {len(caracteristicas)} columnas")
    return pd.DataFrame(filas).set_index("paciente_id")


def titulo(texto):
    print(f"\n{'=' * 78}\n{texto}\n{'=' * 78}")


def comparar(a, b):
    titulo("0. ESTRUCTURA")
    print(f"A: {a.shape[0]} x {a.shape[1]}   B: {b.shape[0]} x {b.shape[1]}   "
          f"mismas columnas: {list(a.columns) == list(b.columns)}   "
          f"mismos pacientes: {list(a.index) == list(b.index)}")

    titulo("1. BINS POR ROI (firstorder_Range / binWidth)")
    rango_a = a.filter(like="_firstorder_Range").to_numpy().ravel()
    rango_b = b.filter(like="_firstorder_Range").to_numpy().ravel()
    for nombre, rango, bw in (("A", rango_a, BIN_A), ("B", rango_b, BIN_B)):
        bins = rango / bw
        print(f"  {nombre} (binWidth {bw:>4}): mediana {np.median(bins):6.1f}   "
              f"min {bins.min():6.1f}   max {bins.max():6.1f}   ({bins.size} ROI)")

    titulo("2. SPEARMAN COLUMNA A COLUMNA (A vs B)")
    rho = pd.Series({c: a[c].corr(b[c], method="spearman") for c in a.columns})
    print(f"  n = {len(a)} pacientes por columna")
    print(f"  rho = 1:           {int((rho >= 1 - 1e-12).sum())}")
    print(f"  rho < {UMBRAL_RHO}:       {int((rho < UMBRAL_RHO).sum())}")
    print(f"  rho indefinido:    {int(rho.isna().sum())}  (columna constante en A o B)")
    clase = pd.Series([c.split("_")[2] for c in a.columns], index=a.columns)
    resumen = pd.DataFrame({"columnas": clase.value_counts(),
                            f"rho<{UMBRAL_RHO}": clase[rho < UMBRAL_RHO].value_counts()})
    print("\n  por clase:")
    print(resumen.fillna(0).astype(int).to_string().replace("\n", "\n    "))

    titulo(f"3. COLUMNAS CON rho < {UMBRAL_RHO}")
    bajas = rho[rho < UMBRAL_RHO].sort_values()
    with pd.option_context("display.max_rows", None):
        print(bajas.round(3).to_string())

    titulo("4. Energy, TotalEnergy y RMS (A | B)")
    for nombre in ("Energy", "TotalEnergy", "RootMeanSquared"):
        cols = [c for c in a.columns if c.endswith(f"_firstorder_{nombre}")]
        filas = []
        for c in cols:
            mod, region = c.split("_")[:2]
            for pid in a.index:
                filas.append({"modalidad": mod, "region": region, "paciente": pid,
                              "A": a.loc[pid, c], "B": b.loc[pid, c],
                              "B/A": b.loc[pid, c] / a.loc[pid, c]})
        rhos = ", ".join(f"{'_'.join(c.split('_')[:2])}={rho[c]:.2f}" for c in cols)
        print(f"\n  firstorder_{nombre}  (spearman A-B por columna: {rhos})")
        with pd.option_context("display.max_rows", None, "display.width", 140,
                               "display.float_format", "{:.4g}".format):
            print(pd.DataFrame(filas).to_string(index=False).replace("\n", "\n    "))


def main():
    radiomics.setVerbosity(logging.ERROR)
    print("Normalizando variante B (z-score x100) ...")
    rutas = normalizar_b()
    print(f"Extrayendo con {PARAMS_B} ...")
    b = extraer_b(rutas)
    b.to_csv(SALIDA / "caracteristicas_b.csv")

    a = pd.read_csv(TABLA_A).set_index("paciente_id")
    comparar(a, b[a.columns].loc[a.index])


if __name__ == "__main__":
    main()
