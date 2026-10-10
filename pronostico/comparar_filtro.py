"""
OBJETIVO 2 - PRONOSTICO: filtro D26 global frente a filtro D26 dentro de cada fold
Reto CIMAT / BraTS 2018

Compara dos corridas de evaluar.py que solo difieren en donde se calcula el
filtro de robustez (Spearman manual vs predicha >= 0.85):
    pred              lista global de 935, calculada con los 163 pacientes
    pred_filtrofold   lista recalculada en cada fold con su entrenamiento externo
Mismas particiones, semillas, brazos e hiperparametros: la diferencia aisla el
efecto de calcular el filtro con los pacientes de prueba.

Salida: resultados/pronostico/pred_filtrofold/comparacion_filtro.txt

Uso (desde la raiz del repo, despues de las dos corridas):
    python pronostico/comparar_filtro.py
"""

import os
import re
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import datos as d     # noqa: E402
import metricas as m  # noqa: E402
from evaluar import t_corregido  # noqa: E402

BASE = Path("resultados/pronostico")
GLOBAL, FOLD = BASE / "pred", BASE / "pred_filtrofold"
BRAZOS = ["edad", "coxnet", "genetico"]


def tiempo_total(etiqueta):
    """Minutos de la ultima corrida registrada en logs/pronostico_<etiqueta>_*.txt."""
    logs = sorted(Path("logs").glob(f"pronostico_{etiqueta}_2*.txt"))
    for log in reversed(logs):
        hit = re.search(r"tiempo total: ([\d.]+) min", log.read_text())
        if hit:
            return float(hit.group(1))
    return np.nan


def main():
    g = pd.read_csv(GLOBAL / "por_fold.csv")
    f = pd.read_csv(FOLD / "por_fold.csv")
    filtro = pd.read_csv(FOLD / "filtro_por_fold.csv")
    rob = pd.read_csv(d.ROBUSTEZ)
    lista_global = set(rob.loc[rob["spearman"] >= d.CORTE_ROBUSTEZ, "columna"])
    universo = rob["columna"].tolist()
    n = 163
    clave = ["repeticion", "fold", "brazo"]
    if not g[clave].sort_values(clave).reset_index(drop=True).equals(
            f[clave].sort_values(clave).reset_index(drop=True)):
        raise SystemExit("las dos corridas no tienen los mismos folds")
    if not (g.set_index(clave).semilla == f.set_index(clave).semilla.reindex(g.set_index(clave).index)).all():
        raise SystemExit("las semillas no coinciden")

    L = []
    L += ["FILTRO D26 GLOBAL (935, 163 pacientes) FRENTE A FILTRO DENTRO DE CADA FOLD EXTERNO",
          f"{len(g) // 3} folds externos (5 x 10), mismas particiones y semillas", ""]

    # 1. c-index externo de cada brazo
    L.append("1. C-INDEX EXTERNO (media y DE sobre 50 folds) Y DIFERENCIA PAREADA fold - global")
    for brazo in BRAZOS:
        a = g[g.brazo == brazo].set_index(["repeticion", "fold"]).cindex
        b = f[f.brazo == brazo].set_index(["repeticion", "fold"]).cindex.reindex(a.index)
        dif = (b - a).to_numpy()
        p = t_corregido(dif, n * 4 / 5, n / 5)
        L.append(f"  {brazo:9s} global {a.mean():.3f} ({a.std():.3f}) | por fold {b.mean():.3f} ({b.std():.3f}) | "
                 f"dif {dif.mean():+.3f} (DE {dif.std(ddof=1):.3f}) | folds identicos {np.mean(dif == 0):.0%} | "
                 f"p corregido = {p:.3f}")
    for met in ("exactitud", "spearman"):
        txt = " | ".join(f"{br} {g[g.brazo == br][met].mean():.3f} -> {f[f.brazo == br][met].mean():.3f}"
                         for br in BRAZOS)
        L.append(f"  {met:9s}: {txt}")
    rmse = " | ".join(f"{br} {np.sqrt(g[g.brazo == br].mse.mean()):.0f} -> {np.sqrt(f[f.brazo == br].mse.mean()):.0f}"
                      for br in BRAZOS)
    L += [f"  rmse dias: {rmse}", ""]

    # 2. Optimismo del genetico
    L.append("2. OPTIMISMO DEL GENETICO (c interna - c externo)")
    for nombre, t in (("global", g), ("por fold", f)):
        ga = t[t.brazo == "genetico"]
        L.append(f"  {nombre:9s}: c interna {ga.c_interna.mean():.3f} | c externo {ga.cindex.mean():.3f} | "
                 f"brecha {ga.brecha.mean():+.3f} (DE {ga.brecha.std():.3f}) | k {ga.k.mean():.1f} | "
                 f"representantes {ga.n_representantes.mean():.0f}")
    L.append("")

    # 3. Diferencias entre brazos en cada corrida
    L.append("3. DIFERENCIAS ENTRE BRAZOS (c-index pareado, t de Nadeau-Bengio)")
    for nombre, t in (("global", g), ("por fold", f)):
        piv = t.pivot_table(index=["repeticion", "fold"], columns="brazo", values="cindex")
        partes = []
        for a, c in (("coxnet", "edad"), ("genetico", "edad"), ("genetico", "coxnet")):
            dif = (piv[a] - piv[c]).to_numpy()
            partes.append(f"{a}-{c} {dif.mean():+.3f} (p={t_corregido(dif, n * 4 / 5, n / 5):.2f})")
        L.append(f"  {nombre:9s}: " + " | ".join(partes))
    L.append("")

    # 4. Estabilidad del filtro
    conj = [set(x.caracteristica) for _, x in filtro.groupby(["repeticion", "fold"])]
    nf = f[f.brazo == "edad"].n_filtro
    cuenta = filtro.caracteristica.value_counts().reindex(universo, fill_value=0)
    M = len(conj)
    L += ["4. ESTABILIDAD DEL FILTRO ENTRE FOLDS",
          f"  caracteristicas retenidas: media {nf.mean():.1f} | DE {nf.std():.1f} | "
          f"min {nf.min()} | max {nf.max()} (global: {len(lista_global)})",
          f"  Jaccard medio entre folds {m.jaccard_medio(conj):.3f} | Nogueira {m.nogueira(conj, universo):.3f}",
          f"  Jaccard medio con la lista global {np.mean([len(s & lista_global) / len(s | lista_global) for s in conj]):.3f}",
          f"  retenidas en los {M} folds: {int((cuenta == M).sum())} | en ninguno: {int((cuenta == 0).sum())} | "
          f"inestables (1 a {M - 1}): {int(((cuenta > 0) & (cuenta < M)).sum())}",
          f"  de la lista global, retenidas en todos los folds: {int((cuenta[list(lista_global)] == M).sum())}; "
          f"fuera de la global pero retenidas en algun fold: {int((cuenta.drop(list(lista_global)) > 0).sum())}"]
    borde = rob[rob.columna.map(cuenta).between(1, M - 1)]
    if len(borde):
        L.append(f"  Spearman global de las inestables: {borde.spearman.min():.3f} a {borde.spearman.max():.3f} "
                 f"(mediana {borde.spearman.median():.3f})")
    L.append("")

    # Estabilidad de la seleccion de los brazos, con el mismo universo (1146) en las dos corridas
    L.append("   Estabilidad de la seleccion de los brazos (Nogueira, universo de 1146 en ambas)")
    for brazo in ("coxnet", "genetico"):
        vals = []
        for t, carpeta in ((g, GLOBAL), (f, FOLD)):
            sel = pd.read_csv(carpeta / "seleccion.csv")
            sel = sel[sel.brazo == brazo]
            cs = [set(x.caracteristica) for _, x in sel.groupby(["repeticion", "fold"])]
            cs += [set()] * (len(t) // 3 - len(cs))
            vals.append((m.nogueira(cs, universo), m.jaccard_medio(cs)))
        L.append(f"   {brazo:9s}: global Nogueira {vals[0][0]:.3f} / Jaccard {vals[0][1]:.3f} | "
                 f"por fold Nogueira {vals[1][0]:.3f} / Jaccard {vals[1][1]:.3f}")
    L.append("")

    # 5. Tiempo
    tg, tf = tiempo_total("pred"), tiempo_total("pred_filtrofold")
    L += ["5. TIEMPO DE EJECUCION",
          f"  corrida completa (10 nucleos): global {tg:.1f} min | por fold {tf:.1f} min",
          f"  filtro por fold: {f[f.brazo == 'edad'].segundos_filtro.mean():.2f} s por fold "
          f"({f[f.brazo == 'edad'].segundos_filtro.sum():.0f} s en total)"]
    for brazo in BRAZOS:
        L.append(f"  {brazo:9s}: {g[g.brazo == brazo].segundos.mean():.1f} s por fold global | "
                 f"{f[f.brazo == brazo].segundos.mean():.1f} s por fold con filtro por fold")

    texto = "\n".join(L)
    print(texto)
    (FOLD / "comparacion_filtro.txt").write_text(texto + "\n")


if __name__ == "__main__":
    os.chdir(Path(__file__).resolve().parent.parent)
    main()
