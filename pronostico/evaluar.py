"""
OBJETIVO 2 - PRONOSTICO: comparacion genetico vs Elastic Net vs edad (D27)
Reto CIMAT / BraTS 2018

Validacion cruzada anidada: 5 folds externos x 10 repeticiones (datos.py). En
cada fold externo los tres brazos (brazos.py) se ajustan solo con el
entrenamiento y se evaluan en la prueba; como comparten los folds, la
comparacion es pareada.

Salida en resultados/pronostico/<etiqueta>/:
    por_fold.csv          una fila por repeticion x fold x brazo: c-index, exactitud
                          3 clases, Spearman, MSE, k y detalles del brazo
    predicciones_oof.csv  riesgo, clase y dias predichos de cada paciente
    seleccion.csv         caracteristicas elegidas por fold y brazo
    curvas_genetico.csv   aptitud por generacion
    resumen.txt           tablas finales (tambien en el log)
    config.json           argumentos, versiones y fecha
Log: logs/pronostico_<etiqueta>_<fecha>.txt

Uso (desde la raiz del repo):
    python pronostico/evaluar.py                          # principal: pred + filtro D26
    python pronostico/evaluar.py --variante manual        # secundario: mascara manual
    python pronostico/evaluar.py --sin-filtro             # secundario: 1146 caracteristicas
    python pronostico/evaluar.py --reseccion              # secundario: + reseccion
    python pronostico/evaluar.py --permutar               # control: supervivencia barajada
    python pronostico/evaluar.py --repeticiones 1 --ga-poblacion 10 --ga-generaciones 3   # prueba rapida
"""

import os

# Un hilo de BLAS por proceso: el paralelismo lo da joblib
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import json
import logging
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

sys.path.insert(0, str(Path(__file__).resolve().parent))
import brazos as b    # noqa: E402
import datos as d     # noqa: E402
import metricas as m  # noqa: E402

LOGS = Path("logs")


def tarea(rep, fold, ids_tr, ids_te, cohorte, ga):
    """Ajusta y evalua los tres brazos en un fold externo."""
    assert not set(ids_tr) & set(ids_te), "entrenamiento y prueba se solapan"
    c = cohorte
    semilla = 1000 * rep + fold
    filas, oof, sel, curvas = [], [], [], []
    for nombre, ajustar in b.BRAZOS.items():
        t0 = time.time()
        kw = {"ga": ga} if nombre == "genetico" else {}
        r = ajustar(c.clin.loc[ids_tr], c.X.loc[ids_tr], c.dias.loc[ids_tr],
                    c.clin.loc[ids_te], c.X.loc[ids_te], semilla, **kw)
        dias_te, clase_te = c.dias.loc[ids_te].to_numpy(), c.clase.loc[ids_te].to_numpy()
        clase_pred = m.clases_desde_riesgo(r["riesgo_train"], c.clase.loc[ids_tr].to_numpy(), r["riesgo_test"])
        extra = {k: v for k, v in r["extra"].items() if k != "curva"}
        fila = dict(repeticion=rep, fold=fold, brazo=nombre, n_test=len(ids_te), semilla=semilla,
                    cindex=m.cindex(dias_te, r["riesgo_test"]),
                    exactitud=m.exactitud(clase_te, clase_pred),
                    spearman=m.spearman(dias_te, r["riesgo_test"]),
                    mse=m.mse(dias_te, r["dias_pred_test"]),
                    k=len(r["seleccion"]), segundos=round(time.time() - t0, 1), **extra)
        if "c_interna" in extra:
            fila["brecha"] = extra["c_interna"] - fila["cindex"]
        filas.append(fila)
        oof += [dict(repeticion=rep, fold=fold, brazo=nombre, paciente_id=pid, riesgo=rr,
                     dias=dd, dias_pred=dp, clase=cl, clase_pred=cp)
                for pid, rr, dd, dp, cl, cp in zip(ids_te, r["riesgo_test"], dias_te,
                                                   r["dias_pred_test"], clase_te, clase_pred)]
        sel += [dict(repeticion=rep, fold=fold, brazo=nombre, caracteristica=s) for s in r["seleccion"]]
        curvas += [dict(repeticion=rep, fold=fold, **g) for g in r["extra"].get("curva", [])]
    return filas, oof, sel, curvas


def t_corregido(dif, n_train, n_test):
    """t de Nadeau-Bengio para CV repetida (corrige el solapamiento de los entrenamientos)."""
    from scipy.stats import t as dist_t
    J = len(dif)
    var = np.var(dif, ddof=1) * (1 / J + n_test / n_train)
    if var == 0:
        return np.nan
    t = np.mean(dif) / np.sqrt(var)
    return float(2 * dist_t.sf(abs(t), J - 1))


def resumir(pf, sel, cohorte, n_rep, log):
    metricas = ["cindex", "exactitud", "spearman", "mse"]
    lineas = []

    tabla = pf.groupby("brazo")[metricas + ["k"]].agg(["mean", "std"])
    tabla.columns = [f"{a}_{b_}" for a, b_ in tabla.columns]
    tabla = tabla.reindex(list(b.BRAZOS))
    tabla["rmse_mean"] = np.sqrt(tabla.pop("mse_mean"))
    tabla = tabla.drop(columns="mse_std")
    lineas += ["METRICAS POR FOLD EXTERNO (media y DE sobre repeticiones x folds)",
               tabla.round(3).to_string(), ""]

    por_rep = pf.groupby(["brazo", "repeticion"])["cindex"].mean().unstack(0)[list(b.BRAZOS)]
    lineas += ["c-index medio por repeticion (DE entre repeticiones = variabilidad de la particion):",
               por_rep.agg(["mean", "std"]).round(3).to_string(), ""]

    n = len(cohorte.dias)
    lineas.append("DIFERENCIAS PAREADAS DE c-index (por fold externo)")
    piv = pf.pivot_table(index=["repeticion", "fold"], columns="brazo", values="cindex")
    for a, c in (("coxnet", "edad"), ("genetico", "edad"), ("genetico", "coxnet")):
        dif = (piv[a] - piv[c]).to_numpy()
        p = t_corregido(dif, n * 4 / 5, n / 5) if len(dif) > 1 else np.nan
        lineas.append(f"  {a:9s} - {c:7s}: {dif.mean():+.3f} (DE {dif.std(ddof=1) if len(dif) > 1 else 0:.3f}) | "
                      f"a favor en {np.mean(dif > 0):.0%} de {len(dif)} folds | p corregido = {p:.3f}")
    lineas.append("")

    ga = pf[pf.brazo == "genetico"]
    if len(ga):
        lineas += ["GENETICO: optimismo de su aptitud",
                   f"  c interna {ga.c_interna.mean():.3f} vs c externo {ga.cindex.mean():.3f} "
                   f"-> brecha {ga.brecha.mean():+.3f} (DE {ga.brecha.std():.3f}) | "
                   f"k medio {ga.k.mean():.1f} | representantes {ga.n_representantes.mean():.0f}", ""]

    lineas.append("ESTABILIDAD DE LA SELECCION entre folds externos")
    for brazo in ("coxnet", "genetico"):
        conj = [set(g.caracteristica) for _, g in sel[sel.brazo == brazo].groupby(["repeticion", "fold"])]
        conj += [set()] * (len(piv) - len(conj))           # folds sin ninguna caracteristica
        lineas.append(f"  {brazo:9s}: Jaccard medio {m.jaccard_medio(conj):.3f} | "
                      f"Nogueira {m.nogueira(conj, cohorte.X.columns):.3f}")
        top = sel[sel.brazo == brazo].caracteristica.value_counts().head(5)
        if len(top):
            lineas.append("     mas elegidas: " + ", ".join(f"{k} ({v}/{len(piv)})" for k, v in top.items()))
    for l_ in lineas:
        log.info(l_)
    return "\n".join(lineas)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variante", choices=["pred", "manual"], default="pred")
    ap.add_argument("--sin-filtro", action="store_true")
    ap.add_argument("--reseccion", action="store_true")
    ap.add_argument("--permutar", action="store_true", help="baraja la supervivencia (control de fugas)")
    ap.add_argument("--repeticiones", type=int, default=d.N_REPETICIONES)
    ap.add_argument("--nucleos", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--ga-poblacion", type=int, default=b.GA["poblacion"])
    ap.add_argument("--ga-generaciones", type=int, default=b.GA["generaciones"])
    ap.add_argument("--etiqueta", help="nombre de la carpeta de salida (por defecto, segun las opciones)")
    args = ap.parse_args()

    cohorte = d.cargar_cohorte(args.variante, filtro=not args.sin_filtro, reseccion=args.reseccion)
    part = d.particiones_externas(cohorte)
    if args.permutar:
        rng = np.random.default_rng(12345)
        cohorte.dias[:] = rng.permutation(cohorte.dias.to_numpy())
        cohorte.clase[:] = d.clase_por_dias(cohorte.dias)

    etiqueta = args.etiqueta or "_".join([args.variante] + (["sinfiltro"] if args.sin_filtro else [])
                                         + (["reseccion"] if args.reseccion else [])
                                         + (["permutado"] if args.permutar else []))
    salida = Path("resultados/pronostico") / etiqueta
    salida.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(exist_ok=True)
    log = logging.getLogger("pronostico")
    log.setLevel(logging.INFO)
    for h in (logging.FileHandler(LOGS / f"pronostico_{etiqueta}_{datetime.now():%Y%m%d_%H%M%S}.txt"),
              logging.StreamHandler()):
        h.setFormatter(logging.Formatter("%(message)s"))
        log.addHandler(h)

    ga = dict(poblacion=args.ga_poblacion, generaciones=args.ga_generaciones)
    log.info(f"{cohorte.descripcion}{' | SUPERVIVENCIA PERMUTADA' if args.permutar else ''}")
    log.info(f"{args.repeticiones} repeticiones x {d.N_FOLDS} folds | genetico {ga} | {args.nucleos} nucleos")

    tareas = []
    for rep in range(args.repeticiones):
        f = part[part.repeticion == rep].set_index("paciente_id")["fold"]
        for k in range(d.N_FOLDS):
            tareas.append((rep, k, f.index[f != k].tolist(), f.index[f == k].tolist()))

    t0 = time.time()
    res = Parallel(n_jobs=args.nucleos, verbose=5)(
        delayed(tarea)(rep, k, tr, te, cohorte, ga) for rep, k, tr, te in tareas)
    log.info(f"tiempo total: {(time.time() - t0) / 60:.1f} min")

    pf = pd.DataFrame([f for r in res for f in r[0]])
    oof = pd.DataFrame([f for r in res for f in r[1]])
    sel = pd.DataFrame([f for r in res for f in r[2]], columns=["repeticion", "fold", "brazo", "caracteristica"])
    curvas = pd.DataFrame([f for r in res for f in r[3]])
    pf.to_csv(salida / "por_fold.csv", index=False)
    oof.to_csv(salida / "predicciones_oof.csv", index=False)
    sel.to_csv(salida / "seleccion.csv", index=False)
    curvas.to_csv(salida / "curvas_genetico.csv", index=False)

    (salida / "resumen.txt").write_text(resumir(pf, sel, cohorte, args.repeticiones, log) + "\n")
    import deap
    import sklearn
    import sksurv
    (salida / "config.json").write_text(json.dumps({
        **vars(args), "descripcion": cohorte.descripcion, "ga": {**b.GA, **ga},
        "fecha": datetime.now().isoformat(timespec="seconds"),
        "versiones": {"sksurv": sksurv.__version__, "deap": deap.__version__,
                      "sklearn": sklearn.__version__, "numpy": np.__version__}}, indent=2))
    log.info(f"\nresultados: {salida}/")


if __name__ == "__main__":
    os.chdir(Path(__file__).resolve().parent.parent)
    main()
