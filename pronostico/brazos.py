"""
OBJETIVO 2 - PRONOSTICO: los tres brazos de la comparacion
Reto CIMAT / BraTS 2018

Todos reciben SOLO el fold de entrenamiento para ajustar (escalado, agrupamiento,
seleccion e hiperparametros) y devuelven el riesgo del fold de prueba (D19).

    edad      Cox sobre las variables clinicas (edad). El suelo honesto.
    coxnet    Cox Elastic Net (l1_ratio 0.5), clinicas sin penalizar; alpha por
              c-index en CV interna (3 folds x 2 repeticiones).
    genetico  DEAP: selecciona <= 10 caracteristicas entre representantes de
              grupos de correlacion (|rho| > 0.9, sin usar la supervivencia).
              Aptitud = c-index medio en la CV interna de un Cox ridge
              (alpha = 1) sobre clinicas + subconjunto, menos 0.002 * k.
              El modelo final del fold es un Cox ridge con alpha por rejilla.

Cada brazo devuelve un dict con: riesgo_train, riesgo_test, dias_pred_test,
seleccion (caracteristicas radiomicas usadas) y extra (detalles del brazo).
"""

import random
import warnings

import numpy as np
import pandas as pd
from deap import algorithms, base, creator, tools
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.preprocessing import StandardScaler
from sksurv.linear_model import CoxnetSurvivalAnalysis, CoxPHSurvivalAnalysis
from sksurv.util import Surv

from metricas import cindex

# DEAP: las clases se crean a nivel de modulo para que joblib pueda serializar individuos
if not hasattr(creator, "FitnessMax"):
    creator.create("FitnessMax", base.Fitness, weights=(1.0,))
if not hasattr(creator, "Individual"):
    creator.create("Individual", list, fitness=creator.FitnessMax)

# Parametros del genetico (R5); se pueden reducir para pruebas rapidas
GA = dict(poblacion=40, generaciones=30, torneo=3, cxpb=0.7, elites=2,
          k_max=10, k_inicial=5, penalizacion=0.002, alpha_aptitud=1.0, umbral_grupo=0.9)
ALPHAS_RIDGE = [0.01, 0.1, 1.0, 10.0, 100.0]


# ---------------------------------------------------------------------------
# Utilidades comunes
# ---------------------------------------------------------------------------

def _surv(dias):
    return Surv.from_arrays(event=np.ones(len(dias), bool), time=np.asarray(dias, float))


def divisiones_internas(dias_train, semilla):
    """CV interna fija por fold externo: 3 folds x 2 repeticiones, estratificada por tercil."""
    tercil = pd.qcut(np.asarray(dias_train), 3, labels=False, duplicates="drop")
    rskf = RepeatedStratifiedKFold(n_splits=3, n_repeats=2, random_state=semilla)
    return list(rskf.split(np.zeros(len(tercil)), tercil))


def escalar(train, test):
    """StandardScaler ajustado solo con el entrenamiento."""
    sc = StandardScaler().fit(train)
    return (pd.DataFrame(sc.transform(train), index=train.index, columns=train.columns),
            pd.DataFrame(sc.transform(test), index=test.index, columns=test.columns))


def mediana_dias(funciones):
    """Mediana de supervivencia de cada funcion S(t); si nunca baja de 0.5, el ultimo tiempo."""
    return np.array([f.x[np.argmax(f.y <= 0.5)] if (f.y <= 0.5).any() else f.x[-1] for f in funciones])


def cox_ridge(Z, dias, alpha):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return CoxPHSurvivalAnalysis(alpha=alpha, ties="breslow").fit(Z, _surv(dias))


def c_interna(Z, dias, divisiones, alpha):
    """c-index medio de un Cox ridge en la CV interna."""
    Z, dias = np.asarray(Z), np.asarray(dias)
    vals = []
    for tr, va in divisiones:
        try:
            m = cox_ridge(Z[tr], dias[tr], alpha)
            vals.append(cindex(dias[va], m.predict(Z[va])))
        except (ArithmeticError, ValueError, np.linalg.LinAlgError):
            vals.append(0.5)
    return float(np.mean(vals))


def _final_ridge(Ztr, dias_tr, Zte, divisiones):
    """Cox ridge final con alpha elegido por rejilla en la CV interna."""
    alpha = max(ALPHAS_RIDGE, key=lambda a: c_interna(Ztr, dias_tr, divisiones, a))
    m = cox_ridge(Ztr.to_numpy(), dias_tr, alpha)
    return m, alpha


# ---------------------------------------------------------------------------
# Brazo 1: edad sola
# ---------------------------------------------------------------------------

def ajustar_edad(clin_tr, X_tr, dias_tr, clin_te, X_te, semilla):
    Ztr, Zte = escalar(clin_tr, clin_te)
    m = cox_ridge(Ztr.to_numpy(), dias_tr, 1e-4)
    return dict(riesgo_train=m.predict(Ztr.to_numpy()), riesgo_test=m.predict(Zte.to_numpy()),
                dias_pred_test=mediana_dias(m.predict_survival_function(Zte.to_numpy())),
                seleccion=[], extra={})


# ---------------------------------------------------------------------------
# Brazo 2: Cox Elastic Net
# ---------------------------------------------------------------------------

def ajustar_modelo_coxnet(clin, X, dias, semilla):
    """Ajusta el Cox Elastic Net completo (escalado + alpha por CV interna) con estos datos.

    Lo usan la validacion cruzada (ajustar_coxnet) y el modelo final (entrenar_final.py).
    Devuelve un dict con escalador, columnas, modelo, alpha, c_interna y seleccion.
    """
    Z = pd.concat([clin, X], axis=1)
    escalador = StandardScaler().fit(Z)
    A, dias = escalador.transform(Z), np.asarray(dias)
    n_clin = clin.shape[1]
    pf = np.r_[np.zeros(n_clin), np.ones(X.shape[1])]
    params = dict(l1_ratio=0.5, penalty_factor=pf, max_iter=100000)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        alphas = CoxnetSurvivalAnalysis(n_alphas=30, alpha_min_ratio=0.01, **params) \
            .fit(A, _surv(dias)).alphas_
        divisiones = divisiones_internas(dias, semilla)
        puntaje = np.zeros(len(alphas))
        for tr, va in divisiones:
            m = CoxnetSurvivalAnalysis(alphas=alphas, **params).fit(A[tr], _surv(dias[tr]))
            for k, a in enumerate(alphas):
                puntaje[k] += cindex(dias[va], m.predict(A[va], alpha=a))
        mejor = int(np.argmax(puntaje))
        # Se ajusta el camino hasta el alpha elegido (arranque en caliente, mas estable)
        modelo = CoxnetSurvivalAnalysis(alphas=alphas[:mejor + 1], fit_baseline_model=True, **params) \
            .fit(A, _surv(dias))
    coef = modelo.coef_[:, -1]
    return dict(escalador=escalador, columnas=list(Z.columns), modelo=modelo, alpha=float(alphas[mejor]),
                c_interna=float(puntaje[mejor] / len(divisiones)),
                seleccion=[c for c, w in zip(Z.columns[n_clin:], coef[n_clin:]) if w != 0])


def predecir_coxnet(ajuste, clin, X):
    """Riesgo y mediana de supervivencia (dias) con un ajuste de ajustar_modelo_coxnet."""
    A = ajuste["escalador"].transform(pd.concat([clin, X], axis=1)[ajuste["columnas"]])
    m, a = ajuste["modelo"], ajuste["alpha"]
    return m.predict(A, alpha=a), mediana_dias(m.predict_survival_function(A, alpha=a))


def ajustar_coxnet(clin_tr, X_tr, dias_tr, clin_te, X_te, semilla):
    aj = ajustar_modelo_coxnet(clin_tr, X_tr, dias_tr, semilla)
    riesgo_tr, _ = predecir_coxnet(aj, clin_tr, X_tr)
    riesgo_te, dias_te = predecir_coxnet(aj, clin_te, X_te)
    return dict(riesgo_train=riesgo_tr, riesgo_test=riesgo_te, dias_pred_test=dias_te,
                seleccion=aj["seleccion"], extra=dict(alpha=aj["alpha"], c_interna=aj["c_interna"]))


# ---------------------------------------------------------------------------
# Brazo 3: algoritmo genetico + Cox ridge
# ---------------------------------------------------------------------------

def representantes(Z, umbral):
    """Un representante (medoide) por grupo de |rho de Spearman| > umbral. No usa la supervivencia."""
    R = np.abs(np.corrcoef(Z.rank().to_numpy(), rowvar=False))
    D = 1 - R
    np.fill_diagonal(D, 0)
    grupos = fcluster(linkage(squareform(D, checks=False), "average"), t=1 - umbral, criterion="distance")
    elegidas = []
    for g in np.unique(grupos):
        idx = np.flatnonzero(grupos == g)
        elegidas.append(Z.columns[idx[np.argmax(R[np.ix_(idx, idx)].mean(axis=1))]])
    return elegidas


def ajustar_genetico(clin_tr, X_tr, dias_tr, clin_te, X_te, semilla, ga=None):
    ga = {**GA, **(ga or {})}
    random.seed(semilla)
    np.random.seed(semilla)

    Ctr, Cte = escalar(clin_tr, clin_te)
    Xtr, Xte = escalar(X_tr, X_te)
    reps = representantes(Xtr, ga["umbral_grupo"])
    n = len(reps)
    dias_tr = np.asarray(dias_tr)
    divisiones = divisiones_internas(dias_tr, semilla)
    C, R = Ctr.to_numpy(), Xtr[reps].to_numpy()

    cache = {}

    def aptitud(ind):
        clave = tuple(ind)
        if clave not in cache:
            on = [i for i, b in enumerate(ind) if b]
            c = c_interna(np.hstack([C, R[:, on]]), dias_tr, divisiones, ga["alpha_aptitud"])
            cache[clave] = (c - ga["penalizacion"] * len(on), c)
        return cache[clave]

    def nuevo():
        ind = [0] * n
        for i in random.sample(range(n), min(ga["k_inicial"], n)):
            ind[i] = 1
        return creator.Individual(ind)

    def reparar(ind):
        on = [i for i, b in enumerate(ind) if b]
        for i in random.sample(on, max(0, len(on) - ga["k_max"])):
            ind[i] = 0

    tb = base.Toolbox()
    tb.register("mate", tools.cxUniform, indpb=0.5)
    tb.register("mutate", tools.mutFlipBit, indpb=1.0 / n)
    tb.register("select", tools.selTournament, tournsize=ga["torneo"])

    pob = [nuevo() for _ in range(ga["poblacion"])]
    for ind in pob:
        ind.fitness.values = (aptitud(ind)[0],)
    curva = []
    for g in range(ga["generaciones"]):
        elites = [tb.clone(e) for e in tools.selBest(pob, ga["elites"])]
        hijos = algorithms.varAnd(tb.select(pob, ga["poblacion"] - ga["elites"]), tb,
                                  cxpb=ga["cxpb"], mutpb=1.0)
        for h in hijos:
            reparar(h)
            h.fitness.values = (aptitud(h)[0],)
        pob = elites + hijos
        apt = [i.fitness.values[0] for i in pob]
        mejor = tools.selBest(pob, 1)[0]
        curva.append(dict(generacion=g + 1, aptitud_max=max(apt), aptitud_media=float(np.mean(apt)),
                          k_mejor=int(sum(mejor))))

    mejor = tools.selBest(pob, 1)[0]
    seleccion = [reps[i] for i, b in enumerate(mejor) if b]
    apt_pen, c_bruta = aptitud(mejor)

    Ztr = pd.concat([Ctr, Xtr[seleccion]], axis=1)
    Zte = pd.concat([Cte, Xte[seleccion]], axis=1)
    m, alpha = _final_ridge(Ztr, dias_tr, Zte, divisiones)
    return dict(riesgo_train=m.predict(Ztr.to_numpy()), riesgo_test=m.predict(Zte.to_numpy()),
                dias_pred_test=mediana_dias(m.predict_survival_function(Zte.to_numpy())),
                seleccion=seleccion,
                extra=dict(aptitud_interna=apt_pen, c_interna=c_bruta, n_representantes=n,
                           alpha=alpha, evaluaciones=len(cache), curva=curva))


BRAZOS = {"edad": ajustar_edad, "coxnet": ajustar_coxnet, "genetico": ajustar_genetico}
