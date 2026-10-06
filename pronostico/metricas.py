"""
OBJETIVO 2 - PRONOSTICO: metricas
Reto CIMAT / BraTS 2018

Principal: c-index de Harrell. Sin censura (BraTS 2018), equivale a la
concordancia entre el riesgo predicho y los dias: un riesgo mayor debe ir con
menos dias.

Secundarias (comparables con BraTS 2018): exactitud en 3 clases, rho de
Spearman y MSE en dias. Ademas, la estabilidad de la seleccion entre folds.
"""

import itertools

import numpy as np
from scipy.stats import spearmanr
from sksurv.metrics import concordance_index_censored


def cindex(dias, riesgo):
    dias = np.asarray(dias, float)
    return concordance_index_censored(np.ones(len(dias), bool), dias, np.asarray(riesgo, float))[0]


def clases_desde_riesgo(riesgo_train, clase_train, riesgo_test):
    """Clase predicha a partir del riesgo, sin entrenar un clasificador.

    Los cortes son cuantiles del riesgo en el entrenamiento, elegidos para que
    alli se reproduzcan las proporciones de clase del mismo entrenamiento.
    Mas riesgo = supervivencia mas corta (clase 0).
    """
    p = np.bincount(np.asarray(clase_train), minlength=3) / len(clase_train)
    corte_corta = np.quantile(riesgo_train, 1 - p[0])
    corte_media = np.quantile(riesgo_train, 1 - p[0] - p[1])
    r = np.asarray(riesgo_test)
    return np.where(r >= corte_corta, 0, np.where(r >= corte_media, 1, 2))


def exactitud(clase, clase_pred):
    return float(np.mean(np.asarray(clase) == np.asarray(clase_pred)))


def spearman(dias, riesgo):
    """rho entre -riesgo y dias: positivo = el modelo ordena bien."""
    return float(spearmanr(-np.asarray(riesgo), dias)[0])


def mse(dias, dias_pred):
    return float(np.mean((np.asarray(dias) - np.asarray(dias_pred)) ** 2))


def jaccard_medio(conjuntos):
    """Jaccard medio entre todos los pares de subconjuntos seleccionados."""
    pares = [len(a & b) / len(a | b) if (a | b) else 1.0
             for a, b in itertools.combinations([set(c) for c in conjuntos], 2)]
    return float(np.mean(pares)) if pares else np.nan


def nogueira(conjuntos, universo):
    """Indice de estabilidad de Nogueira et al. (2018): 1 = identicos, ~0 = al azar."""
    universo = list(universo)
    M, p = len(conjuntos), len(universo)
    Z = np.array([[f in set(c) for f in universo] for c in conjuntos], float)
    k_medio = Z.sum(axis=1).mean()
    if M < 2 or k_medio in (0, p):
        return np.nan
    s2 = M / (M - 1) * Z.mean(axis=0) * (1 - Z.mean(axis=0))
    return float(1 - s2.mean() / ((k_medio / p) * (1 - k_medio / p)))
