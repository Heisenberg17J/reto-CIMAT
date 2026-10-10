"""
OBJETIVO 2 - PRONOSTICO: datos de la cohorte y particiones externas
Reto CIMAT / BraTS 2018

Cohorte: los 163 HGG con supervivencia (los unicos con el dato, D21). En BraTS
2018 todos tienen el evento (muerte): no hay censura.

    cargar_cohorte(variante, filtro, reseccion) -> Cohorte
        X       caracteristicas radiomicas (pred o manual; con o sin el filtro D26)
        clin    variables clinicas que siempre entran sin penalizar: edad
                (y, si reseccion=True, GTR/STR frente a "no reportada")
        dias    supervivencia en dias
        clase   0 corta (< 300), 1 media (300-450), 2 larga (> 450)

    particiones_externas(cohorte) -> DataFrame (paciente_id, repeticion, fold)
        5 folds x 10 repeticiones. La repeticion 0 son los folds de D22; las
        1-9 se estratifican por tercil de supervivencia. Se guardan en
        particiones/folds_pronostico.csv y no se regeneran (como D22).
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold

CLINICA = Path("datos/brats2018/clinica.csv")
FOLDS_D22 = Path("particiones/folds.csv")
FOLDS_PRONOSTICO = Path("particiones/folds_pronostico.csv")
ROBUSTEZ = Path("resultados/pronostico/robustez_segmentacion.csv")
FEATURES = {"pred": Path("resultados/radiomica/pred/features.csv"),
            "manual": Path("resultados/radiomica/manual/features.csv")}

CORTE_ROBUSTEZ = 0.85          # Spearman manual vs predicha (D26)
CORTES_DIAS = (300, 450)       # clases de BraTS 2018: < 10, 10-15, > 15 meses
N_FOLDS, N_REPETICIONES = 5, 10


@dataclass
class Cohorte:
    X: pd.DataFrame
    clin: pd.DataFrame
    dias: pd.Series
    clase: pd.Series
    descripcion: str
    X_manual: pd.DataFrame = None   # solo con filtro="fold": para calcular la robustez en cada fold


def clase_por_dias(dias):
    """0 corta (< 300), 1 media (300-450), 2 larga (> 450)."""
    d = np.asarray(dias)
    return np.where(d < CORTES_DIAS[0], 0, np.where(d <= CORTES_DIAS[1], 1, 2))


def cargar_cohorte(variante="pred", filtro=True, reseccion=False):
    """filtro: True = lista global de D26 (935); False = sin filtro (1146); "fold" = sin
    recortar X, pero con X_manual para que evaluar.py calcule el filtro en cada fold."""
    c = pd.read_csv(CLINICA).set_index("paciente_id")
    c = c[c["tiene_supervivencia"]].sort_index()

    X = pd.read_csv(FEATURES[variante]).set_index("paciente_id").loc[c.index]
    X_manual = None
    if filtro == "fold":
        if variante != "pred":
            raise ValueError("el filtro por fold solo tiene sentido con la variante pred")
        X_manual = pd.read_csv(FEATURES["manual"]).set_index("paciente_id").loc[c.index][X.columns]
        if X_manual.isna().any().any():
            raise ValueError("hay NaN en las caracteristicas manuales de la cohorte")
    elif filtro:
        rob = pd.read_csv(ROBUSTEZ)
        X = X[rob.loc[rob["spearman"] >= CORTE_ROBUSTEZ, "columna"].tolist()]
    if X.isna().any().any():
        raise ValueError("hay NaN en las caracteristicas de la cohorte (ver D25)")

    clin = c[["edad"]].astype(float)
    if reseccion:
        # "no reportada" es la referencia; ojo: identifica casi siempre al centro (D21)
        clin = clin.assign(reseccion_GTR=(c["reseccion"] == "GTR").astype(float),
                           reseccion_STR=(c["reseccion"] == "STR").astype(float))

    dias = c["supervivencia_dias"].astype(float)
    txt_filtro = "por fold" if filtro == "fold" else ("si" if filtro else "no")
    desc = (f"variante={variante} | filtro D26={txt_filtro} | "
            f"reseccion={'si' if reseccion else 'no'} | {len(c)} pacientes | "
            f"{X.shape[1]} caracteristicas + {clin.shape[1]} clinicas")
    return Cohorte(X, clin, dias, pd.Series(clase_por_dias(dias), index=c.index), desc, X_manual)


def particiones_externas(cohorte):
    """5 folds x 10 repeticiones; se crean una vez y despues solo se leen."""
    ids = cohorte.dias.index
    if FOLDS_PRONOSTICO.exists():
        p = pd.read_csv(FOLDS_PRONOSTICO)
        if set(p["paciente_id"]) != set(ids):
            raise SystemExit(f"{FOLDS_PRONOSTICO} no corresponde a la cohorte actual")
        return p

    filas = []
    d22 = pd.read_csv(FOLDS_D22).set_index("paciente_id")["fold"]
    filas += [{"paciente_id": pid, "repeticion": 0, "fold": int(d22[pid])} for pid in ids]

    tercil = pd.qcut(cohorte.dias, 3, labels=False)
    for rep in range(1, N_REPETICIONES):
        skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=rep)
        for k, (_, prueba) in enumerate(skf.split(ids, tercil)):
            filas += [{"paciente_id": ids[i], "repeticion": rep, "fold": k} for i in prueba]

    p = pd.DataFrame(filas).sort_values(["repeticion", "paciente_id"]).reset_index(drop=True)
    FOLDS_PRONOSTICO.parent.mkdir(exist_ok=True)
    p.to_csv(FOLDS_PRONOSTICO, index=False)
    return p
