"""
PARTICION EN FOLDS (compartida por los dos objetivos)
Reto CIMAT / BraTS 2018

Asigna cada paciente a uno de 5 folds. La MISMA particion se usa para:
  1. segmentacion (nnU-Net): el fold k se valida con un modelo entrenado en los
     otros 4, asi cada paciente recibe una mascara predicha "fuera de fold";
  2. pronostico: la validacion cruzada usa los mismos folds, para que ningun
     paciente de prueba haya sido visto por el segmentador que genero su mascara.

Estratifica por grupo: HGG con supervivencia / HGG sin supervivencia / LGG,
para que cada fold tenga la misma proporcion de pacientes utiles para el
pronostico (163 HGG con supervivencia).

Entrada:  data/clinica.csv      (organizar_datos.py)
Salida:   particiones/folds.csv (paciente_id, grado, grupo, fold)

Uso (desde la raiz del repo):
    python scripts/folds.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold

CLINICA = Path("data/clinica.csv")
SALIDA = Path("particiones/folds.csv")
N_FOLDS = 5
SEMILLA = 42


def main():
    if SALIDA.exists():
        # Cambiar los folds despues de entrenar invalida todo lo que se entreno con ellos
        raise SystemExit(f"{SALIDA} ya existe. Los folds no se regeneran: borralo a mano "
                         f"solo si estas seguro de que nada se entreno con ellos.")

    c = pd.read_csv(CLINICA)
    c["grupo"] = np.where(c["grado"] == "LGG", "LGG",
                          np.where(c["tiene_supervivencia"], "HGG_superv", "HGG_sin_superv"))

    skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEMILLA)
    c["fold"] = -1
    for k, (_, idx) in enumerate(skf.split(c, c["grupo"])):
        c.loc[idx, "fold"] = k

    SALIDA.parent.mkdir(exist_ok=True)
    c[["paciente_id", "grado", "grupo", "fold"]].to_csv(SALIDA, index=False)

    print(pd.crosstab(c["grupo"], c["fold"], margins=True).to_string())
    print(f"\nFolds: {SALIDA}")


if __name__ == "__main__":
    main()
