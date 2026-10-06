"""
OBJETIVO 2 - PRONOSTICO: modelo final para pacientes nuevos
Reto CIMAT / BraTS 2018

La validacion cruzada (evaluar.py, D27) estima como funcionara el modelo con
pacientes nuevos, pero descarta los modelos que entrena. Este script ajusta UNA
vez el brazo tradicional (Cox Elastic Net: edad sin penalizar + caracteristicas
de las mascaras predichas con el filtro D26) con los 163 HGG, igual que se
ajusto dentro de cada fold, y lo guarda para la inferencia (D30).

El desempeno esperado es el de la validacion cruzada, NO el c-index aparente en
estos mismos 163 pacientes, que es optimista y solo se imprime como referencia.

Salida (versionada):
    modelos/pronostico_coxnet.joblib   escalador, modelo, alpha, columnas y cortes de clase
    modelos/pronostico_coxnet.json     lo mismo legible: coeficientes, cortes, desempeno CV

Uso (desde la raiz del repo):
    python pronostico/entrenar_final.py
"""

import json
import os
import sys
from datetime import datetime
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import brazos as b    # noqa: E402
import datos as d     # noqa: E402
import metricas as m  # noqa: E402

SALIDA = Path("modelos")
CV = Path("resultados/pronostico/pred/por_fold.csv")
SEMILLA = 2026


def main():
    c = d.cargar_cohorte("pred", filtro=True)
    aj = b.ajustar_modelo_coxnet(c.clin, c.X, c.dias, SEMILLA)
    riesgo, _ = b.predecir_coxnet(aj, c.clin, c.X)

    # Cortes de clase sobre el riesgo de entrenamiento, igual que en la CV (metricas.clases_desde_riesgo)
    p = np.bincount(c.clase.to_numpy(), minlength=3) / len(c.clase)
    cortes = {"corta_si_riesgo_>=": float(np.quantile(riesgo, 1 - p[0])),
              "media_si_riesgo_>=": float(np.quantile(riesgo, 1 - p[0] - p[1]))}

    cv = pd.read_csv(CV).groupby("brazo")[["cindex", "exactitud"]].agg(["mean", "std"]).round(3)
    coef = aj["modelo"].coef_[:, -1]
    coeficientes = {col: float(w) for col, w in zip(aj["columnas"], coef) if w != 0}

    import sklearn
    import sksurv
    meta = {
        "descripcion": "Cox Elastic Net (l1_ratio 0.5), edad sin penalizar + radiomica de mascaras "
                       "predichas con filtro D26. Para HGG; no aplica a tumores sin realce (ET).",
        "fecha": datetime.now().isoformat(timespec="seconds"),
        "cohorte": c.descripcion,
        "alpha": aj["alpha"],
        "c_interna_cv": round(aj["c_interna"], 4),
        "n_caracteristicas_entrada": len(aj["columnas"]),
        "coeficientes_no_cero (escala estandarizada)": coeficientes,
        "cortes_de_clase": cortes,
        "clases": {"0": "corta (< 300 dias)", "1": "media (300-450 dias)", "2": "larga (> 450 dias)"},
        "desempeno_esperado_cv_5x10 (D27)": {
            brazo: {"cindex": f"{cv.loc[brazo, ('cindex', 'mean')]} +/- {cv.loc[brazo, ('cindex', 'std')]}",
                    "exactitud_3_clases": f"{cv.loc[brazo, ('exactitud', 'mean')]} +/- {cv.loc[brazo, ('exactitud', 'std')]}"}
            for brazo in ("edad", "coxnet")},
        "versiones": {"scikit-survival": sksurv.__version__, "scikit-learn": sklearn.__version__,
                      "numpy": np.__version__},
    }

    SALIDA.mkdir(exist_ok=True)
    joblib.dump({**aj, "cortes": cortes, "columnas_clinicas": list(c.clin.columns),
                 "columnas_radiomicas": list(c.X.columns), "meta": meta},
                SALIDA / "pronostico_coxnet.joblib", compress=3)
    (SALIDA / "pronostico_coxnet.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False))

    print(c.descripcion)
    print(f"alpha {aj['alpha']:.4f} | {len(coeficientes) - 1} caracteristicas radiomicas + edad")
    for col, w in sorted(coeficientes.items(), key=lambda kv: -abs(kv[1])):
        print(f"   {w:+.3f}  {col}")
    print(f"c-index APARENTE en los mismos 163 (optimista, solo referencia): {m.cindex(c.dias, riesgo):.3f}")
    print(f"desempeno esperado (CV 5x10, D27): c-index {cv.loc['coxnet', ('cindex', 'mean')]} "
          f"(edad sola {cv.loc['edad', ('cindex', 'mean')]})")
    print(f"\nmodelo: {SALIDA}/pronostico_coxnet.joblib (+ .json)")


if __name__ == "__main__":
    os.chdir(Path(__file__).resolve().parent.parent)
    main()
