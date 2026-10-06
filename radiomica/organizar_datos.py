"""
BLOQUE 1 - ORGANIZACION DE LOS DATOS
Reto CIMAT / BraTS 2018

Recorre los datos originales de BraTS 2018 y genera las tablas de entrada del
pipeline. No mueve ni modifica ningun archivo de imagen.

Espera una estructura:
    datos/brats2018/
      HGG/Brats18_XXXX/Brats18_XXXX_<t1|t1ce|t2|flair|seg>.nii.gz
      LGG/Brats18_XXXX/...
      survival_data.csv        <- BraTS18ID, Age, Survival, ResectionStatus

Salida:
    datos/brats2018/manifest.csv   <- una fila por paciente completo: paciente_id, grado, rutas
    datos/brats2018/clinica.csv    <- una fila por paciente: grado, origen, edad, supervivencia
    logs/organizacion_<fecha>.csv  <- carpetas incompletas o con archivos de mas

survival_data.csv no se modifica: clinica.csv es su version validada y unida
con el grado. Solo los HGG tienen supervivencia en BraTS 2018.

Uso (desde la raiz del repo):
    python radiomica/organizar_datos.py
"""

from datetime import datetime
from pathlib import Path

import pandas as pd

import paciente as p

RAIZ = Path("datos/brats2018")
GRADOS = ["HGG", "LGG"]
SUPERVIVENCIA = RAIZ / "survival_data.csv"
MANIFEST = RAIZ / "manifest.csv"
CLINICA = RAIZ / "clinica.csv"
LOGS = Path("logs")


def escanear():
    """Devuelve (filas_manifest, problemas)."""
    filas, problemas = [], []
    for grado in GRADOS:
        for carpeta in sorted((RAIZ / grado).iterdir()):
            if not carpeta.is_dir():
                problemas.append({"ruta": carpeta.as_posix(), "problema": "no es una carpeta"})
                continue
            pid = carpeta.name
            rutas = {m: carpeta / f"{pid}_{m}.nii.gz" for m in p.MODALIDADES}
            faltan = [m for m, r in rutas.items() if not r.exists()]
            sobran = sorted({f.name for f in carpeta.iterdir()} - {r.name for r in rutas.values()})
            if faltan:
                problemas.append({"ruta": carpeta.as_posix(),
                                  "problema": f"faltan: {', '.join(faltan)}"})
                continue
            if sobran:
                # No impide usar al paciente, pero se deja constancia
                problemas.append({"ruta": carpeta.as_posix(),
                                  "problema": f"archivos extra: {', '.join(sobran)}"})
            filas.append({"paciente_id": pid, "grado": grado,
                          **{m: r.as_posix() for m, r in rutas.items()}})
    return filas, problemas


def leer_supervivencia():
    s = pd.read_csv(SUPERVIVENCIA, dtype=str, keep_default_na=False)
    esperadas = ["BraTS18ID", "Age", "Survival", "ResectionStatus"]
    if list(s.columns) != esperadas:
        raise ValueError(f"columnas inesperadas en {SUPERVIVENCIA}: {list(s.columns)}")
    if s["BraTS18ID"].duplicated().any():
        raise ValueError(f"ids repetidos: {s.loc[s['BraTS18ID'].duplicated(), 'BraTS18ID'].tolist()}")

    # Survival y Age deben ser numericos; un texto (p. ej. "ALIVE") debe fallar, no volverse NaN
    for col in ("Age", "Survival"):
        malos = s.loc[pd.to_numeric(s[col], errors="coerce").isna(), ["BraTS18ID", col]]
        if not malos.empty:
            raise ValueError(f"{col} no numerico:\n{malos.to_string(index=False)}")

    return pd.DataFrame({
        "paciente_id": s["BraTS18ID"],
        "edad": pd.to_numeric(s["Age"]),
        "supervivencia_dias": pd.to_numeric(s["Survival"]).astype(int),
        # "NA" en el original = estado de reseccion no reportado
        "reseccion": s["ResectionStatus"].replace({"NA": pd.NA}),
    })


def main():
    filas, problemas = escanear()
    manifest = pd.DataFrame(filas, columns=["paciente_id", "grado", *p.MODALIDADES])
    if manifest["paciente_id"].duplicated().any():
        raise ValueError("el mismo paciente aparece en HGG y LGG: "
                         f"{manifest.loc[manifest['paciente_id'].duplicated(), 'paciente_id'].tolist()}")
    manifest.to_csv(MANIFEST, index=False)

    superv = leer_supervivencia()
    sin_imagen = sorted(set(superv["paciente_id"]) - set(manifest["paciente_id"]))
    for pid in sin_imagen:
        problemas.append({"ruta": SUPERVIVENCIA.as_posix(),
                          "problema": f"{pid} tiene supervivencia pero no imagenes"})

    clinica = manifest[["paciente_id", "grado"]].copy()
    # Institucion de origen (CBICA, TCIA01, ..., 2013): util para agrupar en la validacion cruzada
    clinica["origen"] = clinica["paciente_id"].str.extract(r"^Brats18_([^_]+)_", expand=False)
    clinica = clinica.merge(superv, on="paciente_id", how="left")
    clinica.insert(3, "tiene_supervivencia", clinica["supervivencia_dias"].notna())
    clinica["supervivencia_dias"] = clinica["supervivencia_dias"].astype("Int64")
    clinica.to_csv(CLINICA, index=False)

    LOGS.mkdir(exist_ok=True)
    ruta_log = LOGS / f"organizacion_{datetime.now():%Y%m%d_%H%M%S}.csv"
    pd.DataFrame(problemas, columns=["ruta", "problema"]).to_csv(ruta_log, index=False)

    print(f"{'=' * 62}\nRESUMEN\n{'=' * 62}")
    print(f"pacientes completos: {len(manifest)}")
    print(manifest["grado"].value_counts().reindex(GRADOS).to_string())
    print(f"\ncon supervivencia: {int(clinica['tiene_supervivencia'].sum())}")
    print(pd.crosstab(clinica["grado"], clinica["tiene_supervivencia"]).to_string())
    print(f"\nreseccion (solo con supervivencia):")
    print(clinica.loc[clinica["tiene_supervivencia"], "reseccion"]
          .value_counts(dropna=False).to_string())
    print(f"\norigen:\n{pd.crosstab(clinica['origen'], clinica['grado']).to_string()}")
    print(f"\nproblemas registrados: {len(problemas)}")
    for pr in problemas:
        print(f"  {pr['ruta']}: {pr['problema']}")
    print(f"\nManifest: {MANIFEST}\nClinica:  {CLINICA}\nLog:      {ruta_log}")


if __name__ == "__main__":
    main()
