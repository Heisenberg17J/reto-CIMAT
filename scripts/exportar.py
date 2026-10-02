"""
BLOQUE 7 - TABLA FINAL
Reto CIMAT / BraTS

Toma la tabla validada por el Bloque 6 y la entrega para el modelado:
    resultados/features.csv   <- paciente_id como primera columna + caracteristicas
    resultados/README.md      <- fecha, version de PyRadiomics, YAML usado,
                                 normalizacion, regiones y pacientes procesados

Solo exporta si control_calidad.py se ejecuto sobre la extraccion actual
(caracteristicas_qc.csv mas reciente que caracteristicas.csv).

Aqui termina el pipeline de extraccion. La seleccion de caracteristicas, el
filtrado por correlacion y la reduccion de dimensionalidad NO van aqui: se
ajustan dentro de la validacion cruzada, con los datos de entrenamiento de cada
fold. Aplicarlos sobre features.csv y guardar el resultado contamina el
experimento.

Uso (desde la raiz del repo):
    python scripts/exportar.py
"""

import hashlib
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

import regiones as r

EXTRACCION = Path("resultados/caracteristicas.csv")
ENTRADA = Path("resultados/caracteristicas_qc.csv")
ELIMINADAS = Path("resultados/columnas_eliminadas.csv")
DIAGNOSTICOS = Path("resultados/diagnosticos.csv")
SALIDA = Path("resultados/features.csv")
README = Path("resultados/README.md")
LOGS = Path("logs")


def verificar_qc():
    """Exige que la tabla validada exista y corresponda a la extraccion actual."""
    if not ENTRADA.exists():
        sys.exit(f"No existe {ENTRADA}: ejecuta primero scripts/control_calidad.py")
    if ENTRADA.stat().st_mtime < EXTRACCION.stat().st_mtime:
        sys.exit(f"{ENTRADA} es anterior a {EXTRACCION}: la extraccion cambio despues "
                 f"del control de calidad. Ejecuta de nuevo scripts/control_calidad.py")


def fecha_extraccion():
    """Fecha de la ultima corrida de extraccion.py, tomada del nombre de su log."""
    logs = sorted(LOGS.glob("extraccion_*.txt"))
    if not logs:
        return "desconocida (no hay logs/extraccion_*.txt)"
    return f"{datetime.strptime(logs[-1].stem, 'extraccion_%Y%m%d_%H%M%S'):%Y-%m-%d %H:%M}"


def escribir_readme(features, eliminadas, n_manifest):
    diag = pd.read_csv(DIAGNOSTICOS)
    versiones = ", ".join(sorted(diag["diagnostics_Versions_PyRadiomics"].unique()))
    sha = hashlib.sha256(r.PARAMS.read_bytes()).hexdigest()[:12]
    regiones = "\n".join(f"| {g} | {' + '.join(map(str, e))} |" for g, e in r.REGIONES.items())
    if eliminadas.empty:
        txt_eliminadas = "ninguna."
    else:
        txt_eliminadas = f"{len(eliminadas)}, listadas con su motivo en `{ELIMINADAS.name}`."

    README.write_text(f"""# features.csv

Tabla de características radiómicas de BraTS 2018: una fila por paciente, y
`paciente_id` como primera columna.

| Campo | Valor |
|-------|-------|
| Fecha de generación | {datetime.now():%Y-%m-%d %H:%M} |
| Fecha de extracción | {fecha_extraccion()} |
| PyRadiomics | {versiones} (según `{DIAGNOSTICOS.name}`) |
| Parámetros | `{r.PARAMS.as_posix()}` (sha256 `{sha}…`) |
| Pacientes procesados | {len(features)} de {n_manifest} en el manifest |
| Columnas | {features.shape[1] - 1} características + `paciente_id` |

## Normalización

Z-score por paciente y por modalidad, calculado **solo sobre los voxeles del cerebro**
(intensidad > 0). El fondo queda en 0. Se hace antes de PyRadiomics
(`scripts/normalizar.py`), con `normalize: false` en el YAML.

## Regiones extraídas

| Región | Etiquetas de seg |
|--------|------------------|
{regiones}

Modalidades: {", ".join(r.MODALIDADES)}.

## Columnas

`<modalidad>_<región>_<clase>_<nombre>` (por ejemplo, `t1ce_ET_glcm_Contrast`). La forma
se extrae una sola vez por región: `mask_<región>_shape_<nombre>`.

Columnas constantes eliminadas en el control de calidad: {txt_eliminadas}

## Uso

Esta tabla **no** tiene selección de características, filtrado por correlación ni
reducción de dimensionalidad. Esos pasos van dentro de la validación cruzada y se ajustan
solo con los datos de entrenamiento de cada fold.

Ver `DECISIONES.md` en la raíz del repositorio.
""", encoding="utf-8")


def main():
    verificar_qc()
    df = pd.read_csv(ENTRADA)
    features = df[["paciente_id"] + [c for c in df.columns if c != "paciente_id"]]
    features.to_csv(SALIDA, index=False)

    n_manifest = len(pd.read_csv(r.MANIFEST_ENTRADA))
    escribir_readme(features, pd.read_csv(ELIMINADAS), n_manifest)

    print(f"features: {SALIDA}  ({features.shape[0]} pacientes x {features.shape[1] - 1} caracteristicas)")
    print(f"readme:   {README}")


if __name__ == "__main__":
    main()
