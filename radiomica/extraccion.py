"""
BLOQUE 5 - EXTRACCION
Reto CIMAT / BraTS

Recorre el plan del Bloque 4 (paciente x region x modalidad) y extrae las
caracteristicas con PyRadiomics usando params_brats2018_v0.yaml.

Convencion de columnas (siempre 4 partes separadas por "_"):
    <modalidad>_<region>_<clase>_<nombre>     ej. t1ce_ET_glcm_Contrast
    mask_<region>_shape_<nombre>              ej. mask_WT_shape_Sphericity
La forma no depende de la modalidad (Bloque 4, D6), por eso su "modalidad" es mask.
Si se activan filtros (LoG/Wavelet), el filtro va dentro de <clase> con guiones:
    t1ce_ET_wavelet-LLH-glcm_Contrast

Entrada:
    datos/normalizada/manifest.csv
    datos/normalizada/manifest_regiones.csv    <- generado por regiones.py
    params_brats2018_v0.yaml

Salida:
    resultados/radiomica/manual/caracteristicas.csv   <- una fila por paciente (para el modelo)
    resultados/radiomica/manual/diagnosticos.csv      <- una fila por llamada a PyRadiomics
                                        (columnas diagnostics_*: versiones,
                                        parametros, hash y tamano de mascara)
    logs/extraccion_<fecha>.txt      <- que paciente fallo y por que

Un fallo en un paciente se registra y el lote continua; ese paciente no entra a
caracteristicas.csv. Una region con estado distinto de "ok" en el Bloque 4 se
omite y sus columnas quedan como NaN.

Uso (desde la raiz del repo):
    python radiomica/extraccion.py
    python radiomica/extraccion.py --pacientes Brats18_2013_2_1 Brats18_2013_3_1
"""

import argparse
import json
import logging
import time
from datetime import datetime
from pathlib import Path

import pandas as pd
import radiomics

import regiones as r
import variante as v

SALIDA = v.RESULTADOS          # resultados/radiomica/manual/ o .../pred/ (variante.py)
LOGS = Path("logs")

PREFIJO_DIAG = "diagnostics_"


def nombre_columna(clave, modalidad, region):
    """'original_glcm_Contrast' -> 't1ce_ET_glcm_Contrast'."""
    filtro, clase, nombre = clave.split("_", 2)
    if filtro != "original":
        clase = f"{filtro}-{clase}"
    return f"{modalidad or 'mask'}_{region}_{clase}_{nombre}"


def a_texto(valor):
    """Los diagnosticos traen dicts y tuplas: se guardan como JSON legible."""
    if isinstance(valor, (dict, list, tuple)):
        return json.dumps(valor, default=str)
    return valor


def extraer_paciente(tareas, extractores):
    """Ejecuta todas las tareas de un paciente.

    Devuelve (fila_de_caracteristicas, filas_de_diagnostico).
    Cualquier excepcion sube con la tarea que fallo en el mensaje.
    """
    caracteristicas, diagnosticos = {}, []
    for _, t in tareas.iterrows():
        modalidad = t["modalidad"] if pd.notna(t["modalidad"]) else None
        try:
            res = extractores[t["tipo"]].execute(t["imagen"], t["mascara"])
        except Exception as e:
            raise RuntimeError(f"{t['tipo']} {modalidad or 'mask'}/{t['region']}: "
                               f"{type(e).__name__}: {e}") from e

        diag = {"paciente_id": t["paciente_id"], "region": t["region"],
                "tipo": t["tipo"], "modalidad": modalidad or "mask"}
        for clave, valor in res.items():
            if clave.startswith(PREFIJO_DIAG):
                diag[clave] = a_texto(valor)
            else:
                caracteristicas[nombre_columna(clave, modalidad, t["region"])] = float(valor)
        diagnosticos.append(diag)

    return caracteristicas, diagnosticos


def configurar_log():
    LOGS.mkdir(exist_ok=True)
    ruta = LOGS / f"extraccion_{v.PREFIJO_LOG}{datetime.now():%Y%m%d_%H%M%S}.txt"
    log = logging.getLogger("extraccion")
    log.setLevel(logging.INFO)
    archivo = logging.FileHandler(ruta, encoding="utf-8")
    archivo.setFormatter(logging.Formatter("%(asctime)s  %(levelname)-7s %(message)s",
                                           "%Y-%m-%d %H:%M:%S"))
    consola = logging.StreamHandler()
    consola.setFormatter(logging.Formatter("%(message)s"))
    log.addHandler(archivo)
    log.addHandler(consola)
    # PyRadiomics escribe mucho por su cuenta: solo errores
    radiomics.setVerbosity(logging.ERROR)
    return log, ruta


def main():
    parser = argparse.ArgumentParser(description="Extraccion de caracteristicas radiomicas")
    parser.add_argument("--pacientes", nargs="+",
                        help="solo estos pacientes (por defecto, todos)")
    args = parser.parse_args()

    log, ruta_log = configurar_log()
    log.info(f"segmentacion: {v.DESCRIPCION}")
    log.info(f"pyradiomics {radiomics.__version__} | params {r.PARAMS}")

    forma, intensidad = r.crear_extractores()
    extractores = {"forma": forma, "intensidad": intensidad}

    pacientes = pd.read_csv(r.MANIFEST_ENTRADA)["paciente_id"].tolist()
    if args.pacientes:
        pacientes = [p for p in pacientes if p in set(args.pacientes)]
    plan = r.plan_extraccion()
    regiones = pd.read_csv(r.MANIFEST_SALIDA)
    log.info(f"pacientes: {len(pacientes)} | tareas en el plan: "
             f"{int(plan['paciente_id'].isin(pacientes).sum())}")

    filas, diagnosticos, fallidos = [], [], []
    for i, pid in enumerate(pacientes, 1):
        omitidas = regiones[(regiones["paciente_id"] == pid) & (regiones["estado"] != "ok")]
        for _, o in omitidas.iterrows():
            log.warning(f"{pid}: region {o['region']} omitida "
                        f"(estado={o['estado']}, {o['n_voxeles']} voxeles) -> NaN")

        inicio = time.time()
        try:
            tareas = plan[plan["paciente_id"] == pid]
            if tareas.empty:
                raise RuntimeError("sin tareas: ninguna region valida o falta en manifest_regiones.csv")
            caracteristicas, diag = extraer_paciente(tareas, extractores)
        except Exception as e:
            fallidos.append(pid)
            log.error(f"[{i}/{len(pacientes)}] {pid} FALLO: {e}")
            continue

        filas.append({"paciente_id": pid, **caracteristicas})
        diagnosticos.extend(diag)
        log.info(f"[{i}/{len(pacientes)}] {pid} ok: {len(tareas)} tareas, "
                 f"{len(caracteristicas)} caracteristicas, {time.time() - inicio:.1f} s")

    SALIDA.mkdir(parents=True, exist_ok=True)
    # pd.DataFrame alinea columnas entre pacientes: regiones omitidas quedan como NaN
    df = pd.DataFrame(filas)
    df.to_csv(SALIDA / "caracteristicas.csv", index=False)
    pd.DataFrame(diagnosticos).to_csv(SALIDA / "diagnosticos.csv", index=False)

    log.info("=" * 62)
    log.info(f"extraidos: {len(filas)}/{len(pacientes)} pacientes, "
             f"{df.shape[1] - 1 if len(df) else 0} columnas")
    if fallidos:
        log.info(f"fallidos:  {', '.join(fallidos)}")
    log.info(f"caracteristicas: {SALIDA / 'caracteristicas.csv'}")
    log.info(f"diagnosticos:    {SALIDA / 'diagnosticos.csv'}")
    log.info(f"log:             {ruta_log}")


if __name__ == "__main__":
    main()
