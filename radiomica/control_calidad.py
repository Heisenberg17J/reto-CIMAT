"""
BLOQUE 6 - CONTROL DE CALIDAD DE LA TABLA
Reto CIMAT / BraTS

Revisa resultados/radiomica/manual/caracteristicas.csv antes de modelar:
  1. cordura    -> volumenes conocidos (WT de Brats18_2013_2_1 ~ 61.000 mm3)
                   y VoxelVolume == voxeles contados en el Bloque 4
  2. columnas   -> el numero y los nombres cuadran con params + plan
  3. NaN / inf  -> por columna y por paciente
  4. constantes -> varianza cero en todos los pacientes: se eliminan
  5. rangos     -> |valor| >= 1e15 sugiere division por cero encubierta
  6. filas      -> una fila por paciente, id unico, todos los del manifest

Las pruebas 1, 2 y 6 son CRITICAS: si fallan no se escribe la tabla limpia y el
script termina con codigo 1. Las demas se reportan para revision humana.

Entrada:
    resultados/radiomica/manual/caracteristicas.csv
    datos/normalizada/manifest.csv, datos/normalizada/manifest_regiones.csv
    params_brats2018_v0.yaml

Salida:
    resultados/radiomica/manual/caracteristicas_qc.csv     <- tabla validada (sin constantes), entrada del Bloque 7
    resultados/radiomica/manual/columnas_eliminadas.csv    <- que columna se quito y por que
    logs/control_calidad_<fecha>.txt

Uso (desde la raiz del repo):
    python radiomica/control_calidad.py
"""

import logging
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from radiomics import getFeatureClasses

import regiones as r
import variante as v

ENTRADA = v.RESULTADOS / "caracteristicas.csv"
SALIDA = v.RESULTADOS / "caracteristicas_qc.csv"
ELIMINADAS = v.RESULTADOS / "columnas_eliminadas.csv"
LOGS = Path("logs")

# Prueba de cordura con un valor conocido (Bloque 6)
PACIENTE_REFERENCIA = "Brats18_2013_2_1"
COLUMNA_REFERENCIA = "mask_WT_shape_MeshVolume"
VALOR_REFERENCIA = 61_000.0   # mm3
TOL_REFERENCIA = 0.05         # 5 %

# MeshVolume (malla) siempre es algo menor que VoxelVolume; mas de esto es sospechoso
TOL_MALLA = 0.10

# |valor| a partir del cual sospechamos division por cero encubierta.
# Referencia: el maximo legitimo observado es ~3e8 (glszm_LargeAreaHighGrayLevelEmphasis).
LIMITE_ABSURDO = 1e15


def nombres_esperados(extractor):
    """Nombres (clase, caracteristica) que el extractor va a calcular.

    Una clase sin lista en el params calcula todas sus caracteristicas no obsoletas.
    """
    clases = getFeatureClasses()
    esperados = []
    for clase, lista in extractor.enabledFeatures.items():
        if not lista:
            lista = [n for n, obsoleta in clases[clase].getFeatureNames().items() if not obsoleta]
        esperados.extend((clase, n) for n in lista)
    return esperados


def columnas_esperadas(regiones_ok):
    """Columnas que deberian existir, con la convencion del Bloque 5 (D13)."""
    forma, intensidad = r.crear_extractores()
    cols = []
    for region in regiones_ok:
        cols += [f"mask_{region}_{c}_{n}" for c, n in nombres_esperados(forma)]
        for mod in r.MODALIDADES:
            cols += [f"{mod}_{region}_{c}_{n}" for c, n in nombres_esperados(intensidad)]
    return cols, len(nombres_esperados(forma)), len(nombres_esperados(intensidad))


class Control:
    """Acumula el resultado de cada prueba y escribe en el log."""

    def __init__(self, log):
        self.log = log
        self.criticos = []

    def ok(self, msg):
        self.log.info(f"  [OK]      {msg}")

    def aviso(self, msg):
        self.log.warning(f"  [REVISAR] {msg}")

    def fallo(self, msg):
        self.criticos.append(msg)
        self.log.error(f"  [FALLO]   {msg}")


def prueba_filas(df, manifest, c):
    c.log.info("\n[6] FILAS E IDENTIFICADOR")
    duplicados = df["paciente_id"][df["paciente_id"].duplicated()].unique().tolist()
    if df["paciente_id"].isna().any():
        c.fallo("hay filas sin paciente_id")
    if duplicados:
        c.fallo(f"paciente_id repetido: {duplicados}")
    else:
        c.ok(f"{len(df)} filas, {df['paciente_id'].nunique()} paciente_id unicos")

    esperados = set(manifest["paciente_id"])
    faltan = sorted(esperados - set(df["paciente_id"]))
    sobran = sorted(set(df["paciente_id"]) - esperados)
    if sobran:
        c.fallo(f"pacientes que no estan en el manifest: {sobran}")
    if faltan:
        # Puede ser un fallo legitimo de extraccion (D14): se reporta, no detiene
        c.aviso(f"{len(faltan)} pacientes del manifest sin fila (ver log de extraccion): {faltan}")
    else:
        c.ok(f"estan los {len(esperados)} pacientes del manifest")


def prueba_cordura(df, regiones, c):
    c.log.info("\n[1] CORDURA (volumenes)")
    # Los id repetidos se reportan en la prueba 6; aqui se usa la primera fila
    x = df.drop_duplicates("paciente_id").set_index("paciente_id")

    if PACIENTE_REFERENCIA not in x.index or COLUMNA_REFERENCIA not in x.columns:
        c.fallo(f"no se puede evaluar {COLUMNA_REFERENCIA} de {PACIENTE_REFERENCIA}")
    else:
        v = float(x.loc[PACIENTE_REFERENCIA, COLUMNA_REFERENCIA])
        dif = abs(v - VALOR_REFERENCIA) / VALOR_REFERENCIA
        msg = (f"{PACIENTE_REFERENCIA} {COLUMNA_REFERENCIA} = {v:,.0f} mm3 "
               f"(esperado ~{VALOR_REFERENCIA:,.0f}, dif {dif:.1%})")
        c.ok(msg) if dif <= TOL_REFERENCIA else c.fallo(msg + " -> error de mascara o spacing")

    # Generalizacion a todos: VoxelVolume debe igualar el conteo independiente del Bloque 4
    # (1 voxel = 1 mm3) y MeshVolume debe quedar cerca de VoxelVolume.
    malos_vox, malos_malla = [], []
    for _, reg in regiones[regiones["estado"] == "ok"].iterrows():
        pid, region = reg["paciente_id"], reg["region"]
        if pid not in x.index:
            continue
        vox = x.loc[pid, f"mask_{region}_shape_VoxelVolume"]
        malla = x.loc[pid, f"mask_{region}_shape_MeshVolume"]
        esperado = reg["volumen_cm3"] * 1000
        if not np.isclose(vox, esperado, rtol=1e-6):
            malos_vox.append(f"{pid}/{region}: {vox:,.0f} vs {esperado:,.0f}")
        if abs(malla - vox) / vox > TOL_MALLA:
            malos_malla.append(f"{pid}/{region}: malla {malla:,.0f} vs voxel {vox:,.0f}")

    if malos_vox:
        c.fallo(f"VoxelVolume no coincide con el conteo del Bloque 4: {malos_vox}")
    else:
        c.ok("VoxelVolume == voxeles contados en el Bloque 4, en todas las regiones")
    if malos_malla:
        c.aviso(f"MeshVolume se aleja mas de {TOL_MALLA:.0%} de VoxelVolume: {malos_malla}")
    else:
        c.ok(f"MeshVolume dentro de {TOL_MALLA:.0%} de VoxelVolume en todas las regiones")


def prueba_columnas(df, regiones, c):
    c.log.info("\n[2] NUMERO Y NOMBRES DE COLUMNAS")
    # Una region entra si es valida en al menos un paciente (si no, nunca genera columnas)
    regiones_ok = [g for g in r.REGIONES if (regiones.query("region == @g")["estado"] == "ok").any()]
    esperadas, n_forma, n_int = columnas_esperadas(regiones_ok)
    reales = [col for col in df.columns if col != "paciente_id"]

    c.log.info(f"  esperado: {len(regiones_ok)} regiones x {n_forma} forma "
               f"+ {len(r.MODALIDADES)} modalidades x {len(regiones_ok)} regiones x {n_int} "
               f"= {len(regiones_ok) * n_forma} + {len(r.MODALIDADES) * len(regiones_ok) * n_int} "
               f"= {len(esperadas)}")

    faltan = sorted(set(esperadas) - set(reales))
    sobran = sorted(set(reales) - set(esperadas))
    if faltan or sobran or len(reales) != len(esperadas):
        c.fallo(f"{len(reales)} columnas vs {len(esperadas)} esperadas; "
                f"faltan {faltan[:5]}{'...' if len(faltan) > 5 else ''} "
                f"sobran {sobran[:5]}{'...' if len(sobran) > 5 else ''}")
    else:
        c.ok(f"{len(reales)} columnas, coinciden una a una con las esperadas")

    mal_formadas = [col for col in reales if len(col.split("_")) != 4]
    if mal_formadas:
        c.fallo(f"nombres que no siguen <modalidad>_<region>_<clase>_<nombre>: {mal_formadas[:5]}")
    else:
        c.ok("todos los nombres tienen 4 partes")


def prueba_nan_inf(x, c):
    c.log.info("\n[3] NaN E INFINITOS")
    nan = x.isna().sum()
    inf = np.isinf(x).sum()
    for nombre, serie in (("NaN", nan), ("inf", inf)):
        malas = serie[serie > 0]
        if malas.empty:
            c.ok(f"sin {nombre}")
            continue
        c.aviso(f"{len(malas)} columnas con {nombre}")
        for col, n in malas.items():
            c.log.info(f"              {col}: {n} pacientes")
    filas_nan = x.index[x.isna().any(axis=1)].tolist()
    if filas_nan:
        c.log.info(f"  pacientes con NaN (regiones omitidas, D4/D11): {filas_nan}")


def prueba_constantes(x, c):
    """Devuelve la tabla de columnas a eliminar con su motivo."""
    c.log.info("\n[4] COLUMNAS CONSTANTES")
    finitos = x.replace([np.inf, -np.inf], np.nan)
    n_valores = finitos.nunique(dropna=True)
    eliminar = []
    for col in x.columns:
        if n_valores[col] == 0:
            eliminar.append((col, "sin valores finitos en ningun paciente"))
        elif n_valores[col] == 1:
            eliminar.append((col, f"varianza cero: vale {finitos[col].dropna().iloc[0]:.6g} "
                                  f"en todos los pacientes"))
    if eliminar:
        c.aviso(f"{len(eliminar)} columnas constantes, se eliminan:")
        for col, motivo in eliminar:
            c.log.info(f"              {col}: {motivo}")
    else:
        c.ok(f"ninguna columna constante en {len(x)} pacientes")
    if len(x) < 30:
        c.aviso(f"solo {len(x)} pacientes: repetir este control con el conjunto completo")
    return pd.DataFrame(eliminar, columns=["columna", "motivo"])


def prueba_rangos(x, c):
    c.log.info(f"\n[5] RANGOS ABSURDOS (|valor| >= {LIMITE_ABSURDO:.0e})")
    maximos = x.replace([np.inf, -np.inf], np.nan).abs().max().sort_values(ascending=False)
    absurdas = maximos[maximos >= LIMITE_ABSURDO]
    if absurdas.empty:
        c.ok("ninguna columna")
    else:
        c.aviso(f"{len(absurdas)} columnas:")
        for col, v in absurdas.items():
            c.log.info(f"              {col}: {v:.3e}")
    c.log.info("  5 mayores |valor| (para comparar):")
    for col, v in maximos.head(5).items():
        c.log.info(f"              {col}: {v:.3e}")


def configurar_log():
    LOGS.mkdir(exist_ok=True)
    ruta = LOGS / f"control_calidad_{v.PREFIJO_LOG}{datetime.now():%Y%m%d_%H%M%S}.txt"
    log = logging.getLogger("control_calidad")
    log.setLevel(logging.INFO)
    for h in (logging.FileHandler(ruta, encoding="utf-8"), logging.StreamHandler()):
        h.setFormatter(logging.Formatter("%(message)s"))
        log.addHandler(h)
    return log, ruta


def main():
    log, ruta_log = configurar_log()
    df = pd.read_csv(ENTRADA)
    manifest = pd.read_csv(r.MANIFEST_ENTRADA)
    regiones = pd.read_csv(r.MANIFEST_SALIDA)
    log.info(f"Control de calidad de {ENTRADA}  ({datetime.now():%Y-%m-%d %H:%M})")
    log.info(f"segmentacion: {v.DESCRIPCION}")
    log.info(f"{df.shape[0]} filas x {df.shape[1] - 1} columnas de caracteristicas")

    c = Control(log)
    prueba_cordura(df, regiones, c)
    prueba_columnas(df, regiones, c)

    x = df.set_index("paciente_id")
    prueba_nan_inf(x, c)
    eliminadas = prueba_constantes(x, c)
    prueba_rangos(x, c)
    prueba_filas(df, manifest, c)

    log.info(f"\n{'=' * 62}\nRESUMEN\n{'=' * 62}")
    if c.criticos:
        log.error(f"{len(c.criticos)} pruebas criticas fallaron; NO se escribe {SALIDA}")
        log.info(f"log: {ruta_log}")
        sys.exit(1)

    limpia = df.drop(columns=eliminadas["columna"])
    limpia.to_csv(SALIDA, index=False)
    eliminadas.to_csv(ELIMINADAS, index=False)
    log.info("pruebas criticas superadas")
    log.info(f"tabla limpia: {SALIDA}  ({limpia.shape[0]} x {limpia.shape[1] - 1})")
    log.info(f"eliminadas:   {ELIMINADAS}  ({len(eliminadas)} columnas)")
    log.info(f"log:          {ruta_log}")


if __name__ == "__main__":
    main()
