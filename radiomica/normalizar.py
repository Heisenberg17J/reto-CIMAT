"""
BLOQUE 3 - NORMALIZACION
Reto CIMAT / BraTS

Normaliza (z-score) cada modalidad de cada paciente usando solo los voxeles
del cerebro (intensidad > 0). El fondo se deja en 0.

Entrada  (no se modifica):
    datos/brats2018/manifest.csv
    datos/brats2018/<paciente>/<paciente>_<modalidad>.nii.gz

Salida:
    datos/normalizada/manifest.csv                      <- para el modelo
    datos/normalizada/<paciente>/<paciente>_<mod>.nii.gz
    logs/normalizacion_<fecha>.csv                     <- una fila por archivo

La segmentacion (seg) se copia tal cual: son etiquetas, no intensidades.

Uso (desde la raiz del repo):
    python radiomica/normalizar.py
    python radiomica/normalizar.py --sobrescribir
"""

import argparse
import shutil
from datetime import datetime
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

import paciente as p

MANIFEST_ENTRADA = Path("datos/brats2018/manifest.csv")
SALIDA = Path("datos/normalizada")
LOGS = Path("logs")
MODALIDADES = ["t1", "t1ce", "t2", "flair"]

# Tolerancia para verificar que el resultado quedo con media 0 y std 1
TOL = 1e-3


def zscore(data):
    """Devuelve (normalizado, estadisticas). Solo usa voxeles > 0."""
    mascara = data > 0
    valores = data[mascara]
    if valores.size == 0:
        raise ValueError("la imagen no tiene voxeles > 0")

    mu = float(valores.mean())
    sigma = float(valores.std())
    if sigma == 0:
        raise ValueError("sigma = 0 dentro del cerebro")

    normalizado = np.zeros_like(data, dtype=np.float32)
    normalizado[mascara] = (valores - mu) / sigma

    nuevos = normalizado[mascara]
    stats = {
        "n_voxeles": int(valores.size),
        "mu": mu,
        "sigma": sigma,
        "min_orig": float(valores.min()),
        "max_orig": float(valores.max()),
        "media_norm": float(nuevos.mean()),
        "std_norm": float(nuevos.std()),
        "min_norm": float(nuevos.min()),
        "max_norm": float(nuevos.max()),
    }
    return normalizado, stats


def procesar_paciente(pac, sobrescribir):
    """Normaliza un paciente. Devuelve (rutas_de_salida, filas_de_log)."""
    destino = SALIDA / pac.paciente_id
    destino.mkdir(parents=True, exist_ok=True)
    rutas, filas = {}, []

    for mod, origen in pac.rutas().items():
        salida = destino / origen.name
        rutas[mod] = salida.as_posix()
        fila = {"paciente_id": pac.paciente_id, "modalidad": mod,
                "origen": origen.as_posix(), "salida": salida.as_posix()}

        try:
            if salida.exists() and not sobrescribir:
                filas.append({**fila, "estado": "omitido", "mensaje": "ya existe"})
            elif mod == "seg":
                shutil.copy2(origen, salida)
                filas.append({**fila, "estado": "copiado", "mensaje": ""})
            else:
                filas.append(normalizar_archivo(origen, salida, fila))
        except Exception as e:
            filas.append({**fila, "estado": "error", "mensaje": str(e)})

        print(f"  {mod:<6} {filas[-1]['estado']}")

    return rutas, filas


def normalizar_archivo(origen, salida, fila):
    """Aplica z-score a un NIfTI, lo guarda en float32 y devuelve su fila de log."""
    img = nib.load(origen)
    normalizado, stats = zscore(img.get_fdata())

    nueva = nib.Nifti1Image(normalizado, img.affine, img.header)
    nueva.set_data_dtype(np.float32)
    nib.save(nueva, salida)

    ok = abs(stats["media_norm"]) < TOL and abs(stats["std_norm"] - 1) < TOL
    return {**fila, **stats,
            "estado": "ok" if ok else "revisar",
            "mensaje": "" if ok else "media/std fuera de tolerancia"}


def main():
    parser = argparse.ArgumentParser(description="Normalizacion z-score de MRI")
    parser.add_argument("--sobrescribir", action="store_true",
                        help="vuelve a generar archivos que ya existen")
    args = parser.parse_args()

    pacientes = p.cargar_pacientes(MANIFEST_ENTRADA)
    print(f"Pacientes en {MANIFEST_ENTRADA}: {len(pacientes)}")

    manifest, log = [], []
    for pac in pacientes:
        print(f"\n{pac.paciente_id}")
        rutas, filas = procesar_paciente(pac, args.sobrescribir)
        log.extend(filas)
        # Solo entran al manifest del modelo los pacientes sin errores
        if not any(f["estado"] == "error" for f in filas):
            manifest.append({"paciente_id": pac.paciente_id, **rutas})

    SALIDA.mkdir(exist_ok=True)
    pd.DataFrame(manifest, columns=["paciente_id", *p.MODALIDADES]).to_csv(
        SALIDA / "manifest.csv", index=False)

    LOGS.mkdir(exist_ok=True)
    ruta_log = LOGS / f"normalizacion_{datetime.now():%Y%m%d_%H%M%S}.csv"
    df_log = pd.DataFrame(log)
    df_log.insert(0, "fecha", datetime.now().isoformat(timespec="seconds"))
    df_log.to_csv(ruta_log, index=False)

    print(f"\n{'=' * 62}\nRESUMEN\n{'=' * 62}")
    print(df_log["estado"].value_counts().to_string())
    print(f"\nPacientes listos para el modelo: {len(manifest)}/{len(pacientes)}")
    print(f"Manifest: {SALIDA / 'manifest.csv'}")
    print(f"Log:      {ruta_log}")


if __name__ == "__main__":
    main()
