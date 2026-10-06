"""
BLOQUE 2 - VERIFICACION
Reto CIMAT / BraTS

Comprueba, para cada paciente, que los 5 archivos son lo que creemos:
  1. shape    -> mismo tamano en los 5
  2. affine   -> mismas matrices (co-registro)
  3. spacing  -> tamano del voxel en mm
  4. seg      -> etiquetas presentes
  5. conteo   -> volumen de cada region

Lee los pacientes de datos/brats2018/manifest.csv (generado por organizar_datos.py,
Bloque 1), asi que no depende de como esten organizadas las carpetas.

Uso (desde la raiz del repo):
    python radiomica/verificacion.py
"""

from pathlib import Path
import numpy as np
import nibabel as nib

import paciente as p
from regiones import leer_umbrales

MANIFEST = Path("datos/brats2018/manifest.csv")
MODALIDADES = ["t1", "t1ce", "t2", "flair"]
ARCHIVOS = MODALIDADES + ["seg"]

# Etiquetas segun la convencion de BraTS. El 3 no existe.
NOMBRES = {1: "necrosis / no realzado", 2: "edema", 4: "realce"}
ETIQUETAS_ESPERADAS = {0, 1, 2, 4}

# Por debajo de este numero de voxeles, forma y textura no son confiables.
# Se lee del params para usar el mismo umbral que PyRadiomics (minimumROISize).
MINIMO_VOXELES, _ = leer_umbrales()

# Tolerancia para comparar affines. Son flotantes: nunca usar ==
TOL = 1e-4


def verificar_paciente(pac):
    """Devuelve (lista_de_alertas, dict_con_volumenes)."""
    alertas = []
    rutas = pac.rutas()
    print(f"\n{'=' * 62}\nPACIENTE: {pac.paciente_id}\n{'=' * 62}")

    # --- existencia de archivos ------------------------------------
    faltantes = [s for s in ARCHIVOS if not rutas[s].exists()]
    if faltantes:
        alertas.append(f"faltan archivos: {', '.join(faltantes)}")
        print(f"  !! FALTAN: {faltantes}")
        return alertas, {}

    # --- cargar solo las cabeceras (no hace falta leer los datos) ---
    imgs = {s: nib.load(rutas[s]) for s in ARCHIVOS}

    # --- 1. shape ---------------------------------------------------
    print("\n[1] SHAPE")
    ref_shape = imgs["t1"].shape
    for s in ARCHIVOS:
        igual = imgs[s].shape == ref_shape
        print(f"    {s:<6} {str(imgs[s].shape):<20} {'ok' if igual else 'DISTINTO'}")
        if not igual:
            alertas.append(f"{s}: shape {imgs[s].shape} != {ref_shape}")

    # --- 2. affine (co-registro) ------------------------------------
    print("\n[2] AFFINE  (comparado contra t1, con tolerancia)")
    ref_aff = imgs["t1"].affine
    for s in ARCHIVOS:
        dif = float(np.abs(imgs[s].affine - ref_aff).max())
        igual = np.allclose(imgs[s].affine, ref_aff, atol=TOL)
        print(f"    {s:<6} dif.max = {dif:.2e}   {'ok' if igual else 'NO COINCIDE'}")
        if not igual:
            alertas.append(f"{s}: affine difiere de t1 en {dif:.2e}")

    # --- 3. spacing -------------------------------------------------
    print("\n[3] SPACING (mm por voxel)")
    for s in ARCHIVOS:
        z = tuple(round(float(v), 4) for v in imgs[s].header.get_zooms()[:3])
        print(f"    {s:<6} {z}")
        if s == "t1" and not np.allclose(z, (1.0, 1.0, 1.0), atol=1e-3):
            alertas.append(f"spacing no isotropico de 1 mm: {z}")
    orient = "".join(nib.aff2axcodes(ref_aff))
    vol_voxel = float(np.prod(imgs["t1"].header.get_zooms()[:3]))
    print(f"    orientacion: {orient}   |   1 voxel = {vol_voxel:.3f} mm3")

    # --- 4. etiquetas de seg ----------------------------------------
    print("\n[4] ETIQUETAS EN seg")
    seg = imgs["seg"].get_fdata().astype(np.int16)
    valores, conteos = np.unique(seg, return_counts=True)
    presentes = set(valores.tolist())
    print(f"    presentes: {sorted(presentes)}")
    if 3 in presentes:
        alertas.append("aparece la etiqueta 3: no sigue la convencion BraTS")
    inesperadas = presentes - ETIQUETAS_ESPERADAS
    if inesperadas:
        alertas.append(f"etiquetas inesperadas: {sorted(inesperadas)}")
    for e in (1, 2, 4):
        if e not in presentes:
            alertas.append(f"no hay voxeles de la etiqueta {e} ({NOMBRES[e]})")

    # --- 5. conteo y volumen ----------------------------------------
    print("\n[5] CONTEO POR REGION")
    vols = {}
    for val, n in zip(valores.tolist(), conteos.tolist()):
        if val == 0:
            continue
        mm3 = n * vol_voxel
        vols[val] = mm3
        marca = "  <-- MUY POCOS" if n < MINIMO_VOXELES else ""
        print(f"    etiqueta {val} ({NOMBRES.get(val, '?'):<22}) "
              f"{n:>8,} vox = {mm3/1000:7.2f} cm3{marca}")
        if n < MINIMO_VOXELES:
            alertas.append(f"etiqueta {val} con solo {n} voxeles")

    # regiones compuestas que usaremos en la extraccion
    wt = float((seg > 0).sum()) * vol_voxel
    tc = float(((seg == 1) | (seg == 4)).sum()) * vol_voxel
    et = float((seg == 4).sum()) * vol_voxel
    print(f"\n    WT (tumor completo) = {wt/1000:7.2f} cm3")
    print(f"    TC (nucleo)         = {tc/1000:7.2f} cm3")
    print(f"    ET (realce)         = {et/1000:7.2f} cm3")

    # coherencia: la segmentacion debe caer dentro del cerebro
    cerebro = imgs["t1"].get_fdata() > 0
    fuera = int(np.logical_and(seg > 0, ~cerebro).sum())
    print(f"    voxeles de tumor fuera del cerebro (t1): {fuera}")
    if fuera > 0:
        alertas.append(f"{fuera} voxeles de tumor caen fuera del cerebro")

    return alertas, {"WT": wt, "TC": tc, "ET": et}


def main():
    pacientes = p.cargar_pacientes(MANIFEST)
    print(f"Pacientes en {MANIFEST}: {len(pacientes)}")
    resumen = {}
    todas = {}

    for pac in pacientes:
        alertas, vols = verificar_paciente(pac)
        resumen[pac.paciente_id] = alertas
        if vols:
            todas[pac.paciente_id] = vols

    # --- tabla comparativa entre pacientes --------------------------
    if todas:
        print(f"\n{'=' * 62}\nVOLUMENES (cm3)\n{'=' * 62}")
        print(f"{'paciente':<24}{'WT':>10}{'TC':>10}{'ET':>10}")
        for nom, v in todas.items():
            print(f"{nom:<24}{v['WT']/1000:>10.2f}{v['TC']/1000:>10.2f}{v['ET']/1000:>10.2f}")

    # --- veredicto --------------------------------------------------
    print(f"\n{'=' * 62}\nRESUMEN\n{'=' * 62}")
    limpios = 0
    for nom, alertas in resumen.items():
        if not alertas:
            limpios += 1
            print(f"  [OK]     {nom}")
        else:
            print(f"  [REVISAR] {nom}")
            for a in alertas:
                print(f"            - {a}")
    print(f"\n{limpios}/{len(resumen)} pacientes sin alertas.")


if __name__ == "__main__":
    main()