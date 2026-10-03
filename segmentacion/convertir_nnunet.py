"""
OBJETIVO 1 - SEGMENTACION: conversion al formato de nnU-Net v2
Reto CIMAT / BraTS 2018

Genera el dataset crudo que espera nnU-Net, igual que su conversor oficial de
BraTS (Dataset137_BraTS2021):

    nnunet_raw/Dataset501_BraTS2018/
      imagesTr/<paciente>_0000.nii.gz   t1     (enlace duro al original: 0 bytes extra)
               <paciente>_0001.nii.gz   t1ce
               <paciente>_0002.nii.gz   t2
               <paciente>_0003.nii.gz   flair
      labelsTr/<paciente>.nii.gz        etiquetas reasignadas (nnU-Net exige 0,1,2,3)
      dataset.json                      entrenamiento por regiones WT / TC / ET
      splits_final.json                 los 5 folds de particiones/folds.csv

Etiquetas: BraTS 0/1/2/4 -> nnU-Net 0/2/1/3
    edema (2) -> 1,   necrosis (1) -> 2,   realce (4) -> 3
Asi las regiones son anidadas y consecutivas:
    WT = {1,2,3}   TC = {2,3}   ET = {3}

Se usan las imagenes ORIGINALES, no data_normalizada: nnU-Net normaliza por su
cuenta (z-score sobre los voxeles no nulos, el mismo criterio que el Bloque 3).

Uso (desde la raiz del repo):
    python segmentacion/convertir_nnunet.py
    python segmentacion/convertir_nnunet.py --zip   # ademas, zip para subir a Drive
"""

import argparse
import json
import os
import shutil
import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import paciente as p  # noqa: E402

MANIFEST = Path("data/manifest.csv")
FOLDS = Path("particiones/folds.csv")
DATASET = Path("nnunet_raw/Dataset501_BraTS2018")

CANALES = {"t1": "0000", "t1ce": "0001", "t2": "0002", "flair": "0003"}
BRATS_A_NNUNET = {0: 0, 2: 1, 1: 2, 4: 3}
NNUNET_A_BRATS = {v: k for k, v in BRATS_A_NNUNET.items()}


def reasignar(seg, mapa):
    """Cambia etiquetas segun mapa; falla con cualquier etiqueta no prevista."""
    inesperadas = set(np.unique(seg).tolist()) - set(mapa)
    if inesperadas:
        raise ValueError(f"etiquetas inesperadas: {sorted(inesperadas)}")
    nueva = np.zeros_like(seg, dtype=np.uint8)
    for origen, destino in mapa.items():
        nueva[seg == origen] = destino
    return nueva


def nnunet_a_brats(ruta_pred, ruta_salida):
    """Para despues: devuelve una prediccion de nnU-Net a etiquetas BraTS 0/1/2/4."""
    img = nib.load(ruta_pred)
    seg = reasignar(np.asanyarray(img.dataobj).astype(np.int16), NNUNET_A_BRATS)
    nib.save(nib.Nifti1Image(seg, img.affine, img.header), ruta_salida)


def enlazar(origen, destino):
    """Enlace duro (no ocupa espacio); si el sistema no lo permite, copia."""
    if destino.exists():
        destino.unlink()
    try:
        os.link(origen, destino)
    except OSError:
        shutil.copy2(origen, destino)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--zip", action="store_true", help="crea nnunet_raw/Dataset501_BraTS2018.zip")
    args = parser.parse_args()

    pacientes = p.cargar_pacientes(MANIFEST)
    folds = pd.read_csv(FOLDS)
    ids = {pac.paciente_id for pac in pacientes}
    if set(folds["paciente_id"]) != ids:
        raise SystemExit("particiones/folds.csv y data/manifest.csv no tienen los mismos pacientes")

    (DATASET / "imagesTr").mkdir(parents=True, exist_ok=True)
    (DATASET / "labelsTr").mkdir(parents=True, exist_ok=True)

    for i, pac in enumerate(pacientes, 1):
        for mod, canal in CANALES.items():
            enlazar(pac.rutas()[mod], DATASET / "imagesTr" / f"{pac.paciente_id}_{canal}.nii.gz")
        img = nib.load(pac.seg)
        seg = reasignar(np.asanyarray(img.dataobj).astype(np.int16), BRATS_A_NNUNET)
        etiqueta = nib.Nifti1Image(seg, img.affine, img.header)
        etiqueta.set_data_dtype(np.uint8)
        nib.save(etiqueta, DATASET / "labelsTr" / f"{pac.paciente_id}.nii.gz")
        if i % 50 == 0 or i == len(pacientes):
            print(f"  {i}/{len(pacientes)}")

    dataset = {
        "channel_names": {str(i): mod.upper() for i, mod in enumerate(CANALES)},
        "labels": {"background": 0, "whole tumor": [1, 2, 3],
                   "tumor core": [2, 3], "enhancing tumor": [3]},
        "regions_class_order": [1, 2, 3],
        "numTraining": len(pacientes),
        "file_ending": ".nii.gz",
        "name": "BraTS2018",
        "reference": "MICCAI BraTS 2018 training data; reto CIMAT",
    }
    (DATASET / "dataset.json").write_text(json.dumps(dataset, indent=2))

    splits = [{"train": sorted(folds.loc[folds["fold"] != k, "paciente_id"]),
               "val": sorted(folds.loc[folds["fold"] == k, "paciente_id"])}
              for k in sorted(folds["fold"].unique())]
    (DATASET / "splits_final.json").write_text(json.dumps(splits, indent=2))

    print(f"\nDataset: {DATASET}  ({len(pacientes)} pacientes, {len(splits)} folds)")

    if args.zip:
        destino = shutil.make_archive(str(DATASET), "zip", root_dir=DATASET.parent,
                                      base_dir=DATASET.name)
        print(f"Zip:     {destino}  ({Path(destino).stat().st_size / 1e9:.2f} GB)")


if __name__ == "__main__":
    main()
