# Inferencia de punta a punta (prototipo de investigación)

> **No es una herramienta clínica.** El pronóstico depende casi solo de la edad: en la validación cruzada, el
> modelo (c-index 0.622) no supera a la edad sola (0.624) (D27, D30).

`predecir.py` procesa un paciente nuevo igual que a los 285 del estudio:

```
4 resonancias (formato BraTS) + edad
 → nnU-Net: los 5 modelos de los folds, promediados
 → posprocesado de ET (umbral de 500 voxeles, D24)
 → z-score → regiones WT/TC/ET → PyRadiomics (config/params_brats2018_v0.yaml)
 → Cox Elastic Net final (modelos/pronostico_coxnet.joblib)
 → volúmenes, riesgo, clase de supervivencia (corta / media / larga) y días estimados
```

## Forma recomendada: segmentar en Colab y lo demás en local

| Paso | Dónde | Cómo |
|---|---|---|
| Segmentación (nnU-Net, 5 folds con TTA) + posprocesado de ET | Colab con GPU | `inferencia/segmentar_colab.ipynb`: lee las resonancias de `MyDrive/reto_cimat/casos/<caso>/` y guarda ahí `segmentacion.nii.gz` |
| Radiómica + pronóstico | Local | `python inferencia/predecir.py --caso CARPETA --edad 62 --segmentacion segmentacion.nii.gz` |

Así nnU-Net usa la GPU y los pesos que ya están en Drive, sin descargarlos, y la radiómica corre en el entorno validado.
PyRadiomics no se instala en Colab a propósito: con otras versiones de numpy podría dar características distintas
sin que nada avise. Con `--segmentacion` basta el entorno `radiomica`; PyTorch no hace falta.

El notebook aplica las mismas reglas que `predecir.py`. Se comprobó con predicciones reales de nnU-Net que encuentra
los mismos archivos y que el posprocesado da máscaras idénticas, tanto si descarta el ET como si lo conserva.

Los pasos 1 y 2 de abajo solo hacen falta para correr **todo** localmente, también nnU-Net.

## 1. Entorno `inferencia`

Es una copia del entorno `radiomica` más PyTorch (CPU) y nnU-Net. numpy queda fijo en 1.26, que PyRadiomics necesita:

```
conda create -y --clone radiomica -n inferencia
conda activate inferencia
pip install -c inferencia/restricciones.txt --extra-index-url https://download.pytorch.org/whl/cpu torch "nnunetv2==2.8.1"
```

`inferencia/restricciones.txt` fija `numpy==1.26.4`, `scikit-learn==1.5.2`, `scipy==1.14.1`, `pandas==2.2.3`, `SimpleITK==2.5.6`,
`nibabel==5.3.2` y `scikit-image==0.24.0` (las versiones de `requirements.txt`). Usar la misma versión de nnU-Net que se
usó para entrenar en Colab: el notebook de entrenamiento la imprime en la sección 2.

## 2. Pesos de nnU-Net (desde Drive)

Solo hacen falta el modelo final de cada fold y los archivos de planes. En Colab, con Drive montado:

```python
R = '/content/drive/MyDrive/reto_cimat/nnUNet_results/Dataset501_BraTS2018/nnUNetTrainer_100epochs_ckpt5__nnUNetPlans__3d_fullres'
!cd "{R}/../.." && zip -r /content/drive/MyDrive/reto_cimat/pesos_nnunet.zip \
    Dataset501_BraTS2018/nnUNetTrainer_100epochs_ckpt5__nnUNetPlans__3d_fullres/{plans.json,dataset.json,dataset_fingerprint.json} \
    Dataset501_BraTS2018/nnUNetTrainer_100epochs_ckpt5__nnUNetPlans__3d_fullres/fold_*/checkpoint_final.pth
```

Descargar `pesos_nnunet.zip` y descomprimirlo en `datos/nnunet_results/`, de modo que quede:

```
datos/nnunet_results/Dataset501_BraTS2018/nnUNetTrainer_100epochs_ckpt5__nnUNetPlans__3d_fullres/
    plans.json  dataset.json  dataset_fingerprint.json  fold_0/checkpoint_final.pth ... fold_4/
```

## 3. Uso

```
conda activate inferencia
python inferencia/predecir.py --caso CARPETA --edad 62
```

- `CARPETA` contiene `t1`, `t1ce`, `t2` y `flair` en NIfTI, con nombre simple (`t1.nii.gz`) o de BraTS (`<id>_t1.nii.gz`).
- Las imágenes deben estar en **formato BraTS**: sin cráneo, co-registradas al atlas SRI24, voxel de 1 mm y tamaño
  240 × 240 × 155. Si no lo están, el script avisa. Una resonancia clínica cruda se prepara antes con BraTS Toolkit o CaPTk.
- **Opciones:**
  - `--segmentacion mascara.nii.gz`: usa una máscara ya hecha y omite nnU-Net (sirve sin los pesos ni PyTorch);
  - `--folds 0`: un solo modelo, más rápido;
  - `--tta`: aumento en prueba con espejos; algo mejor, pero unas 8 veces más lento.
- **Tiempo:** en CPU, nnU-Net tarda varios minutos por paciente con los 5 folds; en GPU, segundos. La radiómica tarda unos 20 s.

**Salida**, en `resultados_inferencia/<caso>/`:
- `segmentacion.nii.gz` (etiquetas BraTS);
- `caracteristicas.csv`;
- `resultado.json`;
- `vista.png`.

El pronóstico solo se calcula si el tumor tiene realce (ET), como todos los HGG con los que se entrenó el modelo. Si no lo
tiene, el script lo dice y no inventa un resultado.

## Comprobaciones hechas

Con la máscara fuera de fold de Brats18_CBICA_AAP_1:
- las 1146 características coinciden con las de `resultados/radiomica/pred/features.csv` (diferencia 2e-16);
- el riesgo coincide con el del modelo final aplicado a esa fila (−0.778728, unos 522 días).

Un caso sin ET predicho devuelve "no aplica".

Ojo: los 285 pacientes de BraTS 2018 se usaron para entrenar los modelos, así que estas pruebas comprueban que el flujo
funciona, **no** su exactitud. Para medirla hacen falta casos externos: BraTS 2018 validación, o los pacientes de BraTS
2019/2020 que no estén en 2018.
