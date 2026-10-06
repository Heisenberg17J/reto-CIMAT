# Reto CIMAT — Gliomas en BraTS 2018

Dos objetivos, con modelos distintos y en este orden:

1. **Segmentación (visión):** un modelo que marca qué parte de la resonancia es tumor y qué parte es cerebro sano.
2. **Pronóstico (predicción):** un modelo que estima la supervivencia a partir de cómo se ve el tumor, usando radiómica.

Están conectados: en uso real, la radiómica del objetivo 2 se calcula sobre la máscara que dibuja el modelo 1. Por eso los dos usan **la misma partición en folds**.

El porqué de cada decisión está en [DECISIONES.md](DECISIONES.md), citado aquí como D*n*.

## Estado (2026-10-06)

| Etapa | Estado |
|---|---|
| Organización y verificación de los datos (285 pacientes) | Terminada |
| Radiómica: `features.csv`, 285 × 1146, con control de calidad | Terminada |
| Folds compartidos (5 × 57 pacientes) | Terminada |
| Segmentación: conversión a nnU-Net y prueba de 5 épocas en Colab | Terminada |
| Segmentación: 5 folds (100 épocas) + posprocesado de ET | Terminada: Dice fuera de fold WT 0.907, TC 0.842, ET 0.769 (D24) |
| Pronóstico | **Siguiente**: radiómica sobre las máscaras predichas y luego el modelo |

## Datos

BraTS 2018, conjunto de entrenamiento: **285 pacientes (210 HGG y 75 LGG)**, cada uno con T1, T1ce, T2, FLAIR y su segmentación. Las etiquetas son 1 necrosis o tumor no realzado, 2 edema y 4 realce.

- La **supervivencia** existe solo para **163 pacientes, todos HGG** (D21).
- **27 LGG no tienen realce** (ET vacío) (D11).
- **Centros:** CBICA aporta el 52 % de los pacientes con supervivencia, y el estado de resección solo está reportado en CBICA y en 2013 (D21).

Los datos **no están en git**. Hay que colocarlos así:

```
data/
  HGG/Brats18_XXXX/Brats18_XXXX_{t1,t1ce,t2,flair,seg}.nii.gz
  LGG/Brats18_XXXX/...
  survival_data.csv
```

## Entorno

```
bash crear_env.sh          # crea el entorno conda "radiomica" (Python 3.11, numpy 1.26, pyradiomics 3.0.1)
conda activate radiomica
```

La segmentación se entrena en Colab, porque la máquina local no tiene GPU (D0, D23).

## Parte 1 · Radiómica (`scripts/`)

Todo se corre desde la raíz del repo, en este orden:

| Bloque | Script | Qué hace | Salida principal |
|---|---|---|---|
| 1 | `organizar_datos.py` | Arma el índice de pacientes y valida la tabla clínica | `data/manifest.csv`, `data/clinica.csv` |
| 2 | `verificacion.py` | Revisa shape, affine, spacing y etiquetas de cada paciente | informe en pantalla |
| 3 | `normalizar.py` | Z-score de cada modalidad sobre los voxeles del cerebro (D2) | `data_normalizada/` |
| 4 | `regiones.py` | Máscaras WT, TC y ET; marca las vacías o pequeñas (D3, D4) | `data_normalizada/manifest_regiones.csv` |
| 5 | `extraccion.py` | PyRadiomics; la forma una vez por región (D6, D13) | `resultados/caracteristicas.csv` |
| 6 | `control_calidad.py` | Pruebas de cordura, columnas, NaN, constantes y filas (D16) | `resultados/caracteristicas_qc.csv` |
| 7 | `exportar.py` | Tabla final más README de procedencia (D18) | `resultados/features.csv` |

```
nohup sh -c "python scripts/organizar_datos.py && python scripts/verificacion.py && \
  python scripts/normalizar.py && python scripts/regiones.py && python scripts/extraccion.py && \
  python scripts/control_calidad.py && python scripts/exportar.py" > logs/flujo_completo.txt 2>&1 &
```

La extracción tarda unos 20 s por paciente, alrededor de 1.5–2 h en total. `nohup` la mantiene corriendo aunque se cierre la terminal. Si se interrumpe, basta con volver a lanzar el flujo: la normalización se salta lo que ya existe.

**Parámetros:** `params_brats2018_v0.yaml`, con binWidth 0.1 sobre z-score y solo la imagen original (D7–D10).

**Columnas de `features.csv`:** `paciente_id` y luego `<modalidad>_<región>_<clase>_<nombre>`, por ejemplo `t1ce_ET_glcm_Contrast`. La forma aparece como `mask_<región>_shape_<nombre>` (D13). Los 27 LGG sin ET tienen NaN en las 382 columnas de ET.

**Regla:** la selección de características, el filtrado por correlación y la reducción de dimensionalidad van **dentro** de la validación cruzada, nunca sobre `features.csv` (D19).

## Folds compartidos (`scripts/folds.py`)

`particiones/folds.csv` asigna cada paciente a 1 de 5 folds de 57, estratificados por grupo (HGG con supervivencia, HGG sin supervivencia y LGG). **No se regenera:** cambiarlo invalida todo lo entrenado (D22).

Se descartó agrupar por centro: CBICA sola es la mitad de los datos. Por eso la validación mide el desempeño con pacientes nuevos de los mismos centros, no la generalización a otro hospital (D22).

## Parte 2 · Segmentación (`segmentacion/`)

1. **Convertir** al formato de nnU-Net (imágenes originales; etiquetas 0/1/2/4 → 0/2/1/3; entrenamiento por regiones WT/TC/ET; nuestros folds) (D23):
   ```
   python segmentacion/convertir_nnunet.py --zip
   ```
   Genera `nnunet_raw/Dataset501_BraTS2018.zip` (2.3 GB, fuera de git).
2. **Subir** el zip a Google Drive, en `MyDrive/reto_cimat/`.
3. **Probar** con `segmentacion/prueba_colab.ipynb`: preprocesa y entrena 5 épocas del fold 0 en una GPU T4.
4. **Entrenar** con `segmentacion/entrenar_colab.ipynb`: 100 épocas por fold, checkpoint cada 5 y reanudación automática. Se ejecuta completo en cada sesión de Colab hasta que el resumen diga que terminó.
5. **Posprocesar** las 285 predicciones fuera de fold (carpeta `predicciones_oof/`, descargada de Drive): `python segmentacion/postproceso_et.py` elige con validación anidada el umbral para descartar ET pequeño y escribe las máscaras finales en `segmentaciones_pred/` (D24).
6. **Ver predicciones** con `segmentacion/visor_colab.ipynb`: visor corte por corte (real frente a predicción), Dice por paciente y análisis del error según el tamaño del tumor.

**Resultado (D23, D24):** 5 folds × 100 épocas en una A100 con RAM amplia (64 s por época, ~1.8 h por fold). Dice fuera de fold de los 285 pacientes: **WT 0.907, TC 0.842 y ET 0.769** tras el posprocesado (0.742 sin él). En los 163 HGG con supervivencia: WT 0.903, TC 0.897 y ET 0.831, y el posprocesado no les cambia nada.

**Siguiente paso:** volver a correr la radiómica (bloques 4–7) sobre `segmentaciones_pred/` para el objetivo 2.

## Decisiones abiertas

- **D10:** si se agregan filtros LoG o Wavelet, lo que obligaría a re-extraer.
- **D11:** cómo tratar el ET vacío de los 27 LGG. No afecta al pronóstico, que usa solo HGG.

## Estructura

```
DECISIONES.md              por qué se hizo cada cosa (D0–D24)
params_brats2018_v0.yaml   parámetros de PyRadiomics
crear_env.sh, requirements.txt
scripts/                   radiómica (bloques 1–7) y folds
segmentacion/              conversión a nnU-Net, notebooks de Colab y posprocesado de ET
particiones/folds.csv      partición compartida (versionada)
predicciones_oof/          predicciones de nnU-Net por fold, descargadas de Drive (fuera de git)
segmentaciones_pred/       máscaras finales fuera de fold, etiquetas BraTS (fuera de git)
data/, data_normalizada/, nnunet_raw/, resultados/, logs/   generados o datos (fuera de git)
```
