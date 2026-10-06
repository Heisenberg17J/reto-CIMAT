# Reto CIMAT — Gliomas en BraTS 2018

Dos objetivos, con modelos distintos y en este orden:

1. **Segmentación (visión):** un modelo que marca qué parte de la resonancia es tumor y qué parte es cerebro sano.
2. **Pronóstico (predicción):** un modelo que estima la supervivencia a partir de cómo se ve el tumor, usando radiómica.

Están conectados: en uso real, la radiómica del objetivo 2 se calcula sobre la máscara que dibuja el modelo 1. Por eso los dos usan **la misma partición en folds**.

El porqué de cada decisión está en [DECISIONES.md](DECISIONES.md), citado aquí como D*n*.

## Estado (2026-10-06): fase experimental terminada

| Etapa | Estado |
|---|---|
| Organización y verificación de los datos (285 pacientes) | Terminada |
| Radiómica: `features.csv`, 285 × 1146, con control de calidad | Terminada |
| Folds compartidos (5 × 57 pacientes) | Terminada |
| Segmentación: conversión a nnU-Net y prueba de 5 épocas en Colab | Terminada |
| Segmentación: 5 folds (100 épocas) + posprocesado de ET | Terminada: Dice fuera de fold WT 0.907, TC 0.842, ET 0.769; HD95 mediano 3.6 / 3.5 / 2.2 mm (D24, D28) |
| Pronóstico: radiómica sobre máscaras predichas | Terminada: `resultados/radiomica/pred/features.csv`, 285 × 1146, sin NaN en los 163 HGG con supervivencia (D25) |
| Pronóstico: robustez frente a la segmentación | Terminada: 70 % con CCC ≥ 0.85, 82 % conservan el orden (D26) |
| Pronóstico: genético frente a Elastic Net (5 × 10 folds) | Terminada: edad 0.624, Elastic Net 0.622, genético 0.605 (c-index); la radiómica no supera a la edad y el genético sobreajusta (+0.09) (D27) |

## Datos

BraTS 2018, conjunto de entrenamiento: **285 pacientes (210 HGG y 75 LGG)**, cada uno con T1, T1ce, T2, FLAIR y su segmentación. Las etiquetas son 1 necrosis o tumor no realzado, 2 edema y 4 realce.

- La **supervivencia** existe solo para **163 pacientes, todos HGG** (D21).
- **27 LGG no tienen realce** (ET vacío) (D11).
- **Centros:** CBICA aporta el 52 % de los pacientes con supervivencia, y el estado de resección solo está reportado en CBICA y en 2013 (D21).

Los datos **no están en git**. Hay que colocarlos así:

```
datos/brats2018/
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

## Parte 1 · Radiómica (`radiomica/`)

Todo se corre desde la raíz del repo, en este orden:

| Bloque | Script | Qué hace | Salida principal |
|---|---|---|---|
| 1 | `organizar_datos.py` | Arma el índice de pacientes y valida la tabla clínica | `datos/brats2018/manifest.csv`, `datos/brats2018/clinica.csv` |
| 2 | `verificacion.py` | Revisa shape, affine, spacing y etiquetas de cada paciente | informe en pantalla |
| 3 | `normalizar.py` | Z-score de cada modalidad sobre los voxeles del cerebro (D2) | `datos/normalizada/` |
| 4 | `regiones.py` | Máscaras WT, TC y ET; marca las vacías o pequeñas (D3, D4) | `datos/normalizada/manifest_regiones.csv` |
| 5 | `extraccion.py` | PyRadiomics; la forma una vez por región (D6, D13) | `resultados/radiomica/manual/caracteristicas.csv` |
| 6 | `control_calidad.py` | Pruebas de cordura, columnas, NaN, constantes y filas (D16) | `resultados/radiomica/manual/caracteristicas_qc.csv` |
| 7 | `exportar.py` | Tabla final más README de procedencia (D18) | `resultados/radiomica/manual/features.csv` |

```
nohup sh -c "python radiomica/organizar_datos.py && python radiomica/verificacion.py && \
  python radiomica/normalizar.py && python radiomica/regiones.py && python radiomica/extraccion.py && \
  python radiomica/control_calidad.py && python radiomica/exportar.py" > logs/flujo_completo.txt 2>&1 &
```

La extracción tarda unos 20 s por paciente, alrededor de 1.5–2 h en total. `nohup` la mantiene corriendo aunque se cierre la terminal. Si se interrumpe, basta con volver a lanzar el flujo: la normalización se salta lo que ya existe.

**Parámetros:** `config/params_brats2018_v0.yaml`, con binWidth 0.1 sobre z-score y solo la imagen original (D7–D10).

**Columnas de `features.csv`:** `paciente_id` y luego `<modalidad>_<región>_<clase>_<nombre>`, por ejemplo `t1ce_ET_glcm_Contrast`. La forma aparece como `mask_<región>_shape_<nombre>` (D13). Los 27 LGG sin ET tienen NaN en las 382 columnas de ET.

**Variante con las máscaras predichas (objetivo 2, D25):** los mismos bloques 4–7 con `SEGMENTACION=pred` usan `datos/segmentaciones_pred/` y escriben en `resultados/radiomica/pred/`, sin tocar la versión manual:
```
SEGMENTACION=pred nohup sh -c "python radiomica/regiones.py && python radiomica/extraccion.py && \
  python radiomica/control_calidad.py && python radiomica/exportar.py" > logs/flujo_pred.txt 2>&1 &
```

**Regla:** la selección de características, el filtrado por correlación y la reducción de dimensionalidad van **dentro** de la validación cruzada, nunca sobre `features.csv` (D19).

## Folds compartidos (`radiomica/folds.py`)

`particiones/folds.csv` asigna cada paciente a 1 de 5 folds de 57, estratificados por grupo (HGG con supervivencia, HGG sin supervivencia y LGG). **No se regenera:** cambiarlo invalida todo lo entrenado (D22).

Se descartó agrupar por centro: CBICA sola es la mitad de los datos. Por eso la validación mide el desempeño con pacientes nuevos de los mismos centros, no la generalización a otro hospital (D22).

## Parte 2 · Segmentación (`segmentacion/`)

1. **Convertir** al formato de nnU-Net (imágenes originales; etiquetas 0/1/2/4 → 0/2/1/3; entrenamiento por regiones WT/TC/ET; nuestros folds) (D23):
   ```
   python segmentacion/convertir_nnunet.py --zip
   ```
   Genera `datos/nnunet_raw/Dataset501_BraTS2018.zip` (2.3 GB, fuera de git).
2. **Subir** el zip a Google Drive, en `MyDrive/reto_cimat/`.
3. **Probar** con `segmentacion/prueba_colab.ipynb`: preprocesa y entrena 5 épocas del fold 0 en una GPU T4.
4. **Entrenar** con `segmentacion/entrenar_colab.ipynb`: 100 épocas por fold, checkpoint cada 5 y reanudación automática. Se ejecuta completo en cada sesión de Colab hasta que el resumen diga que terminó.
5. **Posprocesar** las 285 predicciones fuera de fold (carpeta `datos/predicciones_oof/`, descargada de Drive): `python segmentacion/postproceso_et.py` elige con validación anidada el umbral para descartar ET pequeño y escribe las máscaras finales en `datos/segmentaciones_pred/` (D24).
6. **Métricas oficiales de BraTS** (Dice, Hausdorff 95, sensibilidad y especificidad): `python segmentacion/metricas_segmentacion.py` (D28).
7. **Ver predicciones** con `segmentacion/visor_colab.ipynb`: visor corte por corte (real frente a predicción), Dice por paciente y análisis del error según el tamaño del tumor.

**Resultado (D23, D24):** 5 folds × 100 épocas en una A100 con RAM amplia (64 s por época, ~1.8 h por fold). Dice fuera de fold de los 285 pacientes: **WT 0.907, TC 0.842 y ET 0.769** tras el posprocesado (0.742 sin él). Hausdorff 95 mediano: WT 3.6, TC 3.5 y ET 2.2 mm. En los 163 HGG con supervivencia: Dice WT 0.903, TC 0.897 y ET 0.831, y el posprocesado no les cambia nada.

## Parte 3 · Pronóstico (`pronostico/`)

1. **Robustez** frente a la segmentación: `python pronostico/robustez_segmentacion.py` (D26).
2. **Comparación** de edad sola, Cox Elastic Net y genético + Cox con validación cruzada anidada de 5 × 10 folds, unos 14 minutos en CPU (D27):
   ```
   python pronostico/evaluar.py                      # principal: máscaras predichas + filtro D26
   python pronostico/evaluar.py --variante manual    # secundarios: --sin-filtro, --reseccion
   python pronostico/evaluar.py --permutar           # control: supervivencia barajada
   ```
3. **Análisis:** `pronostico/analisis.ipynb`, con el kernel `radiomica`, local.

**Resultado (c-index fuera de fold):** edad 0.624, Elastic Net 0.622 y genético 0.605, sin diferencias significativas. El genético promete 0.70 en su aptitud interna (brecha de +0.09) y su selección es casi aleatoria entre folds (Nogueira 0.05).

## Decisiones abiertas

Ninguna. D10 (solo imagen original) y D11 (ET vacío) quedaron cerradas el 2026-10-06.

## Referencia

Bakas S. et al. *Identifying the Best Machine Learning Algorithms for Brain Tumor Segmentation, Progression Assessment, and Overall Survival Prediction in the BRATS Challenge.* arXiv:1811.02629 (2018). Define las métricas de segmentación (Dice y Hausdorff 95) y las clases de supervivencia (< 10, 10–15 y > 15 meses).

## Estructura

```
README.md, DECISIONES.md        qué se hizo y por qué (D0–D29)
requirements.txt, crear_env.sh  entorno "radiomica"
config/                         parámetros de PyRadiomics
docs/                           documentos de apoyo al diseño del pronóstico y ejercicios iniciales
radiomica/                      bloques 1–7, variante manual/predicha y folds
segmentacion/                   nnU-Net: conversión, notebooks de Colab, posprocesado y métricas
pronostico/                     robustez, modelos (edad, Elastic Net, genético) y notebook de análisis
particiones/                    folds.csv (D22) y folds_pronostico.csv (5 × 10, D27)
resultados/                     tablas finales, versionadas
  radiomica/manual|pred/          features.csv + README de procedencia
  segmentacion/                   Dice/HD95 por paciente, umbral de ET
  pronostico/                     robustez y los 4 análisis de la comparación
datos/                          fuera de git: BraTS, normalizadas, nnU-Net, predicciones y máscaras
logs/                           fuera de git
```

Las tablas intermedias pesadas (`caracteristicas*.csv`, `diagnosticos.csv`) no se versionan: se regeneran con los bloques 5 y 6.
