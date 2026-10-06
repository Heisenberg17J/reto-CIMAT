# features.csv

Tabla de características radiómicas de BraTS 2018: una fila por paciente, y
`paciente_id` como primera columna.

| Campo | Valor |
|-------|-------|
| Segmentación | manual (BraTS) |
| Fecha de generación | 2026-10-06 14:52 |
| Fecha de extracción | 2026-10-02 08:51 |
| PyRadiomics | v3.0.1 (según `diagnosticos.csv`) |
| Parámetros | `config/params_brats2018_v0.yaml` (sha256 `58f6d14fd794…`) |
| Pacientes procesados | 285 de 285 en el manifest |
| Columnas | 1146 características + `paciente_id` |

## Normalización

Z-score por paciente y por modalidad, calculado **solo sobre los voxeles del cerebro**
(intensidad > 0). El fondo queda en 0. Se hace antes de PyRadiomics
(`radiomica/normalizar.py`), con `normalize: false` en el YAML.

## Regiones extraídas

| Región | Etiquetas de seg |
|--------|------------------|
| WT | 1 + 2 + 4 |
| TC | 1 + 4 |
| ET | 4 |

Modalidades: t1, t1ce, t2, flair.

## Columnas

`<modalidad>_<región>_<clase>_<nombre>` (por ejemplo, `t1ce_ET_glcm_Contrast`). La forma
se extrae una sola vez por región: `mask_<región>_shape_<nombre>`.

Columnas constantes eliminadas en el control de calidad: ninguna.

## Uso

Esta tabla **no** tiene selección de características, filtrado por correlación ni
reducción de dimensionalidad. Esos pasos van dentro de la validación cruzada y se ajustan
solo con los datos de entrenamiento de cada fold.

Ver `DECISIONES.md` en la raíz del repositorio.
