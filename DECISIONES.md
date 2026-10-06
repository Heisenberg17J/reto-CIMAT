# Decisiones — Reto CIMAT / BraTS 2018

Registro de las decisiones metodológicas del pipeline de radiómica. Cada supuesto
de `params_brats2018_v0.yaml` (A1, A2, A3) tiene aquí su entrada.

Estados: **Tomada** (implementada) · **Abierta** (pendiente de discutir o de un análisis de sensibilidad).

---

## Entorno

### D0. Versiones fijadas — Tomada
- Python 3.11, numpy 1.26.4 (pyradiomics no funciona con numpy 2).
- pyradiomics se instala aparte con `--no-build-isolation` (ver `crear_env.sh`).
- Por ahora se mantiene pyradiomics **v3.0.1** (la instalada en el entorno `radiomica`).
  El encabezado del params lo refleja y debe actualizarse si se cambia de versión.

---

## Bloque 1 — Organización de los datos (`scripts/organizar_datos.py`)

### D20. Estructura de entrada y manifest — Tomada
Los datos originales se dejan tal como vienen de BraTS 2018 (`data/HGG/<id>/`,
`data/LGG/<id>/`, `data/survival_data.csv`) y no se mueven ni se modifican.
`organizar_datos.py` genera `data/manifest.csv` (`paciente_id`, `grado` y las 5 rutas) y
todo el pipeline lee los pacientes de ahí, no de las carpetas. Un paciente solo entra
al manifest si tiene sus 5 archivos.

Conjunto completo (2026-10-01): **285 pacientes, 210 HGG y 75 LGG**, todos completos.

### D21. Tabla clínica — Tomada
`data/clinica.csv` contiene una fila por paciente con `grado`, `origen`,
`tiene_supervivencia`, `edad`, `supervivencia_dias` y `reseccion`. Es la versión
validada de `survival_data.csv`, que no se toca: ids únicos, todos con imagen, y edad y
supervivencia numéricas (un texto como "ALIVE" detiene el script en lugar de volverse NaN).

- La supervivencia solo existe para **163 pacientes, todos HGG**. 47 HGG y los 75 LGG no
  la tienen. Ninguno de los 4 pacientes de prueba (`Brats18_2013_{2,3,4,5}_1`) la tiene.
- `reseccion`: GTR 59, STR 24 y sin reportar 80. El `NA` del original se guarda vacío.
  **Cuidado:** el estado de resección solo está reportado en CBICA (81 de 85) y 2013 (2 de
  2); en todos los TCIA falta. Un valor vacío no es un dato faltante al azar: identifica
  al centro.
- `origen` es la institución, tomada del id (CBICA, TCIA01…TCIA13, 2013). Salvo 2013,
  ningún origen tiene a la vez HGG y LGG: en una clasificación HGG/LGG, el centro está
  confundido con la clase. Validación agrupada por centro: descartada (ver D22).

  Los 163 con supervivencia por centro: CBICA 85 (52 %), TCIA02 22, TCIA01 19, TCIA03 12,
  TCIA08 9, TCIA04 7, TCIA06 5, 2013 2 y TCIA05 2. La mediana de supervivencia varía entre
  centros (de 131 días en TCIA05 a 446 en TCIA03), aunque con muy pocos pacientes en
  varios de ellos. Una validación agrupada por centro dejaría folds dominados por CBICA.

---

## Bloque 2 — Verificación (`scripts/verificacion.py`)

### D1. Comprobaciones previas a la extracción — Tomada
Para cada paciente se verifica: mismo shape y affine en las 5 imágenes (tolerancia 1e-4),
spacing isotrópico de 1 mm, etiquetas de `seg` ⊆ {0, 1, 2, 4} y que no haya tumor fuera
del cerebro.

Resultado con los 285 (2026-10-01): 239 pacientes sin alertas. No hay problemas de
shape, affine, spacing ni etiquetas inesperadas. Alertas:
- 27 pacientes **sin etiqueta 4 (ET vacío), todos LGG** (ver D11).
- 1 paciente sin edema: Brats18_TCIA13_615_1 (LGG).
- 21 pacientes con tumor fuera del cerebro según t1 > 0 (mediana 10 voxeles; el peor,
  Brats18_TCIA13_634_1, con 588). Tras la normalización esos voxeles valen 0 dentro de
  la máscara.

  **Decisión (2026-10-02): se documenta y no se re-extrae.** Regla acordada: con 0 o 1
  HGG con supervivencia afectados se documenta; con 5 o más se re-extrae.
  - Según t1 = 0, de los 21 pacientes **3 son HGG con supervivencia**: Brats18_2013_11_1,
    Brats18_2013_27_1 y Brats18_TCIA04_343_1.
  - Contando los voxeles ≤ 0 en cualquiera de las 4 modalidades (los que la
    normalización realmente pone en 0), son 28 pacientes y **8 HGG con supervivencia**.
  - El efecto es mínimo: como máximo el **0.34 %** de los voxeles de una región
    (Brats18_TCIA04_343_1, t1/WT); en los demás HGG con supervivencia, menos de 0.12 %.
    En todo el conjunto, el máximo es 0.66 % (Brats18_2013_0_1, t2/TC, LGG).
  - Re-extraer no arreglaría nada sin cambiar antes las máscaras (por ejemplo,
    intersectarlas con el cerebro común a las 4 modalidades), y eso cambiaría también la forma.

**Por qué:** la deduplicación de forma (D6) y `correctMask: false` (D5) dependen de que
las modalidades estén co-registradas.

---

## Bloque 3 — Normalización (`scripts/normalizar.py`)

### D2 (= A1). Z-score solo sobre voxeles del cerebro — Tomada
Cada modalidad se normaliza con media y desviación estándar calculadas **solo sobre los
voxeles > 0**. El fondo queda en 0 y la imagen se guarda en float32. `seg` se copia sin cambios.

**Por qué:** cerca del 82 % de cada volumen BraTS es fondo en 0. Incluirlo (como hace
`normalize: true` de PyRadiomics, que usa `sitk.Normalize` sobre todo el volumen) sesga
la media y la desviación estándar. Por eso en el params `normalize: false`.

**Consecuencia:** las intensidades están en unidades de desviación estándar, lo que afecta
la elección de `binWidth` (D9).

---

## Bloque 4 — Regiones de interés (`scripts/regiones.py`)

### D3 (= A2). Regiones compuestas, no etiquetas sueltas — Tomada
| Región | Etiquetas | Significado        |
|--------|-----------|--------------------|
| WT     | 1 + 2 + 4 | tumor completo     |
| TC     | 1 + 4     | núcleo tumoral     |
| ET     | 4         | tumor con realce   |

Se guarda una máscara binaria por región (valor 1, coherente con `label: 1`) en
`data_normalizada/<paciente>/<paciente>_mask_<región>.nii.gz`. El script comprueba que
ET ≤ TC ≤ WT.

**Por qué:** son las regiones estándar de evaluación en BraTS y tienen sentido clínico.
Las etiquetas sueltas (sobre todo el edema, 2) no son las que se usan como ROI.

### D4. Las regiones inválidas se registran, no se descartan en silencio — Tomada
`regiones.py` aplica los mismos umbrales que PyRadiomics (leídos del params) y marca cada
región como `ok`, `vacia` o `pequena` en `data_normalizada/manifest_regiones.csv`. Solo
las regiones `ok` entran al plan de extracción; las demás quedarán como NaN.

### D5. Geometría sin corrección ni remuestreo — Tomada
`correctMask: false`, sin `resampledPixelSpacing`, `force2D: false`.

**Por qué:** las imágenes ya vienen co-registradas a 1×1×1 mm (verificado en D1). Si
aparece un error de geometría, el proceso debe fallar y no corregirse en silencio.

### D6. Forma una vez por región — Tomada
Las características de `shape` dependen solo de la máscara, así que serían idénticas en
t1, t1ce, t2 y flair. Se extraen **una vez por región** (sobre t1); `firstorder` y las
texturas se extraen una vez por **modalidad × región**. Por paciente: 3 extracciones de
forma + 12 de intensidad, en lugar de 4 columnas repetidas por cada característica de forma.

Implementación: `crear_extractores()` crea los dos extractores a partir del mismo params
(así no pueden divergir) y `plan_extraccion()` arma la lista de tareas para el Bloque 5.

---

## Parámetros de extracción (`params_brats2018_v0.yaml`)

### D7. Clases de características — Tomada
`shape`, `firstorder`, `glcm`, `glrlm`, `glszm`, `gldm`, `ngtdm`. De GLCM se excluye
`SumAverage` (es exactamente 2 × `JointAverage`) y se conserva `MCC`.

### D8. Ajustes de matrices de textura — Tomada
`distances: [1]`, `symmetricalGLCM: true`, `gldm_a: 0`, `voxelArrayShift: 0` (IBSI).
Sin `weightingNorm`: se calcula por cada una de las 13 direcciones 3D y se promedia.

### D9. Discretización — Tomada
Ancho de bin fijo `binWidth: 0.1` sobre la escala z-score. La rejilla de 5–25 del equipo
supone intensidades escaladas ×100; sobre z-scores sin escalar equivale a 0.05–0.25.

**Decisión (2026-10-01):** se mantiene `binWidth: 0.1` sobre z-score sin escalar, con
`voxelArrayShift: 0`. Con los 4 pacientes actuales da entre 16 y 128 bins por ROI
(mediana 55). El criterio de referencia es que `firstorder_Range / binWidth` caiga entre
~30 y 130 bins en la mayoría de las ROI.

**Revisión con los 285 (2026-10-02):** el **90 %** de las ROI cae entre 30 y 130 bins (7 %
por debajo y 3 % por encima). Medianas: WT 67, TC 62, ET 56. Los casos con pocos bins son
ET y TC pequeños (mínimo 6). Se mantiene 0.1.

**Prueba con `binWidth: 25` sin escalar (2026-10-01, 4 pacientes):** descartado. PyRadiomics
pone los bordes en múltiplos del ancho (−25, 0, 25…) y, como los z-scores van de −4 a 8,
todos los ROI quedan en **2 niveles: por debajo y por encima de la media del cerebro**. Las
texturas se calculan sin error y no salen constantes ni NaN, pero describen una imagen
binarizada en 0. En ET, entre 0 % y 10 % de los voxeles caen bajo 0, así que la región es
casi homogénea. Con 0.1 hay entre 16 y 128 bins (mediana 55). El 25 de la rejilla del
equipo solo tiene sentido con intensidades ×100, y equivale a 0.25 aquí.

**Variante B, z-score ×100 + `binWidth: 25` + `voxelArrayShift: 300` (2026-10-01, rama
`prueba/variante-b`):** descartada. Equivale a binWidth 0.25, así que la mediana baja a 22
bins (mínimo 6). En 260 de 1146 columnas cambia el orden de los pacientes, todas de
textura o de Energy/TotalEnergy/RMS. Con el shift de 300, Energy se vuelve casi
proporcional al volumen, información que ya está en `shape`. Elegir la escala no cambia
nada (×100 con binWidth 10 es lo mismo que A); lo que importa es el ancho de bin.

### D10 (= A3). Tipos de imagen — Abierta
Primera pasada solo con `Original`. LoG (σ = 1–5 mm) y Wavelet (coif1, 8 descomposiciones)
quedan comentados en el params. Si se activan, `padDistance: 10` ya cubre el σ más grande.

**Para decidir:** multiplican el número de características (~×5 con LoG, ×8 con Wavelet),
lo que agrava el problema de dimensionalidad con pocos pacientes.

### D11. Regiones pequeñas o vacías — Abierta
`minimumROISize: 27` y `minimumROIDimensions: 3` son provisionales. En el conjunto
completo, **27 de los 75 LGG (36 %) no tienen ET**; ningún HGG está en ese caso, y como la
supervivencia solo es de HGG, no afecta a ese objetivo. En la clasificación HGG/LGG, en
cambio, la ausencia de ET por sí sola ya predice LGG. Cómo se trate (NaN imputado,
indicadora, exclusión) cambia lo que el modelo puede aprender, así que hay que decidirlo
explícitamente.

ET pequeño en HGG (2026-10-02): ningún HGG tiene ET por debajo de 100 voxeles (mínimo 106,
Brats18_CBICA_BHB_1, con supervivencia). 2 tienen menos de 500 y 3 menos de 1000; solo 1
de ellos con supervivencia. Aun así, 51 HGG (34 con supervivencia) tienen al menos una
modalidad con menos de 30 bins en ET: es la región con texturas menos confiables (D9).

El params es la única fuente del umbral: `verificacion.py` y `regiones.py` lo leen de ahí
(`leer_umbrales()`), así que cambiarlo en el params cambia todo el pipeline.

**Opciones para ET vacío:** dejar NaN e imputar, excluir al paciente de los modelos que
usan ET, o añadir una variable indicadora "sin realce".

### D12. Columnas de diagnóstico — Tomada
`additionalInfo: true` conserva las columnas `diagnostics_*` como trazabilidad. No entran
a la tabla del modelo: se guardan en un archivo aparte (D15).

---

## Bloque 5 — Extracción (`scripts/extraccion.py`)

### D13. Convención de nombres de columnas — Tomada
Cada columna tiene siempre 4 partes separadas por `_`:

| Tipo        | Patrón                                  | Ejemplo                     |
|-------------|-----------------------------------------|-----------------------------|
| Intensidad  | `<modalidad>_<región>_<clase>_<nombre>` | `t1ce_ET_glcm_Contrast`     |
| Forma       | `mask_<región>_shape_<nombre>`          | `mask_WT_shape_Sphericity`  |

La forma usa `mask` como modalidad porque no depende de la imagen (D6). Si se activan
filtros (D10), el filtro va dentro de la clase con guiones (`t1ce_ET_wavelet-LLH-glcm_Contrast`),
para que el nombre siga partiéndose en exactamente 4 piezas.

**Por qué:** con `col.split("_")` se puede filtrar por modalidad, región o clase sin
expresiones regulares, y el nombre se entiende sin consultar otra tabla.

### D14. Un fallo no detiene el lote — Tomada
Las excepciones se capturan por paciente. Si falla cualquier llamada, el paciente
completo queda fuera de `caracteristicas.csv` (no se guardan filas parciales) y el log
`logs/extraccion_<fecha>.txt` registra el paciente, la tarea (tipo, modalidad, región) y
el error. Las regiones omitidas en el Bloque 4 (D4) también quedan en el log, como
advertencia, y sus columnas quedan como NaN.

### D15. Diagnósticos aparte de las características — Tomada
Las columnas `diagnostics_*` (versiones, parámetros, hash de imagen y máscara, número de
voxeles, bounding box) se guardan en `resultados/diagnosticos.csv`, con una fila por llamada
a PyRadiomics. `resultados/caracteristicas.csv` contiene solo características y está lista
para el modelo. Esto concreta D12.

---

## Bloque 6 — Control de calidad (`scripts/control_calidad.py`)

### D16. Pruebas críticas y pruebas de revisión — Tomada
| # | Prueba | Tipo |
|---|--------|------|
| 1 | Cordura: `mask_WT_shape_MeshVolume` de Brats18_2013_2_1 ≈ 61 000 mm³ (±5 %); `VoxelVolume` igual al conteo independiente del Bloque 4 en todas las regiones | Crítica |
| 1b | `MeshVolume` a menos de 10 % de `VoxelVolume` | Revisión |
| 2 | Número y nombres de columnas iguales, uno a uno, a los derivados del params (regiones × 14 de forma + modalidades × regiones × 92) | Crítica |
| 3 | NaN e infinitos por columna y por paciente | Revisión |
| 4 | Columnas constantes | Se eliminan (D17) |
| 5 | \|valor\| ≥ 1e15 (división por cero encubierta) | Revisión |
| 6 | Una fila por paciente, `paciente_id` único y presente en el manifest | Crítica |

**Resultado con los 285 (2026-10-02):** pasan todas las críticas. 382 columnas (todas
las de ET) tienen NaN en los 27 LGG sin ET, como se esperaba (D11). Sin infinitos ni
valores ≥ 1e15 (máximo 1.0e9). En 7 ET pequeños (mediana 126 voxeles) MeshVolume se
aleja más de 10 % de VoxelVolume, lo esperable en regiones chicas o fragmentadas.

Si falla una prueba crítica, no se escribe `caracteristicas_qc.csv` y el script termina
con código 1. Los NaN no se tocan aquí: su tratamiento depende de D11.

**Referencia de rangos:** el mayor valor legítimo observado es ~3e8
(`glszm_LargeAreaHighGrayLevelEmphasis` en WT), muy lejos del umbral de 1e15.

### D17. Eliminación de columnas constantes — Tomada
Se elimina toda columna con un único valor finito (o ninguno) en todos los pacientes. Cada
eliminación queda registrada con su motivo en `resultados/columnas_eliminadas.csv`.

**Por qué:** no aportan información y rompen la estandarización (división por varianza
cero) y algunos modelos.

**Resultado con los 285 (2026-10-02):** ninguna columna constante; no se eliminó ninguna. Como no usa la etiqueta
de clase, aplicarlo sobre todos los pacientes no filtra información del objetivo. Las
columnas *casi* constantes y las redundantes (alta correlación) se tratan en la
selección de características, no aquí.

---

## Bloque 7 — Tabla final (`scripts/exportar.py`)

### D18. Entrega y procedencia — Tomada
`exportar.py` toma `caracteristicas_qc.csv` (la tabla que validó el Bloque 6) y escribe
`resultados/features.csv` (con `paciente_id` como primera columna) y, junto a ella,
`resultados/README.md`. Se niega a exportar si `caracteristicas_qc.csv` no existe o es
más antigua que `caracteristicas.csv`: así nunca sale una tabla que no pasó el control
de calidad. El README registra: fecha, versión de PyRadiomics (tomada de los
diagnósticos, es decir, la que realmente se usó), ruta y sha256 del YAML, método de
normalización, regiones, modalidades y número de pacientes procesados.

### D19. Nada de selección fuera de la validación cruzada — Tomada
`features.csv` es el último artefacto que se guarda antes del modelado. La selección de
características, el filtrado por correlación y la reducción de dimensionalidad se ajustan
**dentro de la validación cruzada, solo con los datos de entrenamiento de cada fold**
(por ejemplo, como pasos de un `sklearn.pipeline.Pipeline`). Aplicarlos antes y guardar
el resultado usaría información de los folds de prueba y contaminaría la evaluación.

La única excepción es D17 (columnas con varianza cero): no usa la etiqueta y su resultado
no depende de cómo se dividan los datos.

---

## Folds compartidos (`scripts/folds.py`)

### D22. Una sola partición para los dos objetivos — Tomada
`particiones/folds.csv` asigna cada paciente a 1 de 5 folds (semilla 42), estratificados
por grupo: HGG con supervivencia, HGG sin supervivencia y LGG. Cada fold tiene 57
pacientes, de los cuales 32 o 33 son HGG con supervivencia. Los centros quedan
repartidos, pero no se estratificó por centro.

**Validación agrupada por centro (dejar un centro fuera cada vez): descartada (2026-10-02).**
Es inviable: entre los 163 con supervivencia, CBICA tiene 85 (52 %) y los demás centros
entre 2 y 22. El fold de CBICA se llevaría la mitad de los datos y el resto serían folds
diminutos, con estimaciones inestables.

**Consecuencia:** la validación mide el desempeño con pacientes nuevos de los mismos
centros, no la generalización a un hospital distinto. Hay que decirlo así al reportar. Como
referencia descriptiva (no como validación), se puede comparar el error fuera de fold de
CBICA frente al del resto.

**Por qué:** en el objetivo 2 la radiómica se calculará sobre máscaras predichas por el
segmentador. Si un paciente de prueba del pronóstico se usó para entrenar el segmentador,
su máscara sería demasiado buena y la evaluación quedaría optimista. Con los mismos folds,
la máscara del paciente del fold k sale de un modelo entrenado con los otros 4.

El archivo se versiona y el script se niega a sobrescribirlo: cambiar los folds invalida
todo lo que se entrenó con ellos.

---

## Objetivo 1 — Segmentación (`segmentacion/`)

### D23. nnU-Net v2, entrenamiento por regiones — Tomada
`convertir_nnunet.py` genera `nnunet_raw/Dataset501_BraTS2018`, siguiendo el conversor
oficial de nnU-Net para BraTS:
- **Imágenes originales**, no `data_normalizada`: nnU-Net hace su propio z-score sobre los
  voxeles no nulos, el mismo criterio que D2. Las imágenes son enlaces duros, así que no
  ocupan espacio extra.
- **Etiquetas** 0/1/2/4 → 0/2/1/3 (edema 1, necrosis 2, realce 3), para que las regiones
  sean anidadas: WT = {1,2,3}, TC = {2,3} y ET = {3}. Se verificó que los volúmenes de las
  tres regiones coinciden con el Bloque 4 en los 285 pacientes.
- `splits_final.json` sale de `particiones/folds.csv` (D22) y reemplaza los folds
  aleatorios de nnU-Net.

**Hardware:** la máquina local no tiene GPU (D0). Prueba inicial en Colab con
`segmentacion/prueba_colab.ipynb`: 5 épocas del fold 0 para medir segundos por época y
VRAM.

**Resultado de la prueba (2026-10-02, Colab gratis, T4 de 15.6 GB):**
- **463 s por época** (444–492) y **6.7 GB de VRAM** como máximo.
- Dice en validación del fold 0 tras solo 5 épocas: WT 0.871, TC 0.752, ET 0.639.
- Tiempo estimado, solo entrenamiento: 100 épocas = 12.9 h por fold (64 h los 5 folds);
  250 épocas = 32 h por fold (161 h); 1000 épocas = 129 h por fold (643 h).

**Resuelto (2026-10-05):** A100 con RAM amplia y 100 épocas, detallado abajo. La sospecha de
que la CPU era el cuello de botella se confirmó a medias: en la A100, el uso mediano de la
GPU fue de 51–60 %.

**Entrenamiento real del fold 0 (2026-10-05):** 100 épocas (unas 13 h en una T4) con
`segmentacion/entrenar_colab.ipynb`:
- **Entrenador `nnUNetTrainer_100epochs_ckpt5`:** hereda de `nnUNetTrainer_100epochs` y solo
  guarda `checkpoint_latest` cada 5 épocas en vez de cada 50. Si se corta la sesión se
  pierden como mucho ~40 min, no ~6 h. No redefine `__init__`, porque nnU-Net guarda los
  argumentos del constructor en el checkpoint.
- **Reanudación:** siempre con `--c`, que retoma o empieza de cero; si ya terminó, solo valida.
- **Preprocesado:** se guarda en Drive como `.tar` si cabe, para no repetir 30–60 min por sesión.
- **Sin `--npz`:** las probabilidades ocuparían GB en Drive y solo sirven para ensamblar
  configuraciones.
- **Uso de la GPU:** se registra cada minuto, para confirmar o descartar el cuello de botella
  de CPU.
- **Hardware:** el plan de Google AI del equipo da acceso a A100, L4 y G4 con RAM amplia.
  Se usa **A100 con RAM amplia** (estimado ~3 h por fold, frente a 13 h en la T4); la RAM
  amplia aporta más núcleos de CPU, el probable cuello de botella. El notebook libera la
  GPU al terminar para no gastar unidades de cómputo. Las TPU no sirven: nnU-Net requiere CUDA.
Con 100 épocas el plan de aprendizaje (poly LR) se ajusta a 100, así que no equivale a
cortar en la época 100 un entrenamiento de 1000.

**Resultado del fold 0, 100 épocas (2026-10-05, A100 con RAM amplia):**
- **64 s por época** (1.8 h por fold), con una mediana de uso de GPU de 52 %: la CPU todavía
  la frena un poco.
- **Dice en validación (57 pacientes): WT 0.903 y TC 0.829, que superan la meta, y ET 0.728,
  un poco por debajo de 0.75.**
- El pseudo-Dice de ET al final del entrenamiento fue 0.869. Hipótesis: el ET falso en los
  LGG sin realce cuenta como Dice 0 y baja el promedio. Se revisa con el análisis 6b del
  visor.
- **Posprocesado candidato:** descartar el ET predicho por debajo de un umbral de voxeles,
  como es estándar en BraTS. **El umbral se elige con las predicciones fuera de fold de los
  5 folds**, no mirando un solo fold.
- **Épocas: 100 en los 5 folds (decidido 2026-10-05, con la curva del fold 0).** El
  pseudo-Dice medio se estabiliza hacia la época 70 en ~0.868, y las últimas 30 épocas no
  aportan. Las pérdidas de entrenamiento y validación bajan juntas, con una separación
  pequeña y estable: no hay sobreajuste. Bajar a 75 épocas obligaría a volver a entrenar el
  fold 0 para que los 5 sean comparables, a cambio de unos 30 min por fold: no compensa.
- **Folds 1 a 4:** con `FOLDS = [1, 2, 3, 4]` el notebook los entrena uno tras otro (unas
  7.5 h en la A100), salta los terminados y libera la GPU cuando terminan todos.

**Resultado de los 5 folds, 100 épocas (2026-10-06, A100):**

| Fold | Dice WT | Dice TC | Dice ET | Pacientes sin ET real |
|---|---|---|---|---|
| 0 | 0.903 | 0.829 | 0.728 | 7 |
| 1 | 0.909 | 0.857 | 0.740 | 6 |
| 2 | 0.914 | 0.838 | 0.826 | 1 |
| 3 | 0.905 | 0.846 | 0.666 | 8 |
| 4 | 0.906 | 0.840 | 0.716 | 5 |
| **Media** | **0.907** | **0.842** | **0.735** | 27 |

- WT y TC son estables entre folds y superan la meta. ET varía mucho (0.666 a 0.826) y sigue
  al número de pacientes sin ET real de cada fold: el fold 2, con 1, tiene el mejor; el 3,
  con 8, el peor. Son 5 puntos (ρ = −0.7, no significativo), pero es consistente con que el
  ET falso en LGG sin realce baja el promedio. Los folds se estratificaron por grupo, no por
  ET vacío, y por eso esos casos quedaron repartidos de forma desigual (de 1 a 8).
- Algunos folds registran más de 100 épocas (102–108): la sesión se cortó y retomó, y las
  épocas posteriores al último checkpoint se repitieron. La reanudación funcionó.
- **Siguiente paso:** descargar las 285 predicciones fuera de fold y elegir el umbral de
  descarte de ET con validación anidada (elegirlo en 4 folds y medirlo en el quinto).

**Visor** (`segmentacion/visor_colab.ipynb`): muestra la resonancia, la segmentación real y
la predicción de los pacientes de validación, con su Dice, y analiza si el error depende
del tamaño del tumor. Con un error de borde constante, el error relativo de volumen ya
crece en los tumores pequeños (comprobado con una simulación), así que ese efecto solo no
prueba que el modelo sea peor con ellos. Para eso hace falta una métrica en mm (Hausdorff
95). Pendiente de repetir con el modelo final.

### D24. Posprocesado de ET y máscaras finales fuera de fold — Tomada
`segmentacion/postproceso_et.py`: si el ET predicho tiene menos de T voxeles, se reetiqueta
como necrosis o tumor no realzado (TC y WT no cambian). T se elige con **validación
anidada**: para cada fold, el mejor umbral en los otros 4, aplicado al fold. Se optimiza el
Dice de ET con la convención oficial de BraTS (sin ET real ni predicho = 1).

Resultado (2026-10-06, 285 pacientes fuera de fold):

| | Antes | Después |
|---|---|---|
| WT | 0.907 | 0.907 |
| TC | 0.842 | 0.842 |
| ET (BraTS) | 0.742 | **0.769** |
| ET (nnU-Net) | 0.735 | **0.751** |

- El modelo predice ET en 20 de los 27 pacientes sin ET real, entre 1 y 3621 voxeles.
- **Es un intercambio, no una mejora gratuita:** elimina 13 ET falsos, pero también pierde 14
  ET reales pequeños (12 LGG y 2 HGG sin supervivencia). El umbral es inestable entre folds
  (de 400 a 1000), y en el fold 2, que tiene un solo paciente sin ET, el posprocesado
  empeora el resultado (0.829 → 0.807).
- **Al pronóstico no le afecta:** en los 163 HGG con supervivencia no se descarta ningún ET y
  todos conservan su ET predicho. Su Dice fuera de fold es WT 0.903, TC 0.897 y ET 0.831.
- Umbral para casos nuevos, elegido con los 5 folds: 500 voxeles
  (`resultados/postproceso_et.json`).

Salida: `segmentaciones_pred/<paciente>_seg.nii.gz` (etiquetas BraTS 0/1/2/4). Cada máscara
sale de un modelo que no vio a ese paciente y de un umbral que no se eligió mirándolo. La
tabla por paciente está en `resultados/segmentacion_oof.csv`.

---

## Objetivo 2 — Pronóstico

### D25. Radiómica sobre las máscaras predichas — Tomada
En uso real no hay segmentación manual: la radiómica saldrá de la máscara del modelo 1. Por
eso el modelo de pronóstico se entrena y evalúa con características extraídas de
`segmentaciones_pred/` (D24), y no con las de la segmentación manual de BraTS.

**Cómo se separan las variantes** (`scripts/variante.py`): la variable de entorno
`SEGMENTACION` (`manual` por defecto, o `pred`) elige de dónde sale la segmentación, dónde
se escriben las máscaras y el manifest de regiones, la carpeta de resultados y el prefijo
de los logs. Las imágenes normalizadas son las mismas en las dos variantes. Así ninguna
variante sobrescribe a la otra.

| | manual | pred |
|---|---|---|
| Máscaras WT/TC/ET | `data_normalizada/<id>/` | `segmentaciones_pred/mascaras/` |
| Manifest de regiones | `data_normalizada/manifest_regiones.csv` | `segmentaciones_pred/manifest_regiones.csv` |
| Tablas | `resultados/` | `resultados/pred/` |
| Logs | `logs/<bloque>_<fecha>` | `logs/<bloque>_pred_<fecha>` |

Se verificó que la variante manual no cambia: `features.csv`, `caracteristicas.csv` y las
máscaras conservan su sha256. En `manifest_regiones.csv` solo cambia la columna `mensaje`
("ya existe"). El README de `features.csv` indica qué segmentación se usó.

**Resultado (2026-10-06, 285 pacientes, ~1.5 h):** `resultados/pred/features.csv`, de 285 × 1146,
con las mismas columnas, en el mismo orden, y los mismos pacientes que la tabla manual.
- Extracción sin fallos; todas las pruebas críticas del control de calidad pasan. La
  cordura dio 61 403 mm³ y VoxelVolume coincide con las máscaras predichas. No hay
  columnas constantes ni infinitos, y el máximo es 1.1e9.
- Regiones vacías en las máscaras predichas: ET en 36 pacientes (34 LGG y 2 HGG sin
  supervivencia) y TC en 2, así que 764 columnas tienen algún NaN. **Los 163 HGG con
  supervivencia no tienen ningún NaN.**
- Siguiente paso (fase B): medir, característica por característica, cuánto cambia entre
  la máscara manual y la predicha en esos 163 pacientes.
