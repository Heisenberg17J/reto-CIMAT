# Metodología del pronóstico: independencia de las máscaras y filtro de robustez

Redactado el 2026-10-10 para la sección de métodos y métricas finales. Cada afirmación se verificó contra el
código, los archivos versionados y el historial de git (las horas son del 2026-10-06, UTC−5).

## 1. ¿Cada máscara salió de un modelo que no había visto a ese paciente?

**Sí, en los 285 pacientes.** Se comprobó así:

| Comprobación | Resultado |
|---|---|
| `splits_final.json` de nnU-Net frente a `particiones/folds.csv` (D22) | Idénticos: 5 folds de 57 y ningún paciente en entrenamiento y validación a la vez |
| Predicciones de `datos/predicciones_oof/fold_k/validation/` | Cada carpeta tiene justo los 57 pacientes del fold k y ninguno de los 228 con los que se entrenó el modelo k |
| Posprocesado de ET (`segmentacion/postproceso_et.py`, D24) | El umbral del fold k se eligió con los otros 4 (400, 1000, 1000, 500 y 500 voxeles). En los 163 HGG con supervivencia no se descartó ningún ET |
| Extracción radiómica (D2, D7–D10) | Cada paciente se procesa por separado: z-score dentro de su propio cerebro y binWidth fijo de 0.1. Su valor no depende de los demás pacientes |

**Texto para métodos:** *"Las características radiómicas de cada paciente se extrajeron de la segmentación
predicha por el modelo nnU-Net del fold en que ese paciente estaba en validación, es decir, de un modelo
entrenado con los otros cuatro folds (228 pacientes). El umbral de posprocesado del realce se eligió con
validación anidada. Ningún paciente contribuyó al entrenamiento del segmentador que generó su propia
máscara."*

Dos matices que conviene declarar:

- **Las particiones del pronóstico se repiten.** La repetición 0 de `folds_pronostico.csv` son los folds de
  D22; las repeticiones 1 a 9 reparten a los pacientes de otra forma. En esas repeticiones, los pacientes de
  prueba del Cox pueden haber estado en el entrenamiento del segmentador de *otro* paciente de prueba. Esto
  no filtra información: la máscara de cada paciente sigue saliendo de un modelo que no lo vio, y el
  segmentador nunca usó la supervivencia (D27).
- **El plan de nnU-Net** (`plan_and_preprocess`) se calculó con los 285 casos. Es estándar y no usa
  etiquetas en la práctica: todos los volúmenes de BraTS tienen 1 mm isotrópico y 240×240×155, y en
  resonancia el z-score se calcula por imagen.

## 2. Filtro de robustez (D26): documentación completa

### 2.1 Definición

Para cada una de las 1146 características se comparó su valor con la máscara manual y con la predicha en los
163 HGG con supervivencia (`pronostico/robustez_segmentacion.py` → `resultados/pronostico/robustez_segmentacion.csv`):

- **CCC de Lin:** acuerdo en valor absoluto. Baja si la máscara predicha desplaza o reescala los valores.
- **ρ de Spearman:** conservación del orden de los pacientes.

**Criterio de inclusión:** ρ de Spearman (manual frente a predicha) ≥ 0.85. Lo aplica
`pronostico/datos.py` (`CORTE_ROBUSTEZ = 0.85`). Quedan 935 de 1146 características: WT 318, TC 326 y
ET 291 de 382 cada una. La edad entra siempre, sin penalizar.

| | Spearman ≥ 0.85 | Spearman < 0.85 | Total |
|---|---|---|---|
| CCC ≥ 0.85 | 765 | 40 | 805 |
| CCC < 0.85 | 170 | 171 | 341 |
| **Total** | **935** | **211** | 1146 |

### 2.2 ¿Se seleccionaron antes de la validación cruzada externa o dentro de cada partición?

**Antes, una sola vez.** La lista de 935 es fija y se aplica al cargar la cohorte (`cargar_cohorte`), antes
de crear las particiones externas. No se recalcula dentro de cada fold.

Todo lo que sí usa la supervivencia ocurre **dentro** de cada fold externo, solo con su entrenamiento (D19):

- el camino de penalización de Elastic Net y la elección de `alpha`;
- los grupos de |ρ| > 0.9 del genético (unos 266 representantes);
- la búsqueda del genético;
- los cortes de las 3 clases.

### 2.3 ¿Se calcularon con los 163 pacientes completos, incluidos los de prueba?

**Sí.** El Spearman y el CCC usan los 163 pacientes, incluidos los que después caen en prueba en cada fold
externo. Hay que decirlo explícitamente. Por qué no invalida la evaluación:

1. **No usa la supervivencia.** Mide el acuerdo entre dos maneras de calcular la misma característica, igual
   que el filtro de columnas constantes (D17). No puede elegir características porque predigan bien en los
   pacientes de prueba.
2. **Sí usa información de los pacientes de prueba** (sus máscaras manuales y predichas). Es una selección
   no supervisada y transductiva. En un caso nuevo no hay máscara manual, así que el filtro solo es aplicable
   como una lista fija definida antes.
3. **El efecto medido es nulo o negativo.** Sin el filtro (1146 características), Elastic Net da 0.628,
   frente a 0.622 con él; el genético, 0.608 frente a 0.605. El filtro no infla los resultados.
4. **Se comprobó con el procedimiento estricto (D31).** Con el filtro calculado solo con el entrenamiento de
   cada fold externo, con las mismas particiones y semillas, Elastic Net da 0.622 (diferencia pareada −0.000,
   p = 0.86) y el genético 0.602 (−0.003, p = 0.90). Cada fold retiene 930 ± 21 características (Jaccard 0.95
   entre folds).

**Texto para métodos:** *"Antes del modelado se excluyeron las características cuyo orden entre pacientes
dependía del contorno (ρ de Spearman < 0.85 entre la máscara manual y la predicha, en los 163 pacientes).
Este filtro no usa la supervivencia y se aplicó una sola vez a toda la cohorte; como análisis de
sensibilidad, el modelo se reentrenó sin él y con el filtro recalculado dentro de cada fold externo, usando
solo los pacientes de entrenamiento; ninguna de las dos variantes cambió las conclusiones."*

### 2.4 ¿Se definieron los criterios y umbrales antes de ver los resultados?

Depende de qué resultados:

| Elemento | ¿Cuándo se fijó? |
|---|---|
| Umbral 0.85 | Antes de cualquier resultado: es el corte habitual de radiómica y figura como `CORTE = 0.85` en el script desde su primera versión |
| Spearman en lugar de CCC como criterio | **Después** de ver la tabla de robustez (CCC ≥ 0.85 en el 70 %, 170 características con desplazamiento sistemático), pero **antes** de cualquier modelo de supervivencia |
| Evidencia temporal | D26 y su filtro se versionaron a las 13:28 (commit `21fd008`). La primera corrida del pronóstico (control permutado) es de las 13:33 y la principal de las 13:34 |
| Análisis con y sin filtro | Declarados en D26 ("se evaluará con y sin él") antes de correrlos |

Redacción honesta: *"El umbral se fijó a priori. La elección de Spearman se hizo tras inspeccionar el acuerdo
entre máscaras, sin ver ningún resultado de supervivencia."*

### 2.5 ¿Por qué Spearman como criterio y CCC como medida complementaria?

- **El modelo se entrena y evalúa siempre con la misma variante, la predicha (D25).** Un desplazamiento o un
  cambio de escala sistemático afecta a todos los pacientes por igual. La estandarización dentro del fold lo
  absorbe, y un modelo de Cox, que solo depende del orden de riesgo, no cambia. Lo que sí daña es que el
  error de contorno reordene a los pacientes, y eso es lo que mide Spearman.
- **El CCC castiga desplazamientos que no afectan al modelo.** Ejemplo: `t1ce_ET_glszm_ZoneVariance` tiene
  CCC 0.08 y Spearman 0.92, porque sale 3.8 veces mayor con la máscara predicha, pero ordena igual a los
  pacientes. Filtrar por CCC habría quitado 170 características válidas para este diseño.
- **El CCC sigue siendo útil para describir y advertir.** Muestra que las máscaras de nnU-Net son más suaves
  (área de superficie −20 %, esfericidad +23 %, volumen −2 %) y prohíbe mezclar variantes: un modelo
  entrenado con características manuales no se puede aplicar a las predichas (D26, consecuencia 1).

### 2.6 Inconsistencia a corregir en el repositorio

En `robustez_segmentacion.py`, la docstring describe el CCC como el criterio y la columna `robusta` del CSV
es `ccc >= 0.85` (805 características). El filtro que usa el modelo es la columna `spearman` (935). Al citar
números hay que usar los de `spearman`. Conviene actualizar la docstring para que no confunda.

## 3. Métricas finales (análisis principal)

Cohorte de 163 HGG, todos con evento. Características predichas con el filtro D26 más la edad. Validación
anidada con 5 folds × 10 repeticiones externas (50 folds) y 3 × 2 internas. Media ± DE sobre los 50 folds
externos (`resultados/pronostico/pred/resumen.txt`).

| Brazo | c-index | Exactitud 3 clases (azar 0.345) | ρ Spearman | RMSE (días) | k |
|---|---|---|---|---|---|
| Edad sola (Cox) | **0.624 ± 0.046** | 0.454 ± 0.071 | 0.358 | 333 | 0 |
| Cox Elastic Net | 0.622 ± 0.046 | **0.477 ± 0.072** | 0.352 | 336 | 9.2 |
| Genético + Cox ridge | 0.605 ± 0.054 | 0.440 ± 0.067 | 0.300 | 354 | 6.9 |

Diferencias pareadas de c-index con la t corregida de Nadeau-Bengio:

- Elastic Net − edad: −0.002 (p = 0.82).
- Genético − edad: −0.020 (p = 0.49).
- Genético − Elastic Net: −0.017 (p = 0.50).

Ninguna diferencia es significativa.

Diagnósticos del genético:

- Su aptitud interna es 0.695, frente a 0.605 fuera de fold: una brecha de +0.090.
- La selección es casi aleatoria entre folds (Nogueira 0.053, frente a 0.193 de Elastic Net).

Análisis de sensibilidad (c-index con Elastic Net):

| Análisis | c-index |
|---|---|
| Máscaras manuales | 0.628 |
| Sin filtro | 0.628 |
| Con resección | 0.617 |
| Supervivencia permutada (control) | Elastic Net 0.525, genético 0.527 |

Segmentación (D24, D28), fuera de fold en los 285 pacientes:

- Dice: WT 0.907, TC 0.842 y ET 0.769.
- HD95 mediano: 3.6, 3.5 y 2.2 mm.
- En los 163 HGG con supervivencia, Dice WT 0.903, TC 0.897 y ET 0.831.

Limitación que debe acompañar a estas cifras: la validación mide el desempeño con pacientes nuevos de los
mismos centros, no con un hospital distinto (D22).
