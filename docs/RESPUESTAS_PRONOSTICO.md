# Modelo de pronóstico: respuesta a las cuatro preguntas

Contexto: README del 2026-10-06 (radiómica sobre máscaras predichas lista, folds compartidos 5 × 57, 163 HGG
con supervivencia). Evidencia:

- **[S]** simulaciones con datos sintéticos del mismo tamaño (n = 163, 1146 características correlacionadas
  en bloques, supervivencia guiada por la edad y opcionalmente por una señal radiómica real). Código en
  `sim_lib.py`, resultados en `sim_metric_noise.csv` y `sim_ga_variants_*.csv`.
- **[D]** lo ya registrado en el README y en DECISIONES.md.

Las simulaciones indican la dirección de los efectos; no predicen sus números reales.

## Resumen

| Pregunta | Recomendación |
|---|---|
| 1. Objetivo | **Modelo de riesgo sobre los días** (Cox). Las 3 clases de BraTS **se derivan** de ese riesgo para reportarlas; no se entrena un clasificador de 3 clases ni una regresión de días |
| 2. Papel del genético | **Solo selección de características.** Los hiperparámetros, con búsqueda en rejilla dentro del fold interno |
| 3. Método tradicional | **Cox con Elastic Net** (`l1_ratio` 0.5, edad sin penalizar). Edad sola como suelo. RSF solo como análisis secundario |
| 4. Implementación | **DEAP 1.4.4** |

---

## 1. Objetivo: entrenar sobre los días, reportar también las 3 clases

Son dos preguntas distintas: qué variable aprende el modelo y con qué métrica se reporta.

**Qué aprende el modelo: un riesgo ordenado (Cox) entrenado sobre los días.**

- Un clasificador de 3 clases trata "corto, medio y largo" como categorías sin orden. Además, un paciente de
  299 días y otro de 301 caen en clases distintas aunque casi no se diferencian. Con 163 pacientes, ese
  desperdicio de información cuesta.
- Una regresión de días (MSE) depende mucho de los pocos supervivientes largos, porque la distribución de los
  días tiene una cola larga a la derecha. Lo que se aprende acaba dominado por ellos.
- El modelo de Cox aprende justamente el orden de riesgo. A partir de ese riesgo se obtienen las tres clases
  con umbrales: los cuantiles del riesgo en el fold de entrenamiento se ajustan a las proporciones de cada
  clase en ese mismo fold. Así se obtiene la métrica de BraTS sin entrenar para ella.

**Con qué métrica: c-index como principal y exactitud de 3 clases como secundaria (para comparar con BraTS).**
Las dos métricas se calcularon sobre el mismo modelo con la misma señal [S]:

| Pacientes evaluados | c-index: (media − 0.5) / DE | Exactitud 3 clases: (media − azar) / DE |
|---|---|---|
| 33 (un fold externo) | 2.0 | 1.2 |
| 163 (todos) | 4.6 | 2.7 |

Con el mismo modelo, la exactitud es unas **1.7 veces más ruidosa** que el c-index respecto de su nivel de
azar (0.345 con proporciones 40/25/35). Si se usa como métrica principal, una diferencia real entre el
modelo con genético y el modelo sin él necesitará casi el triple de pacientes para verse.

Recomendación para DECISIONES:

- Primaria: c-index de Harrell, calculado fuera de fold.
- Secundarias: exactitud de 3 clases, ρ de Spearman y MSE en días. BraTS 2018 reportaba las tres, así que
  todas son comparables con la literatura.
- **CHECK** los cortes oficiales. Se suelen citar como < 10, 10–15 y > 15 meses (unos 300 y 450 días), pero
  no pude verificarlos en la página de evaluación de BraTS 2018. Confirmarlos antes de reportar. Calcular
  también las proporciones reales de cada clase en `clinica.csv`; las de la simulación eran supuestas.
- En BraTS 2018 todos los pacientes tienen el evento (muerte), sin censura. Por eso Cox se ajusta con
  `evento = True` para todos, y el c-index equivale a la concordancia entre el riesgo y los días.

## 2. Papel del genético: solo selección de características

Se comparó, con los mismos folds [S], un genético que solo selecciona características con otro que además
evoluciona un gen real: `log10(alpha)`, la penalización ridge del Cox interno, en el rango [−2, 2].

| Escenario | Solo características | Características + penalización | Cox Elastic Net |
|---|---|---|---|
| Sin señal radiómica | 0.563 | 0.572 | 0.607 |
| Con señal radiómica | 0.598 | 0.596 | 0.633 |

(c-index fuera de fold, media de 4 repeticiones × 5 folds.)

- **No mejora.** Las diferencias (±0.01) están muy por debajo de la variación entre repeticiones (DE de la
  diferencia pareada entre 0.016 y 0.031).
- **El hiperparámetro no queda identificado.** El `log10(alpha)` elegido varía con una DE de 1.2 entre folds,
  más de un orden de magnitud. El genético lo usa como un grado de libertad más para ajustar el ruido del
  fold interno, no como un parámetro que estime bien.
- Cada gen añadido amplía el espacio que el genético puede sobreajustar. El objetivo del proyecto es medir si
  el genético selecciona mejor que una penalización estándar. Si también cambia los hiperparámetros, una
  diferencia ya no se puede atribuir a la selección.

Recomendación:

- El genético **solo** selecciona el subconjunto. El modelo que evalúa cada subconjunto usa un Cox ridge con
  `alpha` fijo (por ejemplo 1.0, decidido antes de ver los resultados).
- El modelo final del fold se ajusta con `alpha` elegido por rejilla en el fold interno.
- Si el reto exige explícitamente optimizar hiperparámetros con evolución, que sea un experimento aparte y
  declarado. El candidato natural son los hiperparámetros de la segmentación (nnU-Net), donde sí hay millones
  de vóxeles, no el pronóstico con 163 pacientes.

## 3. Método tradicional: Cox con Elastic Net

| Opción | Problema con n = 163 y p = 1146 |
|---|---|
| **Cox Elastic Net** (`CoxnetSurvivalAnalysis`, `l1_ratio` 0.5, edad sin penalizar) | Selecciona por sí mismo; es el rival justo del genético. Con señal real fue el mejor de los cinco brazos (0.633) [S] |
| LASSO puro (`l1_ratio` = 1) | Con características muy correlacionadas (las radiómicas lo son) elige una al azar de cada bloque, lo que lo vuelve inestable. Elastic Net reparte el peso dentro del bloque |
| Filtro estadístico + SVM | Obliga a clasificar 3 clases (pregunta 1) o a usar SVM de supervivencia. El filtro univariado ignora la correlación entre características, y hay que tunear C y gamma: más hiperparámetros con pocos pacientes |
| Random Forest / RSF | Muchos hiperparámetros, y se degrada con muchas características de ruido. Usarlo dentro del genético costaría unas 100 veces más cómputo que Cox |

Resultados de los cinco brazos [S] (c-index fuera de fold; diferencia pareada con Elastic Net):

| Brazo | Sin señal | Con señal | Δ frente a Elastic Net (con señal) |
|---|---|---|---|
| Edad sola | 0.621 | 0.618 | −0.014 |
| Cox Elastic Net | 0.607 | **0.633** | — |
| Genético: características | 0.563 | 0.598 | −0.035 |
| Genético: características + penalización | 0.572 | 0.596 | −0.037 |
| Genético: aptitud = exactitud | 0.551 | 0.602 | −0.031 |

Lo que muestran los resultados:

- **Sin señal radiómica**, Elastic Net cae poco por debajo de la edad (0.607 frente a 0.621), mientras que el
  genético cae claramente (0.55–0.57).
- **Con señal**, Elastic Net es el único brazo que supera a la edad.
- **El genético sobreajusta su propia aptitud.** Su aptitud interna (0.73–0.74) supera en unos 0.14 a su
  c-index fuera de fold, y los subconjuntos casi no coinciden entre folds (Jaccard < 0.01).
- **Usar la exactitud como aptitud no ayuda.** Es lo que pasaría si el objetivo fueran las 3 clases. Con señal
  la exactitud fuera de fold sube un poco (0.448 frente a 0.429), pero dentro del ruido, y sin señal es el
  peor brazo. Esto refuerza la respuesta 1: la aptitud del genético también debe ser el c-index.

El resultado honesto que cabe esperar con los datos reales es que el genético quede igual o por debajo de
Elastic Net. Es un resultado válido si se reportan la brecha entre aptitud interna y desempeño fuera de fold,
el número de características y la estabilidad de la selección. Las medidas que le dan una oportunidad real al
genético están en `DECISIONS_RECOMMENDATIONS.md`, R5: reducir el espacio de búsqueda dentro de cada fold,
limitar el subconjunto a unas 10 características, y aptitud con CV interna repetida.

Tres brazos obligatorios, sobre los mismos folds: **edad sola**, **Cox Elastic Net** y **genético + Cox**.

## 4. Implementación: DEAP

- DEAP 1.4.4 (abril de 2026) es Python puro (`py3-none-any`). Se instala sin problemas en el entorno
  `radiomica` (Python 3.11) y no toca numpy 1.26.
- Los operadores ya están probados: `tools.selTournament`, `tools.cxUniform`, `tools.mutFlipBit` y
  `algorithms.eaMuPlusLambda`. Un evaluador reconoce DEAP de inmediato; un genético propio obliga a demostrar
  que no tiene errores, y un error en la selección o en el elitismo es difícil de ver en los resultados.
- Un genético propio sí es factible: el de la simulación ocupa unas 40 líneas de numpy. Pero la ventaja no
  compensa el riesgo, salvo que el equipo quiera presentar el operador como contribución propia.

Precauciones con DEAP:

1. **Semillas.** DEAP usa el módulo `random` de Python, no el de numpy. Hay que fijar `random.seed(s)` y
   `np.random.seed(s)` al inicio de cada fold externo, y guardar la semilla en el log.
2. **`creator` es global.** Hay que definir `creator.create("FitnessMax", ...)` y `creator.create("Individual", ...)`
   en el nivel superior del módulo, no dentro de una función. Si no, la ejecución en paralelo
   (joblib/multiprocessing) falla al serializar los individuos.
3. **Caché de aptitud** por individuo (por ejemplo con `tuple(individuo)` como clave). Las élites y los
   duplicados se reevalúan en cada generación, y con 1146 genes es habitual que se repitan.
4. **Inicialización dispersa.** `tools.initRepeat` con probabilidad 0.5 por gen activa unas 573
   características. Hay que inicializar con probabilidad k₀/p (por ejemplo 5/1146) y penalizar o descartar los
   subconjuntos con más de 15 activas. Sin esto el genético arranca en una zona inválida.
5. **Registrar por fold externo:** aptitud interna, c-index fuera de fold, número de características,
   subconjunto elegido y curva de aptitud por generación (`tools.Logbook`).

---

## Conexión con el README

- **Folds externos = los 5 folds fijos de `particiones/folds.csv`** (D22), unos 33 pacientes con supervivencia
  en cada uno. Como no se pueden regenerar, las repeticiones se hacen sobre la CV interna y las semillas del
  genético, no sobre la partición externa. La comparación entre brazos se reporta como diferencia pareada por
  fold.
- **Características:** en el análisis principal, las de máscaras predichas (`resultados/pred/features.csv`),
  porque es el uso real descrito en el README. En el secundario, las de máscaras manuales. La diferencia entre
  ambos mide cuánto cuesta el error de segmentación (complementa D26).
- **Resección:** el README dice que solo se reporta en CBICA y en 2013. Como covariable, "no reportado"
  equivale casi a "centro distinto de CBICA/2013", así que mezclaría efecto del centro con efecto de la cirugía.
  Hay dos opciones: dejarla fuera del modelo principal, o incluirla y declarar esa confusión.

## Apéndice: simulaciones

- **Ruido de las métricas:** un predictor con c verdadero ≈ 0.61, tiempos de Weibull sin censura,
  3000 réplicas por n. Las clases salen de los cuantiles del tiempo, y la clase predicha de los cuantiles del
  riesgo en un conjunto de entrenamiento independiente.
- **Brazos:** 5 folds externos estratificados × 4 repeticiones × 2 escenarios. Genético con población de 40,
  30 generaciones, torneo de 3, cruce uniforme, mutación 1/p, elitismo de 2 y k ≤ 15. Aptitud = c-index (o
  exactitud) en CV interna de 3 folds de un Cox ridge sobre la edad y el subconjunto, menos 0.002·k.
- **Limitaciones:** datos sintéticos; solo 4 repeticiones, suficiente para ver la dirección pero no la tercera
  cifra decimal; las proporciones 40/25/35 de las clases son supuestas.
