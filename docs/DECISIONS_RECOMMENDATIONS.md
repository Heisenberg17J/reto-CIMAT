# Recommendations for the open decisions (CIMAT / BraTS 2018)

Companion to `DECISIONS.md` (version received 2026-10-02). Every recommendation names the entry it closes,
the evidence behind it, and what to write in `DECISIONS.md` if the team accepts it. Evidence comes from
three sources:

- **[D]** facts already recorded in `DECISIONS.md` (counts from the team's own scripts on the full 285);
- **[S]** two simulations run for this note (`sim_*.csv`, `fig_open_decisions_sim.png`; method in the appendix);
- **[L]** PyRadiomics 3.x source code and the literature collected in `PYRADIOMICS_BRATS2018_NOTES.md`.

Nothing here has been run on the real images or the real `clinica.csv`. Points that need checking on the real
data are marked **CHECK**.

---

## Summary

| # | Open decision | Recommendation | Strength |
|---|---|---|---|
| R1 | HGG only vs both grades (survival) | **HGG only, forced by the data.** Grade cannot be a covariate. | Settled by data |
| R2 | Survival cohort: all 163 vs GTR 59 | **All 163 as primary**, resection status (GTR/STR/NA) as covariate; GTR-only as secondary | Strong |
| R3 | Evaluation protocol | **Repeated nested CV as primary estimate**; no small single hold-out | Strong |
| R4 | Baseline model family | **Penalised Cox (elastic net) with age forced in**; RSF as secondary only | Strong |
| R5 | GA design (follows from R3/R4) | Cheap ridge-Cox wrapper, parsimony cap of about 10 features, reduced search space, report inner vs outer gap | Strong |
| R6 | Wavelet/LoG (D10) | **Original only.** Filters become a separate, pre-declared experiment, if run at all | Strong |
| R7 | binWidth (D9 sensitivity) | Keep 0.1; check bin counts on full 285; sensitivity at {0.05, 0.1, 0.25} judged on stability, not c-index | Moderate |
| R8 | Empty ET / ROI thresholds (D11) | Survival: no empty ET expected in HGG, so the decision is moot (**CHECK** small-ET counts). Classification: "no ET" indicator, never imputation | Strong |
| R9 | Tumour voxels outside the brain (D1) | **Intersect ROI masks with the brain mask.** This is a real bug, not cosmetic (see below) | Strong |
| R10 | CV grouping by site (D21) | Survival: stratify by survival tertile; leave-site-out as sensitivity | Moderate |
| R11 | 2D vs 3D segmentation | **3D patch-based U-Net** (needs a GPU); survival features from expert masks first | Moderate |
| R12 | Metric | Harrell's c-index primary + BraTS 3-class accuracy for comparison with the literature | Strong |

---

## R1. HGG only for survival (closes the grade question for task B)

`DECISIONS.md` D21 records that survival exists for **163 patients, all HGG**. All 75 LGG lack survival.
The worry in the project context, a model that just detects grade, cannot arise in task B, and including
"grade as covariate" is impossible because the covariate would be constant.

- **Task B (survival):** HGG only, n = 163.
- **Task A (segmentation):** all 285. Grade is irrelevant to voxel labels, and LGG adds tumour-appearance
  diversity.
- If the team adds an HGG/LGG classification side-task, the grade-confounding risks from D11 (empty ET) and
  D21 (site confounded with grade) return. Treat that side-task under R8 and R10.

**Write in DECISIONS.md:** "Survival cohort restricted to HGG because no LGG patient has survival data
(D21). The decision follows from the data, not from a methodological choice."

## R2. Primary cohort = all 163, resection as covariate

BraTS 2018 ranked survival only on GTR patients (59 in the training data per D21). That is a ranking rule
of the challenge, not a methodological requirement. The cost of following it is precision:

| Patients evaluated | 95% half-width of c-index [S] |
|---|---|
| 163 (all with survival) | ±0.049 |
| 83 (GTR + STR) | ±0.070 |
| 59 (GTR only) | ±0.084 |

With 59 patients the interval is roughly ±0.08, as wide as the gap between a useless model (0.50) and a
good published one (0.58 to 0.62). Nothing would be distinguishable.

- **Primary analysis:** n = 163, resection status as a 3-level categorical covariate (GTR / STR / not
  reported). It must be one-hot encoded, and "not reported" (80 patients) kept as its own level.
  Dropping those 80 patients would lose half the cohort, and imputing them would invent data.
- **Secondary analysis:** the GTR-only subgroup (59), reported for comparability with the BraTS ranking,
  explicitly labelled as underpowered.
- Resection status is clinical, not radiomic, so it belongs to the "clinical/demographic" input the
  challenge allows, next to age.

## R3. Repeated nested CV, not a small hold-out

The project context plans "a held-out test set untouched until the end". With n = 163 a 20% hold-out is
33 patients, and its c-index has a 95% half-width of **±0.11** [S]; with GTR only it is 12 patients and
**±0.21**. A single hold-out would be a coin toss that cannot confirm or reject anything, and its result
would depend mostly on which patients happened to land in it.

Recommendation:

- **Primary estimate:** outer 5-fold CV, **repeated 10 times** with different seeds, stratified by survival
  tertile. Every patient is scored as a test case in each repeat. Report the mean over repeats and the
  spread across repeats.
- **Inner loop** (inside each outer training fold): all fitting, including feature selection, the GA,
  scaling and penalty tuning. This is exactly what D19 already mandates.
- **Fair comparison:** all arms (age-only, elastic-net Cox, GA + Cox) use the **same outer splits**. The
  comparison is then paired: compute `c(GA) − c(baseline)` per outer fold and repeat, and report that
  difference with its spread. Paired differences are more precise than two separate intervals.
- **If the challenge requires a hold-out:** keep it, stratified, and report it with its interval and the
  sentence "this hold-out cannot discriminate differences below about 0.1".

Note that BraTS 2018 survival data contain only deaths (no censoring), so Harrell's c here is equivalent to
Kendall-type rank concordance between predicted risk and survival time. This does not change the plan,
but it should appear in the methods section.

## R4. Baseline model: elastic-net Cox with age unpenalised

Arguments, ordered by weight:

1. **The baseline must be able to handle p ≫ n on its own.** 1,146 features versus 163 patients rules out
   unpenalised Cox PH. An elastic-net Cox (scikit-survival `CoxnetSurvivalAnalysis`) performs its own
   embedded selection. That makes it the correct "without evolution" control: the question becomes
   whether a GA selects better than a standard penalty does, rather than whether a GA beats no selection
   at all, which would be an unfair, easy win.
2. **Age is forced in** (`penalty_factor = 0` for age, and the same for resection dummies). Age is the
   benchmark the literature says radiomics rarely beat [L: Weninger 2019]. Forcing it in means radiomics
   are always evaluated as added value over age.
3. **Random survival forest is the weaker choice as baseline.** RSF has many hyperparameters and degrades
   with many noise features. With 163 events its advantage (non-linearity) is unlikely to be detectable.
   A GA wrapper around RSF would also cost about 100× the compute of a Cox wrapper. Keep RSF as a
   secondary sensitivity model fitted on the final GA subset, if time allows.
4. **Third arm, mandatory: age only** (plus resection). It costs nothing and is the honest floor.

Simulation result [S], with n = 163, 1,146 correlated features, survival driven by age and optionally by a
real radiomic signal:

| Scenario | Age only | Elastic-net Cox | GA + Cox (outer test) | GA inner fitness |
|---|---|---|---|---|
| Radiomics carry **no** signal | 0.618 | 0.607 | 0.549 | 0.718 |
| Radiomics carry a **real** signal | 0.601 | 0.590 | 0.579 | 0.727 |

(Mean c-index over 4 repeats × 5 outer folds.)

The elastic net degrades gracefully: when radiomics are noise it keeps few features and stays close to age.
The GA, given 1,146 candidates, adds noise features and falls **below** age. This is the overfitting the
challenge brief warns about, now quantified.

## R5. GA design so that it has a chance (follows from R3 and R4)

The simulation shows three failure signatures that the real study should measure and report:

- **Optimism:** inner-CV fitness exceeds the outer test c-index by **0.15 to 0.17** [S]. A GA that reports
  its own fitness (about 0.72) would look like a success while being worse than age. That is the
  "above about 0.75 signals an error" rule from the project context, observed directly.
- **Instability:** the mean Jaccard overlap of selected subsets across outer folds is about **0.007** [S].
  Subsets are practically disjoint, so the GA is fitting chance.
- **Parsimony alone is not enough.** The GA kept about 10 features versus 6 to 54 for elastic net, but
  fewer features did not make it generalise better.

Design recommendations:

1. **Shrink the search space before the GA, inside each fold.** Searching subsets of 1,146 features is
   about 10²⁴ subsets of size ≤ 10. Options, in order of preference:
   - (a) a **pre-declared reduced family** chosen by domain knowledge, not by data: for example shape (WT,
     TC, ET) plus first-order of t1ce and FLAIR in TC/ET, about 100 to 150 features;
   - (b) **unsupervised correlation clustering** fitted on the fold's training data (keep one
     representative per cluster with |ρ| > 0.9). This uses no outcome, so it leaks nothing, but it must
     still be fitted per fold;
   - (c) a univariate pre-filter on the training fold. This is the weakest option, because it already
     uses the outcome.
2. **Cap subset size** at about 10, roughly 163 events / 15 events per variable. Use a hard limit plus a
   small parsimony penalty (`fitness = mean inner c − λ·k`).
3. **Inner model:** ridge-penalised Cox (`alpha ≈ 1`) on the subset, with age and resection always
   included. It is fast (one GA run of 40 × 30 took about 35 s in the simulation on one CPU core) and
   numerically stable with correlated features.
4. **Repeated inner CV for fitness** (e.g. 3-fold × 2) reduces the noise the GA is optimising.
5. **Report per outer fold:** inner fitness, outer c, the gap, k, and selected features. Add a stability
   statistic across folds (Jaccard, or the Nogueira stability index). These are the measures the project
   context already promises ("parsimony and selection stability").
6. **Optional stronger variant:** "stability-GA". Run the GA on B bootstrap samples of the training fold,
   then keep features selected in ≥ 60% of runs. This trades compute for stability.

**Expected honest result:** GA + Cox ≈ elastic-net Cox ≈ age, within ±0.05. That is a publishable,
correct finding for this sample size, provided the inner/outer gap and stability are reported.

## R6. Wavelet/LoG: Original only (closes D10)

- The GA already overfits by about 0.15 with 1,146 features [S]. Wavelet multiplies the feature count to
  9,978 and Wavelet + LoG to 15,498 [L, computed from source], making the search space larger by orders
  of magnitude.
- In the literature, the high-dimensional filtered approach (Original + Wavelet + LoG, all modalities)
  did not beat age in GTR patients [L: Weninger 2019].
- Filtered features are also the least reproducible across scanners. Reproducibility studies in the notes
  report this direction, but it was not tested here.

**Write in DECISIONS.md:** "First and primary pass: Original only. LoG/Wavelet, if run, are a secondary
experiment declared before seeing its results, with its own nested CV and the same reduced-family rule
(R5.1)."

## R7. binWidth sensitivity (D9 is settled; this is how to run the planned sensitivity)

D9 fixed 0.1 on the basis of 4 patients (16 to 128 bins, median 55).

1. **CHECK on the full 285:** compute `firstorder_Range / 0.1` per modality × region. Accept 0.1 if
   about 90% of ROIs fall within 16 to 130 bins. ET in HGG is the likely problem (small and bright).
2. **Sensitivity grid:** {0.05, 0.1, 0.25}. Variant B already showed 0.25 gives a median of about 22 bins
   (min 6), which is the coarse end. Carré 2020 suggests about 32 bins as a compromise, which on the
   z-scale corresponds to about 0.1 to 0.2 for typical ranges.
3. **Judge on feature stability, not on survival c-index:** Spearman rank correlation of each feature
   between bin widths across patients, and the fraction of features with ρ > 0.9. Choosing the bin width
   by c-index would make discretisation part of model selection, so it would need to be nested as well.
4. Only bin-dependent families change (GLCM, GLRLM, GLSZM, GLDM, NGTDM, plus entropy/uniformity-type
   first-order features). Shape and most first-order features do not change, so the sensitivity covers
   less than half the columns.

## R8. Empty ET and ROI thresholds (D11)

- **Survival (task B):** D1 found the 27 empty-ET patients are all LGG, so **no survival patient has an
  empty ET**. The empty-ET decision does not affect task B.
- **CHECK:** small (not empty) ET in HGG. Count HGG patients with an ET volume below 27, 64 and 125 voxels.
  `manifest_regiones.csv` (D4) already flags `pequena`. If the count is 0 to 2, keep 27 and move on. If it
  is larger, texture features of tiny ET are unstable and **minimumROISize = 64** (a 4×4×4 cube) is a
  better floor.
- **NaN handling inside the survival model:** Cox models and DEAP fitness functions do not accept NaN.
  If a few HGG patients have a flagged ET, use within-fold median imputation for their ET features plus an
  indicator column "ET small/absent". Never impute outside the fold (D19).
- **If an HGG/LGG classifier is built:** do **not** impute ET features for LGG. A missing ET is itself the
  strongest grade signal (36% of LGG versus 0% of HGG). Use a binary `has_ET` feature, ET volume = 0, and
  exclude ET texture features from that model. Report that `has_ET` alone gives a non-trivial baseline,
  so radiomics must beat it.

## R9. Tumour voxels outside the brain (D1): intersect masks with the brain mask

D1 states that "after normalisation those voxels are 0 inside the mask". With z-scoring, **0 is not "no
signal"; it is the brain mean intensity.** Those voxels therefore enter first-order and texture
computations as fake mean-intensity voxels. In the worst patient (Brats18_TCIA13_634_1, 588 voxels) this
can shift statistics noticeably, and for the median case (10 voxels) the effect is negligible but not zero.

- **Fix:** `mask_region ∩ brain_mask`, where `brain_mask = (t1>0) | (t1ce>0) | (t2>0) | (flair>0)`,
  computed before normalisation. Apply it in `regiones.py` to all three regions, so that shape and
  intensity use the same voxels (D6 consistency). Log the number of removed voxels per patient.
- **Alternative (if the team prefers untouched expert masks):** `resegmentRange` cannot fix this, because
  0 lies inside the legitimate z-range. Only the mask intersection works.
- Check 1 in D16 compares `VoxelVolume` with the Block 4 count, so it will still pass if the intersection
  is done in Block 4.

## R10. CV grouping by site (D21)

- For **survival (HGG only)**, the site-confounding issue in D21 (site ↔ grade) disappears, because grade
  is constant. Site can still shift intensities and survival (different treatment eras/centres).
- **Primary:** outer folds stratified by survival tertile (R3).
- **Sensitivity:** leave-one-site-group-out (group small TCIA sites together so each test group has ≥ 15
  patients). A drop of more than about 0.05 versus the stratified CV indicates site leakage.
- **CHECK:** tabulate the 163 survival patients by `origen` before deciding the grouping.
- For an HGG/LGG classifier, site **must** be grouped, because every site except 2013 holds a single
  grade. Otherwise the classifier can learn scanner signatures instead of grade.

## R11. 2D vs 3D segmentation

- **3D** is the recommended choice: the data are isotropic 1 mm and co-registered, so 3D context is free.
  BraTS winners from 2018 onward have been 3D networks (Myronenko 2018, nnU-Net 2018–2020). A 2D network
  loses through-plane continuity and usually scores lower Dice on TC/ET.
- **Practical constraint:** 3D U-Net training on 128³ patches needs a GPU with about 8 to 12 GB of memory.
  It cannot be trained on a CPU-only laptop. If no GPU is available, a 2.5D approach (2D network with
  neighbouring slices as channels) is the fallback, documented as compute-driven.
- **Chaining tasks A → B:** for the survival analysis, extract features from **expert masks** first. That
  isolates the prognostic question from segmentation error. Then repeat with **predicted masks** as a
  robustness check. That is what a real deployment would see, and the gap between the two is a result in
  itself.
- A GA can legitimately touch task A too, through hyperparameter search for the network (learning rate,
  loss weights, patch size). That is optional and expensive; the survival GA is the core differentiator.

## R12. Metric

- **Primary:** Harrell's c-index from repeated nested CV, with a paired difference between arms (R3).
- **Secondary, for literature comparability:** BraTS 3-class accuracy (short < 10 months, mid 10–15
  months, long > 15 months), obtained by thresholding predicted survival. 2018 leaders reported about
  0.6 [L]. **CHECK** the exact day cut-offs against the BraTS 2018 summary paper before reporting.
- Also report Spearman ρ and MSE in days, which BraTS also used, so all three numbers are comparable.

---

## Appendix: what the simulations did (and did not) show

**Simulation A, c-index precision.** Cox-type (Weibull) survival times with one predictor of true
c ≈ 0.61, no censoring (as in BraTS 2018), 4,000 replicates per n. The 95% half-width is 1.96 × the SD of
Harrell's c across replicates. Output: `sim_cindex_precision.csv`.

**Simulation B, nested CV.** n = 163, p = 1,146 (same as the real Original-only table), 40 latent factors,
each feature loading 0.6 to 0.95 on one factor, giving blocks of correlated features like real radiomics.
Survival = Weibull with log-hazard 0.45·age (+ 0.35·(F₁ + 0.7·F₂) in the "real signal" scenario).
Outer 5-fold stratified CV × 4 repeats. Three arms on identical splits:

- age-only;
- elastic-net Cox (`l1_ratio` = 0.5, α tuned by inner 3-fold CV, age unpenalised);
- GA (population 40, 30 generations, tournament 3, uniform crossover, bit-flip mutation 1/p, elitism 2;
  fitness = inner 3-fold c of ridge-Cox on age + subset − 0.002·k; k ≤ 15).

Outputs: `sim_nested_cv_ga_vs_coxnet.csv` and `sim_nested_cv_summary.csv`.

**Limitations.**
- The data are synthetic. Real radiomic correlations and real effect sizes are unknown, and the true
  radiomic signal in BraTS could be weaker or stronger than simulated.
- Only 4 repeats were run (20 outer folds per scenario), which is enough for direction but not for
  third-decimal precision.
- The GA was deliberately small. A larger GA (more generations, larger population) would search harder
  and, at this n, would be expected to **overfit more**, not less.
- The direction of every conclusion (large inner/outer gap, near-zero stability, no gain over age at
  n = 163) agrees with the BraTS literature [L], which is the reason they are offered as evidence. They
  are not a forecast of the team's exact numbers.
