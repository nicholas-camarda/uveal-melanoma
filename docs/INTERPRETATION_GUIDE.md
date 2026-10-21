# Interpretation Guide

This guide is for reading the generated outputs after the pipeline has run. It explains how to navigate the folders, read the workbooks, and interpret the tables and figures in plain language.

It is not the canonical source for variable derivations or statistical assumptions:

- Use [CALCULATIONS.md](CALCULATIONS.md) for how values and endpoints are constructed.
- Use [STATISTICAL_METHODS.md](STATISTICAL_METHODS.md) for formal model definitions, thresholds, and validation metrics.
- Use [TECHNICAL.md](TECHNICAL.md) for workflow internals and artifact contracts.

---

## Table of Contents

- [Quick Start](#quick-start)
- [Reading Summary Tables](#reading-summary-tables)
- [Understanding Regression Outputs](#understanding-regression-outputs)
- [Reading the Restricted-Cohort Propensity-Overlap Sensitivity](#reading-the-restricted-cohort-propensity-overlap-sensitivity)
- [Interpreting Survival Curves](#interpreting-survival-curves)
- [Reading Forest Plots](#reading-forest-plots)
- [Understanding RMST Analysis](#understanding-rmst-analysis)
- [Interpreting Proportional Hazards Diagnostics](#interpreting-proportional-hazards-diagnostics)
- [Understanding GEP Analysis](#understanding-gep-analysis)
- [Clinical Significance vs Statistical Significance](#clinical-significance-vs-statistical-significance)
- [Common Pitfalls](#common-pitfalls)
- [Troubleshooting](#troubleshooting)

---

## Quick Start

If you only need the shortest route:

1. Find the runtime analysis folder described in the [README path model](../README.md#path-model).
2. Open the objective folder that matches your question.
3. Start with the consolidated `.xlsx` workbook or summary `.html` artifact in that folder.
4. Use this guide to interpret what you are seeing.
5. If you need to know how a value was calculated, move to [CALCULATIONS.md](CALCULATIONS.md). If you need the formal method, move to [STATISTICAL_METHODS.md](STATISTICAL_METHODS.md).

### Finding Your Analysis

**Step 1:** Identify your cohort
- Full cohort: `uveal_full/`
- Restricted cohort: `uveal_restricted/`
- GKSRS-only cohort: `gksrs/`

Current cohort totals are exported with each run to `Analytic Dataset/cohort_summary_statistics.json`.

**Step 2:** Navigate to objective folder
- Objective 0: `00_General/` (baseline characteristics, patient flow)
- Objective 1: `01_Efficacy/` (survival, recurrence, tumor height)
- Objective 2: `02_Safety/` (vision changes, complications)
- Objective 3: `03_Repeat_Radiation/` (PFS-2 for recurrence patients)
- Objective 4: `04_GEP_Validation/` (GEP validation)

**Step 3:** Open relevant files
- `.xlsx` - Summary statistics, regression tables
- `.html` - Interactive regression tables with details
- `.png` - Survival curves, forest plots, RMST visualizations

### Key Files in Each Analysis

| File Pattern | Content | When to Use |
|--------------|---------|-------------|
| `*_summary.xlsx` | Event rates, means, medians | Quick overview of raw data |
| `*_cox_coxph.html` | Cox regression models | Adjusted treatment effects for survival |
| `*_logistic_glm.html` | Logistic regression | Adjusted odds ratios for binary outcomes |
| `*_linear_lm.html` | Linear regression | Mean differences for continuous outcomes |
| `*_km.png` | Kaplan-Meier plots | Visual survival comparison |
| `*_rmst_pvalue_progression.png` | RMST p-value progression plots | Alternative time-horizon view when PH is weak or RMST is the focus |
| `*_subgroup_forest_plot.png` | Subgroup analysis | Treatment effect heterogeneity |

---

## Reading Summary Tables

### Baseline Characteristics Tables

**File:** `00_General/baseline_characteristics/baseline_characteristics.xlsx`

**Structure:**
- Rows: Patient characteristics (age, sex, tumor features)
- Columns: Treatment groups (PBT, GKSRS) + Overall
- Values: n (%), median [IQR], mean (SD)

**How to Read:**

```
Characteristic          PBT (n=100)    GKSRS (n=60)    p-value
Age at diagnosis        62 [52, 71]    65 [55, 73]     0.23
Sex                                                     0.45
  Male                  55 (55%)       30 (50%)
  Female                45 (45%)       30 (50%)
Initial tumor height    5.2 (2.1)      5.5 (2.3)       0.41
```

**Interpretation:**
- **Age:** Median PBT patient was 62 years old; GKSRS was 65 years
- **Sex:** Similar gender distribution (p=0.45, not significantly different)
- **Tumor height:** Mean PBT height was 5.2 mm; GKSRS was 5.5 mm (p=0.41, similar)
- **P-values:** Test whether groups differ at baseline (ideally p > 0.05, indicating balance)

**Clinical Insight:** Groups well-balanced at baseline → fair comparison

### Event Rate Tables

**File:** `{objective}/{sub_objective}/*_summary.xlsx`

**Example:**
```
Outcome              PBT           GKSRS         p-value
Local recurrence     15/100 (15%)  12/60 (20%)   0.41
```

**Interpretation:**
- 15% of PBT patients experienced recurrence
- 20% of GKSRS patients experienced recurrence
- Difference not statistically significant (p=0.41)

**Clinical Insight:** Crude rates suggest similar local control

---

## Understanding Regression Outputs

### Cox Regression Tables (Survival Outcomes)

**File:** `*_cox_coxph.html`

**Example Output:**

```
Characteristic          HR      95% CI          p-value
Treatment
  PBT                   —       —               
  GKSRS                 0.65    0.45, 0.95      0.025
Age (per year)          1.03    1.01, 1.05      0.008
Sex
  Male                  —       —
  Female                0.82    0.58, 1.15      0.25
```

**How to Read:**

**Hazard Ratio (HR):**
- HR = 1.0: No difference
- HR < 1.0: Lower hazard (protective)
- HR > 1.0: Higher hazard (harmful)

**GKSRS HR = 0.65:**
- GKSRS patients have 65% the hazard of PBT patients
- Equivalently: 35% reduction in hazard (1 - 0.65 = 0.35)
- **Clinical meaning:** GKSRS associated with lower risk of event

**Age HR = 1.03:**
- Each additional year of age increases hazard by 3%
- 10-year difference → 1.03^10 = 1.34 (34% higher hazard)

**95% CI = (0.45, 0.95):**
- We're 95% confident true HR is between 0.45 and 0.95
- Does not include 1.0 → statistically significant
- Width indicates precision (narrow = more precise)

**P-value = 0.025:**
- Probability of observing this difference if no true difference exists
- p < 0.05 → statistically significant
- Likelihood ratio test preferred over Wald test

**Reference Categories (—):**
- Other groups compared to this baseline
- PBT is reference for treatment
- Male is reference for sex

**Clinical Interpretation:**
After adjusting for age, sex, location, and optic nerve involvement, GKSRS treatment is associated with a 35% reduction in hazard compared to PBT (HR=0.65, 95% CI: 0.45-0.95, p=0.025).

### Reading the Restricted-Cohort Propensity-Overlap Sensitivity

This supplemental analysis is available only in the restricted cohort's `01_Efficacy/h_propensity_score_sensitivity/` folder. Start with `restricted_cohort_propensity_overlap_summary.md`, then open `restricted_cohort_propensity_overlap_sensitivity.xlsx` for the prespecified model, included population, covariate balance, weighted Cox estimates, and PH diagnostics.

The companion plots show propensity-score overlap, covariate balance before and after weighting, the four endpoint hazard ratios, and the weighted Schoenfeld residual diagnostics. The latter labels treatment and global PH-test p-values. Treat these artifacts as a sensitivity analysis for the average treatment effect in the overlap population—not as a replacement for the primary adjusted Cox models, and not as evidence that unmeasured confounding has been eliminated.

### Logistic Regression Tables (Binary Outcomes)

**File:** `*_logistic_glm.html`

**Example Output:**

```
Characteristic          OR      95% CI          p-value
Treatment
  PBT                   —       —               
  GKSRS                 1.45    0.85, 2.48      0.17
Age (per year)          1.02    0.98, 1.06      0.35
```

**How to Read:**

**Odds Ratio (OR):**
- OR = 1.0: No difference in odds
- OR < 1.0: Lower odds
- OR > 1.0: Higher odds

**GKSRS OR = 1.45:**
- GKSRS patients have 1.45 times the odds of the outcome
- 45% higher odds compared to PBT
- **Clinical meaning:** GKSRS associated with higher odds of complication

**95% CI = (0.85, 2.48):**
- Includes 1.0 → not statistically significant
- Wide interval → imprecise estimate (small sample or few events)

**Clinical Interpretation:**
GKSRS treatment is associated with 45% higher odds of radiation retinopathy compared to PBT, but this difference is not statistically significant (OR=1.45, 95% CI: 0.85-2.48, p=0.17).

### Linear Regression Tables (Continuous Outcomes)

**File:** `*_linear_lm.html`

**Example Output:**

```
Characteristic          Beta    95% CI          p-value
Treatment
  PBT                   —       —               
  GKSRS                 -0.35   -0.68, -0.02    0.038
Baseline tumor height   -0.42   -0.55, -0.29    <0.001
```

**How to Read:**

**Regression Coefficient (Beta):**
- Units same as outcome variable
- Beta = 0: No difference
- Beta < 0: Decrease in outcome
- Beta > 0: Increase in outcome

**GKSRS Beta = -0.35 mm:**
- GKSRS patients have 0.35 mm greater decrease in tumor height
- **Model meaning:** Within this example model, GKSRS is associated with greater tumor-height decrease

**Baseline tumor height Beta = -0.42:**
- Each 1 mm increase in baseline height → 0.42 mm greater decrease
- Because baseline height is part of the tumor-height change outcome, this coefficient is partly algebraic and should be treated as a diagnostic sensitivity result

**Clinical Interpretation:**
In this example, GKSRS is associated with 0.35 mm greater tumor-height reduction compared to PBT (β=-0.35, 95% CI: -0.68 to -0.02, p=0.038), after model adjustment. For reviewer-response tumor-height analyses, this kind of model output must be interpreted with the tumor-height timing audit. A model that adds baseline tumor height to a tumor-height change score is not appropriate as the reviewer-facing adjustment because baseline height is already part of the outcome definition. If imaging follow-up differs materially by treatment group, the comparative result should be demoted to descriptive context or omitted from the main response.

For Objective 2 vision there are two different linear outputs:
- `*_logmar_vision_change_adjusted_lm.html` models continuous logMAR change.
- `*_snellen_line_change_adjusted_lm.html` models the exact integer `Snellen Line Change` outcome.

The Snellen linear model is the direct adjusted analogue of the descriptive `Snellen Line Change` row, whereas the ordinal Snellen model uses the separate 7-level `Snellen Line Change Distribution`.

### Ordinal Regression Tables (Ordered Outcomes)

**File:** `*_ordinal_polr.html`

**Example Output:**

```
Characteristic          OR      95% CI          p-value
Treatment
  PBT                   —       —               
  GKSRS                 1.40    1.05, 1.88      0.022
Age (per year)          0.98    0.96, 1.00      0.071
```

**How to Read:**

**Ordinal Odds Ratio (OR):**
- OR = 1.0: No difference in distribution across ordered outcome levels
- OR > 1.0: Higher odds of being in an earlier ordered category
- OR < 1.0: Higher odds of being in a later ordered category

**For Snellen Line Change Distribution in Objective 2a:**
- Categories are ordered from `≥3-line improvement` to `≥3-line loss`
- The stable category is `Stable (0-line change)` after nearest-line rounding with halves away from zero
- OR > 1.0 therefore favors better vision outcomes
- OR < 1.0 therefore favors worse vision outcomes

**GKSRS OR = 1.40:**
- GKSRS patients have 40% higher odds of falling into a better Snellen line-change category than PBT

**Effect Summary Workbooks**

- `*_effect_summary.xlsx` files are flat one-sheet workbooks that combine descriptive, unadjusted, and adjusted rows for the corresponding analysis folder.
- In `a_vision_changes`, use `analysis_label` to distinguish `LogMAR Vision Change`, `Latest LogMAR Vision`, `Snellen Line Change`, and `Snellen Line Change Distribution`.
- `Latest LogMAR Vision` is the reviewer-response sensitivity model for final/latest logMAR VA. It adjusts for baseline VA, explicit latest-VA follow-up duration, viable reviewer-requested baseline predictors, and the shared confounder set. Interpret it as sensitivity evidence only; missing macular/foveal proximity and radiation-dose fields still limit causal claims about visual-acuity preservation.
- `model_status = "DESCRIPTIVE"` marks summary rows, `model_status = "FIT"` marks modeled rows, and `model_status = "SKIPPED"` documents analyses that could not be fit.
- Workbook methods are model-family specific and match the reader-facing HTML tables: linear rows use mean differences with Wald CIs/p-values, logistic rows use ORs with model-based Wald CIs and the pipeline's standard term-level p-values, Cox rows use HRs with native Cox CIs/p-values, and ordinal rows use proportional-odds ORs with 95% Wald CIs plus likelihood-ratio-test p-values.
- **Clinical meaning:** the outcome distribution is shifted toward more improvement or less loss

**95% CI = (1.05, 1.88):**
- Excludes 1.0 → statistically significant
- Entire interval above 1.0 → consistent direction toward better ordered outcomes

**Clinical Interpretation:**
After adjustment, GKSRS is associated with higher odds of a better Snellen line-change category than PBT (OR=1.40, 95% CI: 1.05-1.88, p=0.022).

---

## Interpreting Survival Curves

### Kaplan-Meier Plots

**File:** `*_km.png`

**Key Components:**

**1. Survival Curves (Lines)**
- **Y-axis:** Probability of survival (0-100%)
- **X-axis:** Time since treatment (years)
- **Separate lines:** One per treatment group
- **Curve steps down:** Each time an event occurs

**2. Confidence Bands (Shaded Regions)**
- Show uncertainty around survival estimates
- Wider bands = more uncertainty (fewer patients at risk)
- Narrow bands = more precise (more patients)

**3. Numbers at Risk (Below X-axis)**
```
Time:     0      1      2      3      4      5
PBT:    100     95     88     80     72     65
GKSRS:   60     58     54     50     45     40
```
- Shows patients still under observation at each time
- Decreases due to events and censoring
- Smaller numbers → less reliable estimates

**4. Median Survival (Dashed Lines)**
- Time when 50% of patients have experienced event
- Where survival curve crosses 50% line
- "Not reached" if curve stays above 50%

**How to Read:**

**Example:** 5-year overall survival
- PBT curve at year 5: 75%
- GKSRS curve at year 5: 82%
- Interpretation: 75% of PBT patients survived 5 years; 82% of GKSRS patients survived

**Comparing Curves:**
- **Curves separate early:** Treatment effects emerge quickly
- **Curves cross:** Treatment effects change over time (PH violation)
- **Curves parallel:** Consistent treatment effect (PH holds)
- **Curves overlap:** No meaningful difference

**Log-Rank P-value:**
- Tests whether survival distributions differ overall
- p < 0.05: Statistically significant difference
- Doesn't specify direction or magnitude

**Clinical Interpretation:**
The 5-year overall survival was 75% for PBT and 82% for GKSRS (log-rank p=0.12). While GKSRS shows numerically higher survival, the difference is not statistically significant. Larger sample sizes or longer follow-up may be needed.

### Cumulative Incidence Plots (Competing Risks)

**Similar to K-M plots but:**
- Account for competing events (e.g., death before recurrence)
- More accurate event probabilities
- Gray's test instead of log-rank test

Objective 4 draws one cumulative-incidence curve per Class x PRAME group; see [Understanding GEP Analysis](#understanding-gep-analysis).

---

## Understanding RMST Analysis

### What is RMST?

**Restricted Mean Survival Time:** Average time survived up to a specified time point (τ)

**Why Use RMST?**
- Doesn't require proportional hazards assumption
- Clinically interpretable ("months gained/lost")
- Robust alternative when PH violated

### RMST Difference Plots

**File:** `*_rmst_pvalue_progression.png`

**Key Components:**

**1. Survival Curves with Shaded Areas**
- **Shaded regions:** Area under survival curve = RMST
- **Larger area:** Longer average survival
- **Compare shading:** Visual difference in survival time

**2. RMST Estimates (Text)**
```
PBT RMST:    4.2 years (95% CI: 3.9-4.5)
GKSRS RMST:  4.6 years (95% CI: 4.2-4.9)
Difference:  0.4 years (95% CI: 0.1-0.7)
P-value:     0.015
```

**How to Interpret:**

**RMST Values:**
- PBT patients survived an average of 4.2 years up to 5-year follow-up
- GKSRS patients survived an average of 4.6 years
- GKSRS patients lived 0.4 years (4.8 months) longer on average

**Statistical Significance:**
- p = 0.015: Difference is statistically significant
- 95% CI (0.1-0.7): Difference likely between 1.2 and 8.4 months

**Clinical Interpretation:**
Within the first 5 years, GKSRS patients survived an average of 4.6 years compared to 4.2 years for PBT patients, representing a clinically meaningful gain of 4.8 months (RMST difference=0.4 years, 95% CI: 0.1-0.7, p=0.015).

### Comparing RMST to Cox HR

| Scenario | RMST Interpretation | Cox HR Interpretation |
|----------|---------------------|----------------------|
| **RMST diff = +0.5 years, HR = 0.70** | GKSRS patients gained 6 months on average | GKSRS reduced hazard by 30% |
| **RMST diff = 0 years, HR = 1.00** | No survival difference | No hazard difference |
| **PH violated** | RMST still valid | Cox HR may be misleading |

**When Results Disagree:**
- RMST significant, Cox not: Effect real but varying over time
- Cox significant, RMST not: Effect consistent but small in absolute terms
- Trust RMST when PH assumption violated

---

## Reading Forest Plots

### Structure

**File:** `*_subgroup_forest_plot.png`

**Layout:**
```
Subgroup       PBT n/N  GKSRS n/N  HR [95% CI]      Forest Plot    P-interact
Age
  <63 years      8/45      7/45    0.55 [0.30-1.02]     ●——●
  ≥63 years     12/55     11/55    0.72 [0.45-1.15]      ●——●          0.43
Sex
  Male           13/48     12/40    0.60 [0.35-1.03]    ●——●
  Female          7/52      6/60    0.75 [0.42-1.35]      ●——●          0.58
```

**Key Components:**

**1. Subgroup Names (Left)**
- Patient subgroups being compared
- Hierarchical structure (category → levels)

**2. Sample Sizes**
- n = total patients in subgroup
- Events = number who experienced outcome
- Higher numbers = more reliable estimates

**3. Effect Estimates (Center)**
- HR [95% CI] or OR [95% CI]
- Numerical value of treatment effect per subgroup

**4. Forest Plot (Right)**
- **Point estimate (●):** HR or OR value
- **Horizontal lines:** 95% confidence interval
- **Vertical reference line:** HR/OR = 1.0 (no effect)
- **Left of line:** Protective effect (HR/OR < 1)
- **Right of line:** Harmful effect (HR/OR > 1)

**5. P-for-Interaction**
- Tests whether treatment effect differs across subgroup levels
- p < 0.05: Significant heterogeneity (effect modification)
- p ≥ 0.05: Consistent effect across subgroups

**6. Overall Effect (Diamond ♦)**
- Pooled estimate across all patients
- Width = 95% CI
- Reference for comparing subgroups

**How to Read:**

**Age Subgroup Example:**
- **<63 years:** HR=0.55 (45% risk reduction)
- **≥63 years:** HR=0.72 (28% risk reduction)
- **P-interact=0.43:** Difference not statistically significant
- **Interpretation:** Treatment appears more effective in younger patients, but evidence is weak (could be chance)

**Confidence Intervals:**
- **Narrow (●—●):** Precise estimate (adequate sample size)
- **Wide (●————●):** Imprecise (small sample or few events)
- **Crosses 1.0:** Not statistically significant
- **Doesn't cross 1.0:** Statistically significant

**Clinical Interpretation:**
Overall, GKSRS reduces the hazard by 35% (HR=0.65, 95% CI: 0.45-0.95). Treatment effects are consistent across age groups (p-interaction=0.43) and sex (p-interaction=0.58), though point estimates suggest slightly stronger effects in younger patients and males. These subgroup differences are not statistically significant and likely reflect random variation.

---

## Interpreting Proportional Hazards Diagnostics

### What is the PH Assumption?

**Assumption:** Hazard ratio remains constant over time

**Example:**
- PH holds: GKSRS always has 70% the hazard of PBT (HR=0.70 at all times)
- PH violated: GKSRS initially has 50% the hazard, but later 90% the hazard

### Diagnostic Files

**File:** `*_ph_diagnostics.xlsx`

**Example:**

```
Variable            Chi-sq    df    p-value
Treatment            4.52     1      0.034
Age                  0.85     1      0.36
Sex                  1.23     1      0.27
Global test          6.89     3      0.075
```

**How to Interpret:**

**Individual P-values:**
- **p > 0.05:** PH assumption holds for that variable
- **p ≤ 0.05:** PH assumption violated for that variable

**Treatment p = 0.034 (VIOLATED):**
- Treatment effect changes over time
- Cox HR is an average (may not reflect reality at any specific time)
- Use RMST as primary analysis

**Age p = 0.36 (HOLDS):**
- Age effect constant over time
- Cox model appropriate for age

**Global test p = 0.075:**
- Overall model PH assumption borderline
- Some variables violate assumption
- Proceed with caution; consider RMST

### Schoenfeld Residual Plots

**File:** `*_ph_plots_treatment.png`

**Components:**
- **X-axis:** Time
- **Y-axis:** Scaled Schoenfeld residuals
- **Points:** Residual at each event time
- **Smooth line:** Trend in residuals over time
- **Confidence band:** Uncertainty around trend

**How to Interpret:**

**Flat trend (horizontal line):**
- PH assumption holds
- Hazard ratio constant over time
- Cox regression appropriate

**Non-flat trend (sloped or curved line):**
- PH assumption violated
- Hazard ratio changes over time
- Use RMST or time-varying models

**Example Patterns:**

**Downward slope:**
- Treatment effect decreases over time
- Initially strong protection, diminishes later

**Upward slope:**
- Treatment effect increases over time
- Initially weak, strengthens later

**U-shape or inverted U:**
- Complex time-varying effects
- Consider stratified or parametric models

**Clinical Implications:**

**PH Holds:**
- Report Cox HR as primary measure
- RMST as secondary confirmation
- Standard interpretation

**PH Violated:**
- Report RMST as primary measure
- Cox HR as secondary (note violation)
- Describe time-varying pattern
- Consider time-stratified analysis

---

## Understanding GEP Analysis

Objective 4 compares the gene expression profile (GEP) results with outcomes observed in the cohort at 5 years. The imported GEP values are the laboratory's published group-level estimates (four values for the four Class x PRAME groups), so the analysis asks whether each group's published risk fits the patients in that group and how well the four groups order patients by risk. Formal definitions are in [GEP Validation](STATISTICAL_METHODS.md#gep-validation).

### Where to find the files

Each cohort's `04_GEP_Validation/` folder holds one workbook and three figures; the full cohort adds a `d_no_gep/` subfolder.

| File | What it holds |
|------|---------------|
| `<prefix>gep_validation_5yr.xlsx` | Primary workbook: observed vs imported 5-year risk by group (MFS and MSS), discrimination, sensitivity analyses |
| `<prefix>gep_validation_mfs_curves.png`, `<prefix>gep_validation_mss_curves.png` | Cumulative risk curves by group, with the imported 5-year risk marked |
| `<prefix>gep_validation_observed_vs_imported.png` | Observed 5-year risk with its interval beside the imported risk, for each group |
| `d_no_gep/<prefix>no_gep_5yr.xlsx` | Full cohort: outcomes in every GEP group, baseline characteristics by group, baseline-model performance |
| `d_no_gep/<prefix>no_gep_group_risk.png`, `d_no_gep/<prefix>no_gep_model_auc.png` | Observed 5-year risk by group; model and comparator AUC |

A `<prefix>..._diagnostics.xlsx` workbook sits next to each primary workbook and travels with it in published snapshots and the collaborator export. It documents how each result was produced and lists, by patient `id`, the records behind every exclusion or flag. It contains patient identifiers, so share it only with people who already have access to the study data. The three cohorts also share a cross-cohort observed-vs-imported figure in `merged_tables/objective4_cross_cohort_figures/`.

### Reading the group table

`Group_Table_MFS` and `Group_Table_MSS` have one row per Class x PRAME group and an overall row.

| Column | Meaning |
|--------|---------|
| `n` | Patients in the group |
| `events_by_horizon` | Events by 60 months |
| `event_free_through_horizon` | Patients known to be event-free at 60 months |
| `max_followup_months` | Longest follow-up in the group |
| `observed_risk`, `observed_ci_lower`, `observed_ci_upper` | 5-year risk estimated from follow-up (Kaplan-Meier for MFS, cumulative incidence for MSS), with its 95% interval |
| `imported_risk` | 1 minus the imported 5-year survival |
| `difference_pct_points` | Observed minus imported risk, in percentage points |
| `estimable`, `note` | Whether a 5-year value exists, and why not |

**Worked example (generic).** A row with `observed_risk` 0.20 (95% interval 0.08 to 0.45), `imported_risk` 0.15 and `n` 12 says: in this group of 12 patients the estimated 5-year risk is 20%; the imported 15% lies inside the interval, so the data are compatible with it; the interval is wide because few patients were followed. If the row had `event_free_through_horizon` of 2, only two patients were known to be event-free at 5 years and the estimate is imprecise. A row with `estimable` FALSE has no 5-year value because nobody in the group was followed to 60 months. A row with no events shows an observed risk of 0 and no interval.

The overall row averages the group risks weighted by group size and compares them with the group-size-weighted imported risk. It has no interval.

### Reading the figures

- **Curves.** Each step curve is one group's cumulative risk over follow-up, drawn to 60 months. The open diamond at 60 months is that group's imported 5-year risk. A curve that reaches the diamond means the observed and imported risks agree at 5 years. A curve that stops before 60 months marks the end of that group's follow-up. Groups with no events are flat at zero.
- **Observed vs imported.** One panel per endpoint. Each group shows the observed 5-year risk as a point with its interval and the imported risk as a marker; open circles mark groups with no events, whose interval cannot be estimated. A marker inside the interval means the data are compatible with the imported value.
- **No-GEP group risk.** Observed 5-year risk with intervals for all six GEP groups. A group with no point has no patient followed to 60 months.
- **Model AUC.** One panel each for MFS and MSS shows the AUC of the penalized models (filled points) next to the T stage and tumor diameter comparators (open diamonds), each with its interval; a third panel shows the Class 2-like score. The dashed line at 0.5 marks chance.

### Discrimination

`Discrimination` gives Harrell's C for MFS and the competing-risk AUC for MSS with intervals, the number of events, the share of tied pairs, and a `status`.

**Worked example (generic).** A C of 0.80 with interval 0.70 to 0.90, from 25 events and 20% tied pairs, says: in 80% of comparable pairs the patient with the earlier event had the higher imported risk, counting tied pairs as one half; the interval excludes 0.5. Because the imported risk has four values, this describes how well the four groups are ordered; it does not mean an individual's risk is predicted at that accuracy. A `status` of "not estimable" means fewer than 3 events by 60 months or no variation in the imported risk.

### Sensitivity sheets

- `Sensitivity_Late_Biopsy`: each primary quantity next to its value after removing patients sampled more than 180 days after treatment, with the difference.
- `Sensitivity_Cause_of_Death`: each group's risk with unknown-cause deaths counted as melanoma deaths, the number of deaths reclassified, and the difference in percentage points.
- `Sensitivity_Treatment_Era` (no-GEP workbook): each group's risk among patients treated in the testing era; see below.

**Worked example (generic).** A group risk of 12.0% that becomes 12.5% when late samples are removed moved by 0.5 percentage points; the late samples do not drive the result. A risk that moves from 9% to 16% when unknown-cause deaths are counted as melanoma deaths depends on how those deaths are classified, and both values are informative.

### Patients without a definitive call

`Group_Outcomes_MFS` and `Group_Outcomes_MSS` give the same measures as the group table for each of the six GEP groups, with median and maximum follow-up. `Baseline_By_Group` compares baseline characteristics across groups; a group with fewer than 3 patients shows only its size. Compare groups only where both have a 5-year estimate, and read differences between "Not tested" and the Class 1 and Class 2 groups with testing patterns in mind (for example the share treated in or after the year testing appears).

`Sensitivity_Treatment_Era` repeats each group's risk among patients treated in or after that year (`patients_in_era`), next to the primary risk, with the interval and the difference in percentage points. **Worked example (generic).** A group with 100 patients, 30 of them treated before testing began, a primary risk of 10% and an in-era risk of 14% (interval 8% to 24%) shifts by 4 percentage points; if the primary risk lies inside the in-era interval the era does not change the reading, and a group with `patients_before_era` of 0 is unchanged by construction.

### Baseline models without GEP

`Model_Performance` lists, for each endpoint and performance scope, the AUC, its interval, the Brier score, the no-information Brier score and the Brier skill for the penalized models and for the comparators. Rows labeled with the "No GEP" scope are the headline; "All patients" is shown for reference. `AUC_Difference` gives the model AUC minus the comparator AUC with its interval.

**Worked example (generic).** A model AUC of 0.62 (0.48 to 0.76) next to T stage alone at 0.70 gives a difference of -0.08. If the difference interval is -0.22 to 0.04, the data cannot distinguish the model from T stage alone, and the model is not shown to rank better. A Brier skill of 0.00 says the predicted risks are no closer to the outcomes than the overall rate. The intervals come from a bootstrap that does not refit the model, so they understate model-fitting uncertainty; with few events they are wide in any case.

### Diagnostics workbook

Read `Specification_Coding` first: each quantity's method, estimand, units, direction and plain-language reading. `Denominators_Flow` counts patients through each step. `Checks` gives each check a result and a status (no concern, concern, not assessed) with a reason; a check that was not run is "not assessed" and never appears as passed. `Record_Issues` lists, one row per record and issue, each with the patient `id` from the analytic dataset, the patients outside the primary set, the patients flagged for biopsy timing (including a biopsy long before treatment, which stays in the analysis), and the patients removed from a model; its counts reconcile with `Denominators_Flow`.

### Suggested interpretation template

> "Among N patients with a Class 1 or Class 2 result and PRAME reported, the observed 5-year metastasis risk was X% (95% CI a to b) in [group] compared with the imported Y%. Harrell's C for the imported risk was C (95% CI c to d); tied pairs, [Z%] of comparable pairs, count one half."

Take the numbers from the run's workbooks; the run log carries the source revision they came from.

---

## Clinical Significance vs Statistical Significance

### Understanding P-values

**P-value:** Probability of observing this difference if no true difference exists

**Common Misinterpretations:**
- ❌ p < 0.05 means clinically important
- ❌ p > 0.05 means no difference
- ❌ Smaller p-value = larger effect

**Correct Interpretation:**
- ✅ p < 0.05: Evidence against null hypothesis (unlikely due to chance alone)
- ✅ p > 0.05: Insufficient evidence to reject null (doesn't prove no difference)
- ✅ P-value measures statistical evidence, not clinical importance

### Assessing Clinical Significance

**Focus on Effect Size and Confidence Interval, Not Just P-value**

**Example 1: Statistically Significant but Clinically Trivial**
```
Vision change: β = -0.05 logMAR (95% CI: -0.08 to -0.02), p = 0.001
```
- Statistically significant (p < 0.05)
- Effect size tiny (0.05 logMAR ≈ half a line on eye chart)
- **Conclusion:** Not clinically meaningful despite significance

**Example 2: Not Significant but Potentially Important**
```
Overall survival: HR = 0.70 (95% CI: 0.45 to 1.10), p = 0.12
```
- Not statistically significant (p > 0.05)
- Point estimate suggests 30% risk reduction
- Wide CI includes both large benefit (55% reduction) and harm (10% increase)
- **Conclusion:** Possibly clinically important, but underpowered study

**Example 3: Statistically and Clinically Significant**
```
Local recurrence: OR = 0.45 (95% CI: 0.25 to 0.80), p = 0.007
```
- Statistically significant (p < 0.05)
- Large effect (55% reduction in odds)
- Narrow CI (excludes no effect)
- **Conclusion:** Strong evidence for meaningful clinical benefit

### Clinical Significance Thresholds

| Outcome | Minimal Clinically Important Difference (MCID) |
|---------|-----------------------------------------------|
| **Vision change** | 0.2 logMAR (2 lines on eye chart) |
| **Tumor height** | 1.0 mm reduction |
| **Overall survival HR** | 0.70 or 1.43 (30% change) |
| **Complication OR** | 0.50 or 2.0 (50% change) |

**How to Use:**
1. Check if confidence interval includes MCID
2. If yes: Clinically significant difference possible
3. If no: Effect likely too small to matter clinically

---

## Common Pitfalls

### 1. Confusing Association with Causation

**Issue:** Observational data can show associations but cannot prove causation

**Example:**
- GKSRS patients have better survival than PBT patients
- ❌ Conclusion: GKSRS causes better survival
- ✅ Conclusion: GKSRS is associated with better survival, but residual confounding possible

**Why:** Patients weren't randomized → unmeasured differences may explain results

### 2. Ignoring Confidence Intervals

**Issue:** Focusing only on point estimates without considering uncertainty

**Example:**
```
HR = 0.65 (95% CI: 0.45 to 0.95)
```
- Point estimate: 35% risk reduction
- CI range: 5% to 55% reduction
- **Truth likely somewhere in that range, not exactly 35%**

### 3. Over-Interpreting Subgroup Analyses

**Issue:** Finding "significant" subgroup differences when none truly exist

**Multiple Comparisons Problem:**
- Test 20 subgroups → expect 1 false positive by chance (p < 0.05)
- Require p < 0.05 for interaction, not individual subgroup p-values
- Exploratory findings need confirmation

**Example:**
- GKSRS HR = 0.50 in males (p=0.04), HR = 0.80 in females (p=0.35)
- P-for-interaction = 0.25
- **Conclusion:** No strong evidence for sex differences (interaction p > 0.05)

### 4. Ignoring PH Violations

**Issue:** Reporting Cox HR when proportional hazards assumption violated

**Example:**
```
Cox HR = 0.70 (p=0.05)
PH test: p = 0.02 (VIOLATED)
RMST difference = 0.1 years (p=0.30)
```

**Problem:** HR averages time-varying effects; may not reflect reality
**Solution:** Report RMST as primary, note PH violation

### 5. Misinterpreting Non-Significance

**Issue:** Concluding "no difference" when p > 0.05

**Example:**
```
HR = 0.75 (95% CI: 0.50 to 1.12), p = 0.16
```
- ❌ "GKSRS and PBT have equivalent survival"
- ✅ "We cannot rule out a 25% risk reduction or 12% risk increase with GKSRS"

**Reality:** Study may be underpowered; true effect remains uncertain

### 6. Ignoring Baseline Imbalances

**Issue:** Comparing treatment groups with different baseline characteristics

**Example:**
- PBT group: Mean age 55, mean tumor height 4 mm
- GKSRS group: Mean age 68, mean tumor height 7 mm
- Survival differs between groups
- **Problem:** Difference may reflect age/tumor size, not treatment

**Solution:** Check baseline tables; rely on adjusted regression models

---

## Troubleshooting

### "Why are survival curves missing?"

**Possible Reasons:**
1. **Insufficient events:** <5 total events (minimum for survival analysis)
2. **Insufficient groups:** Events in only 1 treatment group
3. **Analysis skipped:** Check logs for explanation

**How to Check:**
- Open `*_summary.xlsx` to see event counts
- Read `logs/txt/run_log_*.txt` for skip reasons
- Look for `*_analysis_not_performed.txt` files

**Example:**
```
GKSRS-only cohort, PFS-2: Only 3 events (need 5)
→ Survival curves not generated
→ Summary tables still available
```

### "Why do forest plots have missing subgroups?"

**Possible Reasons:**
1. **Sample size too small:** <2 patients per treatment arm
2. **No events:** 0 events in subgroup
3. **Statistical instability:** Extreme event imbalances

**How to Check:**
- Open `*_subgroup_diagnostics.xlsx`
- Look for "exclusion_reason" column
- Common reasons: "insufficient_sample", "zero_events", "skipped_non_finite"

**Example:**
```
Age 80+ years subgroup:
- PBT: n=3, events=0
- GKSRS: n=5, events=2
→ Excluded due to zero events in PBT arm
```

### "Why do regression tables show different results than summary tables?"

**Answer:** Adjustment for confounders

**Example:**
```
Summary table (unadjusted):
- PBT recurrence: 15%
- GKSRS recurrence: 20%
- Difference: +5% (GKSRS worse)

Regression table (adjusted):
- GKSRS OR = 0.75 (GKSRS better)

Explanation:
- GKSRS patients older with larger tumors (higher baseline risk)
- After adjusting for age and tumor size, GKSRS actually performs better
```

**Lesson:** Always interpret adjusted estimates, not crude rates

### "Why does RMST show significance but Cox doesn't (or vice versa)?"

**RMST Significant, Cox Not:**
- Treatment effect varies over time (PH violated)
- RMST captures net benefit; Cox averages contradictory effects
- **Trust RMST**

**Cox Significant, RMST Not:**
- Effect is small in absolute terms but consistent relative risk
- RMST may be underpowered for small differences
- Both are valid; report both

**Example:**
```
Cox HR = 0.70 (p=0.04)
RMST difference = 0.2 years (p=0.15)
PH test: p=0.03 (violated)

Interpretation: Treatment reduces hazard by 30% on average, but effect 
varies over time. Absolute survival gain is 2.4 months, which is not 
statistically significant with current sample size. Report RMST as primary 
given PH violation.
```

### "How do I export specific results?"

**Excel Tables:**
- Open in Excel or LibreOffice
- Copy relevant cells
- Paste into manuscript

**HTML Tables:**
- Open in browser
- Right-click → "Save As" → PDF (print to PDF)
- Or copy table and paste into Word

**PNG Figures:**
- Already publication-ready
- Resolution: 300 DPI
- Insert directly into manuscript

**Tips:**
- `merged_tables/` folder has cross-cohort comparisons, including a separately labeled baseline table spanning full, restricted, and GKSRS-only cohorts
- `*_summary.xlsx` has descriptive statistics
- `*.html` tables have full model details with diagnostics

---

## Next Steps

After interpreting results:

1. **Verify Findings:**
   - Check baseline balance
   - Review exclusion reasons
   - Validate event counts

2. **Draft Interpretation:**
   - Focus on adjusted estimates (regression tables)
   - Consider clinical significance (effect sizes, MCID)
   - Note limitations (sample size, PH violations)

3. **Plan Manuscript:**
   - Table 1: Baseline characteristics (`baseline_characteristics.xlsx`)
   - Table 2: Primary outcomes (survival, recurrence rates)
   - Table 3: Secondary outcomes (complications, tumor response)
   - Figure 1: Survival curves
   - Figure 2: Forest plots (subgroup analysis)

4. **Consult Statistical Methods:**
   - See [STATISTICAL_METHODS.md](STATISTICAL_METHODS.md) for methodology details
   - See [TECHNICAL.md](TECHNICAL.md) for implementation specifics
   - See [CALCULATIONS.md](CALCULATIONS.md) for variable definitions

5. **Report Results:**
   - Cite the analysis pipeline (see README.md)
   - Include software versions from logs
   - Describe quality assurance procedures

---

## Getting Help

**Documentation:**
- [README.md](../README.md) - Overview and quick start
- [CALCULATIONS.md](CALCULATIONS.md) - Variable derivations
- [TECHNICAL.md](TECHNICAL.md) - Implementation details
- [STATISTICAL_METHODS.md](STATISTICAL_METHODS.md) - Statistical approaches

**Common Issues:**
- Check logs in `logs/txt/` for error messages
- Review diagnostics workbooks (`*_diagnostics.xlsx`); the sheets and how to read them are in [STATISTICAL_METHODS.md](STATISTICAL_METHODS.md#model-diagnostics-workbooks-objectives-1-3)
- Verify input data quality
- Confirm adequate sample sizes and event counts

**Contact:**
Refer to README.md for author contact information and citation details.
