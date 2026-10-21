# Statistical Methods

This document is the canonical statistical-methods reference for the repository. It defines the model families, assumptions, thresholds, and validation metrics used throughout the analysis pipeline.

For output-reading guidance, use [INTERPRETATION_GUIDE.md](INTERPRETATION_GUIDE.md). For derived-variable definitions and endpoint construction, use [CALCULATIONS.md](CALCULATIONS.md).

---

## Table of Contents

- [Overview](#overview)
- [Survival Analysis](#survival-analysis)
- [Restricted Mean Survival Time (RMST)](#restricted-mean-survival-time-rmst)
- [Proportional Hazards Assumption Testing](#proportional-hazards-assumption-testing)
- [Competing Risks Analysis](#competing-risks-analysis)
- [Binary Outcomes](#binary-outcomes)
- [Continuous Outcomes](#continuous-outcomes)
- [Subgroup Analysis](#subgroup-analysis)
- [GEP Validation](#gep-validation)
- [Statistical Metric Categorization](#statistical-metric-categorization)

---

## Overview

The analysis employs standard statistical methods appropriate for comparative effectiveness research:

| Analysis Type | Method | R Package | Function |
|---------------|--------|-----------|----------|
| **Survival analysis** | Kaplan-Meier, Cox regression, RMST | `survival` | `survfit()`, `coxph()`, `rmst2()` |
| **Binary outcomes** | Logistic regression | Base R | `glm(family="binomial")` |
| **Continuous outcomes** | Linear regression | Base R | `lm()` |
| **Table generation** | Summary tables, regression tables | `gtsummary` | `tbl_summary()`, `tbl_regression()` |
| **Proportional hazards testing** | Scaled Schoenfeld residuals | `survival` | `cox.zph()` |
| **Competing risks** | Cumulative incidence functions | `cmprsk` | `cuminc()` |

---

## Survival Analysis

### Kaplan-Meier Method

**Purpose:** Estimate survival probabilities over time

**Implementation:** `survival::survfit()`

**Outputs:** 
- Survival curves stratified by treatment
- Median survival times with 95% CI
- Survival probability tables

**Key Features:**
- Right-censored survival data
- Non-parametric estimation
- No distributional assumptions

For practical reading of survival curves, see [Kaplan-Meier Plots](INTERPRETATION_GUIDE.md#kaplan-meier-plots).

### Cox Proportional Hazards Regression

**Purpose:** Compare hazard rates between treatment groups while adjusting for confounders

**Implementation:** `survival::coxph()`

**Model Formula:**
```r
Surv(time, event) ~ treatment + age_at_diagnosis + sex + location
```

**Outputs:**
- Hazard ratios with 95% CI
- P-values from likelihood ratio tests
- Model diagnostics

**Assumptions:**
- Proportional hazards (tested with `cox.zph()`)
- Independent censoring
- Multiplicative hazard effects

For practical reading of hazard ratios and confidence intervals, see [Cox Regression Tables (Survival Outcomes)](INTERPRETATION_GUIDE.md#cox-regression-tables-survival-outcomes).

### Restricted-Cohort Propensity-Overlap Sensitivity

The existing adjusted Cox models remain the primary analyses. In the restricted cohort, an exploratory sensitivity uses one shared propensity score to address measured pretreatment differences between PBT and GKSRS for local recurrence, metastatic progression, overall survival, and progression-free survival. The treatment-assignment model estimates the probability of GKSRS:

```r
treatment_group ~ age_at_diagnosis + sex + location +
    initial_tumor_height + initial_tumor_diameter + srf + treatment_year
```

Age, tumor height, tumor diameter, and treatment year are continuous. Objective 0 derives integer `treatment_year` once from `treatment_date`, stores it in each analytic dataset, and validates it for downstream Objective 1 and Objective 2 use; downstream analyses do not derive or center year.

Patients in location levels with fewer than five observations are excluded only from the propensity-model population. For propensity score $e_i=P(\mathrm{GKSRS}\mid X_i)$, overlap weights are

$$
w_i =
\begin{cases}
1-e_i, & \mathrm{GKSRS} \\
e_i, & \mathrm{PBT}.
\end{cases}
$$

These weights target the average treatment effect in the overlap population (ATO). No matching, trimming, weight capping, stabilization, or support-based filtering is applied. The same weights are reused in four treatment-only Cox models fit with robust sandwich variance; effects are reported as GKSRS-versus-PBT hazard ratios with 95% confidence intervals. Local-recurrence and metastatic-progression models estimate cause-specific hazards, with death before the endpoint censored and reported as a competing event. The existing cumulative-incidence and Gray summaries are unweighted descriptive context; Objective 1 does not estimate a propensity-adjusted Fine-Gray or cumulative-incidence treatment effect.

Diagnostics report propensity-score overlap, weight distributions, effective sample size, unweighted and overlap-weighted standardized mean differences against the 0.10 threshold, and treatment/global Schoenfeld results with a combined residual figure that displays the corresponding p-values. The sensitivity runs only for `uveal_melanoma_restricted_cohort` and writes its reader-facing workbook, plots, and narrative summary to that cohort's `01_Efficacy/h_propensity_score_sensitivity/` folder; the `.rds` audit is technical, runtime-only evidence. See the [README path model](../README.md#path-model) for the runtime-root layout. This sensitivity remains associational: it adjusts only measured covariates, cannot address unavailable comparable dose or quantitative tumor-proximity variables, and does not remove residual confounding. Limited outcome events also constrain precision.

### Log-Rank Test

**Purpose:** Compare survival distributions between groups

**Implementation:** Built into `survival::survfit()`

**Interpretation:**
- p < 0.05 indicates significantly different survival distributions
- Non-parametric test (no distributional assumptions)
- Sensitive to differences across entire follow-up period

For practical reading of survival-curve separation and log-rank context, see [Kaplan-Meier Plots](INTERPRETATION_GUIDE.md#kaplan-meier-plots).

---

## Restricted Mean Survival Time (RMST)

### Motivation

**Problem:** Cox regression assumes proportional hazards (constant hazard ratio over time), which is often violated in real data

**Solution:** RMST provides a non-parametric, clinically interpretable alternative that does not require the proportional hazards assumption

### What is RMST?

**Definition:** The area under the survival curve up to a specified time point (τ)

**Clinical Interpretation:** The average survival time (in years) up to time τ

**Example:**
- RMST at 3 years = 2.5 years means patients survived an average of 2.5 years during the first 3 years of follow-up
- RMST difference = 0.3 years means one treatment group survived 0.3 years (3.6 months) longer on average

### RMST Analysis Implementation

**Package:** `survRM2::rmst2()`

**Time Horizon (τ):**
- Objective 1c (Overall Survival): 5 years
- Objective 1d (PFS): 3 years
- Objective 3a (PFS-2): 3 years
- Objective 4 (MFS, MSS): 5 years

**Outputs:**
1. **RMST Estimates:** Mean survival time for each treatment group
2. **RMST Difference:** Treatment effect (years gained/lost)
3. **P-value:** Statistical significance of difference
4. **95% CI:** Uncertainty around difference
5. **Visualization:** RMST curves with shaded areas

### Advantages of RMST

| Feature | RMST | Cox HR |
|---------|------|--------|
| **Assumption-free** | ✓ No proportional hazards assumption | ✗ Requires proportional hazards |
| **Clinical interpretation** | ✓ Years gained/lost | ✗ Multiplicative ratio |
| **Non-parametric** | ✓ Distribution-free | ✗ Assumes exponential-family baseline |
| **Direct comparison** | ✓ Absolute difference | ✗ Relative measure |

### When to Use RMST vs Cox

**Use RMST when:**
- Proportional hazards assumption violated
- Direct clinical interpretation needed ("X months gained")
- Regulatory submissions requiring robust methods
- Communicating results to non-statisticians

**Use Cox when:**
- Proportional hazards assumption holds
- Relative risk measures preferred
- Standard in your clinical field
- Comparing multiple adjusted models

**Best Practice:** Report both RMST and Cox HR when possible

For plain-English interpretation of RMST outputs, see [Understanding RMST Analysis](INTERPRETATION_GUIDE.md#understanding-rmst-analysis) and [Comparing RMST to Cox HR](INTERPRETATION_GUIDE.md#comparing-rmst-to-cox-hr).

---

## Proportional Hazards Assumption Testing

The proportional hazards (PH) assumption states that the hazard ratio between treatment groups remains constant over time.

### Testing Methodology

**Method:** Scaled Schoenfeld residuals test

**Implementation:** `survival::cox.zph()`

**Test Statistic:** Correlation between scaled residuals and time

**Interpretation:**
- **p > 0.05:** PH assumption holds (constant HR over time)
- **p ≤ 0.05:** PH assumption violated (time-varying HR)

### Diagnostic Outputs

For each time-to-event analysis with ≥10 total events, the pipeline generates:

**1. Statistical Tests** (`*_ph_diagnostics.xlsx`)
- Global test p-value (overall model assumption)
- Per-covariate test p-values
- Chi-square statistics
- Degrees of freedom

**2. Schoenfeld Residual Plots** (`*_ph_plots_{variable}.png`)
- Residuals vs time scatter plots
- Smoothed trend line
- 95% confidence band
- Flat trend = PH holds; non-flat = PH violation

### Clinical Implications

**When PH holds:**
- Cox HR provides valid treatment effect estimate
- HR interpretation straightforward (constant over time)
- Standard survival analysis appropriate

**When PH violated:**
- Cox HR may be misleading (average effect that varies over time)
- RMST provides more appropriate analysis
- Consider time-stratified models or parametric alternatives
- Report both Cox and RMST results

### Example Interpretation

```
Treatment PH Test: p = 0.03 (VIOLATION)
```

**Interpretation:** The treatment hazard ratio changes over time. For example, GKSRS may be more protective early but less protective late in follow-up. Use RMST for primary analysis and report Cox HR as secondary.

Model feasibility thresholds are kept separate by analysis family: adjusted adverse-event logistic regression uses a 10-event minimum, PFS-2 survival modeling uses a 10-patient precheck and a 5-event minimum, and PH diagnostics follow a 10-event reporting floor.

### Automation

The pipeline automatically:
1. Tests PH assumption for all Cox models with ≥10 events
2. Generates diagnostic plots for each covariate
3. Documents violations in Excel workbooks
4. Provides RMST as alternative when PH violated
5. Stores diagnostics in `{cohort}/{objective}/h_proportional_hazards_diagnostics/`

**Minimum Event Threshold:** PH testing skipped when <10 total events (insufficient power for reliable testing)

For practical reading of PH tests and Schoenfeld plots, see [Interpreting Proportional Hazards Diagnostics](INTERPRETATION_GUIDE.md#interpreting-proportional-hazards-diagnostics) and [Schoenfeld Residual Plots](INTERPRETATION_GUIDE.md#schoenfeld-residual-plots).

---

## Competing Risks Analysis

### What are Competing Risks?

**Definition:** Events that preclude or alter the probability of the event of interest

**Example:** Death is a competing risk for local recurrence because patients who die cannot subsequently develop recurrence

### Why Standard Kaplan-Meier Fails

**Problem:** Kaplan-Meier treats competing events as censoring, which:
- Overestimates event probabilities
- Assumes patients remain "at risk" after competing event
- Produces biased cumulative incidence estimates

**Solution:** Cumulative Incidence Functions (CIF) account for competing events

### Implementation

**Package:** `cmprsk::cuminc()`

**Event Coding:**
```r
0 = Censored (still at risk)
1 = Event of interest (e.g., local recurrence)
2 = Competing event (e.g., death before recurrence)
```

**Outputs:**
1. **Cumulative Incidence Curves:** Probability of event over time accounting for competing events
2. **Gray's Test:** One global statistical comparison of treatment-group cumulative-incidence curves
3. **Fine-Gray Model:** A possible subdistribution-hazard model for a separately specified estimand; it is not estimated for the Objective 1 propensity sensitivity

Objective 1 recurrence and metastatic-progression cumulative-incidence sheets report this as `gray_test_global_curve_p_value`. The same value can appear on multiple horizon rows because it describes the overall across-group curve comparison, not a separate p-value for each horizon or group.

### When Competing Risks Matter

| Outcome | Competing Event | Impact |
|---------|-----------------|--------|
| **Local Recurrence** | Death before recurrence | High (older patients, comorbidities) |
| **Metastatic Progression** | Death before metastasis | High (aggressive disease) |
| **PFS** | None (composite endpoint) | N/A (already accounts for death) |
| **Overall Survival** | None (death is the outcome) | N/A |
| **Vision Changes** | Death before vision loss | Moderate |

### Interpretation

**Cumulative Incidence at 5 years = 15%**
- 15% of patients developed the event by 5 years, accounting for those who died first
- More realistic than K-M estimate which might overestimate as 20%

**Gray's Test p = 0.03**
- Treatment-group cumulative-incidence curves differ significantly over follow-up
- Accounts for differential competing risk rates
- This is a global curve-comparison p-value, not a timepoint-specific p-value

For practical reading of cumulative-incidence outputs, see [Cumulative Incidence Plots (Competing Risks)](INTERPRETATION_GUIDE.md#cumulative-incidence-plots-competing-risks).

---

## Binary Outcomes

### Logistic Regression

**Purpose:** Compare binary outcome rates (yes/no events) between treatment groups

**Implementation:** `glm(family = binomial(link = "logit"))`

**Model Formula:**
```r
outcome ~ treatment + age_at_diagnosis + sex + location
```

**Outputs:**
- Odds ratios with 95% CI
- P-values from likelihood ratio tests
- Model diagnostics

**Interpretation:**
- OR > 1: Higher odds in treatment group
- OR < 1: Lower odds in treatment group
- OR = 1: No difference

### Applications

| Outcome | Analysis Type | Location |
|---------|---------------|----------|
| **Radiation Retinopathy** | Recorded toxicity burden by available follow-up | Objective 2b |
| **Neovascular Glaucoma** | Recorded toxicity burden by available follow-up | Objective 2c |
| **Serous Retinal Detachment (all recorded causes in the published implementation)** | Recorded toxicity burden by available follow-up | Objective 2d |

Objective 2d includes all recorded SRD causes, including mass-induced SRD when present.

Objective 2 toxicity endpoints are analyzed as recorded burden fields prepared during Objective 0 (`retinopathy_burden_event`, `nvg_burden_event`, `srd_burden_event`). They are descriptive binary burden-by-follow-up endpoints, not time-to-toxicity incidence endpoints. Objective 0 validation requires complete canonical source values and complete binary 0/1 burden fields in included analytic rows before Objective 2 consumes them.

### Assumptions

- Independent observations
- Linearity of log-odds
- No perfect multicollinearity
- Adequate sample size (rule of thumb: 10 events per predictor)

For practical reading of odds ratios and logistic-regression outputs, see [Logistic Regression Tables (Binary Outcomes)](INTERPRETATION_GUIDE.md#logistic-regression-tables-binary-outcomes).

---

## Continuous Outcomes

### Linear Regression

**Purpose:** Compare continuous outcome changes between treatment groups

**Implementation:** `lm()`

**Model Formulas:**

**Unadjusted (Primary):**
```r
outcome_change ~ treatment
```

**Adjusted for Baseline (Sensitivity):**
```r
outcome_change ~ treatment + baseline_value + confounders
```

**Outputs:**
- Regression coefficients (mean difference)
- 95% CI
- P-values
- R² (variance explained)

### Applications

| Outcome | Primary Analysis | Sensitivity Analysis | Location |
|---------|------------------|---------------------|----------|
| **Tumor Height Change** | Linear regression adjusted for the shared covariate set (`age_at_diagnosis`, `sex`, `location`) | Internal diagnostic change-score sensitivity that includes baseline tumor height, which is also part of the outcome definition | Objective 1e, 1f |
| **Vision Change (logMAR)** | Unadjusted | Adjusted linear regression on the implemented change score; baseline vision is not added as a separate covariate. A separate latest-VA sensitivity models final logMAR with baseline VA, explicit latest-VA follow-up duration, viable reviewer-requested baseline predictors, and the shared covariate set | Objective 2a |
| **Snellen Line Change (exact integer lines)** | Descriptive converted summary row | Adjusted linear regression | Objective 2a |

Tumor-height change remains a secondary/descriptive endpoint for reviewer-response interpretation. Because `height_change` subtracts baseline height, adding `initial_tumor_height` to the change-score model creates algebraic coupling between the covariate and outcome. That sensitivity can be used only as an internal diagnostic for dependence on baseline size; it is not ordinary confounder adjustment, it is not a reviewer-facing solution, and it does not address unequal time from treatment to follow-up height measurement. Tumor-height timing summaries must be reviewed alongside any comparative estimate, and current timing evidence supports demoting rather than emphasizing the comparative tumor-height result.

### Interpretation

**Tumor Height Change Coefficient = -0.5 mm**
- Treatment group had 0.5 mm greater decrease (more shrinkage) than reference within the model specification
- Negative change = shrinkage (see [CALCULATIONS.md](CALCULATIONS.md))

**Vision Change Coefficient = -0.2 logMAR**
- Treatment group had 0.2 logMAR less worsening on the implemented visual-acuity change score
- Negative change = worsening vision (see [CALCULATIONS.md](CALCULATIONS.md))

### Assumptions

- Linearity of relationships
- Homoscedasticity (constant variance)
- Normality of residuals
- Independent observations

For practical reading of linear-regression outputs, see [Linear Regression Tables (Continuous Outcomes)](INTERPRETATION_GUIDE.md#linear-regression-tables-continuous-outcomes).

Objective 2 vision change models use the precomputed `vision_change` endpoint (`initial_vision - last_vision`, or `initial_vision - recurrence1_pretreatment_vision` for patients with local recurrence). The adjusted change-score rows include the central confounder set (`age_at_diagnosis`, `sex`, and `location`) but do not include a separate baseline-vision covariate.

Reviewer-response latest-VA sensitivity uses an ANCOVA-style linear model with `last_vision` as the outcome. It adjusts for treatment group, `initial_vision`, explicit treatment-to-`last_followup` timing, viable reviewer-requested baseline predictors (`initial_tumor_height`, `initial_tumor_diameter`, `initial_t_stage_simple`, `srf`, optic-nerve involvement when variable within cohort, and treatment year), plus the central confounder set. The full-cohort sensitivity estimated a GKSRS-vs-PBT mean difference of -0.200 logMAR (95% CI -0.462 to 0.063; p = 0.135; n = 208). The restricted-cohort sensitivity estimated -0.126 logMAR (95% CI -0.429 to 0.177; p = 0.411; n = 126). These sensitivity results should be read as support for cautious interpretation, not causal evidence of superior functional preservation.

Latest-VA minimum-follow-up sensitivity outputs report explicit treatment-to-`last_followup` timing as the primary timing surface and a separately labeled proxy surface that uses derived general `follow_up_months` when explicit timing is missing. The same workbook includes a `reviewer_predictor_availability` sheet documenting requested visual-outcome predictors that were included, excluded for cohort support, or not available as structured baseline/dosimetry fields.

---

## Ordinal Outcomes

### Ordinal Logistic Regression

**Purpose:** Compare ordered categorical outcomes between treatment groups while preserving clinical rank order.

**Implementation:** `MASS::polr()`

**Model Formula:**

```r
ordered_outcome ~ treatment + confounders
```

**Outputs:**
- Odds ratios for a shift toward earlier ordered categories
- 95% Wald CI
- Overall treatment/variable p-values from likelihood ratio tests
- Diagnostics workbook with sample-size and outcome-level summaries

### Applications

| Outcome | Ordered Levels | Location |
|---------|----------------|----------|
| **Snellen Line Change Distribution** | `≥3-line improvement`, `2-line improvement`, `1-line improvement`, `Stable (0-line change)`, `1-line loss`, `2-line loss`, `≥3-line loss` | Objective 2a |

### Interpretation

**Snellen distribution OR = 1.40**
- The comparison treatment has higher odds of falling into a better Snellen Line Change Distribution category
- Because categories are ordered from improvement to loss, OR > 1 favors better vision outcomes and OR < 1 favors worse outcomes

### Assumptions

- Ordered outcome categories are clinically meaningful
- Proportional odds assumption is not formally tested in the pipeline and should be treated as an interpretation assumption
- Independent observations
- Adequate observations across outcome levels

For practical reading of ordered-outcome model tables, see [Ordinal Regression Tables (Ordered Outcomes)](INTERPRETATION_GUIDE.md#ordinal-regression-tables-ordered-outcomes).

For Objective 2 vision, ordinal modeling is reserved for the 7-level `Snellen Line Change Distribution`. The exact integer `Snellen Line Change` outcome is modeled with adjusted linear regression, and a full exact-integer ordinal model is not used because the observed line-count support is very wide and sparse across cohorts. Reader-facing output files therefore separate `Snellen Line Change` from `Snellen Line Change Distribution`, and each Objective 2 subfolder now also includes a flat `*_effect_summary.xlsx` workbook that combines descriptive, unadjusted, and adjusted rows in one sheet.

For avoidance of doubt, the stable distribution level is the exact `0-line` category after nearest-line rounding with halves away from zero. Small non-zero logMAR deltas with absolute magnitude below 0.05 therefore round to `0` and are reported as stable rather than as 1-line gain or loss.

For consistency across the ordinal HTML tables and the effect-summary workbooks, ordinal (`polr`) treatment effects are reported as proportional-odds ORs with 95% Wald confidence intervals and likelihood-ratio-test p-values. Other model families retain their standard reporting conventions: linear models report mean differences with Wald CIs/p-values, logistic models report odds ratios with model-based Wald CIs and the pipeline's standard term-level p-values, and Cox models report hazard ratios with the native Cox confidence intervals and Cox-model p-values.

Objective 2 categorical descriptive summaries that require simulated Fisher p-values use the local seed `20260422`. The seed is restored after each calculation so these p-value notes are reproducible without changing unrelated random-number state.

Model-facing factors use Objective 0 canonical levels or shared level-preserving helpers. Ordered factors are converted to unordered factors before regression when needed so treatment contrasts remain stable while preserving the explicit level order.

---

## Subgroup Analysis

### Purpose

Evaluate whether treatment effects differ across patient subgroups defined by:
- Age at diagnosis (<63 versus ≥63 years)
- Sex
- Tumor location
- Tumor size categories
- Baseline vision
- Disease stage

### Statistical Approach

**Method:** Stratified analysis with interaction testing

The exploratory Objective 1 age subgroup uses the 63-year split. Age remains continuous in the adjusted treatment-effect models outside this subgroup display. The selected age representation is controlled by `OBJECTIVE1_AGE_SUBGROUP_VAR` in `scripts/config/objective1_contracts.R`. Local recurrence, metastatic progression, overall survival, and progression-free survival subgroup effects are all estimated with Cox proportional-hazards models and reported as hazard ratios.

**Steps:**
1. Fit model within each subgroup
2. Estimate treatment effect per subgroup
3. Test treatment × subgroup interaction
4. Generate forest plots visualizing heterogeneity

**Interaction Test:**
```r
outcome ~ treatment * subgroup_variable + confounders
```

**Interpretation:**
- **p < 0.05:** Evidence of interaction (treatment effect may vary by subgroup)
- **p ≥ 0.05:** No strong evidence of interaction (overall treatment effect is the safer summary)

### Filtering Criteria

See [TECHNICAL.md](TECHNICAL.md#subgroup-filtering) for detailed subgroup filtering requirements.

**Summary:**
- Minimum 2 patients per treatment arm
- Minimum 1 event per treatment arm (survival outcomes)
- Automatic exclusion of unstable comparisons

### Outputs

**1. Forest Plots** (`*_forest_plot.png`)
- Treatment effect estimates per subgroup
- 95% confidence intervals
- Overall pooled estimate
- P-value for interaction

**2. Consolidated Diagnostic Workbooks**
- `forest_plots/*_forest_plot_diagnostics.xlsx`
- `tumor_height_primary/*_primary_tumor_height_diagnostics.xlsx`
- `tumor_height_sensitivity/*_sensitivity_tumor_height_diagnostics.xlsx`
- `*_age_decade_subgroup_sensitivity.xlsx`

The event-time sheets identify the outcome, endpoint type, Cox model family, hazard-ratio effect measure, time variable, event variable, and estimand. They also contain per-arm sample/event counts, effect estimates with confidence intervals, exclusions, and non-estimability reasons. When no subgroup effects are supportable, the workbook receives an explicit `Analysis_Status` sheet rather than remaining empty or failing.

### Interpretation Guidelines

**Possible Interaction (p < 0.05):**
- Treatment effect genuinely differs across subgroups
- Clinically meaningful heterogeneity
- Consider subgroup-specific treatment recommendations

**No Strong Interaction Signal (p ≥ 0.05):**
- No strong evidence for differential effects
- Apply overall treatment effect across subgroups
- Observed differences likely due to chance

**Caution:** Subgroup analyses are exploratory and hypothesis-generating. Confirm findings in independent datasets before clinical application.

For practical reading of subgroup outputs, see [Reading Forest Plots](INTERPRETATION_GUIDE.md#reading-forest-plots).

---

## GEP Validation

Objective 4 describes how the laboratory's gene expression profile (GEP) results relate to observed outcomes in this cohort. The claims it supports are descriptive and associational: observed 5-year risk by GEP group, and how well a risk ordering separates patients who had the event from those who did not.

For reading the generated files, see [Understanding GEP Analysis](INTERPRETATION_GUIDE.md#understanding-gep-analysis). Variable definitions are in [CALCULATIONS.md](CALCULATIONS.md#gep-variables).

### Contents

- [Goals and subgoals](#goals-and-subgoals)
- [Imported predictions](#imported-predictions)
- [GEP groups](#gep-groups)
- [Observed 5-year risk](#observed-5-year-risk)
- [Discrimination](#discrimination)
- [Sensitivity analyses](#sensitivity-analyses)
- [Patients without a definitive call](#patients-without-a-definitive-call)
- [Baseline models without GEP](#baseline-models-without-gep)
- [Diagnostics workbooks and record-level issues](#diagnostics-workbooks-and-record-level-issues)

### Goals and subgoals

| Goal | Question | Population | Estimand | Method | Primary output |
|------|----------|------------|----------|--------|----------------|
| **4a** Observed vs imported risk | Does each group's imported 5-year risk lie within the risk observed in that group? | Patients with a Class 1 or Class 2 call and a PRAME result (four Class x PRAME groups) | 5-year probability of metastasis (MFS) and of melanoma death (MSS), per group | Kaplan-Meier (MFS), Aalen-Johansen (MSS) | `*_gep_validation_5yr.xlsx` (`Group_Table_MFS`, `Group_Table_MSS`), observed-vs-imported and curve figures |
| **4b** Discrimination | How well does the imported risk order patients who have the event before those who do not? | Same as 4a | Concordance (MFS); time-dependent AUC with competing risks (MSS) | Tie-aware Harrell's C; cumulative/dynamic AUC; patient-level bootstrap | `Discrimination` sheet |
| **4c** Outcomes by availability group | What are the 5-year outcomes of patients without a definitive call? | All six GEP groups, full cohort | 5-year risk per group | Kaplan-Meier, Aalen-Johansen | `*_no_gep_5yr.xlsx` (`Group_Outcomes_MFS`, `Group_Outcomes_MSS`, `Baseline_By_Group`) |
| **4d** Baseline features | Do baseline clinical features rank 5-year risk better than T stage or tumor diameter alone in patients without a definitive call? | Full cohort, scored in the "No GEP" scope | Out-of-fold AUC and Brier skill; paired AUC difference to a comparator | Repeated nested cross-validation of a ridge model with inverse-probability-of-censoring weights | `Model_Performance`, `AUC_Difference` sheets |

The sensitivity analyses (late biopsy, unknown cause of death, treatment era) belong to 4a-4c and are described below.

### Imported predictions

The analytic dataset carries two imported columns, `biopsy1_gep_mfs` and `biopsy1_gep_mss`: the laboratory's published 5-year metastasis-free and melanoma-specific survival estimates for the patient's Class x PRAME group. They are group-level values, entered by hand from the published estimates, so the four groups of patients with a Class and a PRAME result carry four distinct values and every patient in a group carries the same one. The pipeline copies them into `expected_mfs_5yr` and `expected_mss_5yr` and converts them to the **imported risk** $1 - S(5)$.

Because the imported value is a group property, agreement between a group's observed and imported risk says whether the published estimate fits this cohort's patients in that group. It is not a patient-level calibration.

Patients whose class is known but whose PRAME result is not reported carry, in the source spreadsheet, the value of the matching PRAME-negative group. That is an assumption of PRAME negativity made in the source, and these values are not used: the patients are reported as their own group (see below).

### GEP groups

Each patient belongs to exactly one availability group, assigned once from the laboratory's source label by `classify_gep_group()` in `scripts/config/objective4_policy.R`. The group names, their order and the raw source labels are defined only in that file.

| Group | Meaning |
|-------|---------|
| Class 1 | Definitive Class 1 call with PRAME reported |
| Class 2 | Definitive Class 2 call with PRAME reported |
| Class known, PRAME not reported | Class reported, PRAME not reported; no prediction is imputed |
| Discordant result | The laboratory flagged discordant results |
| Failed or unknown result | Test failed or result unknown |
| Not tested | No GEP result recorded |

Class 1 and Class 2 together are the **primary validation set**. For the four Class x PRAME groups, `biopsy1_gep` carries the label "Class 1, PRAME negative", "Class 1, PRAME positive", "Class 2, PRAME negative" or "Class 2, PRAME positive"; for the other groups it carries the group name. These labels, and the group colors in figures, are the same in every objective.

Terminology follows the current laboratory guideline, which reports Class 1 and Class 2 only. The source spreadsheet still carries the earlier 1A/1B subclass wording in its strings; analysis and outputs use Class 1 and Class 2.

A patient enters the primary validation set for an endpoint when the group is Class 1 or Class 2, the imported value for that endpoint lies in [0, 1], and the endpoint's own follow-up fields are present (`mfs_analysis_eligible`, `mss_analysis_eligible`; metastatic disease at baseline excludes a patient from MFS). Eligibility is defined once, in `derive_gep_analysis_flags()`.

### Observed 5-year risk

**Estimands.** The 5-year probability of metastasis (MFS) and the 5-year probability of melanoma death (MSS). Time starts at treatment.

**Event definitions.**

- MFS: time to documented distant metastasis; death before metastasis is censored.
- MSS: time to death with recorded cause `Metastatic_Uveal_Melanoma`. Every other death, including a death with unknown recorded cause, is a competing event. Patients without death are censored.

**Estimators.** MFS risk is $1 - \hat S(60)$ from the Kaplan-Meier estimator,

$$
\hat S(t) = \prod_{t_i \le t}\left(1 - \frac{d_i}{n_i}\right),
$$

with $d_i$ events and $n_i$ patients at risk at time $t_i$ (months). MSS risk is the Aalen-Johansen cumulative incidence of melanoma death,

$$
\hat F_1(t) = \sum_{t_i \le t} \hat S(t_{i-})\,\frac{d_{1i}}{n_i},
$$

where $\hat S$ is the Kaplan-Meier estimate of remaining free of every event type and $d_{1i}$ counts melanoma deaths. Intervals are 95% log-log intervals. Both are computed within each group.

**Rules.**

- No value is carried past a group's last follow-up. A group whose longest follow-up is under 60 months has no 5-year estimate.
- A group with no event by 60 months has a point estimate of 0 and no interval.
- A group with fewer than 5 patients known event-free at 60 months is flagged as sparse (`Checks`).
- The overall row is the group-size-weighted average of the group risks, compared with the group-size-weighted imported risk. It has no interval.

**Assumption.** Censoring is independent of risk within a group: a patient censored at time $c$ has the same future risk as patients in the same group who are still followed at $c$.

**Rationale.** A patient censored before 5 years is not a known 5-year survivor, so counting them as event-free overstates survival. Kaplan-Meier and Aalen-Johansen use the information from patients followed for part of the interval without assuming their outcome. A single pooled estimate across groups with different follow-up can disagree with the group estimates, which is why the overall row is built from the groups.

**How to read it.** Take a group of 10 patients. One has metastasis at 12 months and nine are followed beyond 60 months: the 5-year risk is 10%. If instead one patient is censored at 12 months and one has metastasis at 60 months, the Kaplan-Meier survival at 60 months is 0 (risk 100%), not 0.5, because the patient censored at 12 months says nothing about months 13 to 60 and the only patient still observed at 60 months had the event. An imported risk inside the observed interval means the data are compatible with it; a wide interval means few patients and events. The assumption fails if patients lost to follow-up are systematically sicker or healthier than those who stay.

### Discrimination

**Harrell's C (MFS).** The proportion of comparable patient pairs in which the patient with the earlier event has the higher imported risk, with pairs tied on imported risk counting one half:

$$
C = \frac{\sum_{i,j} \Delta_{ij}\left[\mathbb{1}(r_i > r_j) + \tfrac{1}{2}\,\mathbb{1}(r_i = r_j)\right]}{\sum_{i,j} \Delta_{ij}},
$$

where $\Delta_{ij} = 1$ when patient $i$ has an observed event before patient $j$'s event or censoring time, and $r$ is the imported risk. Follow-up is truncated at 60 months, and the estimate comes from `survival::concordance()`. The interval is a patient-level percentile bootstrap.

**Competing-risk AUC (MSS).** The cumulative/dynamic time-dependent AUC at 60 months from `timeROC` with inverse-probability-of-censoring weights: the probability that a patient with melanoma death by 60 months has a higher imported risk than a patient without it, with patients who died of another cause counted among those without. The interval is a patient-level percentile bootstrap.

**Rules.** Discrimination is estimated only when at least 3 events occur by 60 months and the imported risk varies. Each bootstrap uses a fixed local seed; the number of valid and failed draws is reported, and the interval is withheld when too many draws fail.

**Assumptions.** Independent censoring (as above); for the AUC, censoring weights estimated from the Kaplan-Meier censoring distribution.

**Rationale.** With four distinct imported values, a large share of patient pairs are tied. Counting ties as one half is the standard convention; dropping tied pairs inflates C. The tied-pair share is reported in `Checks`.

**How to read it.** 0.5 is chance and 1 is perfect ordering. If two groups hold all patients and the group with the higher imported risk contains every earlier event, C is 1 across groups, but pairs within a group are tied and count one half, so C sits below 1 by half the tied share. Because the predictor has four values, C describes how well the four groups are ordered, not individual-level prediction. An interval that includes 0.5 means the data do not separate the ordering from chance; with few events the interval is wide.

**Proportional-hazards check.** A Schoenfeld-residual test (`cox.zph`) of a Cox model with the imported risk as its only term (cause-specific for MSS) checks that the effect of the imported risk is constant over follow-up. A p-value below 0.05 is recorded as a concern in `Checks`. A concern would mean the ordering of risk changes with follow-up time, so a single C or AUC averages different periods.

### Sensitivity analyses

Each sensitivity analysis recomputes the group tables and discrimination under one alternative assumption and reports the difference from the primary value.

- **Late biopsy.** The biopsy date is compared with the treatment date, and patients sampled more than 180 days after treatment are removed. A sample taken after treatment (for example at recurrence) describes a different tumor state, and the interval from the GEP sample to the outcome starts after follow-up began. Biopsy timing is classified as more than 31 days before treatment, within 31 days, 31 to 180 days after, more than 180 days after, or no biopsy date. The windows are defined in `objective4_policy.R`.
- **Unknown cause of death.** Deaths with recorded cause `Unknown` or `Unknown_Cancer` are counted as melanoma deaths instead of competing deaths. The primary analysis counts them as competing, so this is the opposite extreme.
- **Treatment era.** GEP testing appears in this cohort only from a given treatment year on (`GEP_TESTING_ERA_START_YEAR` in `objective4_policy.R`), so patients treated earlier could not have been tested. Each group's 5-year risk is recomputed among patients treated in or after that year, with a log-log 95% interval and the difference from the primary risk, so tested and untested groups cover the same era. Patients treated before it stay in every primary table, are left out of this sensitivity only, and are listed by id in `Record_Issues`.
- **Biopsy date review.** A primary-set patient whose biopsy is more than 365 days before treatment, or has no date, is listed in `Record_Issues` for date review and stays in the analysis. A long interval can reflect an earlier diagnostic sample or a date-entry error; the pipeline does not decide which.

**How to read it.** A small difference from the primary value (for example a group risk moving by 1 percentage point) means the assumption does not drive the result. A large difference means the result depends on that assumption, and both values should be reported.

### Patients without a definitive call

Goal 4c applies the same Kaplan-Meier and Aalen-Johansen estimators to every GEP group in the full cohort, including the groups without a Class 1/Class 2 call. A group with no patient followed to 60 months has no 5-year estimate. Follow-up is summarized from treatment with the reverse Kaplan-Meier method, in which the censoring is treated as the event. The baseline table by group shows only the size of a group with fewer than 3 patients, so a published table cannot describe one or two individuals. It also gives the share treated in or after the testing-era start year, because testing coverage varies with treatment era; the treatment-era sensitivity (above) restricts each group to that era.

**How to read it.** A group's risk is comparable with another group's only where both have estimates and overlapping follow-up. A group followed for a median of 2 years cannot show a 5-year risk. Groups differ in who was tested, so differences between "Not tested" and the Class 1 and Class 2 groups reflect testing patterns as well as biology.

### Baseline models without GEP

Goal 4d asks whether baseline clinical features rank 5-year risk better than a single baseline measurement in patients who have no definitive GEP call.

**Predictors.** Age at diagnosis, sex, tumor height, tumor diameter, location, T stage (T1 to T4), internal reflectivity, subretinal fluid, initial vision, and optic nerve involvement. Internal reflectivity "Unknown" is recoded to missing. A predictor is kept only if at least 90% of values are present in every analysis set and each factor level has at least 5 patients; rows with a sparse level are removed once, so every table for a set uses the same patients. A parsimonious set (age, tumor diameter, location, T stage) is analyzed alongside.

**Model.** A ridge-penalized logistic model of the 5-year event, fit on patients from all GEP groups who are eligible for the endpoint:

$$
\hat{\beta} = \arg\min_{\beta_0,\beta}\left\{-\sum_{i=1}^{n} w_i\left[y_i\eta_i - \log\left(1 + e^{\eta_i}\right)\right] + \lambda\sum_{j=1}^{p}\beta_j^2\right\},\quad \eta_i = \beta_0 + x_i^T\beta,
$$

where $y_i$ is the 5-year event indicator and $w_i$ is an inverse-probability-of-censoring weight. A patient with the event by 60 months is a case; a patient followed past 60 months without the event is a control; for MSS, a patient who died of another cause is a control; a patient censored before 60 months has unknown status and weight 0. Weights are the inverse of the Kaplan-Meier probability of remaining uncensored at the patient's event time (or at 60 months), capped at the 99th percentile. Ridge shrinks coefficients toward zero to reduce overfitting when events are few.

**Out-of-fold evaluation.** Repeated nested cross-validation (5 outer folds, 5 inner folds for choosing $\lambda$, 20 repeats, deterministic fold assignment). Censoring weights are estimated inside each training fold, so no weight uses the patients being scored. Every patient receives an out-of-fold prediction in each repeat. AUC and Brier score are computed on those predictions in the **No GEP scope** (patients without a Class 1/Class 2 call); an All-patients scope is shown for reference. The estimate is the mean over repeats. The interval is a patient-level bootstrap; each draw resamples patients and uses one randomly chosen repeat, so it reflects patient sampling and fold assignment, and it does not refit the model on bootstrap samples.

**Comparators.** T stage (1 to 4) and tumor diameter are used directly as the risk score with no fitting, on the same patients and weights. A model that does not exceed them has not added information beyond one baseline measurement.

**Brier skill.** $1 - B_{\text{model}} / B_{\text{null}}$, where $B_{\text{null}}$ is the weighted Brier score of predicting the overall weighted prevalence for everyone. A value of 0 or below means the predicted risks are no better than that single number; skill under 0.05 is flagged.

**Paired AUC difference.** Model AUC minus the comparator's AUC, with both computed on the same resampled patients and the same partition repeat in each bootstrap draw, so the interval is for the difference itself.

**Class 2-like score.** A ridge logistic model of Class 2 versus Class 1 among patients with a definitive call, scored out of fold, with the AUC and a bootstrap interval. It describes how much baseline features resemble the Class 2 pattern; it is not the probability of a molecular Class 2 result.

**Rules.** Events per design column below 10 are flagged. A failed analysis stops the run.

**How to read it.** An AUC of 0.5 is chance. If the model's AUC is 0.66 and T stage alone gives 0.69, the model does not rank better than T stage; a difference interval that includes 0 says the data cannot distinguish them. Because the bootstrap does not refit the model, the interval understates the uncertainty from model fitting; with few events it is wide regardless. A Brier skill near 0 says the model's predicted risks are no closer to the outcomes than the overall rate.

### Diagnostics workbooks and record-level issues

Each cohort's analysis writes two workbooks. The **primary workbook** holds the tables and figures' source values and contains no patient identifiers. The **diagnostics workbook** holds the evidence for how each result was produced: `Specification_Coding` (method, estimand, units, direction and reading of each quantity), `Denominators_Flow` (patients through each step, reconciling with the cohort), `Checks` (each check with a result and a status of no concern, concern or not assessed, with a reason), `Warnings_Limitations`, and `Record_Issues`. The no-GEP workbook also has `Predictor_Screening`.

`Record_Issues` lists, one row per record and issue, the patients outside the primary set (by group), primary-set patients flagged for biopsy timing, and patients removed from a model, with the step, reason, field and current value. Each row carries the patient `id` used in the analytic dataset, so every left-out or flagged patient can be traced; the counts reconcile with `Denominators_Flow`. The diagnostics workbooks are part of the results package (published snapshot and collaborator export), because the results cannot be read without knowing which patients were left out or flagged. They contain identifiers, so they are never committed to Git or posted publicly and go only to people who already have access to the study data.

A check that cannot be run is recorded as "not assessed" with its reason, so it never appears to have passed.

---

## Statistical Metric Categorization

The analysis pipeline categorizes statistical outputs into three tiers based on clinical interpretability and robustness.

### Primary Metrics (Tier 1)

**Definition:** Model-based estimates with direct clinical interpretation

**Examples:**
- Hazard ratios from Cox regression
- Odds ratios from logistic regression
- Regression coefficients from linear models
- RMST differences
- Risk differences

**Characteristics:**
- Adjusted for confounders
- Direct effect estimates
- Primary basis for clinical conclusions

### Secondary Metrics (Tier 2)

**Definition:** Descriptive statistics and unadjusted comparisons

**Examples:**
- Median survival times
- Event rates by treatment group
- Mean changes in continuous outcomes
- Survival probabilities at specific time points

**Characteristics:**
- Unadjusted for confounders
- Descriptive summaries
- Complement primary metrics

### Tertiary Metrics (Tier 3)

**Definition:** Technical diagnostics and quality metrics

**Examples:**
- P-values from proportional hazards tests
- Model fit statistics (AIC, BIC, R²)
- Residual diagnostics
- Calibration metrics
- Number of observations/events

**Characteristics:**
- Assess model assumptions and quality
- Guide interpretation of primary metrics
- Technical validation

### Reporting Priority

**Clinical Reports:**
1. Primary metrics (adjusted effect estimates)
2. Secondary metrics (descriptive context)
3. Tertiary metrics (footnotes on model quality)

**Statistical Appendices:**
1. All three tiers comprehensively documented
2. Full model diagnostics
3. Sensitivity analyses

**Forest Plots:**
- Display: Primary metrics (effect estimates with CI)
- Annotate: Sample sizes, event counts (secondary metrics)
- Document: Model diagnostics in separate tables (tertiary metrics)

This categorization ensures clinicians focus on actionable findings while maintaining statistical rigor and transparency.

---

## References

Key statistical methods references:

1. **Cox Regression:** Cox DR. Regression models and life-tables. J R Stat Soc Series B. 1972;34(2):187-220.

2. **RMST:** Uno H, et al. Moving beyond the hazard ratio in quantifying the between-group difference in survival analysis. J Clin Oncol. 2014;32(22):2380-2385.

3. **Proportional Hazards Testing:** Grambsch PM, Therneau TM. Proportional hazards tests and diagnostics based on weighted residuals. Biometrika. 1994;81(3):515-526.

4. **Competing Risks:** Fine JP, Gray RJ. A proportional hazards model for the subdistribution of a competing risk. J Am Stat Assoc. 1999;94(446):496-509.

5. **Prediction Model Validation:** Steyerberg EW, et al. Assessing the performance of prediction models. Epidemiology. 2010;21(1):128-138.

See README.md for complete citation of the analysis pipeline and how to cite this work.
