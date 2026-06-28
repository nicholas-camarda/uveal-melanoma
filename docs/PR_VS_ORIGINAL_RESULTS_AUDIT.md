# PR Versus Original Results Audit

This document records the reviewer-facing differences between the original analysis surface and the peer-review statistical-revision surface. It is an internal audit and reporting guide, not manuscript text.

Reviewer-facing scope is limited to the full and restricted cohorts. GKSRS-only runtime outputs may exist for other project purposes, but they are not part of this paper response.

---

## 1. Purpose

The peer-review revision changed several analytic choices in response to reviewer/statistical comments. The main changes were:

- local recurrence and metastatic progression are now led by Cox time-to-event models rather than binary event-rate inference;
- PFS was corrected to include metastatic progression, local recurrence, or death from any cause;
- age is modeled continuously in reviewer-response adjusted Cox models;
- follow-up duration is reported by treatment arm;
- visual-acuity and tumor-height endpoints are caveated for non-standardized follow-up timing;
- reviewer-facing language avoids equivalence, noninferiority, and simulated-randomization claims.

The goal of this audit is to explain why revised numbers differ from older generated results and how to report them without implying a data error.

---

## 2. PFS Definition Change

### Original Definition

The older implementation treated PFS as local recurrence or death. Metastatic progression was unintentionally omitted from the PFS composite.

### Revised Definition

The peer-review revision defines PFS as time from primary treatment to first local recurrence, metastatic progression, or death from any cause:

```r
pfs_event = recurrence_event == 1 | mets_event == 1 | death_event == 1
```

This is implemented in [data_derivation.R](/Users/ncamarda/Projects/uveal_melanoma-peer-review-statistical-revision/scripts/data_helper/data_derivation.R).

### Event-Count Impact

The correction adds metastatic-progression-only PFS events.

| Cohort/layer | Revised PFS events | Older recurrence-or-death PFS events | Difference |
| --- | ---: | ---: | ---: |
| Full cohort, dataset/KM layer | 82 | 78 | +4 |
| Restricted cohort, dataset/KM layer | 52 | 50 | +2 |
| Full cohort, Cox modeled layer | 79 | 75 | +4 |
| Restricted cohort, Cox modeled layer | 51 | 49 | +2 |

The full-cohort added Cox events are patients 22, 92, 215, and 249. The restricted-cohort added Cox events are patients 215 and 249.

### Reporting Interpretation

This is an intentional correction. The response should state the current definition directly and should not describe the old recurrence-or-death-only definition as the current PFS endpoint.

---

## 3. Source Cohort N Versus Cox Model N

The source cohort denominators are:

- full cohort: 260 patients;
- restricted cohort: 167 patients.

The adjusted Cox models use smaller modeled N because sparse `location` levels are excluded before model fitting. This is not ordinary missing-covariate complete-case loss and not a PFS-definition problem.

| Cohort | Source N | Cox modeled N | Excluded patient IDs | Reason |
| --- | ---: | ---: | --- | --- |
| Full | 260 | 254 | 24, 30, 78, 90, 101, 247 | Sparse `location` levels |
| Restricted | 167 | 164 | 24, 30, 78 | Sparse `location` levels |

Full-cohort sparse-location exclusions remove 2 local-recurrence events, 1 metastatic-progression event, 1 death event, and 3 PFS events from the Cox modeled layer. Restricted-cohort sparse-location exclusions remove 1 local-recurrence event and 1 PFS event from the Cox modeled layer.

### Reporting Rule

Use precise denominator labels:

- `source cohort n` for dataset summaries and follow-up summaries;
- `field availability n` for variable-specific data availability;
- `Cox model n` for regression tables and Cox model event counts.

Do not write that the Cox N dropped from 260 to 254 or 167 to 164 because of missing covariates.

---

## 4. Objective 1 Model Differences

### Recurrence And Metastatic Progression

The original result surface emphasized binary/logistic event-rate comparisons. The peer-review response should lead with adjusted Cox models because the reviewer specifically objected that crude event rates ignore unequal follow-up.

The modeled N and event counts are unchanged for the recurrence/metastasis Cox comparisons after the current sparse-location modeling rule:

| Cohort | Endpoint | Cox model N | Cox events |
| --- | --- | ---: | ---: |
| Full | Local recurrence | 254 | 31 |
| Full | Metastatic progression | 254 | 43 |
| Restricted | Local recurrence | 164 | 18 |
| Restricted | Metastatic progression | 164 | 28 |

The response should say these were reanalyzed as time-to-event endpoints using adjusted Cox proportional-hazards models. Binary event-rate results can be descriptive only.

### OS Hazard-Ratio Shift

The reviewer-response adjusted Cox models use continuous age rather than the prior median-split age variable. This changes the OS adjustment surface while keeping the same modeled N and death-event count.

Observed shift from the recovered handoff:

| Cohort | Older OS HR | Revised OS HR | Main reason |
| --- | ---: | ---: | --- |
| Full | 1.12 | 1.20 | Continuous age replaces median-split age |
| Restricted | 1.10 | 1.07 | Continuous age replaces median-split age |

This is a modeling-policy change made to address the reviewer request not to dichotomize continuous variables.

---

## 5. Follow-Up And Timing Findings

### General Follow-Up

Treatment-arm follow-up summaries are available in:

- `<OUTPUT_DIR>/peer_review_revision_audits/full_followup_and_data_availability.xlsx`
- `<OUTPUT_DIR>/peer_review_revision_audits/restricted_followup_and_data_availability.xlsx`

General follow-up is longer after PBT than GKSRS.

| Cohort | Treatment | Source N | Median follow-up, months | Mean follow-up, months |
| --- | --- | ---: | ---: | ---: |
| Full | PBT | 121 | 75.6 | 83.2 |
| Full | GKSRS | 139 | 37.2 | 51.1 |
| Restricted | PBT | 103 | 75.2 | 81.5 |
| Restricted | GKSRS | 64 | 38.0 | 53.7 |

### Latest-VA Timing

Latest-VA timing uses treatment date to `last_followup` as the available timing field associated with `last_vision`.

| Cohort | Treatment | Cohort N | Explicit timing available | Median latest-VA timing, months |
| --- | --- | ---: | ---: | ---: |
| Full | PBT | 121 | 96 | 77.5 |
| Full | GKSRS | 139 | 117 | 32.7 |
| Restricted | PBT | 103 | 80 | 78.3 |
| Restricted | GKSRS | 64 | 53 | 35.3 |

These are not fixed 1-, 3-, or 5-year visual-acuity landmark measurements.

### Latest-VA Predictor-Adjusted Sensitivity

The reviewer asked whether visual-acuity differences could be adjusted for clinically relevant predictors, including baseline VA, tumor geometry/stage, subretinal fluid, optic-nerve proximity, follow-up duration, and radiation dose/proximity fields when available. The implemented sensitivity model uses latest logMAR VA (`last_vision`) as the outcome and adjusts for baseline logMAR VA rather than adding baseline VA to the change-score outcome.

Full-cohort model:

`last_vision ~ treatment_group + initial_vision + last_vision_followup_months_explicit + initial_tumor_height + initial_tumor_diameter + initial_t_stage_simple + srf + optic_nerve + treatment_year_centered + age_at_diagnosis + sex + location`

Restricted-cohort model:

`last_vision ~ treatment_group + initial_vision + last_vision_followup_months_explicit + initial_tumor_height + initial_tumor_diameter + initial_t_stage_simple + srf + treatment_year_centered + age_at_diagnosis + sex + location`

Variable handling:

| Candidate variable | Used in latest-VA sensitivity? | Reason |
| --- | --- | --- |
| `age_at_diagnosis` | Yes, full and restricted | Standard confounder set. |
| `sex` | Yes, full and restricted | Standard confounder set. |
| `location` | Yes, full and restricted | Standard confounder set. |
| `initial_tumor_height` | Yes, full and restricted | Available baseline tumor-size predictor. |
| `initial_tumor_diameter` | Yes, full and restricted | Available baseline basal-diameter predictor. |
| `initial_vision` | Yes, full and restricted | Baseline logMAR VA; adjusted in the latest-VA model rather than in the change-score model. |
| `initial_t_stage_simple` | Yes, full and restricted | Available baseline T-stage grouping with adequate model support. |
| `initial_overall_stage` | No | Not added with `initial_t_stage_simple`; including both would duplicate tumor-stage information and add model burden. |
| `optic_nerve` | Yes in full; no in restricted | Full cohort had usable variation. Restricted cohort had insufficient variation. |
| `srf` | Yes, full and restricted | Available baseline subretinal-fluid field. |
| `visual_field_defect` | No | Not used as an adjustment covariate; the analysis should not be described as a longitudinal visual-field analysis. |
| `vision_loss_blurred_vision` | No | Not used as an adjustment covariate; it is a symptom/presentation descriptor rather than the measured baseline VA adjustment. |

Results:

| Cohort | Adjusted GKSRS-vs-PBT MD, logMAR | 95% CI | p-value | Model n |
| --- | ---: | --- | ---: | ---: |
| Full | -0.200 | -0.462 to 0.063 | 0.135 | 208 |
| Restricted | -0.126 | -0.429 to 0.177 | 0.411 | 126 |

Interpretation: the predictor-adjusted latest-VA sensitivity attenuated the visual-acuity association and was not statistically significant in either cohort. It supports cautious visual-acuity wording and does not support a strong claim that GKSRS provides superior functional preservation after adjustment for available baseline predictors and follow-up duration.

### Tumor-Height Timing

Tumor-height change uses latest available height, not a standardized imaging landmark. Timing differs by treatment group.

| Cohort | Treatment | Cohort N | Height timing available | Median timing, months |
| --- | --- | ---: | ---: | ---: |
| Full | PBT | 121 | 119 | 53.9 |
| Full | GKSRS | 139 | 138 | 26.3 |
| Restricted | PBT | 103 | 101 | 51.8 |
| Restricted | GKSRS | 64 | 64 | 31.9 |

Negative treatment-to-height intervals were identified in the timing workbooks:

- full cohort: patients 1, 2, 80, 95, 166, 174, and 176;
- restricted cohort: patients 1, 80, 166, 174, and 176.

### Reporting Interpretation

Tumor-height change should be descriptive/secondary with a timing caveat. The comparative tumor-height regression should not be presented as a definitive treatment-effect endpoint.

---

## 6. Reporting Guide For Manuscript And Reviewer Response

### PFS Definition

Use:

> PFS was defined as time from primary treatment to first local recurrence, metastatic progression, or death from any cause.

Avoid:

> PFS was recurrence or death.

### Recurrence And Metastasis

Use:

> Local recurrence and metastatic progression were reanalyzed as time-to-event endpoints using adjusted Cox proportional-hazards models to account for unequal follow-up.

Avoid making crude binary event-rate comparisons the lead inference.

### Cox Denominators

Use:

> The source cohort included 260 full-cohort and 167 restricted-cohort patients. Adjusted Cox models included 254 and 164 patients, respectively, because sparse `location` levels were excluded before model fitting.

Avoid:

> Six full-cohort patients and three restricted-cohort patients were dropped for missing data.

### Follow-Up

Use:

> Median follow-up was longer after PBT than GKSRS in both cohorts, supporting the shift from crude event-rate comparisons to time-to-event models.

### Visual Acuity

Use:

> Visual-acuity sensitivity analyses were repeated after restricting to patients with at least 12, 36, and 60 months of latest-VA follow-up. These are minimum-follow-up restrictions, not standardized 1-, 3-, or 5-year landmark VA analyses.

Also state:

> Standardized visual-acuity landmarks and time to clinically meaningful deterioration, such as loss of at least 3 Snellen lines, would require additional longitudinal data collection and were not feasible within the revision timeframe.

### Tumor Height

Use:

> Tumor-height change was based on latest available follow-up height rather than a fixed imaging landmark. Because treatment-to-height-assessment timing differed by treatment arm, tumor-height change is interpreted descriptively/secondarily.

### Toxicity

Use:

> Adverse events were analyzed as recorded toxicity burden by available follow-up. The curated dataset did not contain standardized toxicity grades or dated onset fields for each adverse event.

Do not add unsupported toxicity endpoints to this reviewer response. The checked source/derived fields support SRD.

### Propensity Score

Use:

> Propensity-score feasibility was assessed using available baseline demographic, tumor, ocular, and symptom variables. Full-cohort propensity adjustment was not supportable because of poor overlap/separation. The restricted cohort had acceptable measured-variable overlap, which supports the restricted-cohort comparability rationale but does not simulate randomization.

Do not present the feasibility workbooks as treatment-effect analyses.

### Language Guardrails

Use:

> No statistically significant difference was observed.

Avoid:

> Equivalent, noninferior, randomized, randomized-like, or simulated randomization.

---

## 7. Current Evidence Pointers

- [peer_review_revision_response.md](/Users/ncamarda/Projects/uveal_melanoma-peer-review-statistical-revision/docs/peer_review_revision_response.md)
- [statistical revision plan](/Users/ncamarda/Projects/uveal_melanoma-peer-review-statistical-revision/docs/superpowers/plans/2026-06-26-peer-review-statistical-revision.md)
- [PFS derivation code](/Users/ncamarda/Projects/uveal_melanoma-peer-review-statistical-revision/scripts/data_helper/data_derivation.R)
- [KM display cap config](/Users/ncamarda/Projects/uveal_melanoma-peer-review-statistical-revision/scripts/config/labels_display.R)
- `<OUTPUT_DIR>/peer_review_revision_audits/full_followup_and_data_availability.xlsx`
- `<OUTPUT_DIR>/peer_review_revision_audits/restricted_followup_and_data_availability.xlsx`
- `<OUTPUT_DIR>/peer_review_revision_audits/full_propensity_score_feasibility.xlsx`
- `<OUTPUT_DIR>/peer_review_revision_audits/restricted_propensity_score_feasibility.xlsx`
- `<OUTPUT_DIR>/uveal_full/01_Efficacy/e_tumor_height_primary/full_cohort_tumor_height_timing_summary.xlsx`
- `<OUTPUT_DIR>/uveal_restricted/01_Efficacy/e_tumor_height_primary/restricted_cohort_tumor_height_timing_summary.xlsx`

---

## 8. Unresolved Coauthor Choices

- Tim/Angie choose which VA minimum-follow-up cutoff, if any, to emphasize.
- Tim/Angie choose whether tumor-height change remains in the manuscript as a secondary/descriptive result or is removed from comparative emphasis.
- Tim owns final Table 1 p-value removal/formatting; Nick does not need to change generated working tables.
- Tim/Melhus/Mignano own final clinical technique prose for GK/PBT details.
