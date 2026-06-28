# Peer Review Revision Response And Results

This is the single coauthor-facing response document for Nick's work on Tim-directed reviewer revisions. It is organized around Tim's Google Doc comments and to-dos for Nick, not around internal analysis modules.

Each item below should be updated as analyses are finalized. Runtime artifacts stay in the peer-review runtime root and are cited here with aliases such as `<OUTPUT_DIR>/...` and `<PROCESSED_DATA_DIR>/...`; do not paste machine-specific absolute paths into this committed document.

Each T-NICK item must answer Tim's or the reviewer's actual question in the **Direct answer** field before listing methods, results, or artifacts. Evidence paths support the answer; they are not a substitute for the answer.

Reviewer-facing scope is limited to the full and restricted cohorts. Internal GKSRS-only runtime outputs may exist for other project purposes, but they should not be included in this peer-review response document or in manuscript-facing reviewer answers for this paper.

**Current posture:** This document is the coauthor-facing tracker for Nick's statistical work. Items marked `Complete` are complete for Nick's statistical response unless they explicitly identify a clinical-coauthor or Tim-owned manuscript wording decision.

---

## Source Packet

- Editorial decision email dated June 21, 2026: resubmission due August 20, 2026.
- Tim email dated June 22, 2026: major revisions required; both reviewers concerned about radiation-treatment details and statistical methods; some treatment-detail requests may be infeasible within the response window.
- Angie email dated June 23, 2026: most statistical comments make sense; meet or resolve by email before final response.
- Shared Google Doc dated June 24, 2026: Tim asked Nick to review and address the statistics comments/questions and Angie's suggestions in the comments.
- Statistical planning memo: de-escalate and refocus the analysis around fewer, stronger models; make follow-up central; weaken causal/equivalence language.

---

## Recorded Decisions

| Item | Decision |
| --- | --- |
| T-NICK-15 | Tim removes final Table 1 p-values manually; no Nick-side pipeline change. Generated working tables can keep p-values for now. |
| T-NICK-12 | Nick's statistical dose/proximity limitation wording is complete. Final GK/PBT technique prose remains Tim/Melhus/Mignano-owned. |
| T-NICK-13 | Use recorded-burden toxicity/SRD wording. Do not add unsupported toxicity endpoints to this reviewer response. |
| T-NICK-14 | Nick's statistical language guardrails are documented. Tim owns final manuscript wording and this is not a Nick blocker. |
| RMST | Omit RMST from the reviewer response letter. |
| Propensity score | Treat as feasibility only. Do not add a PS treatment-effect sensitivity unless Tim/Angie explicitly request the extra work. |
| Visual acuity | All 12-, 36-, and 60-month minimum-follow-up sensitivity outputs are available for full/restricted cohorts. Tim/Angie choose manuscript emphasis. |
| Tumor height | Treat tumor-height change as secondary/descriptive with timing caveat. |
| KM x-axis | KM display cap is 216 months (18 years). This is display-only and does not affect Cox models. |

---

## Status Labels

- `Not started`: no implementation or runtime evidence yet.
- `Implemented pending verification`: code/doc work exists but final runtime artifacts or interpretation still need review.
- `Pending runtime result`: method is defined, but result/evidence must be filled from regenerated outputs.
- `Decision needed`: Nick/Tim/Angie need to choose whether to report, demote, or omit the analysis.
- `Clinical coauthor input needed`: statistical response depends on clinical technique or data-collection confirmation.
- `Complete`: result, interpretation, limitation, exact evidence aliases, and final checklist validation are ready for Tim review.

---

## Adjustment Set Definition

The reviewer-response adjusted Cox models use a deliberately parsimonious adjustment set: treatment group plus continuous age at diagnosis, sex, and tumor location when those covariates are estimable. This is the "reviewer-response covariate set" referenced below. It is not a propensity-score model and does not include every available baseline variable; it is the small conventional multivariable Cox adjustment set used for the primary reviewer-response time-to-event analyses.

### Source Cohort N Versus Cox Model N

The source cohort denominators are 260 for the full cohort and 167 for the restricted cohort. Follow-up summaries, data-availability summaries, and dataset-level endpoint checks use those denominators unless a field-specific availability count is stated.

The adjusted Cox model denominators are smaller because the regression-table pipeline excludes sparse categorical levels from model variables before fitting. This is not missing-data complete-case loss: in the current peer-review runtime RDS files, `treatment_group`, `age_at_diagnosis`, `sex`, `location`, and the Objective 1 time/event variables are complete for all full and restricted cohort rows. The reduction comes from sparse `location` levels under the model-layer threshold of fewer than 5 observations.

Full cohort Cox models: 260 source rows enter the model-prep step, and 254 rows are modeled. The 6 excluded rows are the sparse location levels `Ciliary Body` (ids 24 and 101), `Conjunctival` (id 30), `Irido-Ciliary` (ids 78 and 90), and `Iris` (id 247). These exclusions remove 2 local-recurrence events, 1 metastatic-progression event, 1 death event, and 3 PFS events. Therefore full-cohort dataset-level events are recurrence 33, metastasis 44, death 57, and PFS 82, while full-cohort Cox model event counts are recurrence 31, metastasis 43, death 56, and PFS 79. For the 60-month sensitivity models, the same 254 rows are modeled; events within 60 months are OS 38 and PFS 65.

Restricted cohort Cox models: 167 source rows enter the model-prep step, and 164 rows are modeled. The 3 excluded rows are the sparse location levels `Ciliary Body` (id 24), `Conjunctival` (id 30), and `Irido-Ciliary` (id 78). All three are PBT patients. These exclusions remove 1 local-recurrence event, 0 metastatic-progression events, 0 death events, and 1 PFS event. Therefore restricted-cohort dataset-level events are recurrence 19, metastasis 28, death 39, and PFS 52, while restricted-cohort Cox model event counts are recurrence 18, metastasis 28, death 39, and PFS 51. For the 60-month sensitivity models, the same 164 rows are modeled; events within 60 months are OS 27 and PFS 41.

Reviewer-facing text should therefore label denominators precisely as `source cohort n`, `field availability n`, or `Cox model n`. Do not describe the 260-to-254 or 167-to-164 reduction as ordinary missing covariate complete-case filtering.

---

## Tim-Directed Nick To-Do Index

### Status At A Glance

| Task | Topic | Current status |
| --- | --- | --- |
| T-NICK-01 | Cox-only recurrence/metastasis reanalysis | Complete |
| T-NICK-02 | Proportional-hazards diagnostics | Complete |
| T-NICK-03 | Continuous variables and dichotomization | Complete |
| T-NICK-04 | PFS definition and death/progression components | Complete |
| T-NICK-05 | Separate recurrence/metastasis outputs | Complete |
| T-NICK-06 | Five-year limited-follow-up OS/PFS Cox sensitivity | Complete |
| T-NICK-07 | Propensity-score feasibility | Complete: feasibility only |
| T-NICK-08 | Visual-acuity follow-up and minimum-follow-up sensitivity | Complete; Tim/Angie choose emphasis |
| T-NICK-09 | Follow-up duration summaries | Complete |
| T-NICK-10 | Tumor-height timing | Complete; demote/caveat comparative tumor-height change |
| T-NICK-11 | PRAME and T4 subgroup pruning | Complete |
| T-NICK-12 | Dose/proximity/treatment-detail availability | Complete for Nick's statistical limitation; clinical technique prose remains coauthor-owned |
| T-NICK-13 | Adverse-event grading and SRD scope | Complete for checked source/derived fields |
| T-NICK-14 | Statistical language guardrails | Tim final wording pass; not a Nick blocker |
| T-NICK-15 | Table 1 p-value handling | Complete: Tim-owned final table edit; no pipeline change |

---

## Reviewer Statistics Block Coverage

This maps the reviewer/statistical comment block beginning with the concern about binary outcomes and ending with the instruction to leave continuous variables continuous.

| Reviewer issue | Direct answer | Where handled |
| --- | --- | --- |
| Binary recurrence/metastasis outcomes ignore follow-up duration. | Replace crude binary inference with time-to-event Cox models for time to local recurrence and time to metastatic progression; retain event rates only descriptively. | T-NICK-01, T-NICK-05 |
| Cox models require PH assumption checks. | Use scaled Schoenfeld residual tests via `cox.zph()` for every full/restricted Cox model shown or cited in the reviewer-facing response; report treatment-term and global PH results, and qualify any model with a PH concern. | T-NICK-02 |
| Follow-up imbalance affects recurrence, metastasis, PFS, VA, and toxicity interpretation. | Report follow-up explicitly and avoid crude event-rate or late-tail claims without risk-set support. Treatment-arm follow-up summaries are now available for the full and restricted cohorts. | T-NICK-06, T-NICK-08, T-NICK-09, T-NICK-13 |
| Continuous variables should not be dichotomized. | Use continuous age at diagnosis in reviewer-response adjusted Cox models; do not use dichotomized age as the primary adjustment variable. | T-NICK-03 |

---

### T-NICK-01: Cox-Only Local Recurrence And Metastatic Progression Reanalysis

> **STATUS:** Complete.

**Direct answer:** Working answer: yes. The recurrence and metastatic-progression comparisons have been rerun as adjusted Cox time-to-event analyses, which is the appropriate response to the reviewer concern that crude event rates ignore unequal follow-up. These models did not show a statistically significant treatment-group difference for local recurrence or metastatic progression.

**Google Doc context:** Reviewer objected that local recurrence and metastasis were treated as binary outcomes despite unequal follow-up. Tim's comment to Nick, with Angie's agreement, was to rerun these data for Cox regression analyses only.

**Owner:** Nick.

**Action required:** Make Cox proportional-hazards models the lead inferential analyses for time to local recurrence and time to metastatic progression.

**Method / implementation:** Use time from primary treatment to first local recurrence and time from primary treatment to metastatic progression, adjusted with the reviewer-response covariate set defined above. Logistic/Fisher results should not be presented as lead inference.

**Result:** Reviewer-response Cox outputs were regenerated for the full and restricted cohorts, with exact model tables and PH summaries present in the runtime tree. In the full cohort, adjusted Cox models estimated GKSRS versus PBT HR 1.36 for local recurrence (95% CI 0.669-2.78; p=0.394; Cox model n=254, 31 events) and HR 1.25 for metastatic progression (95% CI 0.679-2.30; p=0.472; Cox model n=254, 43 events). In the restricted cohort, adjusted Cox models estimated HR 1.72 for local recurrence (95% CI 0.674-4.41; p=0.256; Cox model n=164, 18 events) and HR 1.00 for metastatic progression (95% CI 0.446-2.24; p=0.999; Cox model n=164, 28 events). The Cox model n values are lower than the source cohort denominators because sparse `location` levels are excluded before modeling, as detailed above.

**Response/manuscript language:** State that local recurrence and metastatic progression were reanalyzed as time-to-event endpoints using adjusted Cox proportional-hazards models to account for unequal follow-up. Do not present binary logistic or Fisher tests as lead inference for these endpoints.

**Limitation / interpretation:** Cox estimates remain observational and limited by event counts, treatment selection, and available covariates. These results support "no statistically significant difference was observed" language, not equivalence or noninferiority.

**Evidence:**

- `<OUTPUT_DIR>/uveal_full/01_Efficacy/a_recurrence/full_cohort_local_recurrence_free_probability_cox_coxph.html`
- `<OUTPUT_DIR>/uveal_full/01_Efficacy/a_recurrence/full_cohort_local_recurrence_free_probability_proportional_hazards_summary.txt`
- `<OUTPUT_DIR>/uveal_full/01_Efficacy/b_metastatic_progression/full_cohort_metastasis_free_survival_probability_cox_coxph.html`
- `<OUTPUT_DIR>/uveal_full/01_Efficacy/b_metastatic_progression/full_cohort_metastasis_free_survival_probability_proportional_hazards_summary.txt`
- `<OUTPUT_DIR>/uveal_restricted/01_Efficacy/a_recurrence/restricted_cohort_local_recurrence_free_probability_cox_coxph.html`
- `<OUTPUT_DIR>/uveal_restricted/01_Efficacy/a_recurrence/restricted_cohort_local_recurrence_free_probability_proportional_hazards_summary.txt`
- `<OUTPUT_DIR>/uveal_restricted/01_Efficacy/b_metastatic_progression/restricted_cohort_metastasis_free_survival_probability_cox_coxph.html`
- `<OUTPUT_DIR>/uveal_restricted/01_Efficacy/b_metastatic_progression/restricted_cohort_metastasis_free_survival_probability_proportional_hazards_summary.txt`

---

### T-NICK-02: Proportional-Hazards Assumption Confirmation

> **STATUS:** Complete.

**Direct answer:** Working answer: yes. The Cox PH assumption was assessed with scaled Schoenfeld residual tests using `cox.zph()` for the full-cohort and restricted-cohort Cox models shown or cited for Objective 1. Objective 2 full/restricted outputs currently use linear or logistic models for vision/toxicity endpoints, not Cox models, so there are no Objective 2 PH assumptions to report. For Objective 1, the treatment-term PH tests did not show evidence against PH for any shown model. Covariate-level PH diagnostics flagged two secondary caveats: full-cohort 60-month OS global PH concern driven by age at diagnosis (treatment term OK), and restricted-cohort full-follow-up OS `location` covariate marginal PH concern (p=0.0485; treatment term and global PH OK).

**Google Doc context:** Reviewer asked whether proportional-hazards assumptions were checked. Tim asked Nick to confirm/report this.

**Owner:** Nick.

**Action required:** Report PH diagnostics for every full/restricted Cox model shown or cited in the reviewer-facing Objective 1 response, with treatment-term p-value/status, global p-value/status, and any caution note.

**Method / implementation:** Use scaled Schoenfeld residual testing through `cox.zph()` and record global and covariate-level diagnostic status or an explicit skip reason when a model is not estimable.

**Result:** The consolidated PH inventory for the full/restricted Objective 1 Cox models is below. It includes the Cox models shown or cited in the reviewer-facing response: time to local recurrence, time to metastatic progression, OS, PFS, and the 60-month OS/PFS sensitivity analyses. Post-baseline recurrence-stratified or metastasis-stratified OS/PFS analyses and GKSRS-only runtime outputs are not part of this paper response.

| Cohort | Cox model | Follow-up horizon | Treatment-term PH p/status | Global PH p/status | Reporting implication |
| --- | --- | --- | --- | --- | --- |
| Full | Local recurrence | Full follow-up | p=0.0967; OK | p=0.111; OK | Reportable with standard PH statement. |
| Full | Metastatic progression | Full follow-up | p=0.510; OK | p=0.964; OK | Reportable with standard PH statement. |
| Full | OS | Full follow-up | p=0.204; OK | p=0.100; OK | Reportable with standard PH statement. |
| Full | PFS | Full follow-up | p=0.123; OK | p=0.104; OK | Reportable with standard PH statement. |
| Full | OS | 60-month limited follow-up | p=0.606; OK | p=0.0464; concern | Report only as secondary sensitivity with PH caveat. |
| Full | PFS | 60-month limited follow-up | p=0.177; OK | p=0.150; OK | Reportable as 60-month sensitivity. |
| Restricted | Local recurrence | Full follow-up | p=0.272; OK | p=0.117; OK | Reportable with standard PH statement. |
| Restricted | Metastatic progression | Full follow-up | p=0.679; OK | p=0.786; OK | Reportable with standard PH statement. |
| Restricted | OS | Full follow-up | p=0.218; OK | p=0.198; OK | Treatment-term PH OK. `location` covariate PH marginal (p=0.0485); report with standard treatment-term PH statement. |
| Restricted | PFS | Full follow-up | p=0.514; OK | p=0.116; OK | Reportable with standard PH statement. |
| Restricted | OS | 60-month limited follow-up | p=0.861; OK | p=0.313; OK | Reportable as 60-month sensitivity. |
| Restricted | PFS | 60-month limited follow-up | p=0.909; OK | p=0.334; OK | Reportable as 60-month sensitivity. |

**Response/manuscript language:** "For Cox models shown in the revised analysis, the proportional-hazards assumption was assessed using scaled Schoenfeld residual tests. The treatment-effect terms did not show evidence against proportional hazards in the reported full-cohort or restricted-cohort treatment-comparison models. The full-cohort 60-month OS sensitivity had a global PH diagnostic concern driven by age at diagnosis, so that result is treated as secondary and caveated rather than unqualified lead evidence. In the restricted-cohort full-follow-up OS model, the `location` covariate had a marginal PH diagnostic (p=0.0485); the treatment-term PH test remained acceptable."

**Limitation / interpretation:** A non-significant PH test does not prove the model is correct; sparse events reduce diagnostic power. A significant PH diagnostic on an adjustment covariate does not automatically invalidate the treatment HR, but it should make the 60-month OS sensitivity secondary and qualified. The restricted-cohort OS `location` covariate PH flag is a covariate-level diagnostic only; the treatment-term PH test remained acceptable.

**Evidence:**

- `<OUTPUT_DIR>/uveal_full/01_Efficacy/a_recurrence/full_cohort_local_recurrence_free_probability_proportional_hazards_summary.txt`
- `<OUTPUT_DIR>/uveal_full/01_Efficacy/b_metastatic_progression/full_cohort_metastasis_free_survival_probability_proportional_hazards_summary.txt`
- `<OUTPUT_DIR>/uveal_full/01_Efficacy/h_proportional_hazards_diagnostics/full_cohort_overall_survival_probability_proportional_hazards_summary.txt`
- `<OUTPUT_DIR>/uveal_full/01_Efficacy/h_proportional_hazards_diagnostics/full_cohort_progression_free_survival_probability_proportional_hazards_summary.txt`
- `<OUTPUT_DIR>/uveal_full/01_Efficacy/c_overall_survival/full_cohort_overall_survival_probability_5yr_capped_proportional_hazards_summary.txt`
- `<OUTPUT_DIR>/uveal_full/01_Efficacy/d_progression_free_survival/full_cohort_progression_free_survival_probability_5yr_capped_proportional_hazards_summary.txt`
- `<OUTPUT_DIR>/uveal_restricted/01_Efficacy/a_recurrence/restricted_cohort_local_recurrence_free_probability_proportional_hazards_summary.txt`
- `<OUTPUT_DIR>/uveal_restricted/01_Efficacy/b_metastatic_progression/restricted_cohort_metastasis_free_survival_probability_proportional_hazards_summary.txt`
- `<OUTPUT_DIR>/uveal_restricted/01_Efficacy/h_proportional_hazards_diagnostics/restricted_cohort_overall_survival_probability_proportional_hazards_summary.txt`
- `<OUTPUT_DIR>/uveal_restricted/01_Efficacy/h_proportional_hazards_diagnostics/restricted_cohort_progression_free_survival_probability_proportional_hazards_summary.txt`
- `<OUTPUT_DIR>/uveal_restricted/01_Efficacy/c_overall_survival/restricted_cohort_overall_survival_probability_5yr_capped_proportional_hazards_summary.txt`
- `<OUTPUT_DIR>/uveal_restricted/01_Efficacy/d_progression_free_survival/restricted_cohort_progression_free_survival_probability_5yr_capped_proportional_hazards_summary.txt`

---

### T-NICK-03: Continuous Variables And Dichotomization

> **STATUS:** Complete.

**Direct answer:** Yes. The reviewer-response adjusted Cox models use age as a continuous covariate rather than dichotomizing it, so this reviewer criticism is addressed for the main time-to-event models.

**Google Doc context:** Reviewer criticized dichotomizing continuous variables. Tim asked Nick to audit and correct this for the statistical response.

**Owner:** Nick.

**Action required:** Use continuous age in reviewer-response adjusted models and avoid unnecessary dichotomization in inferential models.

**Method / implementation:** Reviewer-response Cox models use continuous age at diagnosis with sex and tumor location in the parsimonious adjustment set.

**Result:** Reviewer-response adjusted models use continuous `age_at_diagnosis` with sex and tumor location where estimable. The dichotomized `age_at_diagnosis_general_pop_median` field is not part of the reviewer-response adjusted Cox covariate set.

**Response/manuscript language:** State that adjusted models used continuous age rather than a dichotomized age cutoff.

**Limitation / interpretation:** With limited events, flexible nonlinear terms may not be supportable; linear continuous adjustment is the conservative default.

**Evidence:** `scripts/config/modeling_policy.R`; `tests/testthat/test_peer_review_revision_contract.R`; exact Cox model table aliases listed under T-NICK-01 and T-NICK-06.

---

### T-NICK-04: PFS Definition And Death/Progression Components

> **STATUS:** Complete.

**Direct answer:** Yes. PFS is now defined as time from primary treatment to first local recurrence, metastatic progression, or death from any cause. The prior omission of metastatic progression from the PFS composite was fixed, and the regenerated full and restricted datasets use the corrected definition.

**Google Doc context:** Reviewer asked for a clearer PFS definition and death handling. The internal audit found the prior code did not fully match a conventional PFS label because metastatic progression was omitted from the composite.

**Owner:** Nick.

**Action required:** Ensure code, methods, figures, and response text use one identical PFS definition.

**Method / implementation:** PFS is defined as time from primary treatment to first local recurrence, metastatic progression, or death from any cause. This corrected definition includes metastatic progression; the previous omission of metastatic progression from the composite was an implementation error.

**Result:** The regenerated full and restricted analytic datasets define `pfs_event` exactly as the composite of local recurrence, metastatic progression, or death from any cause. In the full cohort, all 82 dataset-level PFS events correspond to at least one of `recurrence_event`, `mets_event`, or `death_event`; in the restricted cohort, all 52 dataset-level PFS events correspond to at least one of those components. The full-follow-up Cox model event counts reported under T-NICK-06 are lower (79 full-cohort and 51 restricted-cohort PFS events) because the Cox model-prep step excludes sparse `location` levels. That sparse-location exclusion removes 3 full-cohort PFS events and 1 restricted-cohort PFS event; it is not a PFS-definition discrepancy.

**Response/manuscript language:** Define PFS as time from primary treatment to first local recurrence, metastatic progression, or death from any cause. State that death without documented local recurrence or metastatic progression is counted as a PFS event.

**Limitation / interpretation:** The endpoint definition was repaired to include metastatic progression in the composite. This can change PFS results relative to analyses that treated PFS as local recurrence or death only, so manuscript methods and any response text must use the regenerated definition consistently.

**Evidence:** `<PROCESSED_DATA_DIR>/uveal_melanoma_full_cohort.rds`; `<PROCESSED_DATA_DIR>/uveal_melanoma_restricted_cohort.rds`; `docs/CALCULATIONS.md`; `tests/testthat/test_peer_review_revision_contract.R`.

---

### T-NICK-05: Separate Time-To-Local-Recurrence And Time-To-Metastasis Outputs

> **STATUS:** Complete.

**Direct answer:** Working answer: yes. The generated results include separate adjusted Cox models for time to local recurrence and time to metastatic progression. That is the direct answer to Tim's question and the direct statistical response to the reviewer concern about recurrence/metastasis comparisons needing to account for follow-up duration.

**Google Doc context:** The reviewer specifically criticized crude local-recurrence and metastatic-progression event rates because follow-up differed between treatment groups. Tim's follow-up question was whether we could include Cox models for "time to local recurrence" and "time to metastatic progression."

**Owner:** Nick.

**Action required:** Answer yes: local recurrence and metastatic progression can be presented as separate follow-up-aware time-to-event analyses, with adjusted Cox models as the lead inferential results.

**Method / implementation:** The local-recurrence model uses time from primary treatment to first local recurrence. The metastatic-progression model uses time from primary treatment to metastatic progression. Both models compare GKSRS versus PBT and adjust for the reviewer-response covariate set: continuous age at diagnosis, sex, and tumor location. These Cox models answer the reviewer concern more directly than crude event rates because each patient contributes their observed follow-up time, not just a yes/no event indicator.

**Result:** Yes, the generated results include separate adjusted Cox models for time to local recurrence and time to metastatic progression. In the full cohort, GKSRS versus PBT HR was 1.36 for local recurrence (95% CI 0.669-2.78; p=0.394; Cox model n=254, 31 events) and 1.25 for metastatic progression (95% CI 0.679-2.30; p=0.472; Cox model n=254, 43 events). In the restricted cohort, HR was 1.72 for local recurrence (95% CI 0.674-4.41; p=0.256; Cox model n=164, 18 events) and 1.00 for metastatic progression (95% CI 0.446-2.24; p=0.999; Cox model n=164, 28 events). These results do not show a statistically significant treatment-group difference for either endpoint, but the intervals are wide. Cox model n values reflect sparse-location model exclusions, not source cohort denominators.

**Response/manuscript language:** "We reanalyzed local recurrence and metastatic progression as time-to-event endpoints using adjusted Cox proportional-hazards models, rather than relying on crude event-rate comparisons alone. In these models, no statistically significant treatment-group difference was observed for time to local recurrence or time to metastatic progression."

**Limitation / interpretation:** This likely answers what the collaborator wants from the reviewer comment: crude recurrence/metastasis rates should be demoted, and the response should lead with follow-up-aware time-to-event models. It does not by itself answer the separate reviewer request to report median follow-up by treatment arm; that remains tracked under T-NICK-09. It also does not make an equivalence claim, and low event counts limit precision. Competing-risk methods or cumulative-incidence estimates could be described as supportive if the manuscript needs absolute risks at specific time points, but the Cox models are the direct answer to the requested "time to local recurrence" and "time to metastatic progression" analyses.

**Evidence:** Exact Cox model table and PH summary aliases listed under T-NICK-01; `tests/testthat/test_peer_review_revision_contract.R`.

---

### T-NICK-06: Five-Year Limited-Follow-Up OS/PFS HR Sensitivity

> **STATUS:** Complete.

**Direct answer:** Working answer: yes. OS and PFS Cox models were refit after limiting follow-up to the first 60 months, while KM curves can still show observed follow-up beyond 5 years. The 60-month analyses did not materially change the statistical conclusion from the full-follow-up Cox models: OS and PFS treatment comparisons remained statistically non-significant.

**Google Doc context:** Reviewer raised concern about sparse late follow-up/risk sets. Tim asked Nick to address the statistical concern while keeping conclusions conservative.

**Owner:** Nick.

**Action required:** Answer Tim's parallel-analysis request: refit OS and PFS Cox models using only the first 60 months of follow-up for HR estimation, while still allowing the KM curves to display observed follow-up beyond 5 years.

**Method / implementation:** Cox models were refit after limiting follow-up to the first 60 months. The KM curves can still be plotted beyond 5 years to show the observed data, but the 5-year HR is estimated from the portion of follow-up where late sparse risk sets are less dominant. This is separate from KM plotting and does not alter the full-follow-up endpoint definitions.

**Result:** The 60-month Cox sensitivity analyses do not change the statistical conclusion from the full-follow-up Cox models. The n and event counts in this paragraph are Cox model counts after sparse-location exclusions, not source cohort denominators. In the full cohort, full-follow-up OS HR was 1.20 (95% CI 0.694-2.08; p=0.514; Cox model n=254, 56 events) versus 60-month OS HR 1.43 (95% CI 0.751-2.72; p=0.277; Cox model n=254, 38 events). Full-follow-up PFS HR was 1.21 (95% CI 0.771-1.90; p=0.409; Cox model n=254, 79 events) versus 60-month PFS HR 1.35 (95% CI 0.827-2.21; p=0.228; Cox model n=254, 65 events). In the restricted cohort, full-follow-up OS HR was 1.07 (95% CI 0.533-2.15; p=0.847; Cox model n=164, 39 events) versus 60-month OS HR 1.34 (95% CI 0.613-2.92; p=0.465; Cox model n=164, 27 events). Full-follow-up PFS HR was 1.16 (95% CI 0.647-2.07; p=0.623; Cox model n=164, 51 events) versus 60-month PFS HR 1.34 (95% CI 0.713-2.51; p=0.365; Cox model n=164, 41 events). In each cohort, the 60-month estimates remain statistically non-significant with confidence intervals crossing 1.

**Response/manuscript language:** "Because the number at risk becomes limited at later follow-up times, Cox models for OS and PFS were also refit after limiting follow-up to the first 60 months. These 60-month analyses were consistent with the full-follow-up models and did not show a statistically significant treatment-group difference for OS or PFS. KM curves continue to display observed follow-up beyond 5 years, but the parallel HR analysis is explicitly limited to the first 60 months."

**Limitation / interpretation:** This is appropriately described as a sensitivity analysis because it asks whether the OS/PFS treatment HRs change when inference is limited to a clinically interpretable period with stronger risk-set support. The 60-month HRs should be labeled as 60-month limited-follow-up estimates and not mixed with full-follow-up HRs. The full-cohort 60-month OS model has a PH diagnostic concern and should be reported with that caveat.

**Evidence:**

- `<OUTPUT_DIR>/uveal_full/01_Efficacy/c_overall_survival/full_cohort_overall_survival_probability_5yr_capped_coxph.html`
- `<OUTPUT_DIR>/uveal_full/01_Efficacy/c_overall_survival/full_cohort_overall_survival_probability_5yr_capped_proportional_hazards_summary.txt`
- `<OUTPUT_DIR>/uveal_full/01_Efficacy/d_progression_free_survival/full_cohort_progression_free_survival_probability_5yr_capped_coxph.html`
- `<OUTPUT_DIR>/uveal_full/01_Efficacy/d_progression_free_survival/full_cohort_progression_free_survival_probability_5yr_capped_proportional_hazards_summary.txt`
- `<OUTPUT_DIR>/uveal_restricted/01_Efficacy/c_overall_survival/restricted_cohort_overall_survival_probability_5yr_capped_coxph.html`
- `<OUTPUT_DIR>/uveal_restricted/01_Efficacy/c_overall_survival/restricted_cohort_overall_survival_probability_5yr_capped_proportional_hazards_summary.txt`
- `<OUTPUT_DIR>/uveal_restricted/01_Efficacy/d_progression_free_survival/restricted_cohort_progression_free_survival_probability_5yr_capped_coxph.html`
- `<OUTPUT_DIR>/uveal_restricted/01_Efficacy/d_progression_free_survival/restricted_cohort_progression_free_survival_probability_5yr_capped_proportional_hazards_summary.txt`

---

### T-NICK-07: Propensity-Score Feasibility

> **STATUS:** Complete: feasibility only.

**Direct answer:** A full-cohort propensity-score adjusted sensitivity should not be reported because treatment-group overlap is poor and the treatment-assignment model is separation-prone. The restricted cohort has acceptable measured overlap and is the only cohort where a PS-adjusted outcome sensitivity is technically supportable. My recommendation is not to add a full-cohort PS model. If Tim/Angie want a PS-adjusted sensitivity, it should be restricted-cohort only, exploratory, and clearly secondary to the prespecified adjusted Cox models.

**Google Doc context:** Tim asked Nick whether a propensity-score analysis could make the multivariable comparison more robust. This followed reviewer/statistical concerns that treatment comparisons should account for baseline clinical predictors and treatment-selection differences rather than relying on unadjusted or sparsely adjusted comparisons.

**Owner:** Nick with Angie/Tim.

**Action required:** Decide whether to promote a restricted-cohort-only PS-adjusted Cox sensitivity into the analysis workflow. The feasibility audit itself is complete; the PS-adjusted outcome model has not yet been run or reported.

**Method / implementation:** The feasibility tool estimated the probability of receiving GKSRS versus PBT using baseline variables available before treatment assignment. Candidate variables were chosen because they are the available baseline demographic, tumor, ocular, and symptom predictors most directly related to the reviewer concern about clinical differences between treatment groups: age at diagnosis, sex, tumor location, initial tumor height, initial tumor diameter, initial vision, simplified initial T stage, initial overall stage, optic nerve status, subretinal fluid (`srf`), baseline visual-field defect, and baseline vision loss/blurred vision.

Variables were included only if they were already present in the analytic dataset with usable support. Requested or clinically relevant variables that were not available in the curated analytic dataset were not included in the propensity model. These unavailable fields include quantitative macula/fovea/optic-nerve distances, dose to macula/fovea/optic nerve, GK margin dose, GK isodose, GK shots, and GK isocenters. Those would require new chart review or treatment-plan abstraction before they could be used in any adjustment model.

The tool checked whether a propensity-score analysis was technically supportable before estimating any outcome effect. First, it removed variables that could not help distinguish treatment groups because they were completely missing or had only one observed value. Second, it removed variables that were redundant with other variables in the treatment-assignment model. Third, it fit a logistic regression model where the outcome was treatment received (GKSRS versus PBT), using the retained baseline variables as predictors. The fitted value from that model is each patient's estimated probability of receiving GKSRS given the measured baseline variables.

The overlap check asks whether the dataset contains comparable GKSRS and PBT patients after summarizing baseline differences with the propensity score. If one treatment group has patients with propensity scores that the other group does not share, then the model has no real observed comparison group for those patients, so a PS-adjusted treatment comparison would be unreliable.

This feasibility tool did not estimate the treatment effect on recurrence, metastasis, OS, PFS, or vision; it only tested whether a future propensity-adjusted sensitivity analysis would have enough variable support and treatment-group overlap to be worth reporting. Therefore, the conclusion here is about whether a PS-adjusted outcome analysis is defensible to run, not about the treatment effect itself.

**Result:**

The feasibility workbooks were regenerated and inspected in the configured peer-review runtime tree.

- Full cohort: not reportable as a PS-adjusted sensitivity. All 260 patients were complete cases, but fitted propensity scores ranged from 0.000000174 to 0.994. That extreme range indicates separation-prone treatment assignment, and 89.2% of patients were inside common support. The diagnostic workbook sets `separation_prone = 1` and `reportable_sensitivity = 0`.

- Full cohort retained covariates: age at diagnosis, sex, location, initial tumor height, initial tumor diameter, initial vision, simplified initial T stage, optic nerve status, subretinal fluid, baseline visual-field defect, and baseline vision loss/blurred vision. Initial overall stage was excluded as aliased/collinear.

- Restricted cohort: feasible as an exploratory PS sensitivity candidate. All 167 restricted-cohort patients were complete cases, propensity scores ranged from 0.136 to 0.735, there was no separation-prone behavior, and 96.4% of patients were inside common support. The diagnostic workbook sets `separation_prone = 0` and `reportable_sensitivity = 1`.

- Restricted cohort retained covariates: age at diagnosis, sex, initial tumor height, initial tumor diameter, initial vision, simplified initial T stage, subretinal fluid, baseline visual-field defect, and baseline vision loss/blurred vision. Optic nerve status was single-level by restricted-cohort design, and location plus initial overall stage were excluded as aliased/collinear.

- Interpretation of the restricted-cohort finding: the restricted cohort was designed to focus on patients who were more plausibly eligible for both modalities based on tumor size and optic nerve criteria. That design rationale aligns with the PS feasibility result: after restricting to a more clinically comparable cohort, the measured baseline variables show better treatment-group overlap. This does not prove exchangeability or remove unmeasured confounding, but it supports using the restricted cohort, rather than the full cohort, if a PS-adjusted outcome sensitivity is pursued.

- Practical conclusion: the response can say that PS feasibility was assessed and that full-cohort PS adjustment was not supportable. The remaining decision is whether the added complexity of a restricted-only exploratory PS-adjusted Cox sensitivity is worth including. If included, it should be supplemental or response-only sensitivity evidence, not a replacement for the reviewer-response adjusted Cox models.

**Response/manuscript language:**

Do not present a full-cohort PS-adjusted sensitivity. Any propensity-score sensitivity should be limited to the restricted cohort and described as exploratory. Do not describe the current feasibility workbooks as treatment-effect analyses because they do not estimate recurrence, metastasis, OS, PFS, or visual-acuity effects.

Suggested response wording: "We assessed whether a propensity-score sensitivity analysis was supportable using available baseline demographic, tumor, ocular, and symptom variables. The full cohort showed poor overlap/separation and was not suitable for propensity adjustment. The restricted cohort, which was designed to focus on patients more plausibly eligible for both modalities, had acceptable overlap; any propensity-adjusted analysis in that cohort is therefore treated as an exploratory sensitivity rather than evidence of randomized balance. Quantitative dose, macular/foveal distance, optic-nerve distance, and GK planning variables were not available in the curated analytic dataset and therefore could not be included in the propensity model."

**Limitation / interpretation:** Propensity scores are a way to summarize measured baseline treatment-selection factors into the modeled probability of receiving one treatment rather than the other. They do not create randomization, do not adjust for unmeasured confounding, and cannot solve missing radiation dosimetry, quantitative proximity, visual-field detail, or other uncollected treatment-selection factors.

**Evidence:** `<OUTPUT_DIR>/peer_review_revision_audits/full_propensity_score_feasibility.xlsx`; `<OUTPUT_DIR>/peer_review_revision_audits/restricted_propensity_score_feasibility.xlsx`.

---

### T-NICK-08: Visual-Acuity Follow-Up And Minimum-Follow-Up Sensitivity

> **STATUS:** Complete; Tim/Angie choose which minimum-follow-up cutoff, if any, to emphasize.

**Direct answer:** Yes for the statistical/data-availability part, but the three reviewer subrequests need to be kept separate. Latest-VA timing can be reported using treatment date to `last_followup`, which is the available date associated with `last_vision`; there is no separate `last_vision_date` field. Minimum-latest-VA-follow-up sensitivities have now been run using 12-, 36-, and 60-month cutoffs for both full and restricted cohorts. These cutoffs give Tim options aligned with the reviewer's 1-, 3-, and 5-year language, but they are **not** the same as visual acuity measured at standardized 1-, 3-, and 5-year landmarks. Tim's comment that additional data collection is not feasible within the 60-day timeframe applies specifically to standardized visual-acuity landmarks and time-to-clinically-meaningful VA deterioration, such as loss of at least 3 Snellen lines. Those two analyses are not available from the current curated dataset.

**Reviewer subrequest mapping:**

| Reviewer subrequest | Current answer |
| --- | --- |
| a. Visual acuity at standardized landmarks, such as 1, 3, and 5 years | Not feasible from the current curated dataset because this requires dated longitudinal VA measurements near those fixed landmark windows. The current data have baseline VA, latest VA, recurrence pretreatment VA where applicable, and `last_followup`, not a longitudinal VA visit table. |
| b. Time to clinically meaningful deterioration, such as loss of at least 3 Snellen lines | Not feasible as a time-to-event analysis from the current curated dataset because there is no dated series showing when each patient first crossed the deterioration threshold. The existing Objective 2 workflow can summarize latest-available Snellen line change and line-loss distribution, but it cannot estimate time to first 3-line loss. |
| c. Sensitivity analysis restricted to patients with adequate minimum follow-up | Feasible and completed. The adjusted latest-VA model was rerun after requiring at least 12, 36, and 60 months of latest-VA follow-up. These are minimum-follow-up sensitivity analyses, not fixed-landmark VA analyses. |

**Google Doc context:** Reviewer asked for timing of last visual-acuity assessment, median ophthalmologic follow-up by treatment group, standardized landmarks if possible, time to meaningful VA deterioration if possible, minimum-follow-up sensitivity, and adjustment for clinical predictors. Tim directed Nick to address the statistical portion. Tim's note that "These points would require additional data collection which is not feasible within the 60 day timeframe" was scoped to the standardized-landmark VA request and the time-to-clinically-meaningful-deterioration request, not to the latest-VA timing summary or minimum-follow-up sensitivity.

**Owner:** Nick, with clinical coauthor input for unavailable predictors.

**Action required:** Report latest-VA timing by treatment group, state that standardized VA landmarks and time-to-clinically-meaningful VA deterioration would require additional data collection, and have Tim/Angie/Nick decide whether to report the 12-, 36-, or 60-month minimum-follow-up sensitivity. Do not describe the 12/36/60-month cutoff sensitivities as landmark VA analyses.

**Method / implementation:** The Objective 2 VA workflow now derives `last_vision_followup_months` from treatment date to `last_followup`. The primary VA endpoint remains the implemented visual-acuity change score: initial vision minus final or recurrence-pre-treatment vision. Reviewer-response sensitivities then restrict to patients with `last_vision_followup_months >= 12`, `>= 36`, and `>= 60` and rerun the adjusted treatment-effect model for `vision_change ~ treatment_group + age_at_diagnosis + sex + location`.

For the reviewer's separate request to adjust for clinically relevant baseline predictors of visual outcome, the added sensitivity model uses latest logMAR VA (`last_vision`) as the outcome rather than adding baseline VA to the change-score outcome. This directly adjusts for starting VA without putting `initial_vision` on both sides of a change-score equation. The full-cohort model was:

`last_vision ~ treatment_group + initial_vision + last_vision_followup_months_explicit + initial_tumor_height + initial_tumor_diameter + initial_t_stage_simple + srf + optic_nerve + treatment_year_centered + age_at_diagnosis + sex + location`

The restricted-cohort model used the same structure except `optic_nerve` was not included because it had insufficient variation in the restricted cohort:

`last_vision ~ treatment_group + initial_vision + last_vision_followup_months_explicit + initial_tumor_height + initial_tumor_diameter + initial_t_stage_simple + srf + treatment_year_centered + age_at_diagnosis + sex + location`

Requested/available predictor handling:

| Candidate variable | Used in latest-VA sensitivity? | Reason |
| --- | --- | --- |
| `age_at_diagnosis` | Yes, full and restricted | Standard confounder set. |
| `sex` | Yes, full and restricted | Standard confounder set. |
| `location` | Yes, full and restricted | Standard confounder set. |
| `initial_tumor_height` | Yes, full and restricted | Available baseline tumor-size predictor. |
| `initial_tumor_diameter` | Yes, full and restricted | Available baseline basal-diameter predictor. |
| `initial_vision` | Yes, full and restricted | Baseline logMAR VA; adjusted in the latest-VA model rather than in the change-score model. |
| `initial_t_stage_simple` | Yes, full and restricted | Available baseline T-stage grouping with adequate model support. |
| `initial_overall_stage` | No | Not added with `initial_t_stage_simple`; including both would duplicate tumor-stage information and add model burden without answering an additional reviewer-requested VA predictor beyond T stage. |
| `optic_nerve` | Yes in full; no in restricted | Full cohort had usable variation. Restricted cohort had insufficient variation, so it was excluded there. |
| `srf` | Yes, full and restricted | Available baseline subretinal-fluid field. |
| `visual_field_defect` | No | Not used as an adjustment covariate for the VA sensitivity because the reviewer request was for visual-acuity outcome predictors and the response should not imply a longitudinal visual-field analysis. Missing/limited visual-field information remains a limitation. |
| `vision_loss_blurred_vision` | No | Not used as an adjustment covariate; it is a symptom/presentation descriptor rather than the prespecified baseline VA measure used for adjustment. |

**Result:** Latest-VA follow-up timing is available for 213/260 full-cohort patients and 133/167 restricted-cohort patients. Median latest-VA follow-up was longer for PBT than GKSRS: full cohort 77.5 vs 32.7 months; restricted cohort 78.2 vs 35.3 months. The adjusted minimum-follow-up sensitivity results are:

| Cohort | Minimum latest-VA follow-up cutoff | PBT n | GKSRS n | Adjusted MD for GKSRS vs PBT | 95% CI | p-value |
| --- | ---: | ---: | ---: | ---: | --- | ---: |
| Full | 12 months | 86 | 88 | 0.13 logMAR | -0.15 to 0.40 | 0.362 |
| Full | 36 months | 74 | 54 | 0.06 logMAR | -0.27 to 0.40 | 0.704 |
| Full | 60 months | 59 | 39 | -0.05 logMAR | -0.45 to 0.35 | 0.805 |
| Restricted | 12 months | 70 | 43 | 0.34 logMAR | 0.01 to 0.68 | 0.046 |
| Restricted | 36 months | 59 | 25 | 0.31 logMAR | -0.15 to 0.76 | 0.185 |
| Restricted | 60 months | 50 | 18 | 0.25 logMAR | -0.30 to 0.80 | 0.365 |

Only the restricted-cohort 12-month cutoff is nominally statistically significant. The 36- and 60-month restricted sensitivities are not, and the full-cohort sensitivities are not. Because these are optional cutoff sensitivities with shrinking sample size, they should not be overstated as definitive visual-acuity evidence.

The reviewer-predictor latest-VA sensitivity attenuated the apparent visual-acuity association and was not statistically significant. In the full cohort, the adjusted GKSRS-vs-PBT difference was -0.200 logMAR (95% CI -0.462 to 0.063; p=0.135; n=208). In the restricted cohort, the adjusted difference was -0.126 logMAR (95% CI -0.429 to 0.177; p=0.411; n=126). These results do not support a strong claim that GKSRS provides superior visual-acuity preservation after adjustment for available baseline VA, follow-up duration, tumor-size/stage predictors, SRF, treatment year, and the standard confounder set.

**Response/manuscript language:** "Visual-acuity change was analyzed using the latest available visual-acuity value, with `last_followup` used as the available date associated with that value. Median latest-VA follow-up was longer after PBT than GKSRS. We therefore repeated the adjusted visual-acuity model among patients with minimum latest-VA follow-up cutoffs of 12, 36, and 60 months. These are minimum-follow-up sensitivity analyses, not fixed 1-, 3-, or 5-year landmark visual-acuity analyses. Standardized longitudinal VA measurements at fixed 1-, 3-, or 5-year landmarks and time-to-clinically-meaningful deterioration, such as loss of at least 3 Snellen lines, would require additional data collection and were not feasible within the revision timeframe. We also added a latest-VA sensitivity model adjusting for baseline logMAR VA, latest-VA follow-up duration, baseline tumor height, basal diameter, T stage, subretinal fluid, treatment year, and the standard confounders; optic-nerve status was included in the full cohort but not in the restricted cohort because it had insufficient variation. This sensitivity was not statistically significant in either cohort, so visual-acuity findings should be interpreted cautiously rather than as evidence of superior functional preservation."

**Limitation / interpretation:** Minimum-follow-up restrictions answer whether the latest-VA treatment comparison changes after excluding short-follow-up patients; they do not estimate VA at exact 1-, 3-, or 5-year visits. The sensitivity reduces sample size as the cutoff increases, especially in the restricted GKSRS arm, and still cannot adjust for missing dose, foveal/macular/optic-nerve distance, or visual-field variables.

**Evidence:** `<OUTPUT_DIR>/uveal_full/02_Safety/a_vision_changes/full_cohort_vision_followup_sensitivity.xlsx`; `<OUTPUT_DIR>/uveal_restricted/02_Safety/a_vision_changes/restricted_cohort_vision_followup_sensitivity.xlsx`; `<OUTPUT_DIR>/uveal_full/02_Safety/a_vision_changes/full_cohort_minimum_followup_12_months_vision_change_minimum_followup_12_months_diagnostics.xlsx`; `<OUTPUT_DIR>/uveal_full/02_Safety/a_vision_changes/full_cohort_minimum_followup_36_months_vision_change_minimum_followup_36_months_diagnostics.xlsx`; `<OUTPUT_DIR>/uveal_full/02_Safety/a_vision_changes/full_cohort_minimum_followup_60_months_vision_change_minimum_followup_60_months_diagnostics.xlsx`; `<OUTPUT_DIR>/uveal_restricted/02_Safety/a_vision_changes/restricted_cohort_minimum_followup_12_months_vision_change_minimum_followup_12_months_diagnostics.xlsx`; `<OUTPUT_DIR>/uveal_restricted/02_Safety/a_vision_changes/restricted_cohort_minimum_followup_36_months_vision_change_minimum_followup_36_months_diagnostics.xlsx`; `<OUTPUT_DIR>/uveal_restricted/02_Safety/a_vision_changes/restricted_cohort_minimum_followup_60_months_vision_change_minimum_followup_60_months_diagnostics.xlsx`; `<OUTPUT_DIR>/peer_review_revision_audits/full_followup_and_data_availability.xlsx`; `<OUTPUT_DIR>/peer_review_revision_audits/restricted_followup_and_data_availability.xlsx`.

---

### T-NICK-09: Follow-Up Duration Summaries

> **STATUS:** Complete.

**Direct answer:** Working answer: yes. Treatment-arm follow-up summaries are available for the full and restricted cohorts, and the full/restricted follow-up audit workbooks have been regenerated and sheet-checked. In both cohorts, PBT patients had longer median follow-up than GKSRS patients, which supports the reviewer concern that crude recurrence/metastasis event rates are hard to interpret without time-to-event methods. This is why recurrence and metastatic progression should be presented with Cox time-to-event models rather than crude binary comparisons.

**Google Doc context:** Reviewer emphasized unequal follow-up and asked for follow-up timing to interpret recurrence, metastasis, visual acuity, and toxicity comparisons. Tim asked Nick to generate response-ready follow-up summaries.

**Owner:** Nick.

**Action required:** Summarize follow-up by treatment group and endpoint-relevant timing fields.

**Method / implementation:** Use the on-demand follow-up/data-availability audit tool for full and restricted cohorts. The audit workbook sheet `followup_by_treatment_arm` summarizes `follow_up_months` and derived latest-VA timing (`treatment_date` to `last_followup` among rows with both dates) by treatment group with n, min, mean, median, max, and IQR.

**Result:** Treatment-arm follow-up summaries are available for the full and restricted cohorts from the regenerated follow-up/data-availability audit workbooks. The n values below are source cohort denominators for general follow-up (`follow_up_months`). Latest-VA timing n values count rows with non-missing derived latest-VA timing among patients with both `treatment_date` and `last_followup`.

**General follow-up (`follow_up_months`, months)**

| Cohort | Treatment | n | Follow-up available | Min | Mean | Median | Max | IQR |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Full | PBT | 121 | 121 | 3.0 | 83.2 | 75.6 | 376.3 | 73.8 |
| Full | GKSRS | 139 | 139 | 1.1 | 51.1 | 37.2 | 246.1 | 64.5 |
| Restricted | PBT | 103 | 103 | 3.0 | 81.5 | 75.2 | 376.3 | 81.1 |
| Restricted | GKSRS | 64 | 64 | 1.1 | 53.7 | 38.0 | 205.4 | 61.6 |

**Latest-VA timing (`treatment_date` to `last_followup`, months; among rows with both dates)**

| Cohort | Treatment | Cohort n | Timing available | Min | Mean | Median | Max | IQR |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Full | PBT | 121 | 96 | 1.0 | 81.0 | 77.5 | 203.4 | 73.9 |
| Full | GKSRS | 139 | 117 | 0.0 | 47.2 | 32.7 | 194.1 | 66.8 |
| Restricted | PBT | 103 | 80 | 1.0 | 80.1 | 78.3 | 203.4 | 81.7 |
| Restricted | GKSRS | 64 | 53 | 0.0 | 50.6 | 35.3 | 194.1 | 69.1 |

In both cohorts, PBT patients had longer mean and median general follow-up than GKSRS patients. The same pattern holds for latest-VA timing among patients with usable timing fields.

**Response/manuscript language:** "Median and mean follow-up differed by treatment arm. In the full cohort, general follow-up was a median of 75.6 months (mean 83.2; range 3.0-376.3) for PBT and a median of 37.2 months (mean 51.1; range 1.1-246.1) for GKSRS. In the restricted cohort, general follow-up was a median of 75.2 months (mean 81.5; range 3.0-376.3) for PBT and a median of 38.0 months (mean 53.7; range 1.1-205.4) for GKSRS. Because recurrence and metastatic progression are follow-up-dependent endpoints, the revised analysis emphasizes time-to-event Cox models rather than crude event-rate comparisons."

**Limitation / interpretation:** General follow-up duration is not the same as dated endpoint-specific assessment frequency. The latest-VA timing summary uses treatment date to `last_followup` as the available timing field for latest visual acuity; standardized 1-, 3-, and 5-year VA landmarks remain a separate visual-acuity feasibility question under T-NICK-08.

**Evidence:** `<OUTPUT_DIR>/peer_review_revision_audits/full_followup_and_data_availability.xlsx` (sheet `followup_by_treatment_arm`); `<OUTPUT_DIR>/peer_review_revision_audits/restricted_followup_and_data_availability.xlsx` (sheet `followup_by_treatment_arm`).

---

### T-NICK-10: Tumor-Height Timing

> **STATUS:** Complete; recommendation is to demote/caveat comparative tumor-height change.

**Direct answer:** Yes. The timing audit shows that follow-up tumor-height assessment was not standardized and was later for PBT than GKSRS in both cohorts. Because tumor-height change is based on latest observed height rather than a fixed landmark measurement, the comparative regression should not be presented as a strong treatment-effect endpoint. It should be demoted to descriptive/secondary interpretation or reported only with a clear timing limitation.

**Google Doc context:** Reviewer/statistical memo identified latest-minus-baseline tumor height as vulnerable to unequal imaging timing. Tim's response needs a defensible statement about timing.

**Owner:** Nick.

**Action required:** Decide with Tim/Angie whether the manuscript will remove the comparative tumor-height regression or retain it only as a secondary descriptive analysis with the timing limitation stated directly.

**Method / implementation:** The Objective 1 tumor-height workflow now derives `last_height_followup_months` from treatment date to `last_height_date` and writes a timing workbook with treatment-arm summaries and negative-interval details. The timing variable was not added as a regression covariate because it occurs after treatment and may reflect follow-up patterns or outcome processes rather than baseline confounding.

**Result:** Median time from treatment to latest height assessment was longer for PBT than GKSRS. In the full cohort, timing was available for 119/121 PBT and 138/139 GKSRS patients; median timing was 53.9 months for PBT and 26.3 months for GKSRS. In the restricted cohort, timing was available for 101/103 PBT and 64/64 GKSRS patients; median timing was 51.8 months for PBT and 31.9 months for GKSRS. Negative treatment-to-height intervals were also present, indicating source-date issues or pre-treatment measurement carryover that require caution: 7 full-cohort rows and 5 restricted-cohort rows.

**Response/manuscript language:** "Tumor-height change was based on latest available follow-up height rather than a fixed imaging landmark. Time from treatment to latest height assessment differed by treatment arm, with longer follow-up after PBT than GKSRS. We therefore interpret tumor-height change descriptively/secondarily and do not treat the comparative tumor-height regression as a definitive treatment-effect analysis."

**Limitation / interpretation:** Latest observed height is not a standardized landmark or longitudinal model. Conditioning on latest-height timing as a model covariate would not solve the problem because timing is post-treatment; it could introduce bias rather than correct the comparison.

**Evidence:** `<OUTPUT_DIR>/uveal_full/01_Efficacy/e_tumor_height_primary/full_cohort_tumor_height_timing_summary.xlsx`; `<OUTPUT_DIR>/uveal_restricted/01_Efficacy/e_tumor_height_primary/restricted_cohort_tumor_height_timing_summary.xlsx`.

---

### T-NICK-11: PRAME And T4 Subgroup Pruning

> **STATUS:** Complete.

**Direct answer:** Yes. Reviewer-facing Objective 1 subgroup surfaces no longer promote PRAME (`gep12_prame_status`) or T4 (`initial_t_stage_simple`) displays. PRAME was removed from the reviewer-facing subgroup variable list; T4 rows are dropped before subgroup modeling. Forest-plot diagnostics include a `reviewer_pruning_audit` sheet and `reviewer_exclusion_note` metadata so these removals are auditable rather than silent.

**Google Doc context:** Tim/Angie comments indicated reviewer-facing subgroup displays should be pruned, including PRAME and T4 displays where support is sparse or not useful for the response.

**Owner:** Nick.

**Action required:** Verify regenerated full/restricted forest diagnostics and confirm Tim/Angie agree these subgroup surfaces stay demoted/exploratory.

**Method / implementation:** Removed `gep12_prame_status` from `subgroup_vars` and `FOREST_PLOT_VARIABLE_ORDER`. Added `drop_reviewer_excluded_subgroup_levels()` to exclude T4 before subgroup modeling. Forest diagnostics workbooks now include `reviewer_exclusion_note` columns and a `reviewer_pruning_audit` sheet with excluded counts.

**Result:** Code and diagnostics wiring are in place. Regenerated runtime verification should confirm that forest plot diagnostics and rendered forest data contain no displayed T4 rows or PRAME subgroup variable, with audit counts populated from the entering cohort data.

**Response/manuscript language:** Label retained subgroup analyses exploratory and omit unsupported rows from reviewer-facing figures/tables. State that PRAME local-recurrence subgroup displays and T4 subgroup displays were removed from reviewer-facing outputs because of inadequate or sparse support.

**Limitation / interpretation:** Sparse subgroup estimates are unstable and multiplicity is substantial. Remaining subgroup surfaces are exploratory support only.

**Evidence:** `<OUTPUT_DIR>/uveal_full/01_Efficacy/g_subgroup_analysis/forest_plots/full_cohort_forest_plot_diagnostics.xlsx`; `<OUTPUT_DIR>/uveal_restricted/01_Efficacy/g_subgroup_analysis/forest_plots/restricted_cohort_forest_plot_diagnostics.xlsx`; `tests/testthat/test_objective1_primary_outcomes.R`.

---

### T-NICK-12: Dose, Proximity, Optic Nerve, Macula/Fovea, Visual Field, And Treatment-Detail Availability

> **STATUS:** Complete for Nick's statistical limitation; clinical technique prose remains coauthor-owned.

**Direct answer:** For Nick's data/statistical part, the working answer is yes: the curated analytic dataset does not contain quantitative macula/fovea/optic-nerve distance fields or GK/PBT dose-planning fields needed for adjustment. The response can say these factors could not be included in statistical models without new chart or radiation-plan abstraction only after the source-absence audit is verified. Clinical coauthors still need to supply the treatment-technique narrative that is not a dataset question.

**Google Doc context:** Reviewers requested radiation dose, optic-disc/macular/foveal proximity, visual-field, and treatment-detail information. Tim noted some of this may not be feasible in the response window and asked Nick to verify data availability before wording limitations.

**Owner:** Nick plus Tim/Melhus/Mignano/clinical coauthors.

**Action required:** Audit derived analytic data, raw/source columns, and data dictionary before stating whether requested fields are available or require new chart/plan review.

**Method / implementation:** Use the on-demand follow-up/data-availability audit tool; distinguish fields in the dataset from fields that would require new chart abstraction.

**Result:** The working audit answers the data-availability question for the curated analytic dataset and active curated source workbook, but the final follow-up/data-availability workbook aliases still need artifact verification. Available treatment-detail fields include treatment indicators and dates (`initial_gk`, `initial_gk_date`, `initial_plaque`, `initial_plaque_date`), plaque/radionuclide fields for many PBT patients (`radionuclide`, `plaque_size`, `plaque_notch`), and optic nerve status as a categorical proximity/abutment descriptor. In the full cohort, non-missing counts are 139/260 for `initial_gk_date`, 121/260 for `initial_plaque_date`, 119/260 for `radionuclide`, 119/260 for `plaque_size`, 118/260 for `plaque_notch`, and 260/260 for `optic_nerve`. In the restricted cohort, non-missing counts are 64/167 for `initial_gk_date`, 103/167 for `initial_plaque_date`, 101/167 for `radionuclide`, 101/167 for `plaque_size`, 100/167 for `plaque_notch`, and 167/167 for `optic_nerve`.

The requested quantitative dosimetry and proximity fields are not present as derived analytic columns in the curated analytic dataset: `macula_distance`, `fovea_distance`, `optic_nerve_distance`, `dose_to_macula`, `dose_to_fovea`, `dose_to_optic_nerve`, `gk_margin_dose`, `gk_isodose`, `gk_shots`, and `gk_isocenters`. The final response should say these variables require new chart review or radiation-plan abstraction only after the source workbook and data dictionary checks required by the plan's Source-Absence Claim Rule are documented.

**Response/manuscript language:** Report the treatment-detail fields that are available descriptively, and state that quantitative macular/foveal/optic-nerve distances and dose/planning parameters were not available in the curated analytic dataset. Suggested response wording: "The analytic dataset included treatment dates, plaque/radionuclide descriptors, plaque size/notch fields, and categorical optic nerve involvement/abutment status. It did not include quantitative macular/foveal/optic-nerve distances or treatment-plan dosimetry fields, so those factors could not be included as covariates without new treatment-plan abstraction."

**Limitation / interpretation:** This is an answer to Nick's data-availability/statistical question. It does not replace clinical coauthor input for prose details about GK targeting, plaque prescription technique, implant duration, or constraints. Missing dose/proximity/visual-field data limit causal interpretation of visual-acuity and toxicity comparisons.

**Expected evidence after generation/verification:** `<OUTPUT_DIR>/peer_review_revision_audits/full_followup_and_data_availability.xlsx`; `<OUTPUT_DIR>/peer_review_revision_audits/restricted_followup_and_data_availability.xlsx`.

---

### T-NICK-13: Adverse-Event Grading And SRD Scope

> **STATUS:** Complete for checked source/derived fields.

**Direct answer:** Yes for the statistical/data-audit part. The checked derived full/restricted cohorts and current source workbook do not contain standardized toxicity grade fields or dated onset fields for each toxicity. The current Objective 2 toxicity outputs should be described as recorded toxicity burden by available follow-up, not CTCAE-style graded adverse-event incidence and not time-to-toxicity incidence. SRD cause is available for records coded as SRD, but the current Objective 2 SRD endpoint is all recorded SRD burden, not a radiation-induced-only SRD incidence analysis.

**Google Doc context:** Reviewer asked for severity/grading and clearer toxicity reporting. Tim's response needs to avoid overstating adverse-event grading if grading was not collected for each toxicity.

**Owner:** Nick with clinical coauthor confirmation.

**Action required:** Use recorded-burden language for retinopathy, NVG, and SRD unless clinical coauthors identify separate graded/dated toxicity data outside the checked source workbook.

**Method / implementation:** Audited derived cohort columns, the active source workbook, and data-dictionary candidates for toxicity, grading, radiation attribution, SRD, and timing fields. Added a `toxicity_scope` sheet to the Objective 2 reviewer-facing sensitivity workbook so downstream response language states the supported scope.

**Result:** Derived full/restricted toxicity fields include `nvg`, `srd`, `srd_cause`, `nvg_burden_event`, and `srd_burden_event`; no standardized toxicity-grade field was identified. In the full cohort, SRD cause among SRD-positive records was Mass for 5 and Radiation for 48; in the restricted cohort, Mass for 4 and Radiation for 22. Because the current Objective 2 analysis uses `srd_burden_event`, it combines all recorded SRD burden rather than isolating radiation-attributed SRD. Retinopathy and NVG outputs likewise represent recorded burden by available follow-up.

**Response/manuscript language:** "Adverse events were analyzed as recorded toxicity burden by available follow-up. The curated dataset did not contain standardized toxicity grades or dated onset fields for each adverse event, so these analyses should not be interpreted as graded adverse-event incidence or time-to-toxicity analyses. SRD cause was recorded for SRD-positive cases, but the reported SRD endpoint reflects all recorded SRD burden rather than a radiation-induced-only incidence endpoint."

**Limitation / interpretation:** Ever-recorded adverse-event burden can be influenced by follow-up duration and assessment frequency. If the manuscript needs graded toxicity or radiation-induced-only SRD, that would require either new source abstraction or a deliberately separate endpoint using `srd_cause`, not the current all-burden toxicity model.

**Evidence:** `<OUTPUT_DIR>/uveal_full/02_Safety/a_vision_changes/full_cohort_vision_followup_sensitivity.xlsx`; `<OUTPUT_DIR>/uveal_restricted/02_Safety/a_vision_changes/restricted_cohort_vision_followup_sensitivity.xlsx`; `<OUTPUT_DIR>/uveal_full/02_Safety/`; `<OUTPUT_DIR>/uveal_restricted/02_Safety/`; checked source workbook `Ocular Melanoma Master Spreadsheet REVISED FOR STATS (5-10-25, TJM).xlsx`; checked data dictionary candidates.

---

### T-NICK-14: Statistical Language Guardrails

> **STATUS:** Tim final wording pass; not a Nick blocker.

**Direct answer:** The response should not claim equivalence, noninferiority, randomized balance, or broad vision preservation. The statistically honest answer is that adjusted observational analyses did not show statistically significant treatment-group differences for the reported endpoints, with residual confounding, limited events, missing dose/proximity variables, and follow-up limitations acknowledged.

**Google Doc context:** Reviewer/statistical comments objected to language implying equivalence, simulated randomization, or broad visual preservation beyond measured endpoints.

**Owner:** Nick with Tim final manuscript approval.

**Action required:** Align methods, results, discussion, and response language with observational, non-randomized, limited-event evidence.

**Method / implementation:** Remove "same efficacy," "preferred modality," "simulate randomization," and broad "vision preservation" unless directly supported by measured outcomes.

**Result:** Nick's statistical language guardrails are documented. Tim owns the final manuscript wording pass.

**Response/manuscript language:** Use "no statistically significant difference was observed," "restricted cohort intended to improve comparability," and "visual-acuity change" or "visual-acuity preservation" only when that exact endpoint is being discussed.

**Limitation / interpretation:** The study does not establish equivalence, noninferiority, or causal superiority.

**Evidence:** `docs/INTERPRETATION_GUIDE.md`; `docs/STATISTICAL_METHODS.md`; `docs/CALCULATIONS.md`.

---

### T-NICK-15: Table 1 P-Value Handling

> **STATUS:** Complete.

**Direct answer:** No Nick-side table-generation change is needed. Tim will handle final manuscript Table 1 p-value removal or formatting manually, while Nick's generated working tables can keep p-values for internal review. The reviewer-facing interpretation should still treat baseline characteristics as descriptive, not as proof of balance.

**Google Doc context:** Tim's Google Doc comments indicated Table 1 p-values were not Nick's main action item because Tim planned to remove them manually, but the statistical response still needs to remain consistent.

**Owner:** Tim for final manuscript table handling.

**Action required:** No Nick-side generated-table change. Keep p-values in generated analysis tables for now; Tim will remove or edit Table 1 p-values manually in the final manuscript product if needed.

**Method / implementation:** Do not alter the generated table pipeline for this item. The statistical-response principle remains that baseline characteristics are descriptive in this observational cohort; p-values in generated working tables should not be interpreted as evidence of randomization or exchangeability.

**Result:** No runtime or code action planned for Nick. Generated working tables may retain p-values.

**Response/manuscript language:** If this comes up in response wording, state that baseline comparisons are descriptive because treatment assignment was observational and clinically selected. Final table formatting is Tim's manuscript-editing task.

**Limitation / interpretation:** Baseline p-values do not establish exchangeability or absence of confounding.

**Evidence:** `docs/STATISTICAL_METHODS.md`; generated table outputs are intentionally unchanged for now.

---

## Clinical Coauthor Items Nick Should Track But Not Own Alone

| Item | Source context | Primary owner | Nick role |
| --- | --- | --- | --- |
| PBT treatment technique details | Reviewer asked for radionuclide, prescription depth/dose, implant duration, plaque diameter, and notched plaque details. | Tim / Melhus | Provide available dataset fields and statistical limitation wording. |
| GK treatment technique details | Reviewer asked for dose range, isodose approach, target definition, margin/setup, and constraints. | Mignano / Tim | Confirm whether any GK planning fields exist in the dataset; otherwise document as requiring plan review. |
| Radiation dose adjustment | Reviewer asked whether dose should be adjusted for visual/toxicity outcomes. | Tim / Melhus / Mignano | Explain statistical impact of missing/non-comparable dose fields. |
| Reference replacement | Reviewers objected to UpToDate/EyeWiki citations. | Tim / clinical team | None unless statistical-method references are needed. |

---

## Supporting Endpoint And Claim Audit

See also [PR_VS_ORIGINAL_RESULTS_AUDIT.md](/Users/ncamarda/Projects/uveal_melanoma-peer-review-statistical-revision/docs/PR_VS_ORIGINAL_RESULTS_AUDIT.md), especially Section 6, for concise manuscript/reviewer-response reporting guidance.

### PFS

- Outcome label: Progression-Free Survival (PFS)
- Reviewer-response definition: Time from primary treatment to first local recurrence, metastatic progression, or death from any cause.
- Code path: `scripts/data_helper/data_derivation.R`; `scripts/workflow/objective_1_primary_outcomes.R`
- Time variable: `tt_pfs_months`
- Event variable: `pfs_event`
- Event definition: `recurrence_event == 1 | mets_event == 1 | death_event == 1`
- Censoring/competing-event rule: Patients without recurrence, metastasis, or death are censored at last known alive date through the component times.
- Runtime source fields: `recurrence1`, `recurrence1_date`, `mets_progression`, `mets_progression_date`, `dod`, `last_known_alive_date`, `treatment_date`.
- Reviewer question addressed: Clarifies whether death and metastatic progression are included in PFS and repairs the prior metastatic-progression omission.
- Action: rederived and verified in the regenerated full and restricted analytic datasets.

### PFS-2

- Outcome label: PFS-2 Probability / freedom from second recurrence after salvage.
- Reviewer-response definition: Time from salvage treatment after first local recurrence to second local recurrence.
- Code path: `scripts/data_helper/data_derivation.R`; Objective 3 analysis code.
- Time variable: `tt_pfs2_months`
- Event variable: `pfs2_event`
- Event definition: Second local recurrence after salvage treatment.
- Censoring/competing-event rule: Death before second local recurrence and last known alive without second local recurrence are censoring events in the current implementation.
- Runtime source fields: `recurrence1`, `recurrence1_treatment`, `recurrence1_treatment_date`, `recurrence2`, `recurrence2_date`, `dod`, `last_known_alive_date`.
- Reviewer question addressed: Protects the repeat-radiation endpoint from being confused with Objective 1 PFS.
- Action: protect from accidental Objective 1 PFS changes; do not reanalyze Objective 3 unless explicitly requested.

### Local Recurrence

- Outcome label: Local Recurrence-Free Probability.
- Reviewer-response definition: Time from primary treatment to first local recurrence.
- Code path: `scripts/workflow/objective_1_primary_outcomes.R`; `scripts/analysis/survival_outcomes.R`.
- Time variable: `tt_recurrence_months`
- Event variable: `recurrence_event`
- Event definition: `recurrence1 == "Y"` after derived event construction.
- Censoring/competing-event rule: Cox/KM model censors patients without local recurrence at last known alive follow-up; cumulative-incidence support treats death before local recurrence as a competing event.
- Runtime source fields: `recurrence1`, `recurrence1_date`, `last_known_alive_date`, `treatment_date`, `dod`.
- Reviewer question addressed: Replaces unequal-follow-up-sensitive crude/binary inference with time-to-event Cox-led output.
- Action: lead with adjusted Cox model; retain event counts and cumulative-incidence summaries as supportive context.

### Metastatic Progression

- Outcome label: Metastasis-Free Survival Probability.
- Reviewer-response definition: Time from primary treatment to metastatic progression.
- Code path: `scripts/workflow/objective_1_primary_outcomes.R`; `scripts/analysis/survival_outcomes.R`.
- Time variable: `tt_mets_months`
- Event variable: `mets_event`
- Event definition: `mets_progression == "Y"` after derived event construction.
- Censoring/competing-event rule: Cox/KM model censors patients without metastatic progression at last known alive follow-up; cumulative-incidence support treats death before metastasis as a competing event.
- Runtime source fields: `mets_progression`, `mets_progression_date`, `last_known_alive_date`, `treatment_date`, `dod`.
- Reviewer question addressed: Replaces unequal-follow-up-sensitive crude/binary inference with time-to-event Cox-led output.
- Action: lead with adjusted Cox model; retain event counts and cumulative-incidence summaries as supportive context.

### Vision Change

- Outcome label: Latest visual-acuity change.
- Reviewer-response definition: `initial_vision - last_vision`, with recurrence-specific pretreatment vision where applicable.
- Code path: `scripts/data_helper/data_derivation.R`; Objective 2 safety/vision workflow.
- Time variable: `last_followup` is used as the timing field associated with the latest visual-acuity assessment.
- Event variable: Not applicable; continuous visual-acuity change endpoint.
- Event definition: Not applicable; the endpoint is continuous change in visual acuity.
- Censoring/competing-event rule: Not a time-to-event endpoint; minimum-follow-up sensitivity must rerun the treatment-effect model if reported as inferential.
- Runtime source fields: `initial_vision`, `last_vision`, `last_followup`, `recurrence1_pretreatment_vision`, `treatment_date`.
- Reviewer question addressed: Clarifies timing for latest visual-acuity assessment and avoids implying standardized longitudinal visual-field data.
- Action: report available latest-VA timing and complete minimum-follow-up VA sensitivity before making inferential claims.

### Tumor-Height Change

- Outcome label: Tumor-height change.
- Reviewer-response definition: Follow-up height minus baseline height, using recurrence pretreatment height where applicable.
- Code path: `scripts/data_helper/data_derivation.R`; `scripts/analysis/tumor_height_analysis.R`.
- Time variable: `last_height_date`.
- Event variable: Not applicable; continuous height-change endpoint.
- Event definition: Not applicable; the endpoint is continuous change in tumor height.
- Censoring/competing-event rule: Not a time-to-event endpoint; latest observed height is vulnerable to unequal imaging timing.
- Runtime source fields: `initial_tumor_height`, `last_height`, `last_height_date`, `recurrence1_pretreatment_height`, `treatment_date`.
- Reviewer question addressed: Flags unequal assessment timing as a required limitation or sensitivity surface.
- Action: summarize treatment-to-height-assessment timing and demote comparative claims if timing is imbalanced.

### Adverse Events

- Outcome label: Recorded adverse-event burden.
- Reviewer-response definition: Toxicity fields as collected in the analytic dataset.
- Code path: Objective 2 safety/toxicity workflow.
- Time variable: Follow-up timing fields where available; current toxicity outputs are not time-to-toxicity endpoints.
- Event variable: `retinopathy_burden_event`, `nvg_burden_event`, `srd_burden_event`, and related toxicity fields where present.
- Event definition: Recorded burden flags in the analytic dataset.
- Censoring/competing-event rule: Not currently censoring-aware incidence; interpret as recorded burden unless dated toxicity fields are identified.
- Runtime source fields: `retinopathy`, `nvg`, `srd`, adverse-event and follow-up fields in the analytic dataset.
- Reviewer question addressed: Separates recorded adverse-event burden from standardized graded or time-to-toxicity claims.
- Action: clarify whether severity grade, timing, and radiation-induced attribution exist before reviewer-facing toxicity wording.

### Dosimetry/Proximity Availability

- Outcome label: Dose/proximity covariate availability.
- Reviewer-response definition: Data availability audit across derived runtime data, raw/source columns, and data dictionary.
- Code path: `scripts/tools/peer_review_followup_audit.R`.
- Time variable: Not applicable.
- Event variable: Not applicable.
- Event definition: Not applicable.
- Censoring/competing-event rule: Not applicable.
- Runtime source fields: PBT treatment-detail fields, GK plan fields, macula/fovea/proximity fields, source spreadsheet columns checked by the audit tool.
- Reviewer question addressed: Documents whether requested dose/proximity/visual-field covariates exist before claiming they are unavailable.
- Action: distinguish unavailable analytic fields from variables requiring new chart or treatment-plan review; state missing dose/proximity/visual-field data as limitations.

---

## Final Tim-Review Checklist

- [x] Every `T-NICK-XX` item has a non-placeholder status.
- [x] Every completed result has a cohort, estimate, diagnostic status, and evidence alias or explicit coauthor-owned limitation.
- [x] PFS definition is identical across code, methods, figures, tables, and response.
- [x] Recurrence and metastasis reviewer-facing inference is Cox-led.
- [x] Visible log-rank p-values are removed from Objective 1 full/restricted KM plots.
- [x] KM plot display is capped at 216 months (`SURVIVAL_XAXIS_MAX_MONTHS` = 216) while Cox models retain full follow-up.
- [x] RMST omitted from reviewer-facing response emphasis.
- [x] Visual-acuity response states timing and minimum-follow-up sensitivity honestly; Tim/Angie choose which cutoff, if any, to emphasize.
- [x] Propensity-score feasibility documented; no PS treatment-effect sensitivity planned unless Tim/Angie explicitly request the extra work.
- [x] Dose/proximity/visual-field limitations are stated without implying unknown hidden supplemental data.
- [x] Language avoids equivalence, noninferiority, simulated randomization, and unsupported broad vision-preservation claims.
- [x] Tumor-height comparative regression demoted/caveated per T-NICK-10; Tim/Angie choose final manuscript placement.
