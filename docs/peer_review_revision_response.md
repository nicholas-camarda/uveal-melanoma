# Peer-Review Statistical Revision Response Evidence

## Endpoint and Claim Audit

### PFS

- Outcome label: Progression-Free Survival (PFS)
- Manuscript wording: Time from primary treatment to first local recurrence, metastatic progression, or death from any cause.
- Code path: `scripts/data_helper/data_derivation.R`; `scripts/workflow/objective_1_primary_outcomes.R`
- Time variable: `tt_pfs_months`
- Event variable: `pfs_event`
- Event definition: `recurrence_event == 1 | mets_event == 1 | death_event == 1`
- Censoring/competing-event rule: Patients without recurrence, metastasis, or death are censored at last known alive date through the component times.
- Runtime source fields: `recurrence1`, `recurrence1_date`, `mets_progression`, `mets_progression_date`, `dod`, `last_known_alive_date`, `treatment_date`
- Reviewer question addressed: Clarifies whether death without documented progression is included and corrects the internal omission of metastatic progression from the PFS composite.
- Action: rederive

### PFS-2

- Outcome label: PFS-2 Probability (Freedom from 2nd Recurrence)
- Manuscript wording: Time from salvage treatment after first local recurrence to second local recurrence.
- Code path: `scripts/data_helper/data_derivation.R`; `scripts/analysis/survival_outcomes.R`
- Time variable: `tt_pfs2_months`
- Event variable: `pfs2_event`
- Event definition: Second local recurrence after salvage treatment.
- Censoring/competing-event rule: Death before second local recurrence and last known alive without second local recurrence are censoring events.
- Runtime source fields: `recurrence1`, `recurrence1_treatment`, `recurrence1_treatment_date`, `recurrence2`, `recurrence2_date`, `dod`, `last_known_alive_date`
- Reviewer question addressed: Protects Objective 3 from accidental import of the Objective 1 PFS composite while still defining the retained endpoint precisely.
- Action: keep

### Local Recurrence

- Outcome label: Local Recurrence-Free Probability
- Manuscript wording: Time from primary treatment to first local recurrence.
- Code path: `scripts/workflow/objective_1_primary_outcomes.R`; `scripts/analysis/survival_outcomes.R`; `scripts/analysis/binary_outcomes.R`
- Time variable: `tt_recurrence_months`
- Event variable: `recurrence_event`
- Event definition: `recurrence1 == "Y"`
- Censoring/competing-event rule: Cox/KM model censors patients without local recurrence at last known alive follow-up; cumulative-incidence support treats death before local recurrence as a competing event.
- Runtime source fields: `recurrence1`, `recurrence1_date`, `last_known_alive_date`, `treatment_date`, `dod`
- Reviewer question addressed: Replaces unequal-follow-up-sensitive logistic inference with time-to-event Cox-led output.
- Action: rederive

### Metastatic Progression

- Outcome label: Metastasis-Free Survival Probability
- Manuscript wording: Time from primary treatment to metastatic progression.
- Code path: `scripts/workflow/objective_1_primary_outcomes.R`; `scripts/analysis/survival_outcomes.R`; `scripts/analysis/binary_outcomes.R`
- Time variable: `tt_mets_months`
- Event variable: `mets_event`
- Event definition: `mets_progression == "Y"`
- Censoring/competing-event rule: Cox/KM model censors patients without metastatic progression at last known alive follow-up; cumulative-incidence support treats death before metastasis as a competing event.
- Runtime source fields: `mets_progression`, `mets_progression_date`, `last_known_alive_date`, `treatment_date`, `dod`
- Reviewer question addressed: Replaces unequal-follow-up-sensitive logistic inference with time-to-event Cox-led output.
- Action: rederive

### Vision Change

- Outcome label: Vision change / latest visual acuity change
- Manuscript wording: Change from initial visual acuity to latest available visual-acuity assessment.
- Code path: `scripts/data_helper/data_derivation.R`; Objective 2 safety/vision workflow.
- Time variable: `last_followup` as the latest visual-acuity follow-up date context.
- Event variable: Not applicable; continuous visual-acuity change endpoint.
- Event definition: `initial_vision - last_vision`, using `recurrence1_pretreatment_vision` for recurrence cases where applicable.
- Censoring/competing-event rule: Not a time-to-event endpoint; minimum-follow-up sensitivity must rerun the treatment-effect model if reported as inferential.
- Runtime source fields: `initial_vision`, `last_vision`, `last_followup`, `recurrence1_pretreatment_vision`, `treatment_date`
- Reviewer question addressed: Clarifies that `last_vision` is the latest acuity value and `last_followup` is the associated timing field.
- Action: keep

### Tumor-Height Change

- Outcome label: Tumor-height change
- Manuscript wording: Change from initial tumor height to latest observed height or recurrence pretreatment height.
- Code path: `scripts/data_helper/data_derivation.R`; `scripts/analysis/tumor_height_analysis.R`
- Time variable: `last_height_date` for latest height assessment timing.
- Event variable: Not applicable; continuous height-change endpoint.
- Event definition: `last_height - initial_tumor_height`, or `recurrence1_pretreatment_height - initial_tumor_height` for local recurrence cases.
- Censoring/competing-event rule: Not a time-to-event endpoint; assessment timing and negative treatment-to-height intervals require audit before reviewer-facing comparative interpretation.
- Runtime source fields: `initial_tumor_height`, `last_height`, `last_height_date`, `recurrence1_pretreatment_height`, `treatment_date`
- Reviewer question addressed: Flags unequal assessment timing as a required feasibility/sensitivity audit surface.
- Action: keep with audit

### Adverse Events

- Outcome label: Radiation adverse events
- Manuscript wording: Recorded adverse-event burden by available follow-up.
- Code path: Objective 2 safety/toxicity workflow.
- Time variable: Follow-up timing fields where available; not currently a time-to-toxicity endpoint.
- Event variable: `retinopathy_burden_event`, `nvg_burden_event`, `srd_burden_event`, and related toxicity fields.
- Event definition: Recorded burden flags in the analytic dataset.
- Censoring/competing-event rule: Not currently censoring-aware incidence; interpret as recorded burden unless time-to-event toxicity data are identified.
- Runtime source fields: `retinopathy`, `nvg`, `srd`, adverse-event and follow-up fields in the analytic dataset.
- Reviewer question addressed: Separates recorded burden from graded/time-to-toxicity claims.
- Action: demote/clarify

### Dosimetry/Proximity Availability

- Outcome label: Dosimetry and tumor-proximity covariate availability
- Manuscript wording: Treatment-detail and proximity adjustment feasibility.
- Code path: On-demand peer-review follow-up/data-availability audit tool.
- Time variable: Not applicable.
- Event variable: Not applicable.
- Event definition: Not applicable.
- Censoring/competing-event rule: Not applicable.
- Runtime source fields: PBT treatment-detail fields, GK plan fields, macula/fovea/proximity fields, source spreadsheet columns checked by the audit tool.
- Reviewer question addressed: Documents whether requested dosimetry/proximity covariates exist before claiming they are absent or unusable.
- Action: audit before final response
