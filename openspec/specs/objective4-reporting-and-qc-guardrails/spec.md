# objective4-reporting-and-qc-guardrails Specification

## Purpose
Define Objective 4 reporting guardrails that keep estimates tied to the follow-up that supports them and make every check and exclusion visible.

## Requirements
### Requirement: Objective 4 SHALL report only the supported horizon
#### Scenario: No estimate past the last follow-up
- **WHEN** a GEP group has no patient followed to the 5-year horizon
- **THEN** the group has no 5-year estimate and the table records why
- **AND** no estimate is carried forward or extrapolated

#### Scenario: Zero events give no interval
- **WHEN** a group has no event by the horizon
- **THEN** its observed risk is 0 and no interval is reported

### Requirement: Objective 4 SHALL give every check a status
#### Scenario: Checks are visible and never silently skipped
- **WHEN** Objective 4 writes the diagnostics workbook
- **THEN** every check has a result and a status of no concern, concern, or not assessed, with a reason
- **AND** a check that could not run is recorded as not assessed
- **AND** a failed analysis stops the run instead of producing a partial workbook

### Requirement: Objective 4 SHALL make every exclusion traceable by patient id
#### Scenario: Diagnostics workbooks travel with the results
- **WHEN** a snapshot is published or the collaborator export is built
- **THEN** the Objective 4 diagnostics workbooks, including their `Record_Issues` sheets, are included beside the primary workbooks
- **AND** each `Record_Issues` row carries the patient `id` used in the analytic dataset
- **AND** the primary workbooks contain no patient identifiers

#### Scenario: Identifiers stay out of Git and public documents
- **WHEN** documentation, tests or committed files describe or exercise the diagnostics workbooks
- **THEN** they contain no real patient identifiers

#### Scenario: Record_Issues reconciles with the flow
- **WHEN** a diagnostics workbook is written
- **THEN** the number of `Record_Issues` rows equals the patients outside the primary set plus the flagged patients counted in `Denominators_Flow`

### Requirement: Objective 4 SHALL use one set of GEP group labels
#### Scenario: Labels come from the policy file
- **WHEN** any table, figure, palette or contract names a GEP group
- **THEN** the name comes from `GEP_GROUP_LEVELS` or `GEP_DETAIL_LEVELS` in the Objective 4 policy module
- **AND** no pooled failed/indeterminate label and no 1A/1B subclass wording appears in reader-facing output
