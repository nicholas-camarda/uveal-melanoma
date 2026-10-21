# objective4-mss-competing-risk-alignment Specification

## Purpose
Define the active Objective 4 contract for melanoma-specific survival validation: competing-risk cumulative incidence is the estimand for observed risk and for discrimination.

## Requirements
### Requirement: Objective 4 MSS observed risk SHALL be a competing-risk cumulative incidence
Observed 5-year melanoma-specific risk, in every Objective 4 table and figure, MUST be the Aalen-Johansen cumulative incidence of melanoma death with non-melanoma death as a competing event, never a simple censoring-based or count-based quantity.

#### Scenario: Group tables and curves use the same estimator
- **WHEN** Objective 4 reports observed MSS risk for a GEP group, in a table, a curve or the observed-versus-imported figure
- **THEN** the value is the Aalen-Johansen cumulative incidence at the horizon from the same prepared MSS dataset
- **AND** no value is carried past the group's last follow-up

### Requirement: Objective 4 MSS discrimination SHALL be competing-risk consistent
#### Scenario: MSS discrimination treats other-cause death as a competing event
- **WHEN** Objective 4 reports MSS discrimination
- **THEN** it is a cumulative/dynamic time-dependent AUC at the horizon in which patients who died of another cause are counted among patients without melanoma death
- **AND** the estimand, the number of events and the bootstrap interval are reported with it

### Requirement: Objective 4 SHALL define melanoma death by recorded cause only
#### Scenario: Unknown cause of death is competing in the primary analysis
- **WHEN** a death has recorded cause `Unknown` or `Unknown_Cancer`
- **THEN** the primary analysis counts it as a competing death
- **AND** a labeled sensitivity analysis recomputes the group risks counting these deaths as melanoma deaths
