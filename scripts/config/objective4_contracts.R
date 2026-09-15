# =============================================================================
# OBJECTIVE 4 GEP DERIVATION CONTRACTS
# =============================================================================
# This contract protects row-wise consistency for imported probabilities and
# horizon fields. All endpoint times are untruncated months. Event types are
# horizon-independent: MFS uses 0=censored/1=metastasis; MSS uses
# 0=censored/1=melanoma death/2=other-cause death.

OBJECTIVE4_GEP_DERIVATION_CONTRACT <- tibble::tribble(
    ~outcome, ~horizon_years, ~horizon_months, ~source_probability_field, ~expected_survival_field, ~predicted_risk_field, ~event_field, ~event_type_field, ~time_field, ~time_unit, ~eligibility_field,
    "mfs", 5, 60, "biopsy1_gep_mfs", "expected_mfs_5yr", "predicted_mfs_risk_5yr", "metastasis_by_5yr", "mfs_event_type", "tt_mets_months_analysis", "months", "mfs_analysis_eligible",
    "mfs", 7, 84, "biopsy1_gep_mfs", "expected_mfs_7yr", "predicted_mfs_risk_7yr", "metastasis_by_7yr", "mfs_event_type", "tt_mets_months_analysis", "months", "mfs_analysis_eligible",
    "mfs", 10, 120, "biopsy1_gep_mfs", "expected_mfs_10yr", "predicted_mfs_risk_10yr", "metastasis_by_10yr", "mfs_event_type", "tt_mets_months_analysis", "months", "mfs_analysis_eligible",
    "mss", 5, 60, "biopsy1_gep_mss", "expected_mss_5yr", "predicted_mss_risk_5yr", "melanoma_death_by_5yr", "mss_event_type", "tt_death_months", "months", "mss_analysis_eligible",
    "mss", 7, 84, "biopsy1_gep_mss", "expected_mss_7yr", "predicted_mss_risk_7yr", "melanoma_death_by_7yr", "mss_event_type", "tt_death_months", "months", "mss_analysis_eligible",
    "mss", 10, 120, "biopsy1_gep_mss", "expected_mss_10yr", "predicted_mss_risk_10yr", "melanoma_death_by_10yr", "mss_event_type", "tt_death_months", "months", "mss_analysis_eligible"
)
