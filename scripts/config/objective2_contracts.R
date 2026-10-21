# =============================================================================
# OBJECTIVE 2 SAFETY/TOXICITY CONTRACTS
# =============================================================================
# These mappings define the Objective 2 source endpoints and the Objective 0
# burden fields consumed by safety analyses. Keep the mapping separate from the
# downstream input registry because it also drives field resolution and derivation.

OBJECTIVE2_TOXICITY_ENDPOINT_MAP <- tibble::tribble(
    ~source_field, ~analysis_field, ~endpoint_label,
    "retinopathy", "retinopathy_burden_event", "Radiation Retinopathy",
    "nvg", "nvg_burden_event", "Neovascular Glaucoma",
    "srd", "srd_burden_event", "Serous Retinal Detachment"
)
