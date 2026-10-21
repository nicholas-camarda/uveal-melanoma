# Centralized color palette utilities

#' Return a qualitative palette of up to 8 colors with project-consistent hues
#' @param n Number of colors requested
#' @return Character vector of hex colors
get_qualitative_palette <- function(n) {
    base <- c(
        "#0072B5FF", # blue
        "#BC3C29FF", # red
        "#E18727FF", # orange
        "#20854EFF", # green
        "#7876B1FF", # purple
        "#6F99ADFF", # steel
        "#FFDC91FF", # sand
        "#EE4C97FF" # pink
    )
    if (n <= length(base)) {
        return(base[seq_len(n)])
    }
    # Fallback to RColorBrewer Set1 if more are needed
    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
        return(RColorBrewer::brewer.pal(min(n, 12), "Set3")[seq_len(n)])
    }
    return(rep(base, length.out = n))
}

#' Fill Unknown Palette Levels Without Failing on Short Fallbacks
#'
#' @param out Named character vector containing known colors and `NA` for
#'   unknown levels.
#' @param mapping Named character vector of canonical colors already reserved.
#'
#' @return Named character vector with all `NA` entries replaced.
fill_unknown_palette_levels <- function(out, mapping) {
    if (!any(is.na(out))) {
        return(out)
    }

    needed <- sum(is.na(out))
    fill <- setdiff(get_qualitative_palette(max(length(out), needed)), unname(mapping))
    if (length(fill) < needed) {
        fill <- rep(get_qualitative_palette(max(needed, 1L)), length.out = needed)
    }

    out[is.na(out)] <- fill[seq_len(needed)]
    out
}

#' Consistent colors for GEP groups
#'
#' One color per GEP group, shared by every figure: the Class x PRAME detail groups
#' (`GEP_DETAIL_LEVELS`) and the Class 1 / Class 2 availability groups.
#' @param levels Character vector of GEP group or detail labels
#' @return Named vector mapping levels to colors
get_gep_class_palette <- function(levels) {
    mapping <- c(
        "Class 1" = "#0072B5FF", # blue
        "Class 2" = "#BC3C29FF", # red
        "Class 1, PRAME negative" = "#0072B5FF",
        "Class 1, PRAME positive" = "#20854EFF",
        "Class 2, PRAME negative" = "#E18727FF",
        "Class 2, PRAME positive" = "#BC3C29FF",
        "Class known, PRAME not reported" = "#7876B1FF",
        "Discordant result" = "#FFDC91FF",
        "Failed or unknown result" = "#6F99ADFF",
        "Not tested" = "#9E9E9EFF"
    )
    lv <- unique(as.character(levels))
    out <- mapping[lv]
    out <- fill_unknown_palette_levels(out, mapping)
    names(out) <- lv
    out
}

#' Consistent colors for PRAME status
#' @param levels Character vector of PRAME levels (e.g., c("Positive","Negative","Unknown"))
#' @return Named vector mapping levels to colors
get_prame_palette <- function(levels) {
    mapping <- c(
        "Positive" = "#BC3C29FF", # red
        "Negative" = "#20854EFF", # green
        "Unknown"  = "#6F99ADFF" # steel
    )
    lv <- unique(as.character(levels))
    out <- mapping[lv]
    out <- fill_unknown_palette_levels(out, mapping)
    names(out) <- lv
    out
}

#' Consistent colors for treatment group
#' @param levels Character vector of treatment groups (e.g., c("Plaque","GKSRS","Other"))
#' @return Named vector mapping levels to colors
get_treatment_palette <- function(levels) {
    mapping <- c(
        "PBT" = "#0072B5FF", # blue
        "Plaque" = "#0072B5FF", # legacy alias of PBT
        "GKSRS"  = "#BC3C29FF", # red
        "Other"  = "#E18727FF" # orange
    )
    lv <- unique(normalize_treatment_group_values(levels))
    out <- mapping[lv]
    out <- fill_unknown_palette_levels(out, mapping)
    names(out) <- lv
    out
}

#' Consistent colors for recurrence status (Yes/No)
#' Uses colors distinct from treatment palette (blue/red)
#' @param levels Character vector of recurrence levels (e.g., c("Yes","No"))
#' @return Named vector mapping levels to colors
get_recurrence_palette <- function(levels) {
    mapping <- c(
        "No"  = "#20854EFF", # green - no recurrence (good outcome)
        "Yes" = "#7876B1FF"  # purple - recurrence
    )
    lv <- unique(as.character(levels))
    out <- mapping[lv]
    out <- fill_unknown_palette_levels(out, mapping)
    names(out) <- lv
    out
}

#' Dispatcher to choose palette based on variable semantics
#' @param variable_name Name of the variable being visualized (e.g., 'gep_group','prame_status','treatment_group')
#' @param levels Character vector of factor levels in plotting order
#' @return Character vector of colors in the same order as levels
get_palette_by_variable <- function(variable_name, levels) {
    lv <- if (identical(variable_name, "treatment_group")) {
        unique(normalize_treatment_group_values(levels))
    } else {
        unique(as.character(levels))
    }
    pal <- switch(variable_name,
        "gep_group" = get_gep_class_palette(lv),
        "biopsy1_gep" = get_gep_class_palette(lv),
        "prame_status" = get_prame_palette(lv),
        "treatment_group" = get_treatment_palette(lv),
        "recurrence1" = get_recurrence_palette(lv),
        # default fallback
        get_qualitative_palette(length(lv))
    )
    # Ensure order corresponds to levels
    if (!is.null(names(pal))) {
        pal[lv]
    } else {
        pal
    }
}
