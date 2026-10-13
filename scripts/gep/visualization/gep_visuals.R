# GEP Visualization Functions
# Visualization-only: calibration, discrimination, decision curves, integrated visuals

#' Create collapsed MFS survival curves for selected GEP strata
#'
#' Builds reader-facing companion KM plots that collapse PRAME subgroups into
#' Class 1 and Class 2 while optionally retaining one or both no-GEP strata.
#'
#' @param data Data frame with survival data and GEP labels.
#' @param prefix Filename prefix.
#' @param dataset_name Dataset label used in the subtitle.
#' @param km_output_dir Directory used for the saved KM figure.
#' @param include_failed_indeterminate Logical indicating whether
#'   `GEP Failed/Indeterminate` should be retained as a separate stratum.
#' @param subtitle_suffix Additional subtitle text describing the collapsed view.
#' @param output_filename Filename written to `km_output_dir`.
#' @param return_plot When `TRUE`, returns the assembled plot objects for
#'   verification or testing instead of only writing files.
#' @param save_plot When `FALSE`, skips writing the PNG to disk.
#' @param time_var Character name of the endpoint-specific follow-up-time
#'   column. Defaults to the centralized incident-MFS analysis field.
#' @param event_var Character name of the endpoint-specific event column.
#' @param display_group_var Optional precomputed display-group column.
#' @param display_levels Optional display order for the groups.
#' @param show_p_value Logical; when `FALSE`, omit the p-value annotation from
#'   the plot while retaining the fit in the returned contract.
#' @param surv_fit Optional precomputed `survfit` object to reuse for the plot.
#'
#' @return Invisibly returns `NULL` after saving plots, or a list containing the
#'   plot, shared fit, plot data, and p-value-annotation flag when
#'   `return_plot = TRUE`.
create_mfs_collapsed_survival_curves <- function(data, prefix, km_output_dir, dataset_name = "GEP Validation", include_failed_indeterminate = FALSE, subtitle_suffix, output_filename, return_plot = FALSE, save_plot = TRUE, time_var = "tt_mets_months_analysis", event_var = "mets_event_analysis", display_group_var = NULL, display_levels = NULL, show_p_value = TRUE, surv_fit = NULL) {
    target_levels <- display_levels %||% c("Class 1", "Class 2", "GEP Not Tested")
    if (is.null(display_levels) && isTRUE(include_failed_indeterminate)) {
        target_levels <- c(target_levels, "GEP Failed/Indeterminate")
    }

    plot_label <- if (isTRUE(include_failed_indeterminate)) {
        "collapsed four-group"
    } else {
        "simplified three-group"
    }
    logger::log_info(sprintf("Creating %s MFS survival curves", plot_label))

    plot_data <- prepare_incident_mfs_km_data(data)
    if (!identical(time_var, "tt_mets_months_analysis") || !identical(event_var, "mets_event_analysis")) {
        stop(
            "Collapsed Objective 4 MFS KM requires tt_mets_months_analysis and mets_event_analysis.",
            call. = FALSE
        )
    }
    group_var <- "gep_km_simple"
    if (is.null(display_group_var)) {
        plot_data <- plot_data %>%
            dplyr::mutate(
                gep_km_simple = dplyr::case_when(
                    .data$biopsy1_gep == "GEP Not Tested" ~ "GEP Not Tested",
                    isTRUE(include_failed_indeterminate) & .data$biopsy1_gep == "GEP Failed/Indeterminate" ~ "GEP Failed/Indeterminate",
                    .data$gep_class_simple %in% c("Class 1", "Class 2") ~ as.character(.data$gep_class_simple),
                    TRUE ~ NA_character_
                )
            )
    } else {
        if (!display_group_var %in% names(plot_data)) {
            stop(sprintf("MFS display group column is missing: %s", display_group_var), call. = FALSE)
        }
        group_var <- display_group_var
    }

    if (!time_var %in% names(plot_data) || !event_var %in% names(plot_data)) {
        stop(sprintf("MFS analysis columns are missing: %s, %s", time_var, event_var), call. = FALSE)
    }

    plot_data <- plot_data %>%
        dplyr::filter(
            !is.na(.data[[group_var]]),
            !is.na(.data[[time_var]]),
            !is.na(.data[[event_var]]),
            is.finite(suppressWarnings(as.numeric(.data[[time_var]]))),
            suppressWarnings(as.numeric(.data[[time_var]])) >= 0
        ) %>%
        dplyr::mutate(
            !!group_var := factor(as.character(.data[[group_var]]), levels = target_levels)
        ) %>%
        dplyr::filter(!is.na(.data[[group_var]])) %>%
        as.data.frame()

    present_levels <- levels(droplevels(plot_data[[group_var]]))
    if (nrow(plot_data) == 0 || length(unique(stats::na.omit(plot_data[[group_var]]))) < 2) {
        logger::log_warn(sprintf("Insufficient data/groups for %s MFS survival curves", plot_label))
        return(invisible(NULL))
    }

    logger::log_info(sprintf(
        "Collapsed MFS KM groups (%s): %s",
        plot_label,
            paste(
            capture.output(print(table(plot_data[[group_var]], useNA = "no"))),
            collapse = " "
        )
    ))

    surv_formula <- stats::as.formula(sprintf("survival::Surv(%s, %s) ~ %s", time_var, event_var, group_var))
    surv_fit <- surv_fit %||% survival::survfit(surv_formula, data = plot_data)
    surv_fit$call$formula <- surv_formula

    raw_max_time <- max(suppressWarnings(as.numeric(plot_data[[time_var]])), na.rm = TRUE)
    max_time <- min(raw_max_time, SURVIVAL_XAXIS_MAX_MONTHS)
    base_by <- if (max_time <= 60) 6 else 12
    x_breaks <- seq(0, ceiling(max_time / base_by) * base_by, by = base_by)
    plot_scale <- SURVIVAL_PLOT_SCALE
    color_palette <- get_palette_by_variable("biopsy1_gep", present_levels)
    surv_plot <- survminer::ggsurvplot(
        fit = surv_fit,
        data = plot_data,
        palette = unname(color_palette),
        risk.table = TRUE,
        conf.int = FALSE,
        pval = show_p_value,
        pval.size = 6 * plot_scale,
        title = paste("Kaplan-Meier Survival Curves", "Metastasis-Free Survival Probability", sep = "\n"),
        subtitle = paste(
            c(
                if (!is.null(dataset_name)) paste("Cohort:", dataset_name) else NULL,
                subtitle_suffix
            ),
            collapse = "\n"
        ),
        xlab = "Time (months)",
        ylab = "Metastasis-Free Survival Probability",
        risk.table.height = 0.18,
        ggtheme = ggplot2::theme_minimal(),
        break.time.by = base_by,
        xlim = c(0, max(x_breaks)),
        ylim = c(0, 1),
        legend.labs = present_levels,
        risk.table.y.text = TRUE,
        tables.y.text = TRUE,
        risk.table.title = "Number at risk",
        font.x = 14 * plot_scale,
        font.y = 14 * plot_scale,
        font.tickslab = 12 * plot_scale,
        font.legend = 14 * plot_scale,
        censor.size = 7 * plot_scale,
        size = 1.2 * plot_scale
    )

    surv_plot$plot <- remove_plot_scales(surv_plot$plot, aesthetics = c("y"))
    surv_plot$table <- remove_plot_scales(surv_plot$table, aesthetics = c("y"))

    surv_plot$plot <- surv_plot$plot +
        ggplot2::guides(color = ggplot2::guide_legend(ncol = 1, byrow = TRUE)) +
        ggplot2::theme(
            legend.position = "bottom",
            legend.box = "vertical",
            legend.text = ggplot2::element_text(size = 16 * plot_scale, color = "black"),
            legend.title = ggplot2::element_text(size = 16 * plot_scale, color = "black"),
            axis.title = ggplot2::element_text(size = 18 * plot_scale, color = "black"),
            axis.title.x = ggplot2::element_text(size = 18 * plot_scale, color = "black", face = "bold", margin = ggplot2::margin(t = 15, r = 0, b = 0, l = 0)),
            axis.title.y = ggplot2::element_text(size = 18 * plot_scale, color = "black", face = "bold", margin = ggplot2::margin(t = 0, r = 6, b = 0, l = 0)),
            axis.text = ggplot2::element_text(size = 14 * plot_scale, color = "black"),
            axis.text.x = ggplot2::element_text(color = "black"),
            axis.text.y = ggplot2::element_text(color = "black"),
            axis.line = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(size = 18 * plot_scale, face = "bold"),
            plot.subtitle = ggplot2::element_text(size = 14 * plot_scale),
            plot.margin = ggplot2::margin(t = 10, r = 10, b = 6, l = 4)
        ) +
        ggplot2::scale_y_continuous(
            limits = c(0, 1),
            breaks = seq(0, 1, by = 0.1),
            labels = function(x) x * 100,
            name = "Metastasis-Free Survival Probability (%)"
        ) +
        ggplot2::labs(x = "Time (months)") +
        ggplot2::geom_hline(yintercept = 0.5, linetype = "solid", color = "black", linewidth = 0.9, alpha = 0.35)

    surv_plot$table <- surv_plot$table +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.title = ggplot2::element_text(size = 12 * plot_scale, color = "black"),
            axis.text.y = ggplot2::element_text(size = 12 * plot_scale, color = "black"),
            axis.text.x = ggplot2::element_text(size = 12 * plot_scale, color = "black"),
            strip.text = ggplot2::element_text(size = 12 * plot_scale, color = "black"),
            plot.margin = ggplot2::margin(t = 10, r = 10, b = 4, l = 10)
        )

    clean_strata_label <- function(x) {
        x_chr <- as.character(x)
        ifelse(grepl("=", x_chr), sub("^[^=]*=", "", x_chr), x_chr)
    }

    remap_risk_table_rows <- function(table_frame) {
        if (is.null(table_frame) || nrow(table_frame) == 0) {
            return(table_frame)
        }

        mapped_labels <- NULL
        if ("strata" %in% names(table_frame)) {
            candidate_labels <- clean_strata_label(table_frame$strata)
            if (all(stats::na.omit(candidate_labels) %in% present_levels)) {
                mapped_labels <- candidate_labels
            }
        }

        if (is.null(mapped_labels) && "y" %in% names(table_frame)) {
            y_as_integer <- suppressWarnings(as.integer(as.character(table_frame$y)))
            if (length(y_as_integer) > 0 && any(!is.na(y_as_integer))) {
                candidate_labels <- present_levels[y_as_integer]
                if (all(stats::na.omit(candidate_labels) %in% present_levels)) {
                    mapped_labels <- candidate_labels
                }
            }

            if (is.null(mapped_labels)) {
                candidate_labels <- clean_strata_label(table_frame$y)
                if (all(stats::na.omit(candidate_labels) %in% present_levels)) {
                    mapped_labels <- candidate_labels
                }
            }
        }

        if (is.null(mapped_labels)) {
            return(table_frame)
        }

        mapped_factor <- factor(mapped_labels, levels = present_levels)

        if ("strata" %in% names(table_frame)) {
            table_frame$strata <- mapped_factor
        }
        if ("y" %in% names(table_frame)) {
            table_frame$y <- mapped_factor
        }

        ordering_time <- if ("time" %in% names(table_frame)) {
            table_frame$time
        } else if ("x" %in% names(table_frame)) {
            table_frame$x
        } else {
            seq_len(nrow(table_frame))
        }

        table_frame[order(mapped_factor, ordering_time), , drop = FALSE]
    }

    if (length(surv_plot$table$layers) > 0) {
        for (i in seq_along(surv_plot$table$layers)) {
            if ("GeomText" %in% class(surv_plot$table$layers[[i]]$geom)) {
                surv_plot$table$layers[[i]]$aes_params$size <- 3.4 * plot_scale
            }

            if (!is.null(surv_plot$table$layers[[i]]$data)) {
                surv_plot$table$layers[[i]]$data <- remap_risk_table_rows(surv_plot$table$layers[[i]]$data)
            }
        }
    }

    surv_plot$table$data <- remap_risk_table_rows(surv_plot$table$data)
    surv_plot$table$mapping <- ggplot2::aes(
        x = time,
        y = strata,
        label = llabels,
        shape = strata
    )

    surv_plot$table <- surv_plot$table +
        ggplot2::scale_y_discrete(
            limits = rev(present_levels),
            expand = ggplot2::expansion(mult = c(0.18, 0.18))
        )

    combined_km <- cowplot::plot_grid(
        surv_plot$plot,
        surv_plot$table,
        ncol = 1,
        align = "v",
        rel_heights = c(0.78, 0.22)
    )

    simplified_plot_width <- 13
    simplified_plot_height <- min(
        KM_MAX_HEIGHT,
        max(
            SURVIVAL_PLOT_HEIGHT * 1.35,
            KM_BASE_HEIGHT + max(0, length(present_levels) - 2) * KM_HEIGHT_PER_STRATUM + 2.0
        )
    )

    output_path <- NULL
    if (save_plot) {
        target_km_dir <- ensure_output_dir(km_output_dir)
        output_path <- file.path(target_km_dir, paste0(prefix, output_filename))
        ggplot2::ggsave(
            output_path,
            combined_km,
            width = simplified_plot_width,
            height = simplified_plot_height,
            dpi = PLOT_DPI,
            bg = "white"
        )
        logger::log_info(sprintf("Collapsed MFS survival curves saved: %s", output_path))
    }

    if (return_plot) {
        return(list(
            plot = surv_plot,
            combined_plot = combined_km,
            plot_data = plot_data,
            fit = surv_fit,
            p_value_annotation = isTRUE(show_p_value),
            present_levels = present_levels,
            output_path = output_path
        ))
    }

    invisible(NULL)
}

#' Create simplified MFS survival curves for Class 1, Class 2, and GEP Not Tested
#'
#' Adds a reader-facing companion KM plot that collapses PRAME subgroups,
#' retains GEP Not Tested, and excludes GEP Failed/Indeterminate.
#'
#' @param data Data frame with survival data and GEP labels.
#' @param output_dir Output directory for saved plots.
#' @param prefix Filename prefix.
#' @param dataset_name Dataset label used in the subtitle.
#' @param km_output_dir Directory used for the saved simplified KM figure.
#' @param return_plot When `TRUE`, returns the assembled plot objects for
#'   verification or testing instead of only writing files.
#' @param save_plot When `FALSE`, skips writing the PNG to disk.
#'
#' @return Invisibly returns `NULL` after saving plots, or a list of plot
#'   objects when `return_plot = TRUE`.
create_mfs_simplified_survival_curves <- function(data, output_dir, prefix, dataset_name = "GEP Validation", km_output_dir = output_dir, return_plot = FALSE, save_plot = TRUE) {
    create_mfs_collapsed_survival_curves(
        data = data,
        prefix = prefix,
        dataset_name = dataset_name,
        km_output_dir = km_output_dir,
        include_failed_indeterminate = FALSE,
        subtitle_suffix = "Simplified display: Class 1, Class 2, and GEP Not Tested",
        output_filename = "mfs_simplified_gep_km.png",
        return_plot = return_plot,
        save_plot = save_plot
    )
}

#' Build a shared cumulative-incidence curve plot
#'
#' @param ci_obj A `tidycuminc` competing-risk fit.
#' @param outcome Outcome level to display.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param color_title Legend title.
#' @param palette Named color vector.
#' @param xlim Numeric x-axis limits.
#' @param caption Optional initial caption.
#'
#' @return A styled `ggplot` cumulative-incidence plot.
build_cif_curve_plot <- function(ci_obj, outcome, title, subtitle, xlab, ylab, color_title, palette, xlim, caption = NULL) {
    # CIFs are bounded at zero. A small view-only margin keeps a true zero
    # curve visible without changing estimates, limits, or tick labels.
    y_baseline_padding <- 0.02
    p <- ggsurvfit::ggcuminc(ci_obj, outcome = outcome)
    p$layers <- lapply(p$layers, function(layer) {
        layer$aes_params$na.rm <- TRUE
        layer$geom_params$na.rm <- TRUE
        layer$stat_params$na.rm <- TRUE
        layer
    })
    p <- remove_plot_scales(p, aesthetics = c("colour", "color", "x", "y")) +
        ggplot2::labs(
            title = title,
            subtitle = subtitle,
            x = xlab,
            y = ylab,
            color = color_title,
            caption = caption
        ) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            plot.background = ggplot2::element_rect(fill = "white"),
            panel.background = ggplot2::element_rect(fill = "white"),
            plot.title = ggplot2::element_text(size = 16, face = "bold", lineheight = 1.05),
            plot.subtitle = ggplot2::element_text(size = 12, color = "darkgray"),
            plot.caption = ggplot2::element_text(size = 10.5, color = "darkgray", hjust = 0, lineheight = 1.1),
            legend.position = "bottom",
            legend.title = ggplot2::element_text(face = "bold", size = 12),
            legend.text = ggplot2::element_text(size = 11),
            axis.title = ggplot2::element_text(size = 13.5),
            axis.text = ggplot2::element_text(size = 11.5),
            axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
            axis.ticks.x = ggplot2::element_line(color = "black", linewidth = 0.5)
        ) +
        ggplot2::scale_color_manual(values = palette) +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.02))) +
        ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = c(0, 0.04)),
            labels = scales::label_percent(accuracy = 1)
        ) +
        ggplot2::coord_cartesian(
            xlim = xlim,
            ylim = c(-y_baseline_padding, 1)
        )

    p
}
