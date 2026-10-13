# GEP Excel output helpers

#' Write a GEP workbook with automatic column widths
#'
#' Writes one or more data frames to an Excel workbook and applies automatic
#' column sizing so the saved sheets are readable without manual resizing.
#'
#' @param workbook_data A data frame or named list of data frames.
#' @param output_path File path for the workbook.
#' @param freeze_header Logical; if `TRUE`, freeze the header row on each sheet.
#' @param min_width Minimum width used by openxlsx auto-sizing.
#' @param max_width Maximum width used by openxlsx auto-sizing.
#' @return Invisibly returns `output_path`.
write_gep_workbook <- function(workbook_data,
                               output_path,
                               freeze_header = TRUE,
                               min_width = 10,
                               max_width = 60) {
    if (is.data.frame(workbook_data)) {
        workbook_data <- list(Data = workbook_data)
    }

    if (!is.list(workbook_data) || length(workbook_data) == 0) {
        stop("write_gep_workbook() requires a data frame or non-empty named list of data frames")
    }

    if (is.null(names(workbook_data)) || any(!nzchar(names(workbook_data)))) {
        names(workbook_data) <- paste0("Sheet", seq_along(workbook_data))
    }

    write_readable_xlsx(
        workbook_data,
        output_path,
        freeze_header = freeze_header,
        min_width = min_width,
        max_width = max_width
    )
}
