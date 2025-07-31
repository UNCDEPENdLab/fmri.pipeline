# Helper functions for structured logging

#' Create or retrieve a subject/session level logger
#'
#' @param base_logger Name of parent logger (e.g., "glm_pipeline/l1_setup")
#' @param gpa glm_pipeline_arguments object containing logging settings
#' @param id Subject identifier
#' @param session Session identifier
#' @param prefix Prefix used to locate output log files (e.g., "setup_l1")
#' @keywords internal
get_subject_logger <- function(base_logger, gpa, id, session, prefix) {
  lg_name <- glue::glue("{base_logger}/subject_{id}_session_{session}")
  lg <- lgr::get_logger(lg_name)
  lg$set_threshold(gpa$lgr_threshold)
  lg$set_context(id = id, session = session)

  if (isTRUE(gpa$log_json)) {
    json_dir <- dirname(gpa$output_locations[[glue::glue('{prefix}_log_json')]])
    json_file <- file.path(json_dir, glue::glue("subject_{id}_session_{session}.json"))
    if (!"json" %in% names(lg$appenders)) {
      lg$add_appender(lgr::AppenderJson$new(json_file), name = "json")
    }
  }

  if (isTRUE(gpa$log_txt)) {
    txt_dir <- dirname(gpa$output_locations[[glue::glue('{prefix}_log_txt')]])
    txt_file <- file.path(txt_dir, glue::glue("subject_{id}_session_{session}.log"))
    if (!"txt" %in% names(lg$appenders)) {
      lg$add_appender(lgr::AppenderFile$new(txt_file), name = "txt")
    }
  }

  return(lg)
}

#' Log a warning and emit a base R warning
#' @param lg Logger object
#' @param message Message text passed to glue::glue
#' @keywords internal
log_warn <- function(lg, message, ...) {
  msg <- glue::glue(message, ...)
  lg$warn(msg)
  warning(msg, call. = FALSE)
}

#' Log an error and stop execution
#' @param lg Logger object
#' @param message Message text passed to glue::glue
#' @keywords internal
log_error <- function(lg, message, ...) {
  msg <- glue::glue(message, ...)
  lg$error(msg)
  stop(msg, call. = FALSE)
}

