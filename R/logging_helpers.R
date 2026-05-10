# Helper functions for structured logging

#' Create or retrieve a subject level logger
#'
#' @param base_logger Name of parent logger (e.g., "glm_pipeline/l1_setup")
#' @param id Subject identifier
#' @param gpa glm_pipeline_arguments object containing logging settings
#' @param log_prefix Optional filename prefix for subject log files
#' @keywords internal
get_subject_logger <- function(base_logger, id, gpa, log_prefix = "setup_l1_models") {

  # for backwards compatibility, if log directory is not specified, return the base logger
  log_dir <- gpa$output_locations$log_directory
  if (is.null(log_dir)) {
    slg <- lgr::get_logger(base_logger)
    return(slg)
  }

  # create subject log folder if it doesn't exist
  subj_log_folder <- file.path(log_dir, paste0("subj", id))
  if (!dir.exists(subj_log_folder)) dir.create(subj_log_folder, recursive = TRUE)

  # create subject logger
  slg <- lgr::get_logger(glue::glue("{base_logger}/subj{id}"))
  slg$set_threshold(gpa$lgr_threshold)

  # add appenders
  if (isTRUE(gpa$log_json)) {
    appender_name <- glue("{log_prefix}_subj{id}_json")
    if (!appender_name %in% names(slg$appenders)) {
      subj_log_json <- file.path(subj_log_folder, glue("{log_prefix}_subj{id}.json"))
      slg$add_appender(lgr::AppenderJson$new(subj_log_json), name = appender_name)
    }
  }
  if (isTRUE(gpa$log_txt)) {
    appender_name <- glue("{log_prefix}_subj{id}_txt")
    if (!appender_name %in% names(slg$appenders)) {
      subj_log_txt <- file.path(subj_log_folder, glue("{log_prefix}_subj{id}.txt"))
      slg$add_appender(lgr::AppenderFile$new(subj_log_txt), name = appender_name)
    }
  }

  return(slg)

}

#' Add base logger appenders for top-level setup logs
#' This is for backwards compatibility now that logs are moved into a new folder for
#'  new gpa objects
#'
#' @param lg A logger instance
#' @param gpa glm_pipeline_arguments object containing logging settings
#' @param log_txt_path Path to the text log file
#' @param log_json_path Path to the json log file
#' @param txt_appender_name Appender name for text logs
#' @param json_appender_name Appender name for json logs
#' @keywords internal
add_base_logger_appenders <- function(
  lg,
  gpa,
  log_txt_path,
  log_json_path,
  txt_appender_name,
  json_appender_name
) {

  sync_appender <- function(name, enabled, target_path, factory) {
    if (!isTRUE(enabled) || is.null(target_path) || !nzchar(target_path)) {
      return(invisible(NULL))
    }

    existing <- lg$appenders[[name]]
    existing_path <- if (!is.null(existing)) existing$file %||% NULL else NULL
    if (!is.null(existing) && !identical(existing_path, target_path)) {
      lg$remove_appender(name)
      existing <- NULL
    }

    if (is.null(existing)) {
      lg$add_appender(factory(target_path), name = name)
    }

    invisible(NULL)
  }

  sync_appender(
    name = json_appender_name,
    enabled = gpa$log_json,
    target_path = log_json_path,
    factory = function(path) lgr::AppenderJson$new(path)
  )

  sync_appender(
    name = txt_appender_name,
    enabled = gpa$log_txt,
    target_path = log_txt_path,
    factory = function(path) lgr::AppenderFile$new(path)
  )

  invisible(lg)
}
