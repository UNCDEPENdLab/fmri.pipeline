# Helper functions for structured logging

#' Create or retrieve a subject level logger
#'
#' @param base_logger Name of parent logger (e.g., "glm_pipeline/l1_setup")
#' @param gpa glm_pipeline_arguments object containing logging settings
#' @param id Subject identifier
#' @param prefix Prefix used to locate output log files (e.g., "setup_l1")
#' @keywords internal
get_subject_logger <- function(base_logger, gpa, id, prefix) {
  lg_name <- glue::glue("{base_logger}/subject_{id}")
  lg <- lgr::get_logger(lg_name)
  lg$set_threshold(gpa$lgr_threshold)

  log_name <- glue::glue("subject_{id}")
    
  if (isTRUE(gpa$log_json) && prefix == "setup_l1") {
    if (!paste0(log_name, "_json") %in% names(lg$appenders)) {
      folder_name <- glue::glue(gpa$output_locations$feat_sub_directory)
      
      # Insert "logs" before "feat_l1" in the folder path
      folder_parts <- strsplit(folder_name, .Platform$file.sep)[[1]]
      feat_idx <- which(folder_parts == "feat_l1")
      folder_parts <- append(folder_parts, "logs", after = feat_idx - 1)
      folder_name <- do.call(file.path, as.list(folder_parts))
      
      # create folder if it doesn't exist    
      if (!dir.exists(folder_name)) {
        dir.create(folder_name, recursive = TRUE)
      }

      # Create JSON appender
      json_filepath <- file.path(folder_name, paste0(log_name, ".json"))
      lg$add_appender(lgr::AppenderJson$new(json_filepath), name = paste0(log_name, "_json"))
    }
  }

  if (isTRUE(gpa$log_txt) && prefix == "setup_l1") {
    if(!paste0(log_name, "_txt") %in% names(lg$appenders)) {
      folder_name <- glue::glue(gpa$output_locations$feat_sub_directory)
      
      # Insert "logs" before "feat_l1" in the folder path
      folder_parts <- strsplit(folder_name, .Platform$file.sep)[[1]]
      feat_idx <- which(folder_parts == "feat_l1")
      folder_parts <- append(folder_parts, "logs", after = feat_idx - 1)
      folder_name <- do.call(file.path, as.list(folder_parts))
      
      # create folder if it doesn't exist    
      if (!dir.exists(folder_name)) {
        dir.create(folder_name, recursive = TRUE)
      }

      txt_filepath <- file.path(folder_name, paste0(log_name, ".txt"))
      lg$add_appender(lgr::AppenderFile$new(txt_filepath), name = paste0(log_name, "_txt"))
    }
  }

  return(lg)
}

#' Log a warning and emit a base R warning
#' @param lg Logger object
#' @param message Message text passed to glue::glue
#' @keywords internal
log_warn <- function(lg, message, ...) {
  msg <- sprintf(message, ...)
  lg$warn(msg)
  warning(msg, call. = FALSE)
}

#' Log an error and stop execution
#' @param lg Logger object
#' @param message Message text passed to glue::glue
#' @keywords internal
log_error <- function(lg, message, ...) {
  msg <- sprintf(message, ...)
  lg$error(msg)
  stop(msg, call. = FALSE)
}

