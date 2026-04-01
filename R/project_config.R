

#' helper function to create/update project configuration json file in gpa output directory
#'
#' @param gpa a glm_pipeline_arguments object
#' @param job_sequence a list with named logicals indicating whether 
#'                      that part of the pipeline was submitted or not
#'                      (currently finalize, l1, l2, l3, cleanup)
#' @param sequence_id the unique identifer used for this batch sequence submission 
#' @param batch_directory a path to the batch directory
#'
#' @keywords internal
update_project_config <- function(gpa, job_sequence, sequence_id, batch_directory = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_list(job_sequence)
  required_names <- c("finalize", "l1", "l2", "l3", "cleanup")
  checkmate::assert_subset(required_names, names(job_sequence))
  sequence_id <- as.character(sequence_id)
  checkmate::assert_string(sequence_id)
  
  if(is.null(batch_directory)) {
    batch_directory <- file.path(gpa$output_locations$scheduler_scripts, paste0("batch_", sequence_id))
  }
  #checkmate::assert_directory_exists(batch_directory)
  
  # make sure config file has been specified in gpa object
  if(is.null(gpa$output_locations$project_config_json)) return(invisible(NULL))
  config_file <- gpa$output_locations$project_config_json
  
  if(!checkmate::test_file_exists(config_file)) {
    # initialize config file if it doesn't exist yet
    config <- list(project_details = list(
      analysis_name = as.character(gpa$analysis_name),
      output_directory = as.character(gpa$output_directory)
    ), batch_history = list())
    
    
    tryCatch({
      write(jsonlite::toJSON(config, pretty = T), config_file)
    },
    error = function(e) {
      print(e)
    })
    
  }
  
  # read in config file
  config <- jsonlite::read_json(config_file, simplifyVector = T)
  
  # update config file with new batch submission
  config$batch_history[[sequence_id]] <- list(
    time_submitted = Sys.time(),
    user = gpa$sys_info["user"],
    job_sequence = job_sequence,
    batch_directory = batch_directory
  )
  
  # rewrite config file
  tryCatch({
    write(jsonlite::toJSON(config, pretty = T), config_file)
  },
  error = function(e) {
    print(e)
  })
  
}


#' helper to summarize the contents of a project configuration file
#'
#' @param config_file the path to the configuration file
#'
#' @keywords internal
summarize_project_config <- function(config_file = NULL) {

  if (!checkmate::test_file_exists(config_file)) {
    # exit code in the case the config file doesn't exist
    return(invisible("file_dne"))
  }
  
  config <- jsonlite::read_json(config_file, simplifyVector = T)
  
  n_batch <- length(config$batch_history)
  seq_ids <- names(config$batch_history)
  
  if (n_batch == 0) {
    # exit code in the case no submissions have been made yet
    return(invisible("no_history"))
  }
  
  # order config history by time of submission
  submission_times <- sapply(config$batch_history, function(x) return(as.POSIXct(x$time_submitted)))
  ordered_seq_ids <- seq_ids[order(submission_times)]
  ordered_batch_history <- config$batch_history[ordered_seq_ids]
  
  # analysis summary
  cli::cli_h1(
    "Analysis: {config$project_details$analysis_name}"
  )
  cli::cli_ul(id = "analysis_summary")
  cli::cli_li(
    c(
      "Output directory: {config$project_details$output_directory}",
      "Number of runs: {n_batch}",
      "History:"
    )
  )
  
  # run summary
  cli::cli_ol(id = "run_summary")
  
  for (i in 1:n_batch) {
    time <- ordered_batch_history[[i]]$time_submitted
    user <- ordered_batch_history[[i]]$user
    job_seq <- ordered_batch_history[[i]]$job_sequence
    
    # print which jobs were run in the sequence
    seq_str <- paste(
      ifelse(job_seq$finalize, "Finalize Config.", ""),
      ifelse(job_seq$l1, "L1 Mods.", ""),
      ifelse(job_seq$l2, "L2 Mods.", ""),
      ifelse(job_seq$finalize, "L3 Mods.", ""),
      ifelse(job_seq$cleanup, "Clean Up", ""),
      sep = " \u2192 "
    )
    
    if (i == n_batch) {
      tag <- "{.emph (Latest)}"
    } else {
      tag <- ""
    }
    
    # print time of submission & user
    cli::cli_li(
      cli::col_cyan(paste("Run on {format(as.POSIXct(time), '%D')} at {format(as.POSIXct(time), '%I:%M:%S %p')} by {user}", tag)),
    )
    # print submission sequence
    cli::cli_bullets(
      c(" " = "{.emph {seq_str}}")
    )
  }
  
  # end bulleting/numbering
  cli::cli_end(id = "run_summary")
  cli::cli_end(id = "analysis_summary")
  
  # return sequence ids in order of job submission times
  return(invisible(ordered_seq_ids))
  
}


#' helper to retrieve information about a specific sequence from the config file
#' 
#' @param config_file the path to the configuration file
#' @param sequence_id unique identifier of sequence
#' 
#' @keywords internal
#' @noRd
get_sequence_info <- function(config_file = NULL, sequence_id = NULL) {
  if (is.null(config_file) || is.null(sequence_id)) return(invisible(NULL))
  sequence_id <- as.character(sequence_id)
  checkmate::assert_string(sequence_id)
  
  # read config file
  config <- jsonlite::read_json(config_file, simplifyVector = T)
  
  # extract sequence info
  seq_info <- config$batch_history[[sequence_id]]
  
  return(seq_info)
  
}
