

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
  checkmate::assert_set_equal(length(job_sequence), 5)
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

#' helper to retrieve information about a specific sequence from the config file
#' 
#' @param config_file the path to the configuration file
#' @param sequence_id unique identifier of sequence
#' 
#' @keywords internal
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
