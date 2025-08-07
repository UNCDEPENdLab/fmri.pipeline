#' Internal helper function to submit a query to the tracking SQLite database
#'
#' @param str Character string SQLite query
#' @param sqlite_db Path to SQLite database used for tracking
#' @param param List of parameters/arguments to be used in query
#'
#' @keywords internal
submit_tracking_query = function(str, sqlite_db, param = NULL) {
  # previously called submit_sqlite()
  checkmate::assert_string(str)
  
  # check if tracking table exists in sqlite_db; if not, create it
  table_exists <- sqlite_table_exists(sqlite_db, "job_tracking")
  
  if (isFALSE(table_exists)) {
    create_tracking_db(sqlite_db)
  }
  
  # open sqlite connection and execute query
  submit_sqlite_query(str = str, sqlite_db = sqlite_db, param = param)
  
}

#' Internal helper function to reset tracking SQLite database
#'
#' @param sqlite_db Path to SQLite database used for tracking
#'
#' @keywords internal
reset_tracking_sqlite_db = function(sqlite_db) {
  # this file has the SQL syntax to setup (and reset) the database
  # reset_sql <- "
  # SET foreign_key_checks = 0;
  # DROP TABLE IF EXISTS job_tracking;
  # SET foreign_key_checks = 1;
  # "
  reset_sql <- "DELETE FROM job_tracking" # delete all records
  submit_tracking_query(str = reset_sql, sqlite_db = sqlite_db)
}


#' Internal helper function to create the tracking SQLite database
#'
#' @param sqlite_db Path to SQLite database used for tracking
#'
#' @keywords internal
create_tracking_db = function(sqlite_db) {
  # previously called create_sqlite_db()
  job_spec_sql <- "
    CREATE TABLE job_tracking (
      id INTEGER PRIMARY KEY,
      parent_id INTEGER,
      child_level INTEGER DEFAULT 0,
      job_id VARCHAR NOT NULL UNIQUE,
      job_name VARCHAR,
      sequence_id VARCHAR,
      batch_directory VARCHAR,
      batch_file VARCHAR,
      compute_file VARCHAR,
      code_file VARCHAR,
      n_nodes INTEGER CHECK (n_nodes >= 1),
      n_cpus INTEGER CHECK (n_cpus >= 1),
      wall_time VARCHAR,
      mem_per_cpu VARCHAR,
      mem_total VARCHAR,
      scheduler VARCHAR,
      scheduler_options VARCHAR,
      job_obj BLOB,
      time_submitted INTEGER,
      time_started INTEGER,
      time_ended INTEGER,
      status VARCHAR(24),
      FOREIGN KEY (parent_id) REFERENCES job_tracking (id)
    );
    "
  # open sqlite connection
  submit_sqlite_query(str = job_spec_sql, sqlite_db = sqlite_db)
}


#' Internal helper funciton to insert a job into the tracking SQLite database
#'
#' @param sqlite_db Path to SQLite database used for tracking
#' @param job_id Character id of job to insert
#' @param tracking_args List of named tracking arguments
#'
#' @keywords internal
insert_tracked_job = function(sqlite_db, job_id, tracking_args = list()) {
  # previously called sqlite_insert_job()
  if (is.null(sqlite_db) || is.null(job_id)) return(invisible(NULL)) # skip out if not using DB or if job_id is NULL
  if (is.numeric(job_id)) job_id <- as.character(job_id)
  if (is.null(tracking_args$status)) tracking_args$status <- "QUEUED" # default value of first status

  insert_job_sql <- "INSERT INTO job_tracking
    (job_id, job_name, sequence_id, batch_directory,
    batch_file, compute_file, code_file,
    n_nodes, n_cpus, wall_time,
    mem_per_cpu, mem_total,
    scheduler, scheduler_options, job_obj,
    time_submitted, status)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"

  # gather tracking parameters into a list
  param <- list(job_id, tracking_args$job_name, tracking_args$sequence_id, 
                tracking_args$batch_directory, tracking_args$batch_file, 
                tracking_args$compute_file, tracking_args$code_file, 
                tracking_args$n_nodes, tracking_args$n_cpus, tracking_args$wall_time, 
                tracking_args$mem_per_cpu, tracking_args$mem_total, tracking_args$scheduler,
                tracking_args$scheduler_options, tracking_args$job_obj, 
                as.character(Sys.time()), tracking_args$status)
  
  for (i in 1:length(param)) {
    param[[i]] <- ifelse(is.null(param[[i]]), NA, param[[i]]) # convert NULL values to NA for dbExecute
  }
  
  # order the tracking arguments to match the query; status is always 'QUEUED' when first added to the database
  submit_tracking_query(str = insert_job_sql, sqlite_db = sqlite_db, param = param)
}


#' Add parent/child id relationship to tracking database
#'
#' @param sqlite_db Path to SQLite database used for tracking
#' @param job_id Job id of job for which to add a parent
#' @param parent_job_id Job id of the parent job to job_id
#' @param child_level Level of child; currently supports two levels
#'
#' @importFrom DBI dbConnect dbExecute dbDisconnect
#' @export
add_tracked_job_parent = function(sqlite_db = NULL, job_id = NULL, parent_job_id = NULL, child_level = 1) {
  # skip out if not using DB or job_id/parent_id is NULL
  if (is.null(sqlite_db) || is.null(job_id) || is.null(parent_job_id)) return(invisible(NULL))
  if (isFALSE(child_level %in% 1:2)) child_level <- NA # if level not 1 or 2, return NA
  if (is.numeric(job_id)) job_id <- as.character(job_id)
  if (is.numeric(parent_job_id)) parent_job_id <- as.character(parent_job_id)

  # retrieve sequence id from parent
  sequence_id_sql <- "SELECT sequence_id FROM job_tracking WHERE job_id = ?"
  sequence_id <- tryCatch({
    # open sqlite connection and execute query
    id <- submit_sqlite_query(str = sequence_id_sql, sqlite_db = sqlite_db, 
                        param = list(parent_job_id), return_result = TRUE)
    if(!is.null(id[1,1])) { id[1,1] } else { NA }
  }, error = function(e) { print(e); NA})
  
  
  add_parent_sql <- "UPDATE job_tracking
    SET parent_id = (SELECT id FROM job_tracking WHERE job_id = ?), 
      sequence_id = ?,
      child_level = ?
    WHERE job_id = ?"
  tryCatch({
    # open sqlite connection and execute query
    submit_sqlite_query(str = add_parent_sql, sqlite_db = sqlite_db, 
                        param = list(parent_job_id, sequence_id, child_level, job_id))
  }, error = function(e) { print(e); return(NULL)})
  
}


#' Update job status in tracking SQLite database
#' 
#' @param sqlite_db Character string specifying the SQLite database used for job tracking
#' @param job_id Character string specifying the job id to update as failed
#' @param status Character string specifying the job status to set. Must be one of: 
#'   "QUEUED", "STARTED", "FAILED", "COMPLETED", "FAILED_BY_EXT"
#' @param exclude Any job ids to ignore when cascading a status
#' @importFrom glue glue
#' @importFrom DBI dbConnect dbExecute dbDisconnect
#' @export
update_tracked_job_status <- function(sqlite_db = NULL, job_id = NULL, status, cascade = FALSE, exclude = NULL) {
  
  if (!checkmate::test_string(sqlite_db)) return(invisible(NULL))
  if (is.numeric(job_id)) job_id <- as.character(job_id)
  if (!checkmate::test_string(job_id)) return(invisible(NULL)) # quiet failure on invalid job id

  checkmate::assert_string(status)
  status <- toupper(status)
  checkmate::assert_subset(status, c("QUEUED", "STARTED", "FAILED", "COMPLETED", "FAILED_BY_EXT"))
  if (cascade & status %in% c("QUEUED", "STARTED", "COMPLETED")) {
    cascade <- FALSE
    warning("Only status FAILED or FAILED_BY_EXT can cascade in `update_tracked_job_status`")
  }
  
  now <- as.character(Sys.time())
  time_field <- switch(status,
                       QUEUED = "time_submitted",
                       STARTED = "time_started",
                       FAILED = "time_ended",
                       COMPLETED = "time_ended",
                       FAILED_BY_EXT = "time_ended"
  )
  
  tryCatch({
    submit_sqlite_query(str = glue("UPDATE job_tracking SET STATUS = ?, {time_field} = ? WHERE job_id = ?"),
                        sqlite_db = sqlite_db, param = list(status, now, job_id))
  }, error = function(e) { print(e); return(NULL)})
  
  # recursive function for "cascading" failures using status "FAILED_BY_EXT"
  if (cascade) {
    if (is.numeric(exclude)) exclude <- as.character(exclude)
    
    status_tree <- fmri.pipeline::get_tracked_job_status(job_id, return_children = TRUE, sqlite_db = sqlite_db) # retreive current job and children
    job_ids <- status_tree$job_id # get list of job ids
    
    for (child_job in setdiff(job_ids, c(job_id, exclude))) {
      child_status <- with(status_tree, status[which(job_ids == child_job)]) # check status
      if (child_status != "FAILED") fmri.pipeline::update_tracked_job_status(child_job, sqlite_db = sqlite_db, status = "FAILED_BY_EXT", cascade = TRUE)
    }
  }
  
  return(invisible(NULL))
  
}

#' Query job status in tracking SQLite database
#' 
#' @param job_id The job id for which to retreive the status
#' @param sqlite_db Character string of sqlite database
#' @param return_children Return child jobs of this job
#' @param return_parent Return parent jobs of this job
#' @param sequence_id Alternatively return all jobs with sequence ID; cannot specify both job ID and sequence ID
#' 
#' @return An R data.frame version of the tracking database
#' @importFrom DBI dbConnect
#' @importFrom RSQLite SQLite
#' @export
get_tracked_job_status <- function(job_id = NULL, return_children = FALSE, return_parent = FALSE, 
                                   sequence_id = NULL, sqlite_db) {
  
  on.exit(try(dbDisconnect(con)))
  if (!checkmate::test_file_exists(sqlite_db)) {
    warning("Cannot find SQLite database at: ", sqlite_db)
  }
  
  stopifnot("Must specify `job_id` or `sequence_id`, not both" =
              is.null(job_id) || is.null(sequence_id))
  
  if (!is.null(sequence_id)) {
    
    if (is.numeric(sequence_id)) sequence_id <- as.character(sequence_id)
    if (!checkmate::test_string(sequence_id)) return(invisible(NULL))
    
    str <- paste0("SELECT * FROM job_tracking WHERE sequence_id = ?")
    param <- list(sequence_id)
    
  } else {

    if (is.numeric(job_id)) job_id <- as.character(job_id)
    if (!checkmate::test_string(job_id)) return(invisible(NULL))
    checkmate::assert_logical(return_children)
    checkmate::assert_logical(return_parent)
    
    str <- paste0("SELECT * FROM job_tracking WHERE job_id = ?", 
                  ifelse(return_children, " OR parent_id = (SELECT id FROM job_tracking WHERE job_id = ?)", ""),
                  ifelse(return_parent, " OR id = (SELECT parent_id FROM job_tracking WHERE job_id = ?)", ""))
  
    param <- as.list(rep(job_id, 1 + return_children + return_parent))
  }
  
  con <- dbConnect(RSQLite::SQLite(), sqlite_db)
  df <- dbGetQuery(con, str, param = param)
  
  # rehydrate job_obj back into R6 class
  if (nrow(df) > 0L) df$job_obj <- lapply(df$job_obj, function(x) if (!is.null(x)) unserialize(x))
  return(df)
}


#' Internal helper function to update tracker_args object
#'
#' @param list_to_populate The list whose argument will be populated/updated
#' @param arg_name The named list element to update
#' @param new_value The new value to update the element with
#' @param append If TRUE, appends the new value to the current value using the paste function. Default: FALSE
#'
#' @keywords internal
populate_list_arg = function(list_to_populate, arg_name, new_value = NULL, append = FALSE) {
  
  checkmate::assert_list(list_to_populate)
  if (is.null(new_value)) return(list_to_populate) # if the new_value arg is NULL just return the list as is

  if(is.null(list_to_populate[[arg_name]]) || is.na(list_to_populate[[arg_name]])) {
    list_to_populate[[arg_name]] <- new_value # if the current arg is NULL, update with new value
  } else if (append) {
    # if it's not NULL but append is TRUE, appends new value to beginning of old value
    list_to_populate[[arg_name]] <- paste(list_to_populate[[arg_name]], new_value, sep = "\n")  
  }
  return(list_to_populate)
}
