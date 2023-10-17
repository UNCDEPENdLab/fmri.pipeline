#' Function which queries for the status of jobs associated with the given gpa object
#' 
#' @param gpa you glm_pipeline_arguments list
#' @param range how far back in time to search for jobs (in days)
#' @param desc whether to sort the status by descending order (newest first)
#' 
#' @return a data.frame containing the status of jobs
show_job_status <- function(gpa=NULL, range=1, desc=TRUE, latest_only=TRUE) {
    # Get the path to the sqlite database from the gpa object

    # Use read database function to get the status of the jobs
    status_df <- read_df_sqlite(gpa = gpa, table = "job_status")

    # Get the current time
    current_time <- Sys.time()

    # Get the time range to search for jobs
    since_time <- current_time - range*24*60*60

    # Query the status df for the latest value of every job id, within the time period of 'range' days from the current time
    status_df <- status_df %>%
        # Only pick values newer than since_time
        filter(submission_time > since_time) %>%
        # Group by job and pick the latest submission time
        group_by(job_id) %>%
        filter(submission_time == max(submission_time)) %>%
        ungroup()

    # Sort the df
    if(desc) {
        status_df <- status_df %>%
            arrange(desc(submission_time))
    } else {
        status_df <- status_df %>%
            arrange(submission_time)
    }

    # Query for the jobs that still have the status "RUNNING" and get the output as a list
    running_jobs <- status_df %>%
        filter(status == "RUNNING") %>%
        pull(job_id)

    # if running_jobs isn't empty, get the lastest status of those jobs
    if(length(running_jobs) > 0) {
         # Get the latest status of any still-running jobs
        current_status_df <- get_job_status(job_ids = running_jobs, scheduler = gpa$scheduler)
    } else {
        current_status_df <- data.frame()
    }

    # Print the status of the jobs
    
}