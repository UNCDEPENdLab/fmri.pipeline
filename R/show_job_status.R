#' Function which queries for the status of jobs associated with the given gpa object
#' 
#' @param gpa you glm_pipeline_arguments list
#' @param batch_id filter by batch_id from run_glm_pipeline
#' @param days_back how far back in time to search for jobs (in days)
#' @param desc whether to sort the status by descending order (newest first)
#' @param latest_only whether to only show the latest status of each job or all records
#' @param update whether or not to update the status of unresolved jobs
#' 
#' @return a data.frame containing the status of jobs
show_job_status <- function(gpa=NULL, batch_id=NULL, days_back=1, desc=TRUE, latest_only=TRUE, update=TRUE) {
    
    # Use read database function to get the status of the jobs
    status_df <- read_df_sqlite(gpa = gpa, table = "job_status")

    if(update) {
        # Query for the jobs that still have the status "RUNNING" or "SUBMITTED" 
        # and get the output as a list. These jobs do not have a settled state so must be
        # requeried.
        open_jobs <- status_df %>%
            filter(status == "RUNNING" || status == "SUBMITTED") %>%
            pull(job_id)

        # if open_jobs isn't empty, get the lastest status of those jobs and update the DB
        if(length(open_jobs) > 0) {
            # Get the latest status of any still-running jobs
            current_status_df <- get_job_status(job_ids = open_jobs, scheduler = gpa$scheduler)

            # Select only the JobID and State columns from the current_status_df
            current_status_df <- current_status_df %>%
                dplyr::select(JobID, State) %>%
                dplyr::rename(job_id = JobID, status = State)

            # Update state column of open_jobs to the state queried in current_status_df
            open_jobs <- open_jobs %>%
                left_join(current_status_df, by = "job_id") %>%
                #mutate(status = ifelse(is.na(status.y), status.x, status.y)) %>%
                #dplyr::select(-status.x, -status.y) %>%
                #dplyr::rename(status = status)

            # Add records with updated state to status table
            insert_df_sqlite(
                gpa = gpa,
                table = "job_status",
                data = open_jobs,
                append = TRUE
            )

            # Re-pull the status_df with the new records
            # Could also avoid a re-pull by just updating the status_df with the new records
            status_df <- read_df_sqlite(gpa = gpa, table = "job_status")
    }
    
    # Narrow down status id by batch_id if it's present, otherwise use a time range
    if(!is.null(batch_id)) {
        # Only select records from status_df with batch ids matching the batch_id argument
        status_df <- status_df %>%
            filter(batch_id == batch_id)
    } else {
        # Select records by time range

        # Get the current time
        current_time <- Sys.time()

        # Get the time range to search for jobs
        range_seconds <- days_back*24*60*60
        since_time <- current_time - range_seconds

        # Query the status df for the latest value of every job id, within the time period of 'range' days
        # from the current time
        status_df <- status_df %>%
            filter(submission_time > since_time) %>%
        
        if(latest_only) {
            # Group by job and pick the latest submission time
            status_df <- state_df %>%
                group_by(job_id) %>%
                filter(submission_time == max(submission_time)) %>%
                ungroup()
        }
    }
    
    # Sort the status_df
    if(desc) {
        status_df <- status_df %>%
            arrange(desc(submission_time))
    } else {
        status_df <- status_df %>%
            arrange(submission_time)
    }
    
    # Show the status_df
    status_df
}