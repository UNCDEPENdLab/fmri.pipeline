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
#' @export
show_job_status <- function(gpa=NULL, batch_id=NULL, days_back=1, desc=TRUE, latest_only=TRUE, update=TRUE) {

    # Use read database function to get the status of the jobs
    #status_df <- read_df_sqlite(gpa = gpa, table = "job_status", id = batch_id)
    conn <- get_sqlite_conn(gpa$output_locations$sqlite_db)
    status_df <- dbReadTable(conn, "job_status")

    if(update) {
        # Query for the jobs that still have the state "RUNNING" or "SUBMITTED" 
        # and get the output as a list. These jobs do not have a settled state so must be
        # requeried.
        open_jobs <- status_df %>%
            filter(state == "RUNNING" || state == "SUBMITTED") %>%
            pull(id)
        cat("Open jobs:\n")
        print(open_jobs)

        # if open_jobs isn't empty, get the lastest state of those jobs and update the DB
        if(length(open_jobs) > 0) {
            # Get the latest state of any still-running jobs
            current_status_df <- get_job_status(job_ids = open_jobs, scheduler = gpa$scheduler)

            cat("Current status:\n")
            print(current_status_df)

            # Select only the JobID and State columns from the current_status_df
            current_status_df <- current_status_df %>%
                dplyr::select(JobID, State) %>%
                dplyr::rename(id = JobID, state = State)

            cat("Current status transform:\n")
            print(current_status_df)

            # Select only open_jobs from the status_df
            status_df_open <- status_df %>%
                filter(id %in% open_jobs)

            # Update state column of open_jobs to the state queried in current_status_df
            status_df_open <- status_df_open %>%
                dplyr::select(-state) %>%
                left_join(current_status_df, by = "id")
                #mutate(status = ifelse(is.na(status.y), status.x, status.y)) %>%
                #dplyr::rename(state = state.y)

            cat("Status df open:\n")
            print(status_df_open)

            # For each id in status_df_open, add record to the status table
            for(record in status_df_open) {
                # Add records with updated state to status table
                insert_df_sqlite(
                    gpa = gpa,
                    id = status_df_open$id,
                    table = "job_status",
                    data = status_df_open,
                    append = TRUE,
                    delete_extant = FALSE
                )
            }

            # Re-pull the status_df with the new records
            # Could also avoid a re-pull by just updating the status_df with the new records
            status_df <- dbReadTable(conn, "job_status")
        }
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
            filter(timestamp > since_time) %>%
        
        if(latest_only) {
            # Group by job and pick the latest submission time
            status_df <- state_df %>%
                group_by(job_id) %>%
                filter(timestamp == max(timestamp)) %>%
                ungroup()
        }
    }
    
    # Sort the status_df
    if(desc) {
        status_df <- status_df %>%
            arrange(desc(timestamp))
    } else {
        status_df <- status_df %>%
            arrange(timestamp)
    }
    
    # Show the status_df
    cat("Job status:\n")
    print(status_df)

    dbDisconnect(conn)
}