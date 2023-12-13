#' Function which queries for the status of jobs associated with the given gpa object
#' 
#' @param gpa you glm_pipeline_arguments list
#' @param batch_id filter by batch_id from run_glm_pipeline
#' @param days_back how far back in time to search for jobs (in days)
#' @param desc whether to sort the status by descending order (newest first)
#' @param update whether or not to update the status of unresolved jobs
#' 
#' @return a data.frame containing the status of jobs
#' @export
show_job_status <- function(gpa=NULL, batch_id=NULL, days_back=1, desc=TRUE, update=TRUE, debug=FALSE) {

    # Use read database function to get the status of the jobs
    #status_df <- read_df_sqlite(gpa = gpa, table = "job_status", id = batch_id)
    conn <- get_sqlite_conn(gpa$output_locations$sqlite_db)
    status_df <- dbReadTable(conn, "job_status")

    if(update) {
        # Query for the jobs that still have the state "RUNNING" or "SUBMITTED" 
        # and get the output as a list. These jobs do not have a settled state so must be
        # requeried.
        open_jobs <- status_df %>%
            filter(state == "RUNNING" | state == "SUBMITTED") %>%
            pull(id)

        if(debug) {
            cat("Open jobs:\n")
            print(open_jobs)
        }

        # if open_jobs isn't empty, get the lastest state of those jobs and update the DB
        if(length(open_jobs) > 0) {
            # Get the latest state of any still-running jobs
            current_status_df <- get_job_status(job_ids = open_jobs, scheduler = gpa$scheduler)

            if(debug) {
                cat("Current status:\n")
                print(current_status_df)
            }

            # Select only the JobID and State columns from the current_status_df
            current_status_df <- current_status_df %>%
                dplyr::select(JobID, State) %>%
                dplyr::rename(id = JobID, state = State)
            if(debug) {
                cat("Current status transform:\n")
                print(current_status_df)
            }

            # Select only open_jobs from the status_df
            status_df_open <- status_df %>%
                filter(id %in% open_jobs)

            # Update state column of open_jobs to the state queried in current_status_df
            status_df_open <- status_df_open %>%
                dplyr::select(-state) %>%
                left_join(current_status_df, by = "id")
                #mutate(status = ifelse(is.na(status.y), status.x, status.y)) %>%
                #dplyr::rename(state = state.y)

            # Get the current time
            current_time <- Sys.time()
            # Convert current time to ISO
            current_time <- format(current_time, "%Y-%m-%dT%H:%M:%S")

            # Update the timestamp column to the current time
            status_df_open <- status_df_open %>%
                mutate(last_updated = current_time)

            if(debug) {
                cat("Status df open:\n")
                print(status_df_open)
            }

            # Iterate over each row in status_df_open and update the status table
            for(i in 1:nrow(status_df_open)) {
                # Get the row
                row <- status_df_open[i,]

                # Update the status table
                insert_df_sqlite(
                    gpa = gpa,
                    id = row$id,
                    session = row$session,
                    table = "job_status",
                    data = row,
                    append = FALSE,
                    overwrite = TRUE
                )
            }


            # Re-pull the status_df with the new records
            # Could also avoid a re-pull by just updating the status_df with the new records
            status_df <- dbReadTable(conn, "job_status")
        }
    }
    
    # Sort the status_df
    if(desc) { 
        status_df <- status_df %>%
            arrange(desc(submitted))
    } else {
        status_df <- status_df %>%
            arrange(submitted)
    }
    
    # Show the status_df
    cat("Job status:\n")
    status_df <- status_df %>% dplyr::rename(slurm_id=id, batch_id=session)
    print(status_df)

    dbDisconnect(conn)
}

# Query the status df for the latest value of every job id, within the time period of 'range' days
# from the current time
# status_df <- status_df %>%
#     filter(timestamp > since_time)
# Get the current time
    # current_time <- Sys.time()

    # # Get the time range to search for jobs
    # range_seconds <- days_back*24*60*60
    # since_time <- current_time - range_seconds