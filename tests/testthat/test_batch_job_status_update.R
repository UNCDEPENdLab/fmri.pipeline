test_that("Request for job status works, with a query to update status", {
	tmp_dir <- get_temp_dir()
	db_path <- file.path(tmp_dir, "testdb")

	parent_job <- get_simple_batch_job(job_name="parent", tmp_dir=tmp_dir)
	parent_job$batch_id <- 123
	parent_job$sqlite_db <- db_path

	child_job <- get_simple_batch_job(job_name="child1", tmp_dir=tmp_dir)
	child_job$parent_jobs <- parent_job
	child_job$batch_id <- 123
	child_job$sqlite_db <- db_path

	child_job2 <- get_simple_batch_job(job_name="child2", tmp_dir=tmp_dir)
	child_job2$parent_jobs <- parent_job
	child_job2$batch_id <- 123
	child_job2$sqlite_db <- db_path

	jobs <- R_batch_sequence$new(parent_job, child_job, child_job2)
	jobs$submit()
	# Create a gpa list filler object that just has output_locations$sqlite_db
    # populated in order to use insert_df_sqlite
    gpa <- list(output_locations = list(sqlite_db = db_path), id="gpa_id_123", scheduler="slurm")
	class(gpa) <- c("list", "glm_pipeline_arguments")

	# Sleep for 10 seconds to allow Slurm to process the jobs
	cat("Sleeping for 10 seconds to allow Slurm to process the jobs\n")
	Sys.sleep(10)

	show_job_status(gpa=gpa, batch_id="test_batch_id_123", days_back=1, desc=TRUE, update=TRUE)

	cat("Sleeping for another 20 seconds to allow the jobs to complete.\n")
	Sys.sleep(20)

	show_job_status(gpa=gpa, batch_id=123, days_back=1, desc=TRUE, update=TRUE)

	expect_equal(1, 1)
})

# Next steps:
# Change timestamp to time submitted
# Add last updated timestamp
# Add parent job column
# Generate batch id in run_glm_pipeline
# Retrn this id from run_glm_pipeline so it can be used, or attach to GPA