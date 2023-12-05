test_that("Request for job status works without querying to update status", {
	tmp_dir <- get_temp_dir()
	db_path <- file.path(tmp_dir, "testdb")

	parent_job <- get_simple_batch_job(job_name="parent", tmp_dir=tmp_dir)
	parent_job$batch_id <- "test_batch_id_123"
	parent_job$sqlite_db <- db_path

	child_job <- get_simple_batch_job(job_name="child1", tmp_dir=tmp_dir)
	child_job$parent_jobs <- parent_job
	child_job$batch_id <- "test_batch_id_123"
	child_job$sqlite_db <- db_path

	child_job2 <- get_simple_batch_job(job_name="child2", tmp_dir=tmp_dir)
	child_job2$parent_jobs <- parent_job
	child_job2$batch_id <- "test_batch_id_123"
	child_job2$sqlite_db <- db_path

	jobs <- R_batch_sequence$new(parent_job, child_job, child_job2)
	jobs$submit()
	# Create a gpa list filler object that just has output_locations$sqlite_db
    # populated in order to use insert_df_sqlite
    gpa <- list(output_locations = list(sqlite_db = db_path), id="gpa_id_123")
	class(gpa) <- c("list", "glm_pipeline_arguments")

	show_job_status(gpa=gpa, batch_id="test_batch_id_123", days_back=1, desc=TRUE, latest_only=TRUE, update=FALSE)

	expect_equal(1, 1)
})