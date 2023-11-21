get_simple_batch_job <- function(job_name="step1", tmp_dir=get_temp_dir()) {
	return(
		R_batch_job$new(
			batch_directory=tmp_dir,
			job_name=job_name,
			n_nodes=1,
			n_cpus=1,
			wall_time="10:00",
			r_code=c(
					"Sys.sleep(100)",
					"print('hi')",
					"x <- 2+2"
			),
			r_packages=c("lme4"),
			batch_code = c(
					"module use /proj/mnhallqlab/sw/modules",
					"module load r/4.0.3_depend"
			),
			scheduler="slurm",
			repolling_interval=1 #repoll every second for this test (should be slower for real jobs)
		)
	)
}

test_that("Basic job setup works", {
	w <- get_simple_batch_job()

	w$generate()
	print(list.files(w$batch_directory))
	expect_equal(1, 1)
})

test_that("Basic job submission works", {
	tmp_dir <- get_temp_dir()
	
	w <- get_simple_batch_job(tmp_dir=tmp_dir)
	w$batch_id <- "test_batch_id_123"
	w$sqlite_db <- file.path(tmp_dir, "testdb")

	w$submit()
	print(list.files(w$batch_directory))

	preview_db(w$sqlite_db)

	expect_equal(1, 1)
})

test_that("Nested job submission works", {
	tmp_dir <- get_temp_dir()
	db_path <- file.path(tmp_dir, "testdb")

	parent_job <- get_simple_batch_job(job_name="parent", tmp_dir=tmp_dir)
	parent_job$batch_id <- "test_batch_id_123"
	parent_job$sqlite_db <- db_path

	child_job <- get_simple_batch_job(job_name="child", tmp_dir=tmp_dir)
	child_job$parent_jobs <- parent_job
	child_job$batch_id <- "test_batch_id_124"
	child_job$sqlite_db <- db_path

	jobs <- R_batch_sequence$new(parent_job, child_job)
	jobs$submit()

	preview_db(child_job$sqlite_db)

	expect_equal(1, 1)
})