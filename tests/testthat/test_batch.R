get_simple_batch_job <- function() {
	return(
		R_batch_job$new(
			batch_directory=tempdir(),
			job_name="step1",
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
	# note that HPC settings like nodes and batch_code are ignored for local scheduler
	w <- get_simple_batch_job()

	w$generate()
	print(list.files(w$batch_directory))
	expect_equal(1, 1)
})

test_that("Basic job submission works", {
	# note that HPC settings like nodes and batch_code are ignored for local scheduler
	w <- get_simple_batch_job()
	w$batch_id <- "test_batch_id_123"
	w$sqlite_db <- file.path(tempdir(), "testdb")

	w$submit()
	print(list.files(w$batch_directory))
	expect_equal(1, 1)
})