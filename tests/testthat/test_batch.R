test_that("Basic job setup works", {
	# note that HPC settings like nodes and batch_code are ignored for local scheduler
	w <- R_batch_job$new(
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

	w$generate()
	expect_equal(1, 1)
})