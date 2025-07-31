#' test build_l2_models function
test_that("test build_l2_models", {

  analysis_name <- "gpa_tests_build_l1_models"
  test_data_base_dir <- "/proj/mnhallqlab/projects/fmri.pipeline_test_data"
  l1_spec_file <- "sample_2_L1_spec.yaml"
  trial_data_file <- "sample_trial_data.csv.gz"
  run_data_file <- "sample_run_data.csv"
  subject_data_file <- "sample_subject_data.csv"
  cache_file <- "gpa_tests_base_build_l1_models.rds"

  # run set_glm_pipeline()
  gpa <- get_gpa_base(analysis_name = analysis_name, 
    test_data_base_dir = test_data_base_dir, 
    trial_data_file = trial_data_file,
    run_data_file = run_data_file,
    subject_data_file = subject_data_file,
    cache_file = cache_file,
    scheduler = "slurm", drop_volumes = 2,
    exclude_run = "max(FD) > 5 | sum(FD > .9)/length(FD) > .10",
    exclude_subject = "n_good_runs < 4",
    truncate_run = "(FD > 0.9 & time > last_offset) | (time > last_offset + last_isi)",
    spike_volumes = NULL
  )

  # # setup compute enviroment
  gpa <- setup_compute_environment(gpa)

  # Build L1 models
  gpa <- build_l1_models(gpa, from_spec_file = file.path(test_data_base_dir, l1_spec_file))

  # Build L2 models
  gpa <- build_l2_models(gpa)

  # check if the l2_models object is entirely created
  expect_equal(length(unique(gpa$l2_models$models$emotion$metadata$id)), 10)
  expect_equal(length(unique(gpa$l2_models$models$emotion$metadata$run_number)), 8)
  expect_equal(nrow(gpa$l2_models$models$emotion$model_data), 80)

})