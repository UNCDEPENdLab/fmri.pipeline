#' test build_l1_models function
test_that("test build_l1_models", {
  
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

  # MNH mentioned that everytime the gpa should be created afresh instead of loading previously saved
  # gpa <- readRDS(file.path(test_data_base_dir, cache_file))

  # Build L1 models
  gpa <- build_l1_models(gpa, from_spec_file = file.path(test_data_base_dir, l1_spec_file))

  # check if the l1_models object and signals is created
  expect_true(!is.null(gpa$l1_models))
  expect_equal(sum(names(gpa$l1_models$signals) == c("clock", "feedback", "entropy_clock", "pe_feedback")), 4)

  # check if the l1_models object has the correct number of elements
  expect_equal(length(unique(gpa$l1_models$signals$clock$value$id)), 10)
  expect_equal(unique(gpa$l1_models$signals$clock$value$session), 1)
  expect_equal(length(unique(gpa$l1_models$signals$clock$value$run_number)), 8)
  expect_equal(length(unique(gpa$l1_models$signals$clock$value$trial)), 50)
  expect_equal(length(unique(gpa$l1_models$signals$clock$value$value)), 1)

  expect_equal(length(unique(gpa$l1_models$signals$feedback$value$id)), 10)
  expect_equal(unique(gpa$l1_models$signals$feedback$value$session), 1)
  expect_equal(length(unique(gpa$l1_models$signals$feedback$value$run_number)), 8)
  expect_equal(length(unique(gpa$l1_models$signals$feedback$value$trial)), 50)
  expect_equal(length(unique(gpa$l1_models$signals$feedback$value$value)), 1)

  expect_equal(length(unique(gpa$l1_models$signals$entropy_clock$value$id)), 10)
  expect_equal(unique(gpa$l1_models$signals$entropy_clock$value$session), 1)
  expect_equal(length(unique(gpa$l1_models$signals$entropy_clock$value$run_number)), 8)
  expect_equal(length(unique(gpa$l1_models$signals$entropy_clock$value$trial)), 50)

  expect_equal(length(unique(gpa$l1_models$signals$pe_feedback$value$id)), 10)
  expect_equal(unique(gpa$l1_models$signals$pe_feedback$value$session), 1)
  expect_equal(length(unique(gpa$l1_models$signals$pe_feedback$value$run_number)), 8)
  expect_equal(length(unique(gpa$l1_models$signals$pe_feedback$value$trial)), 50)
  
  # check if the l1_models object has the correct values
  expect_equal(gpa$l1_models$signals$clock$value$value[1], 1)
  expect_equal(gpa$l1_models$signals$feedback$value$value[1], 1)
  
})

test_that("test build_l1_models with non-numeric parametric modulator", {
  
  analysis_name <- "gpa_tests_build_l1_models_wrong"
  test_data_base_dir <- "/proj/mnhallqlab/projects/fmri.pipeline_test_data"
  l1_spec_file <- "sample_2_L1_spec_wrong_para_mod.yaml"
  trial_data_file <- "sample_trial_data.csv.gz"
  run_data_file <- "sample_run_data.csv"
  subject_data_file <- "sample_subject_data.csv"
  cache_file <- "gpa_tests_base_build_l1_models_wrong.rds"

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

  # setup compute enviroment
  gpa <- setup_compute_environment(gpa)

  # Build L1 models
  gpa <- build_l1_models(gpa, from_spec_file = file.path(test_data_base_dir, l1_spec_file))

  # check if the l1_models object and signals is created
  expect_true(!is.null(gpa$l1_models))
  expect_equal(sum(names(gpa$l1_models$signals) == c("clock", "feedback", "entropy_clock", "pe_feedback")), 4)

  # check if the l1 model has factor parametric modulator values
  unique(gpa$l1_models$signals$rewFunc_feedback$value$value)

})

