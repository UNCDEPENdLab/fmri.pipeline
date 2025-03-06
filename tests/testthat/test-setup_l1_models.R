#' test finalize_pipeline_configuration function
test_that("test finalize_pipeline_configuration", {
  
  analysis_name <- "gpa_tests_build_l1_models"
  test_data_base_dir <- "/proj/mnhallqlab/projects/fmri.pipeline_test_data"
  l1_spec_file <- "sample_2_L1_spec.yaml"
  trial_data_file <- "sample_trial_data.csv.gz"
  run_data_file <- "sample_run_data.csv"
  subject_data_file <- "sample_subject_data.csv"
  gpa_cache_file <- "gpa_tests_build_l1_models.rds"
  cache_file <- "gpa_tests_base_build_l1_models.rds"

  # run building l1, l2, l3 models
  gpa <- create_gpa(
      analysis_name = analysis_name,
      test_data_base_dir = test_data_base_dir,
      l1_spec_file = l1_spec_file,
      trial_data_file = trial_data_file,
      run_data_file = run_data_file,
      subject_data_file = subject_data_file,
      gpa_cache_file = gpa_cache_file,
      cache_file = cache_file)
 
  # finalize pipeline configuration
  gpa <- finalize_pipeline_configuration(gpa)

  # create all FSF files for level one runs
  gpa <- setup_l1_models(gpa)
  
  # gpa$l1_model_setup is the output created by setup_l1_models()
  nsubjects_clock_testing <- 10
  nruns_clock_testing <- 8

  # test if the class of gpa$l1_model_setup is correct
  expect_equal(class(gpa$l1_model_setup)[[2]], "l1_setup")
  
  # test if the output contains the expected elements
  expect_true("fsl" %in% names(gpa$l1_model_setup))
  expect_true("metadata" %in% names(gpa$l1_model_setup))
  
  # test if the fsl element is a data frame
  expect_true(is.data.frame(gpa$l1_model_setup$fsl))
  
  # test if the metadata element is a list
  expect_true(is.list(gpa$l1_model_setup$metadata))

  # test if fsl table contains the id column and has the correct number of subjects
  expect_true("id" %in% names(gpa$l1_model_setup$fsl))
  if ("id" %in% names(gpa$l1_model_setup$fsl)){
    expect_equal(length(unique(gpa$l1_model_setup$fsl$id)), nsubjects_clock_testing)
  }

  # test if fsl table contains all the run numbers for all subjects
  expect_true("run_number" %in% names(gpa$l1_model_setup$fsl))
  if ("run_number" %in% names(gpa$l1_model_setup$fsl)){
    expect_equal(length(unique(gpa$l1_model_setup$fsl$run_number)), nruns_clock_testing)
    run_numbers_per_subject <- table(gpa$l1_model_setup$fsl$id, gpa$l1_model_setup$fsl$run_number)
    expect_true(all(run_numbers_per_subject > 0))
  }

  # test if fsl table has all the l1_models for all subjects and all runs
  expect_true("l1_model" %in% names(gpa$l1_model_setup$fsl))
  if ("l1_model" %in% names(gpa$l1_model_setup$fsl)){
  l1_models_per_subject_run <- table(gpa$l1_model_setup$fsl$id, gpa$l1_model_setup$fsl$run_number, gpa$l1_model_setup$fsl$l1_model)
  expect_true(all(l1_models_per_subject_run > 0))
  
  # test if the l1 models are "entropy" and "pe"
  unique_l1_models <- unique(gpa$l1_model_setup$fsl$l1_model)
  expect_true(all(c("entropy", "pe") %in% unique_l1_models))
  }

  # test if filepaths in feat_fsf exists
  expect_true("feat_fsf" %in% names(gpa$l1_model_setup$fsl))
  if ("feat_fsf" %in% names(gpa$l1_model_setup$metadata)){
    expect_true(all(file.exists(gpa$l1_model_setup$fsl$feat_fsf)))
  }

  # test if metadata contains the expected elements
  expected_metadata_elements <- c("id", "session", "run_number", "run_nifti", "l1_confound_file", "run_volumes", "exclude_run")
  expect_true(all(expected_metadata_elements %in% names(gpa$l1_model_setup$metadata)))
  expect_equal(length(unique(gpa$l1_model_setup$metadata$id)), nsubjects_clock_testing)
  expect_equal(length(unique(gpa$l1_model_setup$metadata$run_number)), nruns_clock_testing)

})