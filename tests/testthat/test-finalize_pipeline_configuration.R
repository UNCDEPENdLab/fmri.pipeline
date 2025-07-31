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
      spec_file = l1_spec_file,
      trial_data_file = trial_data_file,
      run_data_file = run_data_file,
      subject_data_file = subject_data_file,
      gpa_cache_file = gpa_cache_file,
      cache_file = cache_file)

  # finalize pipeline configuration
  gpa <- finalize_pipeline_configuration(gpa)

  # check if the finalize pipeline configuration worked
  expect_true(gpa$finalize_complete)

  # check if additional run_data from get_l1_confounds output was added to the gpa object
  expect_equal(ncol(gpa$run_data), 33)

  # check if nifti information is added correctly
  expect_equal(sum(gpa$run_data$run_nifti_present), 80) # check if each run_nifti_present is TRUE
  expect_equal(sum(gpa$run_data$run_volumes > 200), 80) # check if each run_volumn is greater than 200
  expect_true(all(gpa$run_data$dim_x > 0)) # check if each element of gpa$run_data$dim_x is greater than 0
  expect_true(all(gpa$run_data$dim_y > 0)) # check if each element of gpa$run_data$dim_y is greater than 0
  expect_true(all(gpa$run_data$dim_z > 0)) # check if each element of gpa$run_data$dim_z is greater than 0

  # check if motion parameter file is present
  expect_equal(sum(gpa$run_data$motion_params_present), 80)

  # is confound settings processed correctly
  expect_equal(sum(gpa$run_data$confound_input_file_present), 80)

  # check if drop volumes is processed correctly
  expect_equal(sum(gpa$run_data$drop_volumes == 2), 80)

  # check if correct number of runs were truncated
  expect_equal(sum(gpa$run_data$truncated_run), 0)

  # check if correct number of runs were excluded
  expect_equal(sum(gpa$run_data$exclude_run), 3)
  expect_equal(sum(gpa$subject_data$n_good_runs == 8), 8) # only 8 out of 10 subjects have all 8 good runs

  # check if correct number of subjects were excluded
  expect_equal(sum(gpa$subject_data$exclude_subject), 0)

  # test with flankers dataset -----------------------------------------------

  analysis_name = "gpa_flanker_tests"
  test_data_base_dir = "/proj/mnhallqlab/no_backup/flanker-fmriprep"
  spec_file = "/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/fmri.pipeline/tests/testthat/testdata/int.yaml"
  trial_data_file = "/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/fmri.pipeline/tests/testthat/testdata/flankers_trial_data.csv"
  run_data_file = "/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/fmri.pipeline/tests/testthat/testdata/flankers_run_data.csv"
  subject_data_file = "/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/fmri.pipeline/tests/testthat/testdata/flankers_subject_data.csv"
  gpa_cache_file = "/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/fmri.pipeline/tests/testthat/testdata/gpa.rds"
  cache_file = "/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/fmri.pipeline/tests/testthat/testdata/gpa_base.rds"

  gpa_flanker <- create_flanker_gpa(
    analysis_name = analysis_name,
    test_data_base_dir = test_data_base_dir,
    spec_file = spec_file,
    trial_data_file = trial_data_file,
    run_data_file = run_data_file,
    subject_data_file = subject_data_file,
    gpa_cache_file = gpa_cache_file,
    cache_file = cache_file
  )

  # finalize pipeline configuration
  gpa_flanker <- finalize_pipeline_configuration(gpa_flanker)

  # check if the finalize pipeline configuration worked
  expect_true(gpa_flanker$finalize_complete)

  gpa_flanker$run_data$exclude_run

  # test that crashes -----------------------------------------------
  
  # ---- point to the wrong confounds file ----
  # all subjects clock8 confound file is replaced with clock1 confound file

  analysis_name <- "gpa_tests_build_l1_models_wrong_confound"
  test_data_base_dir <- "/proj/mnhallqlab/projects/fmri.pipeline_test_data"
  l1_spec_file <- "sample_2_L1_spec.yaml"
  trial_data_file <- "sample_trial_data.csv.gz"
  run_data_file <- "sample_run_data_wrong_confound.csv"
  subject_data_file <- "sample_subject_data.csv"
  gpa_cache_file <- "gpa_tests_build_l1_models_wrong_confound.rds"
  cache_file <- "gpa_tests_base_build_l1_models_wrong_confound.rds"

  # run building l1, l2, l3 models
  gpa_wrong <- create_gpa(
      analysis_name = analysis_name,
      test_data_base_dir = test_data_base_dir,
      l1_spec_file = l1_spec_file,
      trial_data_file = trial_data_file,
      run_data_file = run_data_file,
      subject_data_file = subject_data_file,
      gpa_cache_file = gpa_cache_file,
      cache_file = cache_file)

  # finalize pipeline configuration
  gpa_wrong <- finalize_pipeline_configuration(gpa_wrong)

  # TODO What should I expect different?


})