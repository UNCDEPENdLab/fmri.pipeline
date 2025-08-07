#' test setup_l1_models function
test_that("test setup and run l1 models", {

  # -------------------------------------------------------
  # ---- test_l1 runs smoothly ----------------------------
  # -------------------------------------------------------

  analysis_name <- "gpa_tests_build_l1_models"
  test_data_base_dir <- "/proj/mnhallqlab/projects/fmri.pipeline_test_data"
  spec_file <- "sample_2_L1_spec.yaml"
  trial_data_file <- "sample_trial_data.csv.gz"
  run_data_file <- "sample_run_data.csv"
  subject_data_file <- "sample_subject_data.csv"
  gpa_cache_file <- "gpa_tests_build_l1_models.rds"
  cache_file <- "gpa_tests_base_build_l1_models.rds"

  # run building l1, l2, l3 models
  gpa <- create_gpa(
      analysis_name = analysis_name,
      test_data_base_dir = test_data_base_dir,
      spec_file = spec_file,
      trial_data_file = trial_data_file,
      run_data_file = run_data_file,
      subject_data_file = subject_data_file,
      gpa_cache_file = gpa_cache_file,
      cache_file = cache_file)

  # # finalize pipeline configuration
  # gpa <- finalize_pipeline_configuration(gpa)

  # # create all FSF files for level one runs
  # gpa <- setup_l1_models(gpa)

  lg <- lgr::get_logger("glm_pipeline/run_glm_pipeline")
  lg$set_threshold(gpa$lgr_threshold)

  batch_id <- uuid::UUIDgenerate()
  if (is.null(gpa$batch_run)) gpa$batch_run <- list()
  
  # setup a unique directory for this GLM batch submission
  # this will contain all scripts for this run and contain a snapshot of the gpa object
  # that is amended as each batch job runs and completes.
  batch_directory <- file.path(gpa$output_locations$scheduler_scripts, paste0("batch_", batch_id))
  if (!dir.exists(batch_directory)) dir.create(batch_directory, recursive = TRUE)
  gpa_cache <- file.path(batch_directory, "run_pipeline_cache.RData")
  save(gpa, file=gpa_cache)
  
  # batch job to finalize pipeline configuration
  run_finalize <- TRUE
  f_batch <- R_batch_job$new(
    job_name = "finalize_configuration", batch_directory = batch_directory, scheduler = gpa$scheduler,
    input_rdata_file = gpa_cache, output_rdata_file = gpa_cache,
    n_nodes = 1, n_cpus = 1, wall_time = gpa$parallel$finalize_time,
    mem_total = "16G",
    r_code = "gpa <- finalize_pipeline_configuration(gpa)", r_packages = "fmri.pipeline",
    batch_code = get_compute_environment(gpa),
    scheduler_options = gpa$parallel$sched_args
  )

  # initialize batch sequence elements as NULL so that they are ignored in R_batch_sequence if not used
  l1_setup_batch <- l1_execute_batch <- l2_batch <- l3_batch <- NULL

  model_list <- choose_glm_set(gpa, names(gpa$l1_models$models), names(gpa$l2_models$models), names(gpa$l3_models$models), lg) # need manual input
  # batch job for setting up l1 models -- calls setup_l1_models to create relevant FSFs
  l1_setup_batch <- f_batch$copy(
    job_name = "setup_l1", n_cpus = gpa$parallel$l1_setup_cores,
    wall_time = gpa$parallel$l1_setup_time, mem_total = gpa$parallel$l1_setup_memgb,
    r_code = sprintf(
      "gpa <- setup_l1_models(gpa, l1_model_names=%s)", paste(deparse(model_list$l1_model_name), collapse = "")
    )
  )

  l1_setup_batch$depends_on_parents <- "finalize_configuration"

  # batch job for executing l1 jobs (and waiting) after setup
  l1_execute_batch <- f_batch$copy(
    job_name = "run_l1", n_cpus = 1,
    wall_time = gpa$parallel$fsl$l1_feat_alljobs_time,
    r_code = "child_job_ids <- run_feat_sepjobs(gpa, level = 1L)" # execute l1 jobs
  )

  l1_execute_batch$depends_on_parents <- "setup_l1"
  l1_execute_batch$wait_for_children <- TRUE # need to wait for l1 feat jobs to complete before moving to l2/l3

  if (isTRUE(run_finalize)) {
    glm_batch <- R_batch_sequence$new(f_batch, l1_setup_batch, l1_execute_batch)
  } else {
    glm_batch <- R_batch_sequence$new(l1_setup_batch, l1_execute_batch)
  }
  glm_batch$submit()

  # TODO How do I wait for the batch_l1_run to to finish before running the tests?
  # pull the R_batch job ids for finalize_configuration, setup_l1 and run_l1
  R_batch_jobids <- glm_batch$get_job_ids()

  # wait for the jobs to be over
  message("Waiting for the setup_l1 and run_l1 jobs to be completed.")
  wait_for_job(R_batch_jobids, quiet = FALSE, scheduler = "slurm")

  # gpa_new <- readRDS("/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/gpa_tests_build_l1_models/gpa_tests_build_l1_models.rds")
  load(file.path(batch_directory, "run_pipeline_cache.RData"))
  
  # run one job and refresh false to true feat_complete
  # have a test_l1 that crashes
    # nvolumes time points in regressors
    # nifti with 200 volumes regressors with longer (twice as long) long Tr 2 instead of 1

    # redundant contransts
    # files doesnt exist confound ev file nusance timeseries
    # run some setup and brake it break sff file afte running finalize_config
    # hack run_data dataframe
    # point to the wrong confounds file

  # --- testing when everything runs smoothly ----

  # test if l1 model fsl feat_complete is set to true
  expect_true(all(gpa$l1_model_setup$fsl$feat_complete))
  # subid 10767 run 4 both pe and entropy model fail

  # gpa$l1_model_setup is the output created by setup_l1_models()
  nsubjects_clock_testing <- 10
  nruns_clock_testing <- 8
  number_of_l1_models <- 2

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

  # ---------------------------------------------------------
  # ---- test_l1 that crashes -------------------------------
  # ---------------------------------------------------------
  
  # ---- point to the wrong confounds file ----
  # all subjects clock8 confound file is replaced with clock1 confound file

  analysis_name <- "gpa_tests_build_l1_models_wrong_confound"
  test_data_base_dir <- "/proj/mnhallqlab/projects/fmri.pipeline_test_data"
  spec_file <- "sample_2_L1_spec.yaml"
  trial_data_file <- "sample_trial_data.csv.gz"
  run_data_file <- "sample_run_data_wrong_confound.csv"
  subject_data_file <- "sample_subject_data.csv"
  gpa_cache_file <- "gpa_tests_build_l1_models_wrong_confound.rds"
  cache_file <- "gpa_tests_base_build_l1_models_wrong_confound.rds"

  # run building l1, l2, l3 models
  gpa_wrong <- create_gpa(
      analysis_name = analysis_name,
      test_data_base_dir = test_data_base_dir,
      spec_file = spec_file,
      trial_data_file = trial_data_file,
      run_data_file = run_data_file,
      subject_data_file = subject_data_file,
      gpa_cache_file = gpa_cache_file,
      cache_file = cache_file)

  # # finalize pipeline configuration
  # gpa_wrong <- finalize_pipeline_configuration(gpa_wrong)

  # # create all FSF files for level one runs
  # gpa_wrong <- setup_l1_models(gpa_wrong)
  # save(gpa_wrong, file="/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/fmri.pipeline/tests/testthat/testdata/gpa_wrong_confound.RData")

  lg <- lgr::get_logger("glm_pipeline/run_glm_pipeline")
  lg$set_threshold(gpa_wrong$lgr_threshold)

  batch_id <- uuid::UUIDgenerate()
  if (is.null(gpa_wrong$batch_run)) gpa_wrong$batch_run <- list()
  
  # setup a unique directory for this GLM batch submission
  # this will contain all scripts for this run and contain a snapshot of the gpa object
  # that is amended as each batch job runs and completes.
  batch_directory_wrong <- file.path(gpa_wrong$output_locations$scheduler_scripts, paste0("batch_", batch_id))
  if (!dir.exists(batch_directory_wrong)) dir.create(batch_directory_wrong, recursive = TRUE)
  gpa_cache <- file.path(batch_directory_wrong, "run_pipeline_cache.RData")
  save(gpa_wrong, file=gpa_cache)
  
  # batch job to finalize pipeline configuration
  run_finalize <- TRUE
  f_batch <- R_batch_job$new(
    job_name = "finalize_configuration", batch_directory = batch_directory_wrong, scheduler = gpa_wrong$scheduler,
    input_rdata_file = gpa_cache, output_rdata_file = gpa_cache,
    n_nodes = 1, n_cpus = 1, wall_time = gpa_wrong$parallel$finalize_time,
    mem_total = "16G",
    r_code = "gpa_wrong <- finalize_pipeline_configuration(gpa_wrong)", r_packages = "fmri.pipeline",
    batch_code = get_compute_environment(gpa_wrong),
    scheduler_options = gpa_wrong$parallel$sched_args
  )

  # initialize batch sequence elements as NULL so that they are ignored in R_batch_sequence if not used
  l1_setup_batch <- l1_execute_batch <- l2_batch <- l3_batch <- NULL

  model_list <- choose_glm_set(gpa_wrong, names(gpa_wrong$l1_models$models), names(gpa_wrong$l2_models$models), names(gpa_wrong$l3_models$models), lg) # need manual input
  # batch job for setting up l1 models -- calls setup_l1_models to create relevant FSFs
  l1_setup_batch <- f_batch$copy(
    job_name = "setup_l1", n_cpus = gpa_wrong$parallel$l1_setup_cores,
    wall_time = gpa_wrong$parallel$l1_setup_time, mem_total = gpa_wrong$parallel$l1_setup_memgb,
    r_code = sprintf(
      "gpa_wrong <- setup_l1_models(gpa_wrong, l1_model_names=%s)", paste(deparse(model_list$l1_model_name), collapse = "")
    )
  )

  l1_setup_batch$depends_on_parents <- "finalize_configuration"

  # batch job for executing l1 jobs (and waiting) after setup
  l1_execute_batch <- f_batch$copy(
    job_name = "run_l1", n_cpus = 1,
    wall_time = gpa_wrong$parallel$fsl$l1_feat_alljobs_time,
    r_code = "child_job_ids <- run_feat_sepjobs(gpa_wrong, level = 1L)" # execute l1 jobs
  )

  l1_execute_batch$depends_on_parents <- "setup_l1"
  l1_execute_batch$wait_for_children <- TRUE # need to wait for l1 feat jobs to complete before moving to l2/l3

  if (isTRUE(run_finalize)) {
    glm_batch <- R_batch_sequence$new(f_batch, l1_setup_batch, l1_execute_batch)
  } else {
    glm_batch <- R_batch_sequence$new(l1_setup_batch, l1_execute_batch)
  }
  glm_batch$submit()

  # pull the R_batch job ids
  R_batch_jobids <- glm_batch$get_job_ids()

  # wait for the jobs to be over
  message("Waiting for the setup_l1 and run_l1 jobs to be completed.")
  wait_for_job(R_batch_jobids, quiet = FALSE)

  load(file.path(batch_directory_wrong, "run_pipeline_cache.RData"))
  gpa_wrong <- run_feat_status(gpa_wrong, level = 1L)
  # load("/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/gpa_tests_build_l1_models_wrong_confound/scheduler_scripts/batch_0d6e963f-e60c-4622-8388-99bb35c63c94/run_pipeline_cache.RData")
  
  # testing if things happen as expected when the confound is pointed to the wrong confounds file
  # test if l1 model fsl feat_complete is set to false
  expect_equal(sum(!gpa_wrong$l1_model_setup$fsl$feat_complete), number_of_l1_models*nsubjects_clock_testing)
  
  # l1_confound_file column will be NA for the run with wrong confound file
  expect_true("l1_confound_file" %in% names(gpa_wrong$l1_model_setup$fsl))
  if ("l1_confound_file" %in% names(gpa_wrong$l1_model_setup$fsl)){
    expect_equal(sum(is.na(gpa_wrong$l1_model_setup$fsl$l1_confound_file)), number_of_l1_models*nsubjects_clock_testing) # one run each subject 
  }

  # test that feat_complete is FALSE for all runs with missing confound files
  if ("feat_complete" %in% names(gpa_wrong$l1_model_setup$fsl)) {
    idx <- which(is.na(gpa_wrong$l1_model_setup$metadata$l1_confound_file))
    for (i in idx) {
      subj <- gpa_wrong$l1_model_setup$metadata$id[i]
      run <- gpa_wrong$l1_model_setup$metadata$run_number[i]
      expect_true(all(!gpa_wrong$l1_model_setup$fsl$feat_complete[gpa_wrong$l1_model_setup$fsl$id == subj &
                                                                  gpa_wrong$l1_model_setup$fsl$run_number == run]))
    }
  }
  
  # ---- nvolumes time points in regressors are mismatched ----
  # What are the 4 columns in nuisance_regressor.txt? .csf_ts .csf_ts_deriv .wm_ts .wm_ts_deriv. The rows are the timepoints/volumns




  # ---- nifti with 200 volumes regressors with longer (twice as long) Tr 2 instead of 1 ----
  # set tr to 2 instead of 1 in setup_glm_pipeline() in setup.R
 
  analysis_name <- "gpa_tests_build_l1_models_wrong_tr"
  test_data_base_dir <- "/proj/mnhallqlab/projects/fmri.pipeline_test_data"
  spec_file <- "sample_2_L1_spec.yaml"
  trial_data_file <- "sample_trial_data.csv.gz"
  run_data_file <- "sample_run_data_wrong_confound.csv"
  subject_data_file <- "sample_subject_data.csv"
  gpa_cache_file <- "gpa_tests_build_l1_models_wrong_tr.rds"
  cache_file <- "gpa_tests_base_build_l1_models_wrong_tr.rds"

  # run building l1, l2, l3 models
  gpa_wrong_tr <- create_gpa_wrong_tr(
      analysis_name = analysis_name,
      test_data_base_dir = test_data_base_dir,
      spec_file = spec_file,
      trial_data_file = trial_data_file,
      run_data_file = run_data_file,
      subject_data_file = subject_data_file,
      gpa_cache_file = gpa_cache_file,
      cache_file = cache_file)

  lg <- lgr::get_logger("glm_pipeline/run_glm_pipeline")
  lg$set_threshold(gpa_wrong_tr$lgr_threshold)

  batch_id <- uuid::UUIDgenerate()
  if (is.null(gpa_wrong_tr$batch_run)) gpa_wrong_tr$batch_run <- list()
  
  # setup a unique directory for this GLM batch submission
  # this will contain all scripts for this run and contain a snapshot of the gpa object
  # that is amended as each batch job runs and completes.
  batch_directory_wrong <- file.path(gpa_wrong_tr$output_locations$scheduler_scripts, paste0("batch_", batch_id))
  if (!dir.exists(batch_directory_wrong)) dir.create(batch_directory_wrong, recursive = TRUE)
  gpa_cache <- file.path(batch_directory_wrong, "run_pipeline_cache.RData")
  save(gpa_wrong_tr, file=gpa_cache)
  
  # batch job to finalize pipeline configuration
  run_finalize <- TRUE
  f_batch <- R_batch_job$new(
    job_name = "finalize_configuration", batch_directory = batch_directory_wrong, scheduler = gpa_wrong_tr$scheduler,
    input_rdata_file = gpa_cache, output_rdata_file = gpa_cache,
    n_nodes = 1, n_cpus = 1, wall_time = gpa_wrong_tr$parallel$finalize_time,
    mem_total = "16G",
    r_code = "gpa_wrong_tr <- finalize_pipeline_configuration(gpa_wrong_tr)", r_packages = "fmri.pipeline",
    batch_code = get_compute_environment(gpa_wrong_tr),
    scheduler_options = gpa_wrong_tr$parallel$sched_args
  )

  # initialize batch sequence elements as NULL so that they are ignored in R_batch_sequence if not used
  l1_setup_batch <- l1_execute_batch <- l2_batch <- l3_batch <- NULL

  model_list <- choose_glm_set(gpa_wrong_tr, names(gpa_wrong_tr$l1_models$models), names(gpa_wrong_tr$l2_models$models), names(gpa_wrong_tr$l3_models$models), lg) # need manual input
  # batch job for setting up l1 models -- calls setup_l1_models to create relevant FSFs
  l1_setup_batch <- f_batch$copy(
    job_name = "setup_l1", n_cpus = gpa_wrong_tr$parallel$l1_setup_cores,
    wall_time = gpa_wrong_tr$parallel$l1_setup_time, mem_total = gpa_wrong_tr$parallel$l1_setup_memgb,
    r_code = sprintf(
      "gpa_wrong_tr <- setup_l1_models(gpa_wrong_tr, l1_model_names=%s)", paste(deparse(model_list$l1_model_name), collapse = "")
    )
  )

  l1_setup_batch$depends_on_parents <- "finalize_configuration"

  # batch job for executing l1 jobs (and waiting) after setup
  l1_execute_batch <- f_batch$copy(
    job_name = "run_l1", n_cpus = 1,
    wall_time = gpa_wrong_tr$parallel$fsl$l1_feat_alljobs_time,
    r_code = "child_job_ids <- run_feat_sepjobs(gpa_wrong_tr, level = 1L)" # execute l1 jobs
  )

  l1_execute_batch$depends_on_parents <- "setup_l1"
  l1_execute_batch$wait_for_children <- TRUE # need to wait for l1 feat jobs to complete before moving to l2/l3

  if (isTRUE(run_finalize)) {
    glm_batch <- R_batch_sequence$new(f_batch, l1_setup_batch, l1_execute_batch)
  } else {
    glm_batch <- R_batch_sequence$new(l1_setup_batch, l1_execute_batch)
  }
  glm_batch$submit()

  # pull the R_batch job ids
  R_batch_jobids <- glm_batch$get_job_ids()

  # wait for the jobs to be over
  message("Waiting for the setup_l1 and run_l1 jobs to be completed.")
  wait_for_job(R_batch_jobids, quiet = FALSE)

  load(file.path(batch_directory_wrong, "run_pipeline_cache.RData"))

  # TODO NEXT Add tests for different tr


  # ---- break fsf file by pointing to a non-existing nifti file ----
  # run setup and then break fsf file after finalize_configuration

  analysis_name <- "gpa_tests_build_l1_models_break_fsf"
  test_data_base_dir <- "/proj/mnhallqlab/projects/fmri.pipeline_test_data"
  spec_file <- "sample_2_L1_spec.yaml"
  trial_data_file <- "sample_trial_data.csv.gz"
  run_data_file <- "sample_run_data_wrong_confound.csv"
  subject_data_file <- "sample_subject_data.csv"
  gpa_cache_file <- "gpa_tests_build_l1_models_break_fsf.rds"
  cache_file <- "gpa_tests_base_build_l1_models_break_fsf.rds"

  # run building l1, l2, l3 models
  gpa_breakfsf <- create_gpa( # Using create_gpa() which creates the correct gpa object
      analysis_name = analysis_name,
      test_data_base_dir = test_data_base_dir,
      spec_file = spec_file,
      trial_data_file = trial_data_file,
      run_data_file = run_data_file,
      subject_data_file = subject_data_file,
      gpa_cache_file = gpa_cache_file,
      cache_file = cache_file)

  lg <- lgr::get_logger("glm_pipeline/run_glm_pipeline")
  lg$set_threshold(gpa_breakfsf$lgr_threshold)

  batch_id <- uuid::UUIDgenerate()
  if (is.null(gpa_breakfsf$batch_run)) gpa_breakfsf$batch_run <- list()
  
  # setup a unique directory for this GLM batch submission
  # this will contain all scripts for this run and contain a snapshot of the gpa object
  # that is amended as each batch job runs and completes.
  batch_directory <- file.path(gpa_breakfsf$output_locations$scheduler_scripts, paste0("batch_", batch_id))
  if (!dir.exists(batch_directory)) dir.create(batch_directory, recursive = TRUE)
  gpa_cache <- file.path(batch_directory, "run_pipeline_cache.RData")
  save(gpa_breakfsf, file=gpa_cache)
  
  # finalize pipeline configuration
  gpa_breakfsf <- finalize_pipeline_configuration(gpa_breakfsf)

  # load fsf file for each run and each subject and change the run_nifti to a non-existing file
  for (i in 1:nrow(gpa_breakfsf$l1_model_setup$fsl)) {
    gpa_breakfsf$l1_model_setup$fsl$run_nifti[i] <- "/proj/mnhallqlab/projects/fmri.pipeline_test_data/clock8/func/sub-10001_ses-1_task-clock8_run-1_bold.nii.gz"
  }
  # save the gpa object with the broken fsf file
  save(gpa_breakfsf, file=gpa_cache)

  



  # ---- files doesnt exist confound ev file nuisance timeseries ----
  # confound regressors are white matter and csf timeseries, nuisance timeseries are the same thing




})
