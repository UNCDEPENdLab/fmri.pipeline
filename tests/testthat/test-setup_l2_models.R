#' test setup_l2_models function
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
  # gpa <- finalize_pipeline_configuration(gpa)

  # # create all FSF files for level one runs
  # gpa <- setup_l1_models(gpa)

  # # running setup_l2_models
  # gpa <- setup_l2_models(gpa)


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

  model_list <- choose_glm_set(gpa, names(gpa$l1_models$models), names(gpa$l2_models$models), names(gpa$l3_models$models), lg)
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

  # only run level 2 if this is a multi-run dataset and user requests l2 model estimation
  # setup of l2 models (should follow l1)
  l2_batch <- f_batch$copy(
    job_name = "setup_run_l2", n_cpus = gpa$parallel$l2_setup_cores,
    wall_time = gpa$parallel$l2_setup_run_time,
    r_code = c(
      "gpa <- setup_l2_models(gpa)",
      "child_job_ids <- run_feat_sepjobs(gpa, level = 2L)"
    )
  )
  l2_batch$depends_on_parents <- "run_l1"
  l2_batch$wait_for_children <- TRUE # need to wait for l2 feat jobs to complete before moving to l3
  
  if (isTRUE(run_finalize)) {
    glm_batch <- R_batch_sequence$new(f_batch, l1_setup_batch, l1_execute_batch, l2_batch)
  } else {
    glm_batch <- R_batch_sequence$new(l1_setup_batch, l1_execute_batch, l2_batch)
  }
  glm_batch$submit()

  load("/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/gpa_tests_build_l1_models/scheduler_scripts/batch_6e2f513e-e0a3-4a76-a267-0f89314d787b/run_pipeline_cache.RData")
  gpa$l1_model_setup$fsl$feat_complete
  gpa$l2_model_setup$fsl$feat_complete


  # ---- testing ----
  # test if l1 model fsl feat_complete is now set to true
  expect_true(all(gpa$l1_model_setup$fsl$feat_complete))

  # gpa$l2_model_setup is the output created by setup_l2_models()
  nsubjects_clock_testing <- 10
  nruns_clock_testing <- 8

  # test if the class of gpa$l2_model_setup is correct
  expect_equal(class(gpa$l2_model_setup)[[2]], "l2_setup")

  # test if the output contains the expected elements
  expect_true("fsl" %in% names(gpa$l2_model_setup))

  # test if the fsl element is a data frame
  expect_true(is.data.frame(gpa$l2_model_setup$fsl))

  # test if fsl table contains the id column and has the correct number of subjects
  expect_true("id" %in% names(gpa$l2_model_setup$fsl))
  if ("id" %in% names(gpa$l2_model_setup$fsl)){
    expect_equal(length(unique(gpa$l2_model_setup$fsl$id)), nsubjects_clock_testing)
  }

  # test if l1_model column has all the l1 models
  expect_true("l1_model" %in% names(gpa$l2_model_setup$fsl))
  if ("l1_model" %in% names(gpa$l2_model_setup$fsl)){
    expect_true(all(c("entropy", "pe") %in% unique(gpa$l2_model_setup$fsl$l1_model)))
  }

  # test if l2_model column has all the l2 models
  expect_true("l2_model" %in% names(gpa$l2_model_setup$fsl))
  if ("l2_model" %in% names(gpa$l2_model_setup$fsl)){
    expect_equal(all(unique(gpa$l2_model_setup$fsl$l2_model), "emotion"))
  }
  
  # test if filepaths in feat_fsf exists
  expect_true("feat_fsf" %in% names(gpa$l2_model_setup$fsl))
  if ("feat_fsf" %in% names(gpa$l2_model_setup$fsl)){
    expect_true(all(file.exists(gpa$l2_model_setup$fsl$feat_fsf)))
  }

})
