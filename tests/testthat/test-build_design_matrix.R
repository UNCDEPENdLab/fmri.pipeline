#' test build_design_matrix function
setup_l1_models_partial <- function(gpa, l1_model_names = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_character(l1_model_names, null.ok = TRUE)
  checkmate::assert_subset(l1_model_names, names(gpa$l1_models$models))
  checkmate::assert_data_frame(gpa$subject_data)
  checkmate::assert_data_frame(gpa$run_data)
  checkmate::assert_subset(c("id", "session", "run_number", "run_nifti", "exclude_run"), names(gpa$run_data)) # required columns

  #if no model subset is requested, output all models
  if (is.null(l1_model_names)) l1_model_names <- names(gpa$l1_models$models)

  lg <- lgr::get_logger("glm_pipeline/l1_setup")
  lg$set_threshold(gpa$lgr_threshold)

  if (isTRUE(gpa$log_json) && !"setup_l1_log_json" %in% names(lg$appenders)) {
    lg$add_appender(lgr::AppenderJson$new(gpa$output_locations$setup_l1_log_json), name = "setup_l1_log_json")
  }

  if (isTRUE(gpa$log_txt) && !"setup_l1_log_txt" %in% names(lg$appenders)) {
    lg$add_appender(lgr::AppenderFile$new(gpa$output_locations$setup_l1_log_txt), name = "setup_l1_log_txt")
  }

  # TODO: This is a mess if use the l1_model_names since we will always be overwriting what's in the l1_model_setup
  # data.frame. We need more of an append/update approach, perhaps like the sqlite setup for the overall pipeline.

  # refresh status of feat inputs and outputs at level 1 before continuing
  gpa <- refresh_feat_status(gpa, level=1L, lg=lg)
  gpa$parallel$l1_setup_cores <- 2
  #setup parallel worker pool, if requested
  if (!is.null(gpa$parallel$l1_setup_cores) && gpa$parallel$l1_setup_cores[1L] > 1L) {
    lg$info("Initializing l1 setup cluster with %d cores", gpa$parallel$l1_setup_cores[1L])
    cl <- parallel::makeCluster(gpa$parallel$l1_setup_cores[1L])
    doParallel::registerDoParallel(cl)
    on.exit(try(parallel::stopCluster(cl))) #cleanup pool upon exit of this function
  } else {
    lg$info("Initializing l1 setup with serial execution")
    foreach::registerDoSEQ() #formally register a sequential 'pool' so that dopar is okay
  }

  # loop over each subject, identify relevant fMRI data, and setup level 1 analysis files
  # all_subj_l1_list <- lapply(seq_len(nrow(gpa$subject_data)), function(ii) {
  #   subj_df <- gpa$subject_data[ii, ]

  all_subj_l1_list <- foreach(
    subj_df = iter(gpa$subject_data, by = "row"), .inorder = FALSE, 
    .packages = c("dplyr", "fmri.pipeline", "glue"), 
    .errorhandling = "stop",
    .export = c(
      "truncate_runs", "fsl_l1_model", "get_mr_abspath",
      "get_output_directory", "run_fsl_command", "get_feat_status", "add_custom_feat_syntax"
    ), .combine = rbind
  ) %dopar% {
  # all_subj_l1_list <- data.frame(subj_id = 0, subj_session = 0, bdm_out_files = "")
  # for (i in 1:nrow(gpa$subject_data)) {
    # print(toString(i))
    # subj_df <- gpa$subject_data[i, , drop = FALSE]
    subj_df <- subj_df # avoid complaints about visible global binding in R CMD check
    subj_id <- subj_df$id
    subj_session <- subj_df$session
    all_subj_l1_list_sub <- data.frame(subj_id = as.numeric(), subj_session = as.numeric(), bdm_out_files = as.character())  # Store bdm_out_file paths for this subject

    # find the run data for analysis
    rdata <- gpa$run_data %>% dplyr::filter(id == !!subj_id & session == !!subj_session & run_nifti_present == TRUE)
    run_nifti <- get_mr_abspath(rdata, "run_nifti")
    mr_run_nums <- rdata$run_number
    run_volumes <- rdata$run_volumes
    nvoxels <- rdata$nvoxels
    exclude_run <- rdata$exclude_run

    # ensure that we have a TR column and that TR does not vary by run
    if (!"tr" %in% names(rdata)) {
      lg$error("No tr column in run data for subject: %s, session: %d", subj_id, subj_session)
      return(NULL)
    }

    if (length(unique(rdata$tr)) > 1L) {
      lg$error("More than one TR value for runs within a subject. This is not currently supported! subject: %s, session: %d", subj_id, subj_session)
      return(NULL)
    }

    if ("l1_confound_file" %in% names(rdata)) {
      l1_confound_files <- rdata$l1_confound_file
    } else {
      l1_confound_files <- rep("", length(run_nifti))
    }

    subj_mr_dir <- subj_df$mr_dir

    # process files
    if (length(run_nifti) == 0L) {
      lg$warn("Unable to find any preprocessed fMRI files in dir: %s", subj_mr_dir)
      return(NULL)
    } else if (any(!file.exists(run_nifti))) {
      nii_present <- file.exists(run_nifti)
      lg$warn("Could not find some of the expected preprocessed fMRI files. These will be dropped.")
      lg$warn("Missing: %s", run_nifti[!nii_present])
      run_nifti <- run_nifti[nii_present]
      mr_run_nums <- mr_run_nums[nii_present]
      l1_confound_files <- l1_confound_files[nii_present]
      run_volumes <- run_volumes[nii_present]
      exclude_run <- exclude_run[nii_present]
    } else {
      lg$debug(paste("MR files to analyze:", run_nifti)) # log files that were found
    }

    lg$debug("Volumes in run_nifti: %s", paste(run_volumes, collapse = ", "))

    # get all events that pertain to this participant
    m_events <- data.table::rbindlist(
      lapply(gpa$l1_models$events, function(this_event) {
        this_event$data %>% dplyr::filter(id == !!subj_id & session == !!subj_session)
      })
    )

    # lookup subject directory for placing truncated files (now handled upstream)
    # subj_output_directory <- get_output_directory(id = subj_id, session = subj_session, gpa = gpa, create_if_missing = FALSE, what = "sub")

    # initialize mr data frame with elements of $run_data
    mr_df <- data.frame(
      id = subj_id, session = subj_session, run_number = mr_run_nums, run_nifti, l1_confound_file = l1_confound_files,
      run_volumes, exclude_run, row.names = NULL
    )

    # Tracking list containing data.frames for each software, where we expect one row per run-level model (FSL)
    # or subject-level model (AFNI). The structure varies because FSL estimates per-run GLMs, while AFNI concatenates.
    l1_file_setup <- list(fsl = list(), spm = list(), afni = list(), metadata = mr_df)

    # loop over models to output
    for (ii in seq_along(l1_model_names)) {
      this_model <- l1_model_names[ii]

      # setup design matrix for any given software package
      m_signals <- lapply(gpa$l1_models$signals[gpa$l1_models$models[[this_model]]$signals], function(this_signal) {
        # filter down to this id if the signal is a data.frame
        if (inherits(this_signal$value, "data.frame")) {
          this_signal$value <- this_signal$value %>%
            dplyr::filter(id == !!subj_id & session == !!subj_session)

          if (nrow(this_signal$value) == 0L) {
            msg <- glue(
              "In L1 model setup, failed to find any rows in in gpa$l1_models$signals${this_signal$name}$value",
              " for id: {subj_id}, session: {subj_session}, model: {this_model}.\n  We will create an empty regressor.",
              .trim = FALSE
            )
            lg$warn(msg)
            warning(msg)
          }

          # refit wi model if needed
          if (!is.null(this_signal$wi_model) && nrow(this_signal$value) > 0L) {
            this_signal <- fit_wi_model(this_signal)
          }
        } else {
          msg <- "Unable to sort out how to refit wi_model with signal that doesn't have a data.frame in the $value slot"
          lg$error(msg)
          stop(msg)
        }
        return(this_signal)
      })

      l1_output_dir <- get_output_directory(
        id = subj_id, session = subj_session,
        l1_model = this_model, gpa = gpa, what = "l1"
      )

      if (!dir.exists(l1_output_dir)) {
        lg$info("Creating subject output directory: %s", l1_output_dir)
        dir.create(l1_output_dir, showWarnings = FALSE, recursive = TRUE)
      }

      bdm_out_file <- file.path(l1_output_dir, paste0(gpa$l1_models$models[[this_model]]$name, "_bdm_setup.RData"))
      # Save the bdm_out_file path
      # bdm_files[[this_model]] <- bdm_out_file

      run_bdm <- TRUE
      if (run_bdm) {
        t_out <- gpa$glm_software
        if (isTRUE(gpa$use_preconvolve)) t_out <- c("convolved", t_out) # compute preconvolved regressors
        bdm_args <- gpa$additional$bdm_args
        bdm_args$events <- m_events
        bdm_args$signals <- m_signals
        bdm_args$tr <- rdata$tr[1L]
        bdm_args$write_timing_files <- t_out
        bdm_args$drop_volumes <- gpa$drop_volumes
        bdm_args$run_data <- mr_df
        bdm_args$runs_to_output <- mr_run_nums
        bdm_args$output_directory <- file.path(l1_output_dir, "timing_files")
        bdm_args$lg <- lg

        d_obj <- tryCatch(do.call(build_design_matrix, bdm_args), error = function(e) {
          lg$error("Failed build_design_matrix for id: %s, session: %s, model: %s", subj_id, subj_session, this_model)
          lg$error("Error message: %s", as.character(e))
          return(NULL)
        })

        save(d_obj, bdm_args, mr_df, mr_run_nums, subj_mr_dir, run_nifti, run_volumes, exclude_run, nvoxels,
          l1_confound_files, subj_id, subj_session, this_model,
          file = bdm_out_file
        )
        all_subj_l1_list_sub <- all_subj_l1_list_sub %>% rbind(data.frame(subj_id = subj_id, subj_session = subj_session, bdm_out_files = bdm_out_file))
      }
    }
    # return(list(subj_id = subj_id, subj_session = subj_session, bdm_out_files = bdm_files))
    return(all_subj_l1_list_sub)
  }
  
  # bdm_paths <- lapply(all_subj_l1_list, function(x) x$bdm_out_files)

  return(all_subj_l1_list)
}

test_that("test build design matrix", {
  
  analysis_name <- "gpa_tests_build_design_matrix"
  test_data_base_dir <- "/proj/mnhallqlab/projects/fmri.pipeline_test_data"
  l1_spec_file <- "sample_2_L1_spec.yaml"
  trial_data_file <- "sample_trial_data.csv.gz"
  run_data_file <- "sample_run_data.csv"
  subject_data_file <- "sample_subject_data.csv"
  gpa_cache_file <- "gpa_tests_build_design_matrix.rds"
  cache_file <- "gpa_tests_base_build_design_matrix.rds"

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

  all_subj_l1_list <- setup_l1_models_partial(gpa)

  # temporary for testing code
  # save.image("/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/fmri.pipeline/tests/testthat/testdata/temp-test-build_design_matrix.RData")
  
  # test if build_design_matrix did what it was supposed to do
  # check if d_obj is created properly in bdm_files
  nsignal_testing <- 3 # for both entropy and pe l1_models we have 3 signals
  nruns_testing <- 8

  for (s in nrow(all_subj_l1_list)){
    load(all_subj_l1_list[s, "bdm_out_files"])

    # testing design element of d_obj
    expect_true("design" %in% names(d_obj))
    if ("design" %in% names(d_obj)){
      expect_equal(ncol(d_obj$design), nsignal_testing)
      expect_equal(nrow(d_obj$design), nruns_testing)
      for (j in 1:nsignal_testing){
        for (i in 1:nruns_testing){
          expect_equal(ncol(d_obj$design[i,j][[1]]), 4)
        }
      }
    }

    # testing design_convolved element of d_obj
    expect_true("design_convolved" %in% names(d_obj))
    if ("design_convolved" %in% names(d_obj)){
      for (r in 1:nruns_testing){
        expect_equal(ncol(d_obj$design_convolved[[paste0("run",toString(r))]]), 6)
        expect_equal(sum(grepl("base", colnames(d_obj$design_convolved[[paste0("run",toString(r))]]))), 3) # testing if columns for polynomial baseline regressors were added
      }
    }

    # testing design_unconvolved element of d_obj
    expect_true("design_unconvolved" %in% names(d_obj))
    if ("design_unconvolved" %in% names(d_obj)){
      for (r in 1:nruns_testing){
        expect_equal(ncol(d_obj$design_unconvolved[[paste0("run",toString(r))]]), 3)
      }
    }

    # testing collin_events (collinearity of regressors before convolution) element of d_obj
    expect_true("collin_events" %in% names(d_obj))
    if ("collin_events" %in% names(d_obj)){
      for(r in 1:nruns_testing){
        # correlation matrix of the regressors
        expect_equal(nrow(d_obj$collin_events[[paste0("run_number",toString(r))]][["r"]]), nsignal_testing)
        expect_equal(ncol(d_obj$collin_events[[paste0("run_number",toString(r))]][["r"]]), nsignal_testing)
        expect_true(is.na(d_obj$collin_events[[paste0("run_number",toString(r))]][["r"]][1,1]))
        expect_equal(d_obj$collin_events[[paste0("run_number",toString(r))]][["r"]][nsignal_testing,nsignal_testing], 1) # only entropy_clock or pe_clock will have the correlation as 1, others will be NA
      }
    }

    # testing collin_convolved (collinearity of the convolved regressors) element of d_obj
    expect_true("collin_convolved" %in% names(d_obj))
    if ("collin_convolved" %in% names(d_obj)){
      for(r in 1:nruns_testing){
        expect_equal(nrow(d_obj$collin_convolved[[paste0("run",toString(r))]][["r"]]), nsignal_testing)
        expect_equal(ncol(d_obj$collin_convolved[[paste0("run",toString(r))]][["r"]]), nsignal_testing)
        expect_true(all(diag(d_obj$collin_convolved[[paste0("run",toString(r))]][["r"]]) == 1))
        expect_equal(length(d_obj$collin_convolved[[paste0("run",toString(r))]][["vif"]]), nsignal_testing)
      }
    }

    # testing concat_onsets (event onset times for each signal) element of d_obj
    expect_true("concat_onsets" %in% names(d_obj))
    if ("concat_onsets" %in% names(d_obj)){
      for(r in 1:nruns_testing){
        expect_equal(length(names(d_obj$concat_onsets)), nsignal_testing)
      }
    }
 
    # testing run_volumes element of d_obj
    expect_true("run_volumes" %in% names(d_obj))
    if ("run_volumes" %in% names(d_obj)){
      for(r in 1:nruns_testing){
        expect_equal(length(d_obj$run_volumes), nruns_testing)
      }
    }

    # testing design_plot element of d_obj
    expect_true("design_plot" %in% names(d_obj))
    if ("design_plot" %in% names(d_obj)){
      expect_equal(class(d_obj$design_plot)[[2]], "ggplot")
    }
  }

})


