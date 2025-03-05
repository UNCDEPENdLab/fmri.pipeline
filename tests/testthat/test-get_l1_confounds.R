#' sub-functions to test get_l1_confounds
finalize_pipeline_partial <- function(gpa, refinalize = FALSE) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_logical(refinalize, len = 1L)
  lg <- lgr::get_logger("glm_pipeline/setup_glm_pipeline")
  lg$set_threshold(gpa$lgr_threshold)

  if (isTRUE(gpa$finalize_complete) && isFALSE(refinalize)) {
    lg$debug("In finalize_pipeline_configuration, finalization of gpa already complete. Returning object unchanged.")
    return(gpa)
  }

  if (is.null(gpa$sqlite_con) || !DBI::dbIsValid(gpa$sqlite_con)) {
    lg$info("Opening SQLite connection to: %s", gpa$output_locations$sqlite_db)
    gpa$sqlite_con <- DBI::dbConnect(RSQLite::SQLite(), gpa$output_locations$sqlite_db)
  }

  # final checks on compute environment now that we're running inside the compute environment
  test_compute_environment(gpa, stop_on_fail=TRUE)

  # l1 models must be specified to get started (hard enforcement)
  if (!checkmate::test_class(gpa$l1_models, "l1_model_set")) {
    msg <- "Could not find valid $l1_models specification in gpa. Be sure to run build_l1_models() before proceeding!"
    lg$error(msg)
    stop(msg)
  }

  # l2 models are optional and depend on whether it's a multi-run setup
  if (isTRUE(gpa$multi_run) && !checkmate::test_class(gpa$l2_models, "hi_model_set")) {
    msg <- "Could not find valid $l2_models specification in gpa. We will proceed, but weird things may happen! I suggest running build_l2_models()."
    lg$warn(msg)
    warning(msg)
  }

  # l3 models are optional but almost always should be in place
  if (!checkmate::test_class(gpa$l3_models, "hi_model_set")) {
    msg <- "Could not find valid $l3_models specification in gpa. We will proceed, but weird things may happen! I suggest running build_l3_models()."
    lg$warn(msg)
    warning(msg)
  }

  # new approach: use internal model names for creating output directories at subject level
  # default to <analysis_name>/<l1_model>
  # add suffix if using preconvolution approach
  gpa$l1_models$models <- lapply(gpa$l1_models$models, function(mm) {
    mm$outdir <- file.path(gpa$analysis_name, paste0(mm$name, ifelse(gpa$use_preconvolve, "_preconvolve", "")))
    return(mm)
  })

  if (!is.null(gpa$run_number_regex)) {
    if (stringr::str_count(gpa$run_number_regex, stringr::fixed("(")) != 1L) {
      stop(
        "run_number_regex: ", gpa$run_number_regex,
        " must have exactly one opening parenthesis, denoting start of run number capture"
      )
    }

    if (stringr::str_count(gpa$run_number_regex, stringr::fixed(")")) != 1L) {
      stop(
        "run_number_regex: ", gpa$run_number_regex,
        " must have exactly one closing parenthesis, denoting end of run number capture"
      )
    }
  }

  # setup l1 copes, cope names, and contrasts.
  gpa$l1_cope_names <- lapply(gpa$l1_models$models, function(mm) {
    rownames(mm$contrasts)
  }) # names of level 1 copes for each model


  if (is.null(gpa$center_l3_predictors)) gpa$center_l3_predictors <- TRUE
  if (is.null(gpa$bad_ids)) gpa$bad_ids <- c()
  if (is.null(gpa$scheduler)) gpa$scheduler <- "slurm" # HPC batch system
  gpa$scheduler <- tolower(gpa$scheduler)
  checkmate::assert_subset(gpa$scheduler, c("slurm", "sbatch", "torque", "qsub", "local", "sh"))

  if (is.null(gpa$zthresh)) gpa$zthresh <- 3.09 # 1-tailed p=.001 for z stat
  if (is.null(gpa$clustsize)) gpa$clustsize <- 50 # arbitrary reasonable lower bound on cluster size
  if (is.null(gpa$glm_software)) gpa$glm_software <- "fsl" # default to FSL FEAT

  if (is.null(gpa$log_json)) gpa$log_json <- TRUE # whether to write JSON log files
  if (is.null(gpa$log_txt)) gpa$log_txt <- TRUE # whether to write text log files

  if (is.null(gpa$l1_setup_log)) {
    l1_setup_log <- paste0(names(gpa$l1_models$models), "_l1setup") %>% setNames(names(gpa$l1_models$models))
    lg$debug("l1_setup_log is %s", l1_setup_log)
    gpa$l1_setup_log <- l1_setup_log
  }
  if (is.null(gpa$l1_execution_log)) {
    l1_execution_log <- paste0(names(gpa$l1_models$models), "_l1execution") %>% setNames(names(gpa$l1_models$models))
    lg$debug("l1_execution_log is %s", l1_execution_log)
    gpa$l1_execution_log <- l1_execution_log
  }

  if (is.null(gpa$n_expected_runs)) gpa$n_expected_runs <- 1 # assume single run case

  # remove bad ids before running anything further
  if (!is.null(gpa$bad_ids) && length(gpa$bad_ids) > 0L) {
    lg$info("Removing the following IDs from data structure before beginning analysis: %s", paste(gpa$bad_ids, collapse = ", "))
    gpa$subject_data <- gpa$subject_data %>% dplyr::filter(!id %in% gpa$bad_ids) # remove bad ids
    gpa$run_data <- gpa$run_data %>% dplyr::filter(!id %in% gpa$bad_ids) # remove bad ids
    gpa$trial_data <- gpa$trial_data %>% dplyr::filter(!id %in% gpa$bad_ids) # remove bad ids
  }

  # build design matrix default arguments
  if (is.null(gpa$additional$bdm_args)) {
    gpa$additional$bdm_args <- list(
      baseline_coef_order = 2, center_values = TRUE,
      plot = FALSE, convolve_wi_run = TRUE, output_directory = "run_timing"
    )
  }

  # default settings for feat l1
  if (is.null(gpa$additional$feat_l1_args)) {
    gpa$additional$feat_l1_args <- list(z_thresh = 1.96, prob_thresh = .05, paradigm_hp = 120)
  }

  # identify and validate niftis for each run
  gpa <- lookup_nifti_inputs(gpa)

  #####
  # handle GLM settings and defaults
  if (is.null(gpa$glm_settings) || gpa$glm_settings[1L] == "default") {
    lg$info("Using default settings for GLM implementation")
    gpa$glm_settings <- list(
      fsl = list(),
      afni = list(),
      spm = list()
    )
  }

  fsl_defaults <- list(
    force_l1_creation = FALSE, # whether to overwrite existing level 1 setup files (e.g., .fsf)
    failed_l1_folder_action = "delete", # whether to 'delete', 'archive', or 'skip' failed folders before running feat l1 fsf jobs
    incomplete_l1_folder_action = "delete", # whether to 'delete', 'archive', or 'skip' incomplete folders before running feat l1 fsf jobs
    force_l2_creation = FALSE, # whether to overwrite existing level 2 setup files (e.g., .fsf)
    failed_l2_folder_action = "delete", # whether to 'delete', 'archive', or 'skip' failed folders before running feat l1 fsf jobs
    incomplete_l2_folder_action = "delete", # whether to 'delete', 'archive', or 'skip' incomplete folders before running feat l1 fsf jobs
    force_l3_creation = FALSE, # whether to overwrite existing level 3 setup files (e.g., .fsf)
    failed_l3_folder_action = "delete", # whether to 'delete', 'archive', or 'skip' failed folders before running feat l1 fsf jobs
    incomplete_l3_folder_action = "archive", # whether to 'delete', 'archive', or 'skip' incomplete folders before running feat l1 fsf jobs
    replace_l1_nifti_symlink = TRUE # whether to replace filtered_func_data with symlink to the (same) input data
  )

  gpa$glm_settings$fsl <- populate_defaults(gpa$glm_settings$fsl, fsl_defaults)

  output <- list()
  output$gpa <- gpa
  output$lg <- lg

  return(output)
}

finalize_confound_settings_partial <- function(gpa, lg) {
  checkmate::assert_class(lg, "Logger")
  checkmate::assert_class(gpa, "glm_pipeline_arguments")

  # validate confound settings
  confound_defaults <- list(
    motion_params_file = NULL,
    motion_params_colnames = NULL,
    confound_input_file = "confounds.tsv",
    l1_confound_regressors = NULL, # column names in motion_params_file and/or confound_input_file
    exclude_run = "mean(FD) > 0.5 | max(FD) > 6",
    truncate_run = NULL, # example: FD > 1 & time > last_onset
    spike_volumes = "FD > 0.9",
    na_strings = getOption("datatable.na.strings", "NA") # default na.strings argument for data.table::fread calls
  )

  if (is.null(gpa$confound_settings)) {
    lg$info("Using default settings for confounds and exclusions")
    lg$info("Look for confounds in confounds.tsv")
    lg$info("Exclude run if mean(FD) > 0.5 or max(FD) > 6")
  } else {
    checkmate::assert_string(gpa$confound_settings$exclude_run, null.ok = TRUE)
    checkmate::assert_string(gpa$confound_settings$exclude_subject, null.ok = TRUE)
  }

  gpa$confound_settings <- populate_defaults(gpa$confound_settings, confound_defaults)

  #####
  # populate confounds
  if ("motion_params_file" %in% names(gpa$run_data)) {
    lg$info("motion_params_file column already in run_data. Not using motion_params_file confounds specification.")
  } else if (!is.null(gpa$confound_settings$motion_params_file)) {
    checkmate::assert_string(gpa$confound_settings$motion_params_file)
    gpa$run_data$motion_params_file <- gpa$confound_settings$motion_params_file # this gets expanded by get_mr_abspath in get_l1_confounds
  }

  if (!"motion_params_file" %in% names(gpa$run_data)) {
    gpa$run_data$motion_params_file <- gpa$run_data$motion_params_present <- NA_character_
  } else {
    gpa$run_data$motion_params_present <- file.exists(get_mr_abspath(gpa$run_data, "motion_params_file"))
  }

  if ("confound_input_file" %in% names(gpa$run_data)) {
    lg$info("confound_input_file column already in run_data. Not using confound_input_file confounds specification.")
  } else if (!is.null(gpa$confound_settings$confound_input_file)) {
    checkmate::assert_string(gpa$confound_settings$confound_input_file)
    gpa$run_data$confound_input_file <- gpa$confound_settings$confound_input_file # this gets expanded by get_mr_abspath in get_l1_confounds
  }

  # determine whether confound input files are present
  if (!"confound_input_file" %in% names(gpa$run_data)) {
    gpa$run_data$confound_input_file <- gpa$run_data$confound_input_file_present <- NA_character_
  } else {
    gpa$run_data$confound_input_file_present <- file.exists(get_mr_abspath(gpa$run_data, "confound_input_file"))
  }

  rhs_to_vars <- function(str) {
    if (is.null(str)) {
      NULL # return NULL if input is NULL
    } else if (checkmate::test_string(str)) {
      if (!grepl("^\\s*~", str)) str <- paste("~", str)
      all.vars(as.formula(paste("~", str)))
    } else if (checkmate::test_formula(str)) {
      all.vars(str)
    } else {
      stop("rhs_to_vars input is not a string or formula")
    }
  }

  # figure out all confound columns that will be used in the pipeline
  gpa$confound_settings$run_exclusion_columns <- rhs_to_vars(gpa$confound_settings$exclude_run)

  # figure out all confound columns that will be used in the pipeline
  gpa$confound_settings$run_truncation_columns <- rhs_to_vars(gpa$confound_settings$truncate_run)

  # TODO: Should this become 'id_exclusion_columns' and should we support session versus subject exclusion
  # (E.g., in longitudinal analysis)
  gpa$confound_settings$subject_exclusion_columns <- rhs_to_vars(gpa$confound_settings$exclude_subject)

  gpa$confound_settings$all_confound_columns <- unique(c(
    gpa$confound_settings$l1_confound_regressors,
    gpa$confound_settings$run_exclusion_columns,
    gpa$confound_settings$run_truncation_columns,
    gpa$confound_settings$subject_exclusion_columns
  ))

  # handle lookup and creation of all confound files

  # numeric row number of each input to aid in tracking
  gpa$run_data$input_number <- seq_len(nrow(gpa$run_data))

  # populate confounds in SQLite database and calculate run exclusions
  # TODO: allow external $exclude_run from user, add internal calculated $calc_exclude run

  # confound_info <- parallel::mclapply(seq_len(nrow(gpa$run_data)), function(ii) {
  #   # this should add rows to the SQLite data for a subject if not yet present, or just return those rows if they exist
  #   l1_info <- get_l1_confounds(
  #     id = gpa$run_data$id[ii], session = gpa$run_data$session[ii], run_number = gpa$run_data$run_number[ii],
  #     gpa = gpa, drop_volumes = gpa$drop_volumes
  #   )[c("l1_confound_file", "exclude_run")]
  #   return(l1_info)
  # }, mc.cores = gpa$parallel$finalize_cores)

  # add onset + offset + isi data to run_data so that get_l1_confounds can handle run truncation appropriately
  gpa <- populate_last_events(gpa, lg)

  return(gpa)
}


#' test get_l1_confounds function
test_that("test get_l1_confounds", {

  # gpa <- create_gpa() # 2: Cell and 1 for others
  gpa <- create_gpa(
      analysis_name = "gpa_tests_build_l1_models",
      test_data_base_dir = "/proj/mnhallqlab/projects/fmri.pipeline_test_data",
      l1_spec_file = "sample_2_L1_spec.yaml",
      trial_data_file = "sample_trial_data.csv.gz",
      run_data_file = "sample_run_data.csv",
      subject_data_file = "sample_subject_data.csv",
      gpa_cache_file = "gpa_tests_build_l1_models.rds",
      cache_file = "gpa_tests_base_build_l1_models.rds") 
  # gpa <- readRDS("/proj/mnhallqlab/projects/fmri.pipeline_test_data/gpa_l2l3.rds")
  
  output <- finalize_pipeline_partial(gpa)

  gpa <- finalize_confound_settings_partial(output$gpa, output$lg)

  # for each run, calculate confounds, exclusions, and truncation
  get_l1_confounds_output <- lapply(seq_len(nrow(gpa$run_data)), function(ii) {
    get_l1_confounds(run_df = gpa$run_data[ii, , drop = FALSE], gpa = gpa)
  }) # get_l1_confounds_output should contain a list of list 'confounds' of locations of confound text files

  # save.image("test-get_l1_confound_runned_l2l3_local.RData")
  
  # test if confound text files were created in file.path(analysis_outdir, paste0("run", run_number, "_l1_confounds.txt"))
  # for each subject for each run and that they are not empty. 
  # I dont know how the nrows in this confound file is calculated, 4 columns should exists.
  id_list <- unique(gpa$run_data$id)
  run_list <- c(1:8)
  lapply(id_list, function(id) {
    analysis_outdir <- paste0("/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/gpa_tests_l2l3_local/feat_l1/sub-", toString(id))
    lapply(run_list, function(run) {
      confound_filepath <- file.path(analysis_outdir, paste0("run", run, "_l1_confounds.txt"))
      # check if file exists
      expect_true(file.exists(confound_filepath))
      # check if file has 4 columns
      expect_equal(ncol(read.table(confound_filepath)), 4)
      # check if file has more than 0 rows
      expect_gt(nrow(read.table(confound_filepath)), 0)
      # check for missing elements
      expect_false(any(is.na(read.table(confound_filepath))))
      # check if the elements are numeric
      expect_true(all(sapply(read.table(confound_filepath), is.numeric)))
    })
  })
  
  # testing if get_l1_confounds added dataframes to the sqlite database

  # testing if sqlite database has the motion parameters dataframe
  lapply(id_list, function(id) {
    lapply(run_list, function(run) {
      saved_df <- read_df_sqlite(gpa, id = id, session = 1, run_number = run, table = "l1_motion_parameters")
      # check if the dataframe is not empty
      expect_gt(nrow(saved_df), 0)
      # check for missing elements
      expect_false(any(is.na(saved_df)))
      # check if the elements are numeric
      expect_true(all(sapply(saved_df, is.numeric)))
    })
  })

  # testing if sqlite database has the truncation data dataframe
  lapply(id_list, function(id) {
    lapply(run_list, function(run) {
      saved_df <- read_df_sqlite(gpa, id = id, session = 1, run_number = run, table = "l1_truncation_data")
      # check if the dataframe is not empty
      expect_gt(nrow(saved_df), 0)
      # check for missing elements
      expect_false(any(is.na(saved_df)))
      # check if the elements are numeric
      expect_true(all(sapply(saved_df, is.numeric)))
    })
  })

  # testing if sqlite database has the exclusion data dataframe
  lapply(id_list, function(id) {
    lapply(run_list, function(run) {
      saved_df <- read_df_sqlite(gpa, id = id, session = 1, run_number = run, table = "l1_exclusion_data")
      # check if the dataframe is not empty
      expect_gt(nrow(saved_df), 0)
      # check for missing elements
      expect_false(any(is.na(saved_df)))
      # check if the elements are numeric
      expect_true(all(sapply(saved_df, is.numeric)))
    })
  })

  # saved_df <- read_df_sqlite(gpa, id = 10638, session = 1, run_number = 1, table = "l1_confound_inputs")
  # saved_df <- read_df_sqlite(gpa, id = 10638, session = 1, run_number = 1, table = "l1_confounds") # this is what is saved in the confounds text file
  # saved_df <- read_df_sqlite(gpa, id = 10638, session = 1, run_number = 1, table = "l1_run_calculations") # nifti location

  # NEXT TODO screw up the data (how?) and test again




})
