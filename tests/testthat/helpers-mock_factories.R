# Mock Factory Functions for Unit Testing
# =========================================
# This file provides lightweight mock object creators for testing individual
# functions without requiring external data files or the full pipeline setup.
#
# Design principles:
# - Each factory creates the minimal valid object needed for testing
# - Factories use default arguments that can be overridden
# - Objects should validate against the same checks as real objects
# - No external file dependencies

# ==============================================================================
# GPA Mock Factories
# ==============================================================================

#' Create a minimal mock glm_pipeline_arguments (gpa) object
#'
#' @param n_subjects Number of subjects
#' @param n_runs Number of runs per subject
#' @param n_trials Number of trials per run
#' @param tr Repetition time in seconds
#' @param include_l1_models Whether to include mock L1 models
#' @param include_l2_models Whether to include mock L2 models
#' @param include_l3_models Whether to include mock L3 models
#' @return A mock gpa object with class "glm_pipeline_arguments"
create_mock_gpa <- function(
    n_subjects = 3,
    n_runs = 2,
    n_trials = 10,
    tr = 1.0,
    analysis_name = "mock_analysis",
    include_l1_models = FALSE,
    include_l2_models = FALSE,
    include_l3_models = FALSE,
    output_directory = tempdir()
) {
  # Create subject data
  subject_data <- data.frame(
    id = paste0("sub", seq_len(n_subjects)),
    session = rep(1L, n_subjects),
    age = runif(n_subjects, 18, 65),
    sex = sample(c("M", "F"), n_subjects, replace = TRUE),
    stringsAsFactors = FALSE
  )
  
  # Create run data
  run_data <- expand.grid(
    id = subject_data$id,
    run_number = seq_len(n_runs),
    stringsAsFactors = FALSE
  )
  run_data$session <- 1L
  run_data$run_volumes <- 100L
  run_data$drop_volumes <- 2L
  run_data$mr_dir <- file.path(output_directory, run_data$id, paste0("run", run_data$run_number))
  run_data$run_nifti <- file.path(run_data$mr_dir, "bold.nii.gz")
  run_data$exclude_run <- FALSE
  run_data$has_l1_setup <- FALSE
  run_data$has_l1_complete <- FALSE
  
  # Create trial data
  trial_data <- expand.grid(
    id = subject_data$id,
    run_number = seq_len(n_runs),
    trial = seq_len(n_trials),
    stringsAsFactors = FALSE
  )
  trial_data$session <- 1L
  trial_data$onset <- (trial_data$trial - 1) * 2.5 + runif(nrow(trial_data), 0, 0.5)
  trial_data$duration <- 1.0
  trial_data$event <- "stimulus"
  trial_data$rt <- runif(nrow(trial_data), 0.3, 1.5)
  trial_data$value <- rnorm(nrow(trial_data))
  
  gpa <- list(
    analysis_name = analysis_name,
    tr = tr,
    drop_volumes = 2L,
    subject_data = subject_data,
    run_data = run_data,
    trial_data = trial_data,
    n_expected_runs = n_runs,
    vm = c(id = "id", session = "session", run_number = "run_number"),
    output_directory = output_directory,
    output_locations = list(
      root = output_directory,
      feat_l1 = file.path(output_directory, "feat_l1"),
      feat_l2 = file.path(output_directory, "feat_l2"),
      feat_l3 = file.path(output_directory, "feat_l3"),
      scheduler_scripts = file.path(output_directory, "scripts")
    ),
    scheduler = "local",
    glm_software = "fsl",
    parallel = list(
      l1_setup_cores = 1L,
      l1_setup_time = "1:00:00",
      l1_setup_memgb = "8G"
    ),
    confound_settings = list(
      motion_params_file = "motion.par",
      confound_input_file = "nuisance.txt",
      l1_confound_regressors = c("csf", "wm")
    ),
    l1_models = NULL,
    l2_models = NULL,
    l3_models = NULL,
    lgr_threshold = "info"
  )
  
  if (include_l1_models) {
    gpa$l1_models <- create_mock_l1_models(n_subjects, n_runs, n_trials)
  }
  
  if (include_l2_models) {
    gpa$l2_models <- create_mock_l2_models()
  }
  
  if (include_l3_models) {
    gpa$l3_models <- create_mock_l3_models()
  }
  
  class(gpa) <- c("glm_pipeline_arguments", "list")
  return(gpa)
}

#' Create a minimal mock gpa for validation testing only
#'
#' @return A bare-minimum gpa object
create_mock_gpa_minimal <- function() {
  gpa <- list(
    analysis_name = "minimal_test",
    tr = 1.0,
    subject_data = data.frame(id = c("sub1", "sub2"), session = c(1L, 1L)),
    trial_data = data.frame(id = c("sub1", "sub2"), trial = c(1, 1)),
    l1_models = NULL,
    l2_models = NULL,
    l3_models = NULL
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")
  return(gpa)
}

# ==============================================================================
# L1 Model Mock Factories
# ==============================================================================

#' Create mock L1 models structure
#'
#' @param n_subjects Number of subjects
#' @param n_runs Number of runs per subject
#' @param n_trials Number of trials per run
#' @return A mock l1_models list structure
create_mock_l1_models <- function(n_subjects = 3, n_runs = 2, n_trials = 10) {
  # Create signals
  signals <- create_mock_signals(n_subjects, n_runs, n_trials)
  
  # Create a simple model specification
  models <- list(
    model1 = list(
      name = "model1",
      signals = c("event_signal", "parametric_signal"),
      contrasts = list(
        event_vs_baseline = list(
          name = "event_vs_baseline",
          weights = c(event_signal = 1)
        )
      )
    )
  )
  
  list(
    signals = signals,
    models = models,
    wi_factors = NULL,
    n_l1_copes = 1L
  )
}

#' Create mock signals structure
#'
#' @param n_subjects Number of subjects
#' @param n_runs Number of runs
#' @param n_trials Number of trials per run
#' @return A list of mock signal specifications
create_mock_signals <- function(n_subjects = 3, n_runs = 2, n_trials = 10) {
  # Create value data.frame for parametric signal
  value_df <- expand.grid(
    id = paste0("sub", seq_len(n_subjects)),
    session = 1L,
    run_number = seq_len(n_runs),
    trial = seq_len(n_trials),
    stringsAsFactors = FALSE
  )
  value_df$value <- rnorm(nrow(value_df))
  
  list(
    event_signal = list(
      name = "event_signal",
      event = "stimulus",
      value = 1,
      normalization = "none",
      demean_convolved = FALSE,
      add_deriv = FALSE
    ),
    parametric_signal = list(
      name = "parametric_signal",
      event = "stimulus",
      value = value_df,
      normalization = "evtmax_1",
      demean_convolved = TRUE,
      add_deriv = FALSE
    )
  )
}

# ==============================================================================
# L2/L3 Model Mock Factories
# ==============================================================================

#' Create mock L2 models structure
#'
#' @return A mock l2_models list structure
create_mock_l2_models <- function() {
  list(
    models = list(
      l2_model1 = list(
        name = "l2_model1",
        l1_model = "model1",
        l1_cope_names = "event_vs_baseline",
        regressors = list(
          intercept = list(name = "intercept", value = 1)
        )
      )
    )
  )
}

#' Create mock L3 models structure
#'
#' @return A mock l3_models list structure
create_mock_l3_models <- function() {
  list(
    models = list(
      l3_model1 = list(
        name = "l3_model1",
        l2_model = "l2_model1",
        regressors = list(
          intercept = list(name = "intercept", value = 1)
        )
      )
    )
  )
}

# ==============================================================================
# Design Matrix Mock Factories
# ==============================================================================

#' Create mock events data.frame for testing design matrix functions
#'
#' @param n_runs Number of runs
#' @param n_trials Number of trials per run
#' @param tr Repetition time
#' @param event_names Character vector of event type names
#' @return A data.frame of mock events
create_mock_events <- function(n_runs = 2, n_trials = 5, tr = 1.0, event_names = "stimulus") {
  do.call(rbind, lapply(seq_len(n_runs), function(r) {
    do.call(rbind, lapply(event_names, function(evt) {
      data.frame(
        event = evt,
        run_number = r,
        trial = seq_len(n_trials),
        onset = cumsum(runif(n_trials, min = 2 * tr, max = 5 * tr)),
        duration = rep(1, n_trials),
        stringsAsFactors = FALSE
      )
    }))
  }))
}

#' Create mock run_data for testing build_design_matrix
#'
#' @param n_runs Number of runs
#' @param run_volumes Number of volumes per run
#' @param drop_volumes Number of dropped volumes
#' @return A data.frame of mock run data
create_mock_run_data <- function(n_runs = 2, run_volumes = 50, drop_volumes = 0) {
  data.frame(
    run_number = seq_len(n_runs),
    run_volumes = rep(run_volumes, n_runs),
    drop_volumes = rep(drop_volumes, n_runs),
    stringsAsFactors = FALSE
  )
}

#' Create mock signals for build_design_matrix testing
#'
#' @param n_runs Number of runs
#' @param n_trials Number of trials per run
#' @param include_parametric Whether to include a parametric modulator
#' @return A list of mock signal specifications
create_mock_bdm_signals <- function(n_runs = 2, n_trials = 5, include_parametric = TRUE) {
  signals <- list(
    event_sig = list(
      name = "event_sig",
      event = "stimulus",
      value = 1,
      normalization = "none"
    )
  )
  
  if (include_parametric) {
    value_df <- expand.grid(
      run_number = seq_len(n_runs),
      trial = seq_len(n_trials),
      stringsAsFactors = FALSE
    )
    value_df$value <- rnorm(nrow(value_df), mean = 5, sd = 2)
    
    signals$param_sig <- list(
      name = "param_sig",
      event = "stimulus",
      value = value_df,
      normalization = "evtmax_1"
    )
  }
  
  return(signals)
}

# ==============================================================================
# Collinearity Diagnostics Mock Factories
# ==============================================================================

#' Create a mock l1_collinearity_summary object
#'
#' @param n_subjects Number of subjects
#' @param n_runs Number of runs per subject
#' @param n_regressors Number of regressors
#' @param vif_threshold VIF threshold for flagging
#' @param cor_threshold Correlation threshold for flagging
#' @param include_high_vif Whether to include high VIF values
#' @param include_high_cor Whether to include high correlations
#' @return A mock l1_collinearity_summary object
create_mock_collinearity_summary <- function(
    n_subjects = 3,
    n_runs = 2,
    n_regressors = 4,
    vif_threshold = 5,
    cor_threshold = 0.8,
    include_high_vif = TRUE,
    include_high_cor = TRUE) {
  
  regressors <- paste0("reg", seq_len(n_regressors))
  
  # Build VIF summary data
  vif_records <- list()
  for (subj in seq_len(n_subjects)) {
    for (run in seq_len(n_runs)) {
      vifs <- runif(n_regressors, min = 1.1, max = 3.5)
      if (include_high_vif && subj == 1 && run == 1) {
        vifs[1] <- vif_threshold + 2
      }
      
      vif_df <- data.frame(
        id = paste0("sub", subj),
        session = 1,
        run_number = run,
        l1_model = "test_model",
        regressor = regressors,
        vif = vifs,
        stringsAsFactors = FALSE
      )
      vif_records[[length(vif_records) + 1]] <- vif_df
    }
  }
  vif_summary <- do.call(rbind, vif_records)
  vif_summary$high_vif <- vif_summary$vif > vif_threshold
  
  # Build correlation summary data
  cor_records <- list()
  if (n_regressors >= 2) {
    pairs <- combn(regressors, 2)
    for (subj in seq_len(n_subjects)) {
      for (run in seq_len(n_runs)) {
        for (p in seq_len(ncol(pairs))) {
          cor_val <- runif(1, min = -0.5, max = 0.5)
          if (include_high_cor && subj == 1 && run == 1 && p == 1) {
            cor_val <- cor_threshold + 0.1
          }
          
          cor_df <- data.frame(
            id = paste0("sub", subj),
            session = 1,
            run_number = run,
            l1_model = "test_model",
            regressor1 = pairs[1, p],
            regressor2 = pairs[2, p],
            correlation = cor_val,
            stringsAsFactors = FALSE
          )
          cor_records[[length(cor_records) + 1]] <- cor_df
        }
      }
    }
  }
  
  correlation_summary <- if (length(cor_records) > 0) {
    do.call(rbind, cor_records)
  } else {
    data.frame()
  }
  
  if (nrow(correlation_summary) > 0) {
    correlation_summary$abs_correlation <- abs(correlation_summary$correlation)
    correlation_summary$high_correlation <- correlation_summary$abs_correlation > cor_threshold
  }
  
  high_vif <- vif_summary[vif_summary$high_vif, , drop = FALSE]
  high_correlation <- if (nrow(correlation_summary) > 0) {
    correlation_summary[correlation_summary$high_correlation, , drop = FALSE]
  } else {
    data.frame()
  }
  
  summary_stats <- list(
    n_subjects = n_subjects,
    n_runs = n_subjects * n_runs,
    n_models = 1,
    n_regressors = n_regressors,
    vif_mean = mean(vif_summary$vif, na.rm = TRUE),
    vif_max = max(vif_summary$vif, na.rm = TRUE),
    vif_n_flagged = sum(vif_summary$high_vif, na.rm = TRUE),
    cor_mean_abs = if (nrow(correlation_summary) > 0) mean(correlation_summary$abs_correlation, na.rm = TRUE) else NA,
    cor_max_abs = if (nrow(correlation_summary) > 0) max(correlation_summary$abs_correlation, na.rm = TRUE) else NA,
    cor_n_flagged = if (nrow(correlation_summary) > 0) sum(correlation_summary$high_correlation, na.rm = TRUE) else 0
  )
  
  result <- list(
    vif_summary = vif_summary,
    correlation_summary = correlation_summary,
    high_vif = high_vif,
    high_correlation = high_correlation,
    summary_stats = summary_stats,
    thresholds = list(vif = vif_threshold, correlation = cor_threshold)
  )
  
  class(result) <- c("l1_collinearity_summary", "list")
  return(result)
}

#' Create an empty collinearity summary for edge case testing
#'
#' @param vif_threshold VIF threshold
#' @param cor_threshold Correlation threshold
#' @return An empty mock l1_collinearity_summary object
create_empty_collinearity_summary <- function(vif_threshold = 5, cor_threshold = 0.8) {
  result <- list(
    vif_summary = data.frame(),
    correlation_summary = data.frame(),
    high_vif = data.frame(),
    high_correlation = data.frame(),
    summary_stats = list(
      n_subjects = 0, n_runs = 0, n_models = 0, n_regressors = 0,
      vif_mean = NA, vif_max = NA, vif_n_flagged = 0,
      cor_mean_abs = NA, cor_max_abs = NA, cor_n_flagged = 0
    ),
    thresholds = list(vif = vif_threshold, correlation = cor_threshold)
  )
  class(result) <- c("l1_collinearity_summary", "list")
  return(result)
}

# ==============================================================================
# Batch Job Mock Factories
# ==============================================================================

#' Create a mock R_batch_job for testing
#'
#' @param job_name Name of the job
#' @param scheduler Scheduler type ("local", "slurm", "torque")
#' @return A mock job specification list
create_mock_batch_job <- function(job_name = "test_job", scheduler = "local") {
  list(
    job_name = job_name,
    scheduler = scheduler,
    n_cpus = 1L,
    mem_total = "4G",
    wall_time = "1:00:00",
    r_code = "print('test')",
    submitted = FALSE,
    job_id = NULL
  )
}

# ==============================================================================
# Confound Mock Factories
# ==============================================================================

#' Create mock confound data for testing
#'
#' @param n_volumes Number of volumes
#' @param n_confounds Number of confound regressors
#' @return A data.frame of mock confounds
create_mock_confounds <- function(n_volumes = 100, n_confounds = 4) {
  confound_names <- c("csf", "dcsf", "wm", "dwm", "motion1", "motion2")
  n_confounds <- min(n_confounds, length(confound_names))
  
  confounds <- as.data.frame(
    matrix(rnorm(n_volumes * n_confounds), 
           nrow = n_volumes, 
           ncol = n_confounds)
  )
  names(confounds) <- confound_names[seq_len(n_confounds)]
  
  # Add framewise displacement
  confounds$framewise_displacement <- abs(rnorm(n_volumes, mean = 0.2, sd = 0.15))
  
  return(confounds)
}

#' Create mock motion parameters
#'
#' @param n_volumes Number of volumes
#' @return A matrix of mock motion parameters (6 columns)
create_mock_motion_params <- function(n_volumes = 100) {
  motion <- matrix(
    rnorm(n_volumes * 6, mean = 0, sd = 0.5),
    nrow = n_volumes,
    ncol = 6
  )
  colnames(motion) <- c("trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z")
  return(motion)
}
