proj_dir <- "/proj/mnhallqlab/projects/fmri.pipeline_test_data"
test_data_dir <- "testdata"

get_trial_df <- function() {
  trial_df <- data.table::fread(file.path(test_data_dir, "sample_trial_data.csv"))
  return(trial_df)
}

get_subj_df <- function() {
  #subj_df <- data.table::fread(file.path(test_dir, "mmy3_demographics.tsv"), data.table = FALSE)
  subj_df <- data.table::fread(file.path(test_data_dir, "sample_subject_data.tsv"))
  return(subj_df)
}

#' Provide a minimal instantiation of the gpa list
#' 
#' @return a minimal gpa list
get_gpa_minimal <- function() {
  setup_glm_pipeline(
    analysis_name = "gpa_tests",
    output_directory = tempdir(),
    subject_data = data.frame(id=c(1, 2, 3)),
    trial_data = data.frame(id=c(1, 2, 3)),
    tr = 1.0,
    l1_models=NULL
  )
}

#' Provide a gpa list with subject/trial data,
#' but no built models.
#' 
#' @return a minimal gpa list
get_gpa_no_models <- function() {
  setup_glm_pipeline(
    analysis_name = "gpa_tests",
    output_directory = tempdir(),
    subject_data = get_subj_df(),
    trial_data = get_trial_df(),
    tr = 1.0,
    l1_models=NULL
  )
}

get_gpa <- function(
    scheduler = "slurm", drop_volumes = 2,
    exclude_run = "max(FD) > 5 | sum(FD > .9)/length(FD) > .10",
    exclude_subject = "n_good_runs < 4",
    truncate_run = "(FD > 0.9 & time > last_offset) | (time > last_offset + last_isi)",
    spike_volumes = NULL) {

  setup_glm_pipeline(
    analysis_name = "gpa_tests",
    scheduler = scheduler,
    subject_data = get_subj_df(proj_dir),
    trial_data = get_trial_df(proj_dir),
    output_directory = tempdir(),
    n_expected_runs = 8,
    tr = 1.0,
    drop_volumes = drop_volumes,
    l1_models=NULL, l2_models=NULL, l3_models=NULL,
    confound_settings = list(
      motion_params_file = "motion.par",
      confound_input_file = "nuisance_regressors.txt",
      confound_input_colnames = c("csf", "dcsf", "wm", "dwm"),
      l1_confound_regressors = c("csf", "dcsf", "wm", "dwm"),
      exclude_run = exclude_run,
      exclude_subject = exclude_subject,
      truncate_run = truncate_run,
      spike_volumes = NULL
    )
  )
}