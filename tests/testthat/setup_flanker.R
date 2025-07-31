# setwd("/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/fmri.pipeline")
# devtools::load_all()

#' @title Build the trial dataframe from test data. This file is already saved in the repo, so should only be used to adjust the run data file.
#' @author Nidhi Desai
#' @return trial_df
build_flanker_trial_data_file <- function(fmriprep_path = "/proj/mnhallqlab/no_backup/flanker-fmriprep", task_name = "flanker", save_filepath = "/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/fmri.pipeline/tests/testthat/testdata/flankers_trial_data.csv"){
  trial_df <- generate_trial_data_from_bids(fmriprep_path, task_name = task_name) %>%
    dplyr::rename(id = subject, run_number = run) %>%
    mutate(run_number=as.integer(run_number))

  trial_df <- trial_df %>%
    group_by(id, run_number) %>%
    arrange(onset, .by_group = TRUE) %>%
    mutate(trial = 1:n()) %>%
    ungroup() %>%
    filter(correctness=="correct") %>% # only a tiny number of incorrect trials, drop for now
    mutate(trial_type = stringr::str_extract(trial_type, "[^_]+"))

  # to test within-subject factors, create a second fake factor for which finger was used
  trial_df <- trial_df %>%
      group_by(id, run_number, StimVar) %>%
      mutate(wifake = sample(c("left", "right"), n(), replace = TRUE)) %>%
      ungroup()  

  trial_df %>%
    group_by(StimVar) %>%
    summarize(mean(response_time), sd(response_time))

  write.csv(trial_df, save_filepath, row.names = FALSE)

  return(trial_df)
}
# trial_df <- build_flanker_trial_data_file(fmriprep_path = "/proj/mnhallqlab/no_backup/flanker-fmriprep", task_name = "flanker", save_filepath = "/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/fmri.pipeline/tests/testthat/testdata/flankers_trial_data.csv")

#' @title Build the run dataframe from test data. This file is already saved in the repo, so should only be used to adjust the run data file.
#' @author Nidhi Desai
#' @return run_df
build_flanker_run_data_file <- function(fmriprep_path = "/proj/mnhallqlab/no_backup/flanker-fmriprep", task_name = "flanker", desc = "postproc", space = "MNI152NLin2009cAsym", save_filepath = "/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/fmri.pipeline/tests/testthat/testdata/flankers_run_data.csv"){
  run_df <- generate_run_data_from_bids(fmriprep_path, task_name = task_name, desc = desc, space = space)
  write.csv(run_df, save_filepath, row.names = FALSE)
  return(run_df)
}
# run_df <- build_flanker_run_data_file(fmriprep_path = "/proj/mnhallqlab/no_backup/flanker-fmriprep", task_name = "flanker", desc = "postproc", space = "MNI152NLin2009cAsym", save_filepath = "/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/fmri.pipeline/tests/testthat/testdata/flankers_run_data.csv")

#' @title Build the subject dataframe from test data. This file is already saved in the repo, so should only be used to adjust the subject data file.
#' @author Nidhi Desai
#' @return subject_df
build_flanker_subject_data_file <- function(fmriprep_path = "/proj/mnhallqlab/no_backup/flanker-fmriprep", save_filepath = "/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/fmri.pipeline/tests/testthat/testdata/flankers_subject_data.csv"){
  subject_df <- generate_subject_data_from_bids(fmriprep_path) 
  write.csv(subject_df, save_filepath, row.names = FALSE)
  return(subject_df)
}
# subject_df <- build_flanker_subject_data_file(fmriprep_path = "/proj/mnhallqlab/no_backup/flanker-fmriprep", save_filepath = "/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/fmri.pipeline/tests/testthat/testdata/flankers_subject_data.csv")


#' Build a gpa object with all information that doesn't require the CLI, and export to RDS.
#' 
#' @param test_data_base_dir the base directory of the test data.
build_flanker_gpa_base <- function(
    analysis_name = "gpa_flanker_tests",
    test_data_base_dir = "/proj/mnhallqlab/no_backup/flanker-fmriprep",
    trial_data_file = "/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/fmri.pipeline/tests/testthat/testdata/flankers_trial_data.csv",
    run_data_file = "/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/fmri.pipeline/tests/testthat/testdata/flankers_run_data.csv",
    subject_data_file = "/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/fmri.pipeline/tests/testthat/testdata/flankers_subject_data.csv",
    cache_file ="/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/fmri.pipeline/tests/testthat/testdata/gpa_base.rds",
    scheduler = "slurm",
    n_expected_runs = 2,
    tr = 2.0) {

  # Read in trial, run, and subject dataframes
  trial_df <- read.csv(trial_data_file)
  run_df <- read.csv(run_data_file)
  subj_df <- read.csv(subject_data_file)

  gpa <- setup_glm_pipeline(
    analysis_name = analysis_name,
    scheduler = scheduler,
    trial_data = trial_df,
    run_data = run_df,
    subject_data = subj_df, # output_directory = tempdir(), # vm = c(id = "id", run_number = "run"),
    n_expected_runs = n_expected_runs,
    tr = tr,
    confound_settings = list(
      exclude_run = "mean(framewise_displacement) > 0.5 | max(framewise_displacement) > 6",
      truncate_run = NULL,
      exclude_subject = NULL,
      spike_volumes = "framewise_displacement > 0.9"
    ),
    l1_models=NULL, l2_models=NULL, l3_models=NULL
  )

  saveRDS(gpa, file = cache_file)
  return(gpa)
}

#' Populate base GPA with models. Save models & final GPA to cache files.
#' Currently interactively-only until we can build all models from spec YAML,
#' so must be called manually, before get_gpa to rebuild models.
#' This COULD be split up so that each model gets their own cache file, but that feels
#' excessive for now. Just assume all models need to be rebuilt if gpa needs rebuilding.
get_flanker_gpa <- function( 
  analysis_name = "gpa_flanker_tests",
  test_data_base_dir = "/proj/mnhallqlab/no_backup/flanker-fmriprep",
  gpa_cache_file = "/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/fmri.pipeline/tests/testthat/testdata/gpa.rds",
  spec_file = "/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/fmri.pipeline/tests/testthat/testdata/int.yaml",
  ...
) {
  # run set_gm_pipeline()
  gpa <- build_flanker_gpa_base(analysis_name = analysis_name, test_data_base_dir = test_data_base_dir, ...)

  # setup compute enviroment
  gpa <- setup_compute_environment(gpa, preselect_action = 5)

  # Build L1 models
  gpa <- build_l1_models(gpa, from_spec_file = spec_file)

  # Build L2 models
  gpa <- build_l2_models(gpa, from_spec_file = spec_file)

  # Build L3 models
  gpa <- build_l3_models(gpa, from_spec_file = spec_file)

  # Save final RDS object
  saveRDS(gpa, file = gpa_cache_file)
  return(gpa)
}

#' Create a populated GPA object using all the inputs
create_flanker_gpa <- function(
  analysis_name = "gpa_flanker_tests",
  test_data_base_dir = "/proj/mnhallqlab/no_backup/flanker-fmriprep",
  spec_file = "/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/fmri.pipeline/tests/testthat/testdata/int.yaml",
  trial_data_file = "/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/fmri.pipeline/tests/testthat/testdata/flankers_trial_data.csv",
  run_data_file = "/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/fmri.pipeline/tests/testthat/testdata/flankers_run_data.csv",
  subject_data_file = "/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/fmri.pipeline/tests/testthat/testdata/flankers_subject_data.csv",
  gpa_cache_file = "/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/fmri.pipeline/tests/testthat/testdata/gpa.rds",
  cache_file = "/proj/mnhallqlab/users/nidhi/fmri.pipeline_testing/fmri.pipeline/tests/testthat/testdata/gpa_base.rds"
) {
  gpa <- get_flanker_gpa(
    analysis_name = analysis_name,
    test_data_base_dir = test_data_base_dir,
    gpa_cache_file = gpa_cache_file,
    spec_file = spec_file,
    
    trial_data_file = trial_data_file,
    run_data_file = run_data_file,
    subject_data_file = subject_data_file,

    cache_file = cache_file,
    scheduler = "slurm",
    n_expected_runs = 2,
    tr = 2.0
  )
  # TODO check if we need to run get_gpa each time or we could add a if statment to not run this function if the gpa object already exists in a saved file?
  return(gpa)
}
