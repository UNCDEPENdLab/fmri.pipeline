library(data.table)
library(dplyr)
library(fmri.pipeline)

#' Build the trial dataframe from test data. This file is already saved in the repo, so should only be used to
#' adjust the run data file.
#' Currently this is just copying the trial data from the MRI test folder - do not have script used to build it here yet.
build_trial_data_file <- function(test_data_base_dir = "tests/testthat/testdata", mri_data_folder = "mri", cache_file = "mmy3_trial_df_selective_groupfixed.csv.gz", trial_data_file = "sample_trial_data.csv.gz") {
  # Copy the cached trial data file in the MRI folder to the test data folder
  file.copy(file.path(test_data_base_dir, mri_data_folder, cache_file), file.path(test_data_base_dir, trial_data_file))
}

#' Build the run dataframe from test data. This file is already saved in the repo, so should only be used to
#' adjust the run data file.
build_run_data_file <- function(test_data_base_dir = "tests/testthat/testdata", mri_data_folder = "mri", run_data_file = "sample_run_data.csv") {
    # Initialize empty lists to store data
    id_list <- vector("character")
    task_name_list <- vector("character")
    mr_dir_list <- vector("character")
    run_number_list <- vector("integer")
    run_nifti_list <- vector("character")
    confound_input_file_list <- vector("character")

    # Iterate over folders in base directory
    subject_folders <- list.dirs(file.path(test_data_base_dir, mri_data_folder), recursive = FALSE)
    for (subject_folder in subject_folders) {
        subject_dir_name <- basename(subject_folder)

        # Get subject id from folder name
        parts <- strsplit(subject_dir_name, "_")[[1]]
        id <- parts[1]  # Subject ID

        # Change level to mni_5mm_aroma
        subject_folder <- file.path(subject_folder, "mni_5mm_aroma")
        
        task_name <- "clock"

        run_folders <- list.dirs(subject_folder, recursive = FALSE)
        for (run_folder in run_folders) {
            run_dir_name <- basename(run_folder)

            run_number <- as.integer(gsub("clock", "", run_dir_name))  
            # mr_dir <- file.path(subject_dir_name, paste0("clock", run_number), "ica_aroma")
            mr_dir <- file.path(subject_folder, paste0("clock", run_number)) #, "ica_aroma")
            run_nifti <- file.path(mr_dir, "ica_aroma", "melodic_IC_thr_MNI2mm.nii")
            confound_input_file <- file.path(mr_dir, "nuisance_regressors.txt")

            # Appending all once per run
            id_list <- append(id_list, id)
            task_name_list <- append(task_name_list, task_name)
            mr_dir_list <- append(mr_dir_list, mr_dir)
            run_number_list <- append(run_number_list, run_number)
            run_nifti_list <- append(run_nifti_list, run_nifti)
            confound_input_file_list <- append(confound_input_file_list, confound_input_file)
        }
    }

    # Create a DataFrame
    df <- data.frame(id = id_list, task_name = task_name_list, 
                    mr_dir = mr_dir_list, run_number = run_number_list, 
                    run_nifti = run_nifti_list, confound_input_file = confound_input_file_list)

    # write.csv(run_df, file.path(test_data_base_dir, run_data_file), row.names = FALSE)
    write.csv(df, file.path(test_data_base_dir, run_data_file), row.names = FALSE)
}

# The following changes made by NVD in Jan/Feb 2025
# build_run_data_file(test_data_base_dir = "/proj/mnhallqlab/projects/fmri.pipeline_test_data", mri_data_folder = "")
# add emotion as a factor to the run data so we can add it as a l2 covariate
# test_data_base_dir <- "/proj/mnhallqlab/projects/fmri.pipeline_test_data/"
# trial_data_file <- "sample_trial_data.csv.gz"
# trial_df <- read.csv(file.path(test_data_base_dir, trial_data_file))
# # rt_csv column in currently in ms, converting to s
# trial_df <- trial_df %>% mutate(rt_csv = rt_csv / 1000)
# # save the updated trial data
# write.csv(trial_df, file.path(test_data_base_dir, trial_data_file), row.names = FALSE)
# trial_df <- trial_df %>% mutate(run_number = run)
# run_data_file <- "sample_run_data.csv"
# run_df <- read.csv(file.path(test_data_base_dir, run_data_file))
# reduce duplicate rows for same run number
# run_df <- run_df %>% distinct(id, run_number, .keep_all = TRUE)
# change nifti filepath to the current file named nfaswuktm_clock<num>_5.nii
# run_df <- run_df %>% rowwise() %>% mutate(run_nifti = gsub("ica_aroma/melodic_IC_thr_MNI2mm.nii", paste0("nfaswuktm_clock", toString(run_number), "_5.nii"), run_nifti)) %>% ungroup()
# run_df <- run_df %>% left_join(trial_df %>% select(id, run_number, emotion), by = c("id", "run_number"))
# write.csv(run_df, file.path(test_data_base_dir, run_data_file), row.names = FALSE)

# Correc the wrong trial number sequence
# test_data_base_dir <- "/proj/mnhallqlab/projects/fmri.pipeline_test_data/"
# trial_data_file <- "sample_trial_data.csv.gz"
# trial_df <- read.csv(file.path(test_data_base_dir, trial_data_file))
# # for run 2 trial number starts from 51, it should start from 1 again. Apply this to all run numbers
# trial_df <- trial_df %>% group_by(run) %>% mutate(trial = trial - min(trial) + 1) %>% ungroup()
# trial_df <- trial_df %>% rename(run_number = run) # rename column run to run_number
# # save the updated trial data
# write.csv(trial_df, file.path(test_data_base_dir, trial_data_file), row.names = FALSE)

#' Build the subject dataframe from test data. This file is already saved in the repo, so should only be used to
#' adjust the run data file.
#' Currently this is just copying the subject data from the MRI test folder - do not have script used to build it here yet.
build_subject_data_file <- function(test_data_base_dir = "tests/testthat/testdata", mri_data_folder = "mri", cache_file = "mmy3_demographics.tsv", trial_data_file = "sample_subject_data.csv") {
  # Copy the cached trial data file in the MRI folder to the test data folder
  subj_df <- read.csv(file.path(test_data_base_dir, mri_data_folder, cache_file))

  # Rename lunaid column to id
  subj_df <- dplyr::rename(subj_df, id = lunaid)

  write.csv(subj_df, file.path(test_data_base_dir, trial_data_file), row.names = FALSE)
}


#' Provide a minimal instantiation of the gpa list
#' 
#' @return a minimal gpa list
get_gpa_minimal <- function() {
  setup_glm_pipeline(
    analysis_name = "gpa_tests_l2l3", # output_directory = tempdir(),
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
    analysis_name = "gpa_tests_nol2l3",# output_directory = tempdir(),
    subject_data = get_subj_df(),
    trial_data = get_trial_df(),
    tr = 1.0,
    l1_models=NULL
  )
}

#' Build a gpa object with all information that doesn't require the CLI, and export to RDS.
#' 
#' @param test_data_base_dir the base directory of the test data.
build_gpa_base <- function(
    analysis_name = "gpa_tests",
    test_data_base_dir = "tests/testthat/testdata",
    trial_data_file = "sample_trial_data.csv.gz",
    run_data_file = "sample_run_data.csv",
    subject_data_file = "sample_subject_data.csv",
    cache_file = "gpa_base.rds",
    scheduler = "slurm", drop_volumes = 2,
    exclude_run = "max(FD) > 5 | sum(FD > .9)/length(FD) > .10",
    exclude_subject = "n_good_runs < 4",
    truncate_run = "(FD > 0.9 & time > last_offset) | (time > last_offset + last_isi)",
    spike_volumes = NULL) {


  # Read in trial, run, and subject dataframes
  trial_df <- read.csv(file.path(test_data_base_dir, trial_data_file))
  run_df <- read.csv(file.path(test_data_base_dir, run_data_file))
  subj_df <- read.csv(file.path(test_data_base_dir, subject_data_file))

  gpa <- setup_glm_pipeline(
    analysis_name = analysis_name,
    scheduler = scheduler,
    trial_data = trial_df,
    run_data = run_df,
    subject_data = subj_df, # output_directory = tempdir(), # vm = c(id = "id", run_number = "run"),
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

  saveRDS(gpa, file = file.path(test_data_base_dir, cache_file))
  return(gpa)
}

#' Get a populated GPA object. Pull the base gpa from cache if present; if not, run build_gpa_base.
#' Model objects only loaded in from file for now until they can be programatically created from YAML.
get_gpa_base <- function(
  analysis_name = "gpa_tests",
  test_data_base_dir = "tests/testthat/testdata",
  cache_file = "gpa_base.rds",
  ...
) {
  # Check if gpa object is already cached, if not, build it
  # if (file.exists(file.path(test_data_base_dir, cache_file))) {
  #   gpa <- readRDS(file.path(test_data_base_dir, cache_file))
  # } else {
    gpa <- build_gpa_base(analysis_name = analysis_name, test_data_base_dir = test_data_base_dir, cache_file = cache_file, ...)
  # }

  return(gpa)
}

#' Populate base GPA with models. Save models & final GPA to cache files.
#' Currently interactively-only until we can build all models from spec YAML,
#' so must be called manually, before get_gpa to rebuild models.
#' This COULD be split up so that each model gets their own cache file, but that feels
#' excessive for now. Just assume all models need to be rebuilt if gpa needs rebuilding.
build_gpa <- function(
  analysis_name = "gpa_tests",
  test_data_base_dir = "tests/testthat/testdata",
  gpa_cache_file = "gpa.rds",
  l1_spec_file = "sample_spec.yaml", #"sample_L1_spec.yaml",
  ...
) {
  # run set_gm_pipeline()
  gpa <- get_gpa_base(analysis_name = analysis_name, test_data_base_dir = test_data_base_dir, ...)

  # setup compute enviroment
  gpa <- setup_compute_environment(gpa)

  # Build L1 models
  gpa <- build_l1_models(gpa, from_spec_file = file.path(test_data_base_dir, l1_spec_file))

  # Build L2 models
  gpa <- build_l2_models(gpa)

  # Build L3 models
  gpa <- build_l3_models(gpa)

  # Save final RDS object
  saveRDS(gpa, file = file.path(test_data_base_dir, gpa_cache_file))
  return(gpa)
}

#' Get a populated GPA object. Pull the base gpa and models from cache if present; if not, run build_gpa_base.
get_gpa <- function(
  analysis_name = "gpa_tests",
  test_data_base_dir = "tests/testthat/testdata",
  gpa_cache_file = "gpa.rds",
  ...
) {
  # Check if gpa object is already cached, if not, build it
  # if (file.exists(file.path(test_data_base_dir, gpa_cache_file))) {
  #   gpa <- readRDS(file.path(test_data_base_dir, gpa_cache_file))
  # } else {
    gpa <- build_gpa(analysis_name = analysis_name, test_data_base_dir = test_data_base_dir, gpa_cache_file = gpa_cache_file, ...)
  # }
  
  return(gpa)
}

#' Create a populated GPA object using all the inputs
create_gpa <- function(
  analysis_name = "gpa_tests",
  test_data_base_dir = "/proj/mnhallqlab/projects/fmri.pipeline_test_data",
  l1_spec_file = "sample_2_L1_spec.yaml",
  trial_data_file = "sample_trial_data.csv.gz",
  run_data_file = "sample_run_data.csv",
  subject_data_file = "sample_subject_data.csv",
  gpa_cache_file = "gpa.rds",
  cache_file = "gpa_base.rds"
) {
  gpa <- get_gpa(
    analysis_name = analysis_name,
    test_data_base_dir = test_data_base_dir,
    gpa_cache_file = gpa_cache_file,
    l1_spec_file = l1_spec_file,
    
    # need to figure out what to pass here
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
  # TODO check if we need to run get_gpa each time or we could add a if statment to not run this function if the gpa object already exists in a saved file?
  return(gpa)
}
