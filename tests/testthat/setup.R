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
            mr_dir <- file.path(subject_dir_name, paste0("clock", run_number), "ica_aroma")
            run_nifti <- file.path(mr_dir, "melodic_IC_thr_MNI2mm.nii")
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

    write.csv(run_df, file.path(test_data_base_dir, run_data_file), row.names = FALSE)
}


#' Build the subject dataframe from test data. This file is already saved in the repo, so should only be used to
#' adjust the run data file.
#' Currently this is just copying the subject data from the MRI test folder - do not have script used to build it here yet.
build_trial_data_file <- function(test_data_base_dir = "tests/testthat/testdata", mri_data_folder = "mri", cache_file = "mmy3_demographics.tsv", trial_data_file = "sample_subject_data.csv") {
  # Copy the cached trial data file in the MRI folder to the test data folder
  subj_df <- data.table::fread(file.path(test_data_base_dir, mri_data_folder, cache_file))

  # Rename lunaid column to id
  subj_df <- dplyr::rename(subj_df, id = lunaid)

  write.csv(subj_df, file.path(test_data_base_dir, trial_data_file), row.names = FALSE)
}

build_trial_data_file()


#' Load in subject dataframe from test data.
get_subj_df <- function(test_data_base_dir = "local/test_data", demographics_df_file_name = "mmy3_demographics.tsv") {
  subj_df <- data.table::fread(file.path(test_data_base_dir, demographics_df_file_name))

  # Rename lunaid column to id
  subj_df <- dplyr::rename(subj_df, id = lunaid)

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

#' Provide a gpa list populated from test data.
#' 
#' @param test_data_base_dir the base directory of the test data.
get_gpa <- function(
    test_data_base_dir = "local/test_data",
    scheduler = "slurm", drop_volumes = 2,
    exclude_run = "max(FD) > 5 | sum(FD > .9)/length(FD) > .10",
    exclude_subject = "n_good_runs < 4",
    truncate_run = "(FD > 0.9 & time > last_offset) | (time > last_offset + last_isi)",
    spike_volumes = NULL) {

  setup_glm_pipeline(
    analysis_name = "gpa_tests",
    scheduler = scheduler,
    trial_data = get_trial_df(test_data_base_dir),
    run_data = get_run_df(test_data_base_dir),
    subject_data = get_subj_df(test_data_base_dir),
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

gpa <- get_gpa()

saveRDS(gpa, file = "tests/testthat/testdata/gpa_no_models.rds")

run_df.to_csv("tests/testthat/testdata/sample_trial_data.csv")

# Save run_df to csv at location "tests/testthat/testdata/sample_trial_data.csv"

