get_trial_df <- function(test_data_base_dir) {
  trial_df <- data.table::fread(file.path(base_dir, "sample_trial_data.csv"))
  return(trial_df)
}

get_run_df <- function(test_data_base_dir) {
    # Initialize empty lists to store data
    id_list <- vector("character")
    task_name_list <- vector("character")
    mr_dir_list <- vector("character")
    run_number_list <- vector("integer")
    run_nifti_list <- vector("character")
    confound_input_file_list <- vector("character")

    # Iterate over folders in base directory
    subject_folders <- list.dirs(test_data_base_dir, recursive = FALSE)
    for (subject_folder in subject_folders) {
        subject_dir_name <- basename(subject_folder)

        # Get subject id from folder name
        parts <- strsplit(subject_dir_name, "_")[[1]]
        id <- parts[1]  # Subject ID
        print(id)
        print(subject_folder)

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

    # Print the DataFrame
    return(df)
}

get_subj_df <- function(test_data_base_dir) {
  subj_df <- data.table::fread(file.path(test_data_base_dir, "sample_subject_data.tsv"))
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