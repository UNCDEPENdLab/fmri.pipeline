gpa <- NULL

trial_data <- data.frame(matrix(ncol=12, nrow=5))
trial_data[is.na(trial_data)] <- 0
trial_data

colnames(trial_data) <- c(
  "id", "run_number", "trial", "trial_type", "outcome_fac", "choice_onset",
  "choice_time", "reaction_time", "feedback_onset", "feedback_isi", "iti_actual", "feedback_duration")

trial_data
gpa <- setup_glm_pipeline(analysis_name="basic_apr262023", scheduler="slurm",
    output_directory = "/theryansmith/Downloads", trial_data=trial_df,
    subject_data=subject_df, run_data=run_df,
    tr=0.635, drop_volumes=2,
    l1_modes=NULL, l2_models=NULL, l3_models=NULL,
    n_expected_runs=4,
    confound_settings=list(
        #confound_input_colnames = c("csf", "csf_derivative1", "white_matter", "white_matter_derivative1"), #assumption
        l1_confound_regressors = c("csf")
    ))
gpa <- build_l1_models(gpa)

summarize_l1_models <- function(gpa, what='all')


