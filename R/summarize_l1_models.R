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

summarize_l1_models <- function(gpa) {
    # if there is a gpa object, set that to l1_model_set
    if (!is.null(gpa)) {
        lg$info("In build_l1_models, using existing gpa object to build l1 models (ignoring trial_data argument etc.)")
        use_gpa <- TRUE
        trial_data <- gpa$trial_data
        l1_model_set <- gpa$l1_models
    } else {
        lg$debug("In build_l1_models, using trial_data passed in, rather than gpa object")
        use_gpa <- FALSE
    }

    ## l1_model_set creation if there is no gpa object

    checkmate::assert_data_frame(trial_data, null.ok = FALSE)
    checkmate::assert_class(l1_model_set, "l1_model_set", null.ok = TRUE)
    checkmate::assert_subset(onset_cols, names(trial_data)) #make sure that all onset columns are in the data frame
    checkmate::assert_string(onset_regex, null.ok = TRUE)
    checkmate::assert_subset(duration_cols, names(trial_data)) # make sure that all duration columns are in the data frame
    checkmate::assert_string(duration_regex, null.ok = TRUE)
    checkmate::assert_subset(value_cols, names(trial_data)) #make sure all parametric regressor columns are in the data frame
    checkmate::assert_string(value_regex, null.ok = TRUE)
    checkmate::assert_subset(c("id", "session", "run_number", "trial"), names(trial_data)) # required metadata in trial_data

  if (is.null(l1_model_set)) {
    ## initialize overall l1 design object (holds events, signals, and models)
    l1_model_set <- list(onsets=NULL, durations=NULL, values=NULL, wi_factors=NULL, events=NULL, signals=NULL, models=NULL)
    class(l1_model_set) <- c("list", "l1_model_set")
    new_l1 <- TRUE
  } else {
    new_l1 <- FALSE
  }

    summarize_events <- function(l1_model_set) {
        cat("Summary of events available in l1 models:\n--------\n")
        if (!is.null(l1_model_set$events)) {
        lapply(l1_model_set$events, function(x) {
            cat("Event name:     ", x$name, "\n")
            cat("  Event onset:    ", x$onset, "\n")
            cat(
            "    mean [min -- max]: ",
            round(mean(x$data$onset, na.rm = TRUE), 2), "[",
            round(min(x$data$onset, na.rm = TRUE), 2), "--",
            round(max(x$data$onset, na.rm = TRUE), 2), "]\n"
            )
            cat("    First 6 values:", paste(head(x$data$onset), collapse = ", "), "\n\n")
            cat("  Event duration: ", x$duration, "\n")
            cat(
            "    mean [min -- max]: ",
            round(mean(x$data$duration, na.rm = TRUE), 2), "[",
            round(min(x$data$duration, na.rm = TRUE), 2), "--",
            round(max(x$data$duration, na.rm = TRUE), 2), "]\n"
            )
            cat("    First 6 values:", paste(round(head(x$data$duration), 2), collapse = ", "), "\n\n")
            if (!is.null(x$isi)) {
            cat("  Event ISI/ITI:    ", x$isi, "\n")
            cat(
                "    mean [min -- max]: ",
                round(mean(x$data$isi, na.rm = TRUE), 2), "[",
                round(min(x$data$isi, na.rm = TRUE), 2), "--",
                round(max(x$data$isi, na.rm = TRUE), 2), "]\n"
            )
            }
        })
        } else {
        cat("No events in l1 models\n")
        }
    }

    summarize_events(l1_model_set)
}
