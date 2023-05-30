
summarize_l1_models <- function(gpa) {

    checkmate::assert_list(gpa, null.ok = FALSE)

    l1_model_set <- gpa$l1_models

   summarize_events(l1_model_set)
    print("hi")
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

  summarize_pipeline <- function(gpa) {
        cat("\nSummary of GLM Pipeline Analysis: \n\n")
        if (!is.null(gpa)) {
            cat("Analysis Name: ", gpa$analysis_name[1])
          # GLM and Confound Settings
            cat("\n--------\n Settings: \n")
            cat(" -  Scheduler Name: ", gpa$scheduler, "\n")
            cat(" -  TR: ", gpa$tr, "\n")
            cat(" -  Multi-Level Run: ", gpa$multi_run[1], "\n")
            cat(" -  GLM Settings: ", gpa$glm_settings, "\n")
            cat(" -  NA Strings:", paste(gpa$confound_settings$na_strings, collapse = ", "), "\n")
            if (!is.null(gpa$confound_settings$exclude_run)) {
              cat(" -  Exclude Run: ", gpa$confound_settings$exclude_run, "\n")
            } else {
              cat(" -  Exclude Run: NA\n")
            }
            if (!is.null(gpa$confound_settings$exclude_subject)) {
              cat(" -  Exclude Subject: ", gpa$confound_settings$exclude_subject, "\n")
            } else {
              cat(" -  Exclude Subject: NA\n")
            }
            if (!is.null(gpa$confound_settings$truncate_run)) {
              cat(" -  Truncate Run: ", gpa$confound_settings$truncate_run, "\n")
            } else {
              cat(" -  Truncate Run: NA\n")
            }
            if (!is.null(gpa$confound_settings$spike_volumes)) {
              cat(" -  Spike Volumes: ", gpa$confound_settings$spike_volumes, "\n")
            } else {
              cat(" -  Spike Volumes: NA\n")
            }
            # VM
            cat("--------\n VM: \n")
            cat("  ", paste(gpa$vm, collapse = ", "), "\n")
            # Output Directory and System Info
            cat("--------\n System Info: \n")
            cat(" -  Output Directory: ", gpa$output_directory[1], "\n")
            cat(" -  System Info: ", gpa$sys_info, "\n")
        } else {
        cat("No info in gpa object\n")
        }
    }

    # summarize_pipeline(gpa)
    # summarize_l1_models(gpa)