#' helper function to locate nifti inputs for each run of data
#'
#' @return a modified version of the \code{gpa} object that contains a $run_nifti field in
#'   the $run_data data.frame.
#' @param gpa a \code{glm_pipeline_arguments} object created by setup_glm_pipeline
#' @keywords internal
#' @author Michael Hallquist
#' @importFrom dplyr count group_by left_join filter
#' @importFrom magrittr %>%
lookup_nifti_inputs <- function(gpa, add_run_volumes = TRUE, add_nvoxels = TRUE) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_logical(add_run_volumes, len=1L)

  # look at whether we have the number of expected runs for each subject
  n_subj_runs <- gpa$run_data %>%
    group_by(id, session) %>%
    count()

  lg <- lgr::get_logger("glm_pipeline/pipeline_setup")
  lg$set_threshold(gpa$lgr_threshold)

  not_expected <- n_subj_runs %>% dplyr::filter(n != gpa$n_expected_runs)

  if (nrow(not_expected) > 1L) {
    lg$warn("Found an unexpected number of runs for some subjects.")
    lg$warn("Subject %s, Session %s, Number of runs %d.", not_expected$id, not_expected$session, not_expected$n)
  }

   #use specific run NIfTIs included in run_data, rather finding these by regex
  if ("run_nifti" %in% names(gpa$run_data)) {
    lg$info("Using run_data to identify NIfTI files for analysis")
  } else {
    # find run nifti files based on directory and regular expression settings
    lg$info("Using regex-based find approach to identify run NIfTIs (this may take a few minutes)")

    mr_list <- list()

    for (ii in seq_len(nrow(gpa$subject_data))) {
      subj_mr_dir <- gpa$subject_data$mr_dir[ii]
      subj_id <- gpa$subject_data$id[ii]
      subj_session <- gpa$subject_data$session[ii]

      if (!dir.exists(file.path(subj_mr_dir))) {
        lg$warn("Unable to find subject data directory: %s for id: %s, session: %s", subj_mr_dir, subj_id, subj_session)
      }

      ## Find processed fMRI run-level data for this subject
      # run_nifti <- list.files(subj_mr_dir, pattern=gpa$fmri_file_regex, full.names=TRUE, recursive=TRUE)
      ## cat(paste0("command: find ", subj_mr_dir, " -iname '", expectfile, "' -ipath '*", expectdir, "*' -type f\n"))

      # -ipath '*", expectdir, "*' -type f | sort -n"), intern=TRUE)
      if (!is.null(gpa$fmri_path_regex)) {
        addon <- paste0(" -ipath '*/", gpa$fmri_path_regex, "/*'")
      } else {
        addon <- ""
      }
      find_string <- paste0(
        "find ", subj_mr_dir, " -regextype posix-egrep -iregex '.*",
        gpa$fmri_file_regex, "'", addon, " -type f | sort -n"
      )
      lg$debug("run_nifti find syntax: %s", find_string)
      run_nifti <- system(find_string, intern = TRUE)

      # extract run number from file name
      mr_run_nums <- as.integer(sub(paste0(gpa$run_number_regex), "\\1", run_nifti, perl = TRUE))

      mr_list[[ii]] <- data.frame(
        id = subj_id, session = subj_session, run_number = mr_run_nums,
        run_nifti = basename(run_nifti), mr_dir = dirname(run_nifti)
      )
    }

    mr_df <- dplyr::bind_rows(mr_list)

    gpa$run_data <- dplyr::left_join(gpa$run_data, mr_df, by = c("id", "session", "run_number"))
  }

  system.time(run_nifti <- get_mr_abspath(gpa$run_data, "run_nifti"))
  mr_found <- file.exists(run_nifti)

  if (any(mr_found != TRUE & !is.na(run_nifti))) {
    lg$warn(
      "Could not find the following run files: %s. Dropping from analysis",
      paste(run_nifti[!mr_found], collapse = ", ")
    )
  }

  # populate field in run_data used to determine availability of data
  gpa$run_data$run_nifti_present <- mr_found

  # add number of volumes for each run
  if (isTRUE(add_run_volumes)) {
    lg$info("Lookup up number of volumes from NIfTI headers")
    gpa$run_data$run_volumes <- sapply(run_nifti, lookup_run_volumes, USE.NAMES=FALSE)
  }

  # add number of voxels for each run (used in FSL FEAT)
  if (isTRUE(add_nvoxels)) {
    lg$info("Lookup up number of voxels (x * y * z * t) from NIfTI headers")
    gpa$run_data$nvoxels <- sapply(run_nifti, lookup_nvoxels, USE.NAMES=FALSE)
  }

  return(gpa)

}