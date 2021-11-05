#' Identify the last valid volume acquired in a given run.
#' Subjects often exhibit head movement after run ends (MATLAB closes), but scan hasn't stopped
#' This occurs because the MB raw transfer of the prior run is occurring, but does not finish before the current run
#' Thus, truncate mr files to be 12 seconds after final feedback presentation, which is how the paradigm timing files are setup
#' 
#' @keywords internal
#' @importFrom dplyr filter
#' @importFrom checkmate assert_data_frame test_directory_exists test_number test_logical
#' @importFrom data.table fread fwrite
#' @importFrom glue glue_data
truncate_runs <- function(mr_df, subj_outdir=NULL, lg=NULL) {
  checkmate::assert_data_frame(mr_df)
  if (!is.null(subj_outdir) && !checkmate::test_directory_exists(subj_outdir)) {
    lg$warn("In truncate_runs, cannot locate intended output directory: %s", subj_outdir)
    subj_outdir <- NULL # fall back to location of run_nifti
  }
  checkmate::assert_class(lg, "Logger")

  mr_df <- do.call(rbind, lapply(seq_len(nrow(mr_df)), function(r) {
    # if no truncation, default to using all volumes
    this_run <- mr_df %>% dplyr::slice(r)
    run_volumes <- this_run$run_volumes
    this_run$orig_volumes <- run_volumes # always keep record of original length
    last_volume <- run_volumes # default to all volumes (no tail truncation)
    has_cfile <- checkmate::test_file_exists(this_run$l1_confound_file)

    if (!is.null(gpa$confound_settings$truncate_run)) {
      truncation_data <- read_df_sqlite(
        gpa = gpa, id = this_run$id, session = this_run$session,
        run_number = this_run$run_number, table = "l1_truncation_data"
      )

      if (is.null(truncation_data)) {
        lg$warn(
          "Unable to find l1_truncation data in SQLite database for id: %s, session: %d, run_number: %d. Defaulting to no truncation.",
          this_run$id, this_run$session, this_run$run_number
        )
      } else {
        # add last onset and offset to truncation_data for calculating expression
        truncation_data <- truncation_data %>% bind_cols(this_run %>% dplyr::select(last_onset, last_offset))
        last_volume <- tryCatch(with(truncation_data, eval(parse(text = gpa$confound_settings$truncate_run))),
          error = function(e) {
            lg$error(
              "Problem evaluating truncation for subject: %s, session: %s, run_number: %s, expr: %s",
              this_run$id, this_run$session, this_run$run_number,
              gpa$confound_settings$truncate_run
            )
            lg$error("Defaulting to last volume (no truncation).")
            return(run_volumes)
          }
        )

        if (checkmate::test_number(last_volume)) {
          # I guess nothing to do here, right? -- the expression returns a number of the last good volume
        } else if (checkmate::test_logical(last_volume)) {
          # if we have a vector of logicals (more common, probably), find the volume before the first TRUE occurs
          # only truncate if any of the volumes evaluates to TRUE in the truncation expression (otherwise stick with default of last volume)
          if (any(last_volume == TRUE)) {
            last_volume <- min(which(last_volume == TRUE)) - 1
          } else {
            last_volume <- run_volumes
          }
        }

        if (last_volume > mr_df$run_volumes[r]) {
          lg$warn(
            "truncate_run expression evaluated to %d, but run length is %d. Will use %d",
            last_volume, run_volumes, run_volumes
          )
          last_volume <- run_volumes
        } else if (last_volume < 3) {
          lg$warn("truncate_run expression evaluated to %d, which is awfully low! Falling back to %d", last_volume, this_run$run_volumes)
          last_volume <- run_volumes
        }

      }
    }

    # first volume to use for analysis
    first_volume <- this_run$drop_volumes + 1

    # length of truncated file
    final_volumes <- last_volume - first_volume + 1

    # update run_volumes to reflect any drops and truncation
    this_run$run_volumes <- final_volumes

    # number of volumes that were truncated at the end
    this_run$truncate_volumes <- run_volumes - last_volume

    # see whether truncation is needed/specified
    if (final_volumes < run_volumes) {
      # put truncated file in subject output directory, if provided
      if (has_cfile) {
        conext <- file_ext(this_run$l1_confound_file)
        cname <- glue_data("sub-{id}_ses-{session}_run-{run_number}_drop-{drop_volumes}_trunc-{truncate_volumes}_confounds{conext}", .x = this_run)
        c_odir <- ifelse(is.null(subj_outdir), dirname(this_run$l1_confound_file), subj_outdir)
        trunc_confounds <- file.path(c_odir, cname)
        if (!checkmate::test_file_exists(trunc_confounds)) {
          cdata <- data.table::fread(this_run$l1_confound_file, na.strings = gpa$confound_settings$na_strings, data.table = FALSE, header = FALSE)
          cdata <- cdata[first_volume:last_volume, , drop = FALSE]
          lg$debug("Creating truncated confounds file: %s", trunc_confounds)
          data.table::fwrite(cdata, file = trunc_confounds, row.names = FALSE, col.names = FALSE)
        }
        this_run$l1_confound_file <- trunc_confounds
      }

      imgext <- file_ext(this_run$run_nifti)
      fname <- glue_data("sub-{id}_ses-{session}_run-{run_number}_drop-{drop_volumes}_trunc-{truncate_volumes}{imgext}", .x = this_run)
      img_odir <- ifelse(is.null(subj_outdir), dirname(this_run$run_nifti), subj_outdir)
      trunc_file <- file.path(img_odir, fname)
      if (!file.exists(trunc_file)) {
        lg$debug("Creating truncated file: %s", trunc_file)
        # create truncated image with fslroi, which uses 0-based indexing
        runFSLCommand(paste("fslroi", this_run$run_nifti, trunc_file, first_volume - 1, final_volumes))
      }
      this_run$run_nifti <- trunc_file
    }

    return(this_run)
  }))

  mr_df

}
