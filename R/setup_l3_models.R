#' This function generates the inputs for level 3 analyses, where multi-subject data are analyzed
#'   in group analyses.
#'
#' @param gpa a \code{glm_pipeline_arguments} object containing analysis speecification
#' @param l3_model_names a subset of L3 models to be setup by this function. If not specified,
#'   all models in gpa$l2_models will be included
#' @param l2_model_names a subset of L2 models to be setup by this function. If not specified,
#'   all models in gpa$l2_models will be included
#' @param l1_model_names a subset of L1 models to be passed to L2 by this function. If not
#'   specified, all models in gpa$l1_models will be included
#'
#' @details
#'   This function will setup FSL level 2 (subject) .fsf files for all combinations of
#'   \code{l2_model_names} and \code{l1_model_names}.
#'
#' @author Michael Hallquist
#' @importFrom checkmate assert_class assert_character assert_data_frame
#' @importFrom lgr get_logger
#' @export
setup_l2_models <- function(gpa, l3_model_names=NULL, l2_model_names=NULL, l1_model_names=NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_character(l3_model_names, null.ok = TRUE)
  checkmate::assert_character(l2_model_names, null.ok = TRUE)
  checkmate::assert_character(l1_model_names, null.ok = TRUE)
  checkmate::assert_data_frame(gpa$run_data)

  # if no l3 model subset is requested, output all models
  if (is.null(l3_model_names)) l3_model_names <- names(gpa$l3_models$models)

  # if no l2 model subset is requested, output all models
  if (is.null(l2_model_names)) l2_model_names <- names(gpa$l2_models$models)

  # if no l1 model subset is requested, output all models
  if (is.null(l1_model_names)) l1_model_names <- names(gpa$l1_models$models)

  lg <- lgr::get_logger("glm_pipeline/l3_setup")
  lg$debug("In setup_l3_models, setting up the following L3 models:")
  lg$debug("L3 model: %s", l3_model_names)
  lg$debug("In setup_l3_models, passing the following L2 models to L3:")
  lg$debug("L2 model: %s", l2_model_names)
  lg$debug("In setup_l3_models, passing the following L1 models to L3:")
  lg$debug("L1 model: %s", l1_model_names)

  l1_meta <- gpa$l1_model_setup$metadata %>%
    dplyr::select(id, session, run_number, exclude_run, exclude_subject)

  excluded_runs <- l1_meta %>%
    dplyr::filter(exclude_run == TRUE | exclude_subject == TRUE)

  if (nrow(excluded_runs) > 1L) {
    lg$info("In setup_l2_models, the following runs will be excluded from L2 modeling: ")
    lg$info(
      "  subject: %s, session: %s, run_number: %s",
      excluded_runs$id, excluded_runs$session, excluded_runs$run_number
    )
  }

  # only retain good runs and subjects
  l1_meta <- l1_meta %>%
    dplyr::filter(exclude_run == FALSE & exclude_subject == FALSE)

  run_data <- gpa$run_data %>%
    right_join(l1_meta, by=c("id", "session", "run_number"))

  if (nrow(l1_meta) == 0L) {
    msg <- "In setup_l2_models, no runs survived the exclude_subject and exclude_run step."
    lg$warn(msg)
    warning(msg)
    return(NULL)
  }

  if (isTRUE(gpa$log_txt)) {
    # TODO: abstract the log file name to finalize_pipeline_configuration function
    lg$add_appender(lgr::AppenderFile$new("setup_l2_models.txt"), name = "txt")
  }

  if (is.null(gpa$l1_model_setup) || !inherits(gpa$l1_model_setup, "l1_setup")) {
    lg$error("No l1_model_setup found in the glm pipeline object.")
    lg$error("You must run setup_l1_models before running setup_l2_models.")
    stop("No l1_model_setup found in the glm pipeline object.",
    "You must run setup_l1_models before running setup_l2_models.")
  }

  # respecify L2 models for each subject based on available runs
  for (mname in l2_model_names) {
    lg$info("Recalculating per-subject L2 models based on available runs for model: %s", mname)
    gpa$l2_models$models[[mname]] <- respecify_l2_models_by_subject(gpa$l2_models$models[[mname]], run_data)
  }

  # loop over and setup all requested combinations of L1, L2, and L3 models
  feat_l3_df <- list()
  ff <- 1
  for (ii in seq_along(l1_model_names)) {
    this_l1_model <- l1_model_names[ii]
    for (jj in seq_along(l2_model_names)) {
      this_l2_model <- l2_model_names[jj]
      for (kk in seq_along(l3_model_names)) {
        this_l3_model <- l3_model_names[kk]

      if ("fsl" %in% gpa$glm_software) {
        # get list of runs to examine/include
        to_run <- gpa$l1_model_setup$fsl %>%
          dplyr::filter(l1_model == !!this_l1_model) %>%
          dplyr::select(id, session, run_number, l1_model, l1_feat_fsf, l1_feat_dir)

        # handle run and subject exclusions (exclude_run should be FALSE in l1_meta, per filter above)
        to_run <- dplyr::left_join(l1_meta, to_run, by = c("id", "session", "run_number"))
        data.table::setDT(to_run) #convert to data.table for split

        by_subj_session <- split(to_run, by=c("id", "session"))

        # setup Feat L2 files for each id and session
        for (l1_df in by_subj_session) {
          subj_id <- l1_df$id[1L]
          subj_session <- l1_df$session[1L]
          feat_l3_df[[ff]] <- tryCatch({
              fsl_l2_model(
                l1_df = l1_df,
                l2_model_name = this_l2_model, gpa = gpa
              )
            },
            error = function(e) {
              lg$error(
                "Problem with fsl_l2_model. L1 Model: %s, L2 Model: %s, Subject: %s, Session: %s",
                this_l1_model, this_l2_model, subj_id, subj_session
              )
              lg$error("Error message: %s", as.character(e))
              return(NULL)
            }
          )

          if (!is.null(feat_l3_df[1L])) ff <- ff + 1 #increment position in multi-subject list
        }

      }

      if ("spm" %in% gpa$glm_software) {
        lg$warn("spm not supported in setup_l2_models")
      }

      if ("afni" %in% gpa$glm_software) {
        lg$warn("afni not supported in setup_l2_models")
      }
      }
    }
  }

  all_subj_l2_combined <- list(
    fsl=rbindlist(feat_l3_df)
  )

  class(all_subj_l2_combined) <- c("list", "l2_setup")

  # append l2 setup to gpa
  gpa$l2_model_setup <- all_subj_l2_combined

  return(gpa)
}
