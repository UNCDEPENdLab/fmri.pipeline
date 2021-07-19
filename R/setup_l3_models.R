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
setup_l3_models <- function(gpa, l3_model_names = NULL, l2_model_names = NULL, l1_model_names = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_character(l3_model_names, null.ok = TRUE)
  checkmate::assert_character(l2_model_names, null.ok = TRUE)
  checkmate::assert_character(l1_model_names, null.ok = TRUE)
  checkmate::assert_data_frame(gpa$run_data)
  checkmate::assert_data_frame(gpa$subject_data)
  checkmate::assert_class(gpa$l3_models, "hi_model_set")
  checkmate::assert_class(gpa$l1_models, "l1_model_set")

  # full 3-level analysis (runs, subject, sample)
  if (isTRUE(gpa$multi_run)) {
    checkmate::assert_class(gpa$l2_models, "hi_model_set")

    # if no l2 model subset is requested, output all models
    if (is.null(l2_model_names)) l2_model_names <- names(gpa$l2_models$models)

  }

  # if no l3 model subset is requested, output all models
  if (is.null(l3_model_names)) l3_model_names <- names(gpa$l3_models$models)

  # if no l1 model subset is requested, output all models
  if (is.null(l1_model_names)) l1_model_names <- names(gpa$l1_models$models)

  lg <- lgr::get_logger("glm_pipeline/l3_setup")
  lg$debug("In setup_l3_models, setting up the following L3 models:")
  lg$debug("L3 model: %s", l3_model_names)
  if (isTRUE(gpa$multi_run)) {
    lg$debug("In setup_l3_models, passing the following L2 models to L3:")
    lg$debug("L2 model: %s", l2_model_names)
  }  
  lg$debug("In setup_l3_models, passing the following L1 models to L3:")
  lg$debug("L1 model: %s", l1_model_names)

  excluded_runs <- gpa$run_data %>%
    dplyr::select(id, session, run_number, exclude_run, exclude_subject) %>%
    dplyr::filter(exclude_run == TRUE | exclude_subject == TRUE)

  if (nrow(excluded_runs) > 1L) {
    lg$info("In setup_l3_models, the following runs will be excluded from L3 modeling: ")
    lg$info(
      "  subject: %s, session: %s, run_number: %s",
      excluded_runs$id, excluded_runs$session, excluded_runs$run_number
    )
  }

  # only retain good runs and subjects
  run_data <- gpa$run_data %>%
    dplyr::filter(exclude_run == FALSE & exclude_subject == FALSE)

  if (nrow(run_data) == 0L) {
    msg <- "In setup_l3_models, no runs survived the exclude_subject and exclude_run step."
    lg$warn(msg)
    warning(msg)
    return(NULL)
  }

  if (isTRUE(gpa$log_txt)) {
    # TODO: abstract the log file name to finalize_pipeline_configuration function
    lg$add_appender(lgr::AppenderFile$new("setup_l3_models.txt"), name = "txt")
  }

  if (is.null(gpa$l1_model_setup) || !inherits(gpa$l1_model_setup, "l1_setup")) {
    lg$error("No l1_model_setup found in the glm pipeline object.")
    lg$error("You must run setup_l1_models before running setup_l3_models.")
    stop(
      "No l1_model_setup found in the glm pipeline object.",
      "You must run setup_l1_models before running setup_l3_models."
    )
  }

  if (is.null(gpa$l2_model_setup) || !inherits(gpa$l2_model_setup, "l2_setup")) {
    lg$error("No l2_model_setup found in the glm pipeline object.")
    lg$error("You must run setup_l2_models before running setup_l3_models.")
    stop(
      "No l2_model_setup found in the glm pipeline object.",
      "You must run setup_l2_models before running setup_l3_models."
    )
  }

  # loop over and setup all requested combinations of L1, L2, and L3 models
  feat_l3_df <- list()
  model_set <- expand.grid(
    l1_model = l1_model_names, l2_model = l2_model_names,
    l3_model = l3_model_name, stringsAsFactors = FALSE
  )

  ff <- 1
  all_l3_list <- foreach(
    model_info = iter(model_set, by = "row"), .inorder = FALSE, 
    .packages = c("dependlab", "dplyr", "data.table"), .export = c("lg", "gpa", "fsl_l3_model")
  ) %dopar% {
    this_l1_model <- model_info$l1_model
    this_l2_model <- model_info$l2_model
    this_l3_model <- model_info$l3_model

    if ("fsl" %in% gpa$glm_software) {
      # get list of runs to examine/include
      to_run <- gpa$l1_model_setup$fsl %>%
        dplyr::filter(l1_model == !!this_l1_model) %>%
        dplyr::select(id, session, run_number, l1_model, feat_fsf, feat_dir)

      # handle run and subject exclusions (exclude_run should be FALSE in l1_meta, per filter above)
      to_run <- dplyr::left_join(l1_meta, to_run, by = c("id", "session", "run_number"))
      data.table::setDT(to_run) # convert to data.table for split

      by_subj_session <- split(to_run, by = c("id", "session"))

      # setup Feat L3 files for each combination of lower-level models

      for (l1_df in by_subj_session) {
        subj_id <- l1_df$id[1L]
        subj_session <- l1_df$session[1L]
        feat_l3_df[[ff]] <- tryCatch(
          {
            fsl_l3_model(
              l1_df = l1_df,
              l2_model_name = this_l2_model, gpa = gpa
            )
          },
          error = function(e) {
            lg$error(
              "Problem with fsl_l3_model. L1 Model: %s, L2 Model: %s, Subject: %s, Session: %s",
              this_l1_model, this_l2_model, subj_id, subj_session
            )
            lg$error("Error message: %s", as.character(e))
            return(NULL)
          }
        )

        if (!is.null(feat_l3_df[1L])) ff <- ff + 1 # increment position in multi-subject list
      }
    }

    if ("spm" %in% gpa$glm_software) {
      lg$warn("spm not supported in setup_l3_models")
    }

    if ("afni" %in% gpa$glm_software) {
      lg$warn("afni not supported in setup_l3_models")
    }
  }

  all_subj_l3_combined <- list(
    fsl = rbindlist(feat_l3_df)
  )

  class(all_subj_l3_combined) <- c("list", "l3_setup")

  # append l3 setup to gpa
  gpa$l3_model_setup <- all_subj_l3_combined

  return(gpa)
}
