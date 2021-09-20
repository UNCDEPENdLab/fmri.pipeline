#' This function generates the inputs for an FSL level 2 analysis, where multiple runs for a subject are combined using
#' fixed effects estimation.
#'
#' @param gpa a \code{glm_pipeline_arguments} object containing analysis speecification
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
#' @importFrom iterators iter
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach registerDoSEQ
#' @export
setup_l2_models <- function(gpa, l1_model_names=NULL, l2_model_names=NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_data_frame(gpa$run_data)
  checkmate::assert_class(gpa$l2_models, "hi_model_set")
  checkmate::assert_class(gpa$l1_models, "l1_model_set")
  checkmate::assert_character(l1_model_names, null.ok = TRUE)
  checkmate::assert_subset(l1_model_names, names(gpa$l1_models$models))
  checkmate::assert_character(l2_model_names, null.ok = TRUE)
  checkmate::assert_subset(l2_model_names, names(gpa$l2_models$models))

  # if no l2 model subset is requested, output all models
  if (is.null(l2_model_names)) l2_model_names <- names(gpa$l2_models$models)

  # if no l1 model subset is requested, output all models
  if (is.null(l1_model_names)) l1_model_names <- names(gpa$l1_models$models)

  lg <- lgr::get_logger("glm_pipeline/l2_setup")
  if (isTRUE(gpa$log_txt) && !"setup_l2_log_txt" %in% names(lg$appenders)) {
    lg$add_appender(lgr::AppenderFile$new(gpa$output_locations$setup_l2_log_txt), name = "setup_l2_log_txt")
  }

  if (isTRUE(gpa$log_json) && !"setup_l2_log_json" %in% names(lg$appenders)) {
    lg$add_appender(lgr::AppenderJson$new(gpa$output_locations$setup_l2_log_json), name = "setup_l2_log_json")
  }

  lg$debug("In setup_l2_models, setting up the following L2 models:")
  lg$debug("L2 model: %s", l2_model_names)
  lg$debug("In setup_l2_models, passing the following L1 models to L2:")
  lg$debug("L1 model: %s", l1_model_names)

  # setup parallel worker pool, if requested
  if (!is.null(gpa$parallel$l2_setup_cores) && gpa$parallel$l2_setup_cores > 1L) {
    lg$info("Initializing l2 setup cluster with %d cores", gpa$parallel$l2_setup_cores)
    cl <- parallel::makeCluster(gpa$parallel$l2_setup_cores)
    doParallel::registerDoParallel(cl)
    on.exit(try(parallel::stopCluster(cl))) # cleanup pool upon exit of this function
  } else {
    lg$info("Initializing l2 setup with serial execution")
    foreach::registerDoSEQ() # formally register a sequential 'pool' so that dopar is okay
  }

  excluded_runs <- gpa$run_data %>%
    dplyr::select(id, session, run_number, exclude_run, exclude_subject) %>%
    dplyr::filter(exclude_run == TRUE | exclude_subject == TRUE)

  if (nrow(excluded_runs) > 0L) {
    lg$info("In setup_l2_models, the following runs will be excluded from L2 modeling: ")
    lg$info(
      "  subject: %s, session: %s, run_number: %s",
      excluded_runs$id, excluded_runs$session, excluded_runs$run_number
    )
  }

  # only retain good runs and subjects
  run_data <- gpa$run_data %>%
    dplyr::filter(exclude_run == FALSE & exclude_subject == FALSE)

  # subset basic metadata to merge against a given l1 model to enforce run/subject exclusions
  good_runs <- run_data %>%
    dplyr::select(id, session, run_number, exclude_run, exclude_subject)

  if (nrow(run_data) == 0L) {
    msg <- "In setup_l2_models, no runs survived the exclude_subject and exclude_run step."
    lg$warn(msg)
    warning(msg)
    return(NULL)
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

  # refresh l1 model status in $l1_model_setup
  gpa <- refresh_feat_status(gpa, level=1L, lg=lg)

  # loop over and setup all requested combinations of L1 and L2 models
  model_set <- expand.grid(l1_model = l1_model_names, l2_model = l2_model_names, stringsAsFactors = FALSE)
  all_l2_list <- foreach(
    model_info = iter(model_set, by = "row"), .inorder = FALSE,
    .packages = c("dependlab", "dplyr", "data.table"),
    .export = c("lg", "gpa", "fsl_l2_model")
  ) %dopar% {

    model_info <- model_info # avoid complaints about visible global binding in R CMD check
    this_l1_model <- model_info$l1_model
    this_l2_model <- model_info$l2_model

    l2_file_setup <- list(fsl = list(), spm = list(), afni = list())

    if ("fsl" %in% gpa$glm_software) {
      # get list of runs to examine/include
      to_run <- gpa$l1_model_setup$fsl %>%
        dplyr::filter(l1_model == !!this_l1_model) %>%
        dplyr::select(id, session, run_number, l1_model, feat_fsf, feat_dir)

      # handle run and subject exclusions by joining against good runs
      to_run <- dplyr::inner_join(good_runs, to_run, by = c("id", "session", "run_number"))
      data.table::setDT(to_run) # convert to data.table for split

      by_subj_session <- split(to_run, by = c("id", "session"))

      # setup Feat L2 files for each id and session
      for (l1_df in by_subj_session) {
        subj_id <- l1_df$id[1L]
        subj_session <- l1_df$session[1L]
        feat_l2_df <- tryCatch({
          fsl_l2_model(
            l1_df = l1_df,
            l2_model = this_l2_model, gpa = gpa
          )},
          error = function(e) {
            lg$error(
              "Problem with fsl_l2_model. L1 Model: %s, L2 Model: %s, Subject: %s, Session: %s",
              this_l1_model, this_l2_model, subj_id, subj_session
            )
            lg$error("Error message: %s", as.character(e))
            return(NULL)
          }
        )

        if (!is.null(feat_l2_df)) {
          # add to tracking data.frame
          l2_file_setup$fsl <- rbind(l2_file_setup$fsl, feat_l2_df)
        }
      }
    }

    if ("spm" %in% gpa$glm_software) {
      lg$warn("spm not supported in setup_l2_models")
    }

    if ("afni" %in% gpa$glm_software) {
      lg$warn("afni not supported in setup_l2_models")
    }

    return(l2_file_setup)
  }

  all_subj_l2_combined <- list(
    fsl = rbindlist(lapply(all_l2_list, "[[", "fsl"))
    # spm = rbindlist(lapply(all_l2_list, "[[", "spm"))
    # afni = rbindlist(lapply(all_l2_list, "[[", "afni"))
  )

  class(all_subj_l2_combined) <- c("list", "l2_setup")

  # append l2 setup to gpa
  gpa$l2_model_setup <- all_subj_l2_combined

  return(gpa)
}
