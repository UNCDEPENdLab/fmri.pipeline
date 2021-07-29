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
#' FSL 2-level versus 3-level setup
#'
#' 2-level setup (one run per subject)
#'   - Pass L1 .feat folders as input to L3 .fsf setup
#'   - In this approach, the copes in the .fsf pertain to the L1 cope numbers
#'   - Requires one .fsf per L3 model
#'
#' 3-level setup (multiple runs per subject, combined at L2)
#'   - Pass individual cope*.feat folders within subject .gfeat folders
#'   - The folder cope numbers pertain to L1 copes
#'   - The cope*.nii.gz in the cope*.feat subfolders pertain to the L2 contrasts
#'   - Requires one .fsf per L1 cope x L3 model combination
#'   - Example: FSL_L2.gfeat/cope3.feat/stats/cope1.nii.gz
#'      ==> cope3 is the third contrast in the L1 feat model
#'      ==> cope1 is the first contrast in the L2 feat model
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

  # refresh l2 model status in $l2_model_setup
  gpa <- refresh_feat_status(gpa, level = 2L, lg = lg)

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

  if (is.null(gpa$l1_model_setup) || !inherits(gpa$l1_model_setup, "l1_setup")) {
    lg$error("No l1_model_setup found in the glm pipeline object.")
    lg$error("You must run setup_l1_models before running setup_l3_models.")
    stop(
      "No l1_model_setup found in the glm pipeline object.",
      "You must run setup_l1_models before running setup_l3_models."
    )
  }

  if (isTRUE(gpa$multi_run)) {
    lg$info("In setup_l3_models, using a multi-run 3-level setup with runs (l1), subjects (l2), sample (l3)")

    # in multi-run setup, an l2_model_setup must be present
    if (is.null(gpa$l2_model_setup) || !inherits(gpa$l2_model_setup, "l2_setup")) {
      lg$error("No l2_model_setup found in the glm pipeline object.")
      lg$error("You must run setup_l2_models before running setup_l3_models.")
      stop(
        "No l2_model_setup found in the glm pipeline object.",
        "You must run setup_l2_models before running setup_l3_models."
      )
    }

  } else {
    lg$info("In setup_l3_models, using a single run 2-level setup with subjects (l1), sample (l3)")
      # only retain good runs and subjects
      run_data <- gpa$run_data %>%
        dplyr::filter(exclude_run == FALSE & exclude_subject == FALSE)

      if (nrow(run_data) == 0L) {
        msg <- "In setup_l3_models, no runs survived the exclude_subject and exclude_run step."
        lg$warn(msg)
        warning(msg)
        return(NULL)
      }
  }

  # subjects and sessions to run at l3
  subj_df <- gpa$subject_data %>%
    dplyr::filter(exclude_subject==FALSE) %>%
    select(id, session)

  if (isTRUE(gpa$log_txt)) {
    # TODO: abstract the log file name to finalize_pipeline_configuration function
    lg$add_appender(lgr::AppenderFile$new("setup_l3_models.txt"), name = "txt")
  }

  # loop over and setup all requested combinations of L1, L2, and L3 models
  feat_l3_df <- list()
  if (isTRUE(gpa$multi_run)) {
    model_set <- expand.grid(
      l1_model = l1_model_names, l2_model = l2_model_names,
      l3_model = l3_model_names, stringsAsFactors = FALSE
    )
  } else {
    model_set <- expand.grid(l1_model = l1_model_names, l3_model = l3_model_names, stringsAsFactors = FALSE)
  }

  # Get the copes for each contrast, respecting differences in l2 copes across subjects
  # This has the full combination of l1, l2, and l3 copes
  l3_cope_config <- get_fsl_l3_model_df(gpa, model_set, subj_df)

  # For obtaining inputs to l3 models, we don't need the multiple rows for each third level cope
  # since these are part of the l3 model (contrasts to be specified). Thus, just get the first row of
  # each l3 model for calculating inputs to l3.
  l3_cope_input_df <- l3_cope_config %>%
    dplyr::filter(l3_cope_number == 1L) %>%
    dplyr::select(-l3_cope_number, -l3_cope_name)

  to_run <- get_feat_l3_inputs(gpa, l3_cope_input_df, lg)

  all_l3_list <- foreach(
    model_info = iter(to_run), .inorder = FALSE,
    .packages = c("dependlab", "dplyr", "data.table"), .export = c("lg", "gpa", "fsl_l3_model")
  ) %dopar% {
    model_info <- model_info # to avoid complaints about global variable binding in R CMD check

    l3_file_setup <- list(fsl = list(), spm = list(), afni = list())

    if ("fsl" %in% gpa$glm_software) {

      if (nrow(model_info) <= 3) {
        lg$warn(
          "Fewer than 4 complete feat input directories for l1 model %s, l2 model %s, l3 model %s",
          model_info$l1_model, model_info$l2_model, model_info$l3_model
        )
        l3_file_setup$fsl <- NULL
      }

      # setup Feat L3 files for each combination of lower-level models
      l3_file_setup$fsl <- tryCatch(fsl_l3_model(model_info, gpa = gpa),
        error = function(e) {
          lg$error(
            "Problem with fsl_l3_model. L1 Model: %s, L2 Model: %s, L3 model %s",
            model_info$l1_model, model_info$l2_model, model_info$l3_model
          )
          lg$error("Error message: %s", as.character(e))
          return(NULL)
        }
      )
    }

    if ("spm" %in% gpa$glm_software) {
      lg$warn("spm not supported in setup_l3_models")
    }

    if ("afni" %in% gpa$glm_software) {
      lg$warn("afni not supported in setup_l3_models")
    }

    return(l3_file_setup)
  }

  all_subj_l3_combined <- list(
    metadata = l3_cope_config,
    fsl = rbindlist(lapply(all_l3_list, "[[", "fsl"))
  )

  class(all_subj_l3_combined) <- c("list", "l3_setup")

  # append l3 setup to gpa
  gpa$l3_model_setup <- all_subj_l3_combined

  return(gpa)
}


############
# helper function to get a cope data.frame for all level 1 models
get_l1_cope_df <- function(gpa, model_set, subj_df) {
  checkmate::assert_data_frame(model_set)
  dt <- data.table::rbindlist(
    lapply(unique(model_set$l1_model), function(mm) {
      data.frame(
        l1_model = mm,
        l1_cope_number = seq_along(gpa$l1_cope_names[[mm]]),
        l1_cope_name = gpa$l1_cope_names[[mm]]
      )
    })
  )

  #return subject-specific rows for each cope
  dt <- dt %>% tidyr::crossing(subj_df)
  return(dt)
}

# helper function to get a cope data.frame for all level 3 models
# handles the per-subject cope numbering problem
get_l2_cope_df <- function(gpa, model_set, subj_df) {
  checkmate::assert_data_frame(model_set)
  data.table::rbindlist(
    lapply(unique(model_set$l2_model), function(mm) {
      if (!is.null(gpa$l2_models$models[[mm]]$by_subject)) {
        # combine as single data frame from nested list columns
        l2_df <- rbindlist(gpa$l2_models$models[[mm]]$by_subject$cope_list)
        l2_df$l2_model <- mm #retain model name
      } else {
        cope_names <- rownames(gpa$l2_models$models[[mm]]$contrasts)
        l2_df <- data.frame(
          l2_model = mm, l2_cope_number = seq_along(cope_names),
          l2_cope_name = cope_names
        )
        l2_df <- l2_df %>% tidyr::crossing(subj_df)
      }
      return(l2_df)

    })
  )
}

# helper function to get a cope data.frame for all level 3 models
get_l3_cope_df <- function(gpa, model_set, subj_df) {
  checkmate::assert_data_frame(model_set)
  dt <- data.table::rbindlist(
    lapply(unique(model_set$l3_model), function(mm) {
      cope_names <- rownames(gpa$l3_models$models[[mm]]$contrasts)
      data.frame(
        l3_model = mm,
        l3_cope_number = seq_along(cope_names),
        l3_cope_name = cope_names
      )
    })
  )

  # return subject-specific rows for each cope
  dt <- dt %>% tidyr::crossing(subj_df)
  return(dt)

}


# data.table cross-join (tidyr::crossing is a bit faster and already exists)
# https://stackoverflow.com/questions/10600060/how-to-do-cross-join-in-r
# CJ.table <- function(X, Y) {
#   setkey(X[, c(k = 1, .SD)], k)[Y[, c(k = 1, .SD)], allow.cartesian = TRUE][, k := NULL]
# }

get_fsl_l3_model_df <- function(gpa, model_df, subj_df) {
  model_df$model_id <- seq_len(nrow(model_df))

  l1_df <- get_l1_cope_df(gpa, model_df, subj_df)

  l3_df <- get_l3_cope_df(gpa, model_df, subj_df)

  if (isTRUE(gpa$multi_run)) {
    # model_df has l1_model, l2_model, l3_model
    l2_df <- get_l2_cope_df(gpa, model_df)

    combined <- model_df %>%
      tidyr::crossing(subj_df) %>%
      left_join(l1_df, by = c("id", "session", "l1_model")) %>%
      left_join(l2_df, by = c("id", "session", "l2_model")) %>%
      left_join(l3_df, by = c("id", "session", "l3_model"))
  } else {
    combined <- model_df %>%
      tidyr::crossing(subj_df) %>%
      left_join(l1_df, by = c("id", "session", "l1_model")) %>%
      left_join(l3_df, by = c("id", "session", "l3_model"))
  }

  return(combined)
}

get_feat_l3_inputs <- function(gpa, l3_cope_config, lg=NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_data_frame(l3_cope_config)

  if (is.null(lg)) { lg <- lgr::get_logger() }
  if (isTRUE(gpa$multi_run)) {
    # feat directories in $l2_model_setup
    feat_inputs <- gpa$l2_model_setup$fsl %>%
      dplyr::filter(feat_complete == TRUE)

    #join up combination of all models with cope directories
    feat_inputs <- feat_inputs %>%
      dplyr::inner_join(l3_cope_config, by=c("id", "session", "l1_model", "l2_model"))

    # sort out expected cope files for each model combination
    feat_inputs <- feat_inputs %>%
      dplyr::mutate(cope_file = file.path(
        feat_dir,
        paste0("cope", l1_cope_number, ".feat"),
        "stats",
        paste0("cope", l2_cope_number, ".nii.gz")
      )) %>%
      dplyr::select(
        id, session, l1_model, l2_model, l3_model, l1_cope_name, l2_cope_name, feat_dir, cope_file
      )

    split_on <- c("l1_cope_name", "l2_cope_name", "l1_model", "l2_model", "l3_model")
  } else {
    # feat directories in $l1_model_setup
    feat_inputs <- gpa$l1_model_setup$fsl %>%
      dplyr::filter(feat_complete == TRUE)

    feat_inputs <- feat_inputs %>%
      dplyr::left_join(l3_cope_config, by = c("id", "session", "l1_model"))

    feat_inputs <- feat_inputs %>%
      dplyr::mutate(cope_file = file.path(
        feat_dir,
        "stats",
        paste0("cope", l1_cope_number, ".nii.gz")
      )) %>%
      dplyr::select(
        id, session, l1_model, l3_model, l1_cope_name, feat_dir, cope_file
      )

    split_on <- c("l1_cope_name", "l1_model", "l3_model")

  }

  if (!is.data.table(feat_inputs)) data.table::setDT(feat_inputs)
  feat_inputs <- split(feat_inputs, by = split_on)

  return(feat_inputs)
}