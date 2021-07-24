#' Estimate a level 2 (subject) model using FSL FEAT with fixed effects integration of runs
#'
#' @param l3_df a data.frame containing cope inputs for a given l3 model, as well as metadata
#'   that identify the model.
#' @param gpa a \code{glm_pipeline_arguments} object containing model specification
#'
#' @importFrom dplyr mutate filter select right_join pull
#' @author Michael Hallquist
#' @export
fsl_l3_model <- function(l3_df=NULL, gpa) {
  checkmate::assert_data_frame(l3_df)
  checkmate::assert_subset(c("id", "session", "l3_model"), names(l3_df))
  checkmate::assert_class(gpa, "glm_pipeline_arguments")

  lg <- lgr::get_logger("glm_pipeline/l3_setup")

  if (length(unique(l3_df$session)) > 1L) {
    msg <- "fsl_l3_model is designed for execution on a single session data.frame"
    lg$error(msg)
    stop(msg)
  }

  if (length(unique(l3_df$l3_model)) > 1L) {
    msg <- "fsl_l2_model is designed for execution on a single l3 model"
    lg$error(msg)
    stop(msg)
  }

  # elements of metadata for l2
  id <- l3_df$id[1L]
  session <- l3_df$session[1L]
  l3_model <- l3_df$l3_model[1L]
  n_l3_copes <- gpa$n_l3_copes[l3_model] # number of lvl3 copes to combine for this model
  l3_feat_dirs <- l3_df$feat_dir

  # tracking data frame for this model
  feat_l2_df <- data.frame(
    id = id, session = session,
    l3_model = l3_model, l2_model = l2_model_name
  )

    #find rows in run_data that match this subject
    dmat_rows <- gpa$run_data %>%
      dplyr::mutate(rownum = 1:n()) %>%
      dplyr::filter(id == !!id & session == !!session) %>%
      dplyr::filter(exclude_run == FALSE) %>%
      dplyr::pull(rownum)

    if (length(dmat_rows) != nrow(l3_df)) {
      msg <- "Number of rows in gpa$run_data does not match l3_df in fsl_l2_model"
      lg$error(msg)
      stop(msg)
    }

    # should never happen, but sanity check the model matrix against the run data
    stopifnot(nrow(gpa$run_data) == nrow(gpa$l2_models$models[[l2_model_name]]$model_matrix))

    # obtain rows of design for this subject
    dmat <- gpa$l2_models$models[[l2_model_name]]$model_matrix[dmat_rows, , drop=FALSE]

    # obtain contrasts for this L2 model
    cmat <- gpa$l2_models$models[[l2_model_name]]$contrasts # l2 model contrasts

  # generate FSL EV syntax for these regressors
  ev_syntax <- generate_fsf_ev_syntax(inputs = l1_feat_dirs, dmat = dmat)

  # generate FSF contrast syntax for this setup
  contrast_syntax <- generate_fsf_contrast_syntax(cmat)

  l2_fsf_syntax <- readLines(system.file("feat_lvl2_nparam_template.fsf", package = "fmri.pipeline"))

  # Add EVs and contrasts into FSF
  l2_fsf_syntax <- c(l2_fsf_syntax, ev_syntax, contrast_syntax)

  # need to determine number of copes (contrasts) at level 1, which depends on the model being fit
  # FSL usually reads this from the .feat directories itself, but for batch processing, better to insert into the FSF ourselves
  # Need to put this just after the high pass filter cutoff line for Feat to digest it happily

  l2_fsf_syntax <- c(
    l2_fsf_syntax,
    "# Number of lower-level copes feeding into higher-level analysis",
    paste0("set fmri(ncopeinputs) ", n_l1_copes)
  )

  # tell FSL to analyze all lower-level copes in LVL2
  for (n in seq_len(n_l1_copes)) {
    l2_fsf_syntax <- c(
      l2_fsf_syntax,
      paste0("# Use lower-level cope ", n, " for higher-level analysis"),
      paste0("set fmri(copeinput.", n, ") 1"), ""
    )
  }

  # TODO: make output location more flexible and not always relative to L1 outputs
  l2_feat_dir <- file.path(dirname(l1_feat_dirs[1L]), paste0("FEAT_LVL2_", l2_model_name, ".gfeat"))
  l2_feat_fsf <- file.path(dirname(l1_feat_dirs[1L]), paste0("FEAT_LVL2_", l2_model_name, ".fsf"))

  # add columns regarding whether inputs already exist and FEAT is already complete
  feat_l2_df <- feat_l2_df %>%
    dplyr::bind_cols(get_feat_status(feat_dir = l2_feat_dir, feat_fsf = l2_feat_fsf, lg = lg))

  feat_l2_df$to_run <- !feat_l2_df$feat_complete

  lg$debug("Expected L2 feat directory is: %s", l2_feat_dir)
  lg$debug("Expected L2 feat fsf is: %s", l2_feat_fsf)

  # specify output directory (removing .gfeat suffix)
  # .OUTPUTDIR. : the feat output location
  l2_fsf_syntax <- gsub(".OUTPUTDIR.", sub("\\.gfeat$", "", l2_feat_dir), l2_fsf_syntax, fixed = TRUE)

  # handle custom L1 FSF syntax
  l2_fsf_syntax <- add_custom_feat_syntax(l2_fsf_syntax, gpa$additional$feat_l2_args, lg)

  # skip re-creation of FSF and do not run below unless force==TRUE
  if (!file.exists(l2_feat_fsf) || isTRUE(gpa$glm_settings$fsl$force_l2_creation)) {
    lg$info("Writing L2 FSF syntax to: %s", l2_feat_fsf)
    cat(l2_fsf_syntax, file = l2_feat_fsf, sep = "\n")
  } else {
    lg$info("Skipping existing L2 FSF syntax: %s", l2_feat_fsf)
  }

  # not currently supporting l2 execution here
  # if (isTRUE(execute_feat)) {
  #   nnodes <- min(length(all_l1_feat_fsfs), parallel::detectCores())
  #   lg$info("Starting fork cluster with %d workers", nnodes)

  #   cl_fork <- parallel::makeForkCluster(nnodes=ncpus)
  #   runfeat <- function(fsf) {
  #     runname <- basename(fsf)
  #     runFSLCommand(paste("feat", fsf),
  #       stdout = file.path(dirname(fsf), paste0("feat_stdout_", runname)),
  #       stderr = file.path(dirname(fsf), paste0("feat_stderr_", runname))
  #     )
  #     system(paste0("feat_lvl2_to_afni.R --gfeat_dir ", sub(".fsf", ".gfeat", fsf, fixed=TRUE), " --no_subjstats --no_varcope --stat_outfile ", sub(".fsf", "_gfeat_stats", fsf, fixed=TRUE))) #aggregate FEAT statistics into a single file
  #   }
  #   parallel::clusterApply(cl_fork, allFeatRuns, runfeat)
  #   parallel::stopCluster(cl_fork)
  # }

  return(feat_l2_df)

}
