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

  if (length(unique(l3_df$l1_model)) > 1L) {
    msg <- "fsl_l3_model is designed for execution on a single l1 model"
    lg$error(msg)
    stop(msg)
  }

  if (isTRUE(gpa$multi_run)) {
    if (length(unique(l3_df$l2_model)) > 1L) {
      msg <- "fsl_l3_model is designed for execution on a single l2 model"
      lg$error(msg)
      stop(msg)
    }
  }

  if (length(unique(l3_df$l3_model)) > 1L) {
    msg <- "fsl_l3_model is designed for execution on a single l3 model"
    lg$error(msg)
    stop(msg)
  }

  # elements of metadata for l3
  id <- l3_df$id[1L]
  session <- l3_df$session[1L]
  l1_model <- l3_df$l1_model[1L]
  l2_model <- l3_df$l2_model[1L]
  l3_model <- l3_df$l3_model[1L]
  l1_cope_name <- l3_df$l1_cope_name[1L]
  l2_cope_name <- l3_df$l2_cope_name[1L]

  # tracking data frame for this model
  feat_l3_df <- data.frame(
     id = id, session = session,
     l1_model=l1_model, l3_model = l3_model, l2_model = l2_model
  )

  # we need to regenerate the l3 model for the inputs provided
  # l3_df should contain FEAT copes that have been vetted in setup_l3_models.R to exist and be complete
  # merge input data against data used to fit the l3 model at the time of initial specification
  to_model <- l3_df %>%
    dplyr::select(id, session) %>%
    left_join(gpa$l3_models$models[[l3_model]]$model_data, by = c("id", "session"))

  mobj <- respecify_l3_model(gpa$l3_models$models[[l3_model]], to_model)

  # generate FSL EV syntax for these regressors
  ev_syntax <- generate_fsf_ev_syntax(inputs = l3_df$cope_file, dmat = mobj$model_matrix)

  # generate FSF contrast syntax for this setup
  contrast_syntax <- generate_fsf_contrast_syntax(mobj$contrasts)

  l3_fsf_syntax <- readLines(system.file("feat_lvl3_copefiles_template.fsf", package = "fmri.pipeline"))

  # Add EVs and contrasts into FSF
  l3_fsf_syntax <- c(l3_fsf_syntax, ev_syntax, contrast_syntax)

  # need to determine number of copes (contrasts) at level 1, which depends on the model being fit
  # FSL usually reads this from the .feat directories itself, but for batch processing, better to insert into the FSF ourselves
  # Need to put this just after the high pass filter cutoff line for Feat to digest it happily

  n_l3_models <- length(gpa$l3_models$models)
  n_l1_models <- length(gpa$l1_models$models)

  # TODO: shift away from group_output_directory in favor of broader localization of outputs
  if (isTRUE(gpa$multi_run)) {
    n_l2_models <- length(gpa$l2_models$models)

    if (n_l1_models > 1L && n_l2_models > 1L) {
      # structure as L1cope/L1model/L2cope/L2model/L3fsf
      l3_outdir <- file.path(
        gpa$group_output_directory, make.names(l1_cope_name), make.names(l1_model),
        make.names(l2_cope_name), make.names(l2_model)
      )
    } else if (n_l1_models > 1L) {
      # structure as L1cope/L1model/L2cope/L3fsf
      l3_outdir <- file.path(
        gpa$group_output_directory, make.names(l1_cope_name), make.names(l1_model),
        make.names(l2_cope_name)
      )
    } else if (n_l2_models > 1L) {
      # structure as L1cope/L2cope/L2model/L3fsf
      l3_outdir <- file.path(
        gpa$group_output_directory, make.names(l1_cope_name),
        make.names(l2_cope_name), make.names(l2_model)
      )
    } else {
      #structure as L1cope/L2cope/L3fsf
      l3_outdir <- file.path(gpa$group_output_directory, make.names(l1_cope_name), make.names(l2_cope_name))
    }
    
    l3_feat_dir <- file.path(
      l3_outdir,
      paste0(
        "FEAT_L3-", make.names(l3_model), "_L1COPE-", make.names(l1_cope_name),
        "_L2COPE-", make.names(l2_cope_name), ".gfeat"
      )
    )
  } else {
    l3_outdir <- file.path(gpa$group_output_directory, l1_model)
    if (n_l1_models > 1L) {
      # structure as L1cope/L1model/L3fsf
      l3_outdir <- file.path(gpa$group_output_directory, make.names(l1_cope_name), make.names(l1_model))
    } else {
      # structure as L1cope/L3fsf
      l3_outdir <- file.path(gpa$group_output_directory, make.names(l1_cope_name))
    }

    l3_feat_dir <- file.path(
      l3_outdir,
      paste0("FEAT_L3-", make.names(l3_model), "_L1COPE-", make.names(l1_cope_name), ".gfeat")
    )
  }

  l3_feat_fsf <- sub(".gfeat$", ".fsf", l3_feat_dir)

  if (!dir.exists(l3_outdir)) {
    lg$debug("Creating L3 output directory: %s", l3_outdir)
    dir.create(l3_outdir, recursive = TRUE)
  }

  # add columns regarding whether inputs already exist and FEAT is already complete
  feat_l3_df <- feat_l3_df %>%
    dplyr::bind_cols(get_feat_status(feat_dir = l3_feat_dir, feat_fsf = l3_feat_fsf, lg = lg))

  feat_l3_df$to_run <- !feat_l3_df$feat_complete

  lg$debug("Expected L3 feat directory is: %s", l3_feat_dir)
  lg$debug("Expected L3 feat fsf is: %s", l3_feat_fsf)

  # specify output directory (removing .gfeat suffix)
  # .OUTPUTDIR. : the feat output location
  l3_fsf_syntax <- gsub(".OUTPUTDIR.", sub("\\.gfeat$", "", l3_feat_dir), l3_fsf_syntax, fixed = TRUE)

  # handle custom L3 FSF syntax
  l3_fsf_syntax <- add_custom_feat_syntax(l3_fsf_syntax, gpa$additional$feat_l3_args, lg)

  # skip re-creation of FSF and do not run below unless force==TRUE
  if (!file.exists(l3_feat_fsf) || isTRUE(gpa$glm_settings$fsl$force_l3_creation)) {
    lg$info("Writing L3 FSF syntax to: %s", l3_feat_fsf)
    cat(l3_fsf_syntax, file = l3_feat_fsf, sep = "\n")
  } else {
    lg$info("Skipping existing L3 FSF syntax: %s", l3_feat_fsf)
  }

  # not currently supporting l3 execution here
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

  return(feat_l3_df)

}
