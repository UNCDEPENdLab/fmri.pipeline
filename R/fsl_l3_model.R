#' Estimate a level 2 (subject) model using FSL FEAT with fixed effects integration of runs
#'
#' @param l3_df a data.frame containing cope inputs for a given l3 model, as well as metadata
#'   that identify the model.
#' @param gpa a \code{glm_pipeline_arguments} object containing model specification
#'
#' @importFrom dplyr mutate filter select left_join inner_join pull
#' @author Michael Hallquist
#' @keywords internal
fsl_l3_model <- function(l3_df=NULL, gpa) {
  checkmate::assert_data_frame(l3_df)
  checkmate::assert_subset(c("id", "session", "l3_model"), names(l3_df))
  checkmate::assert_class(gpa, "glm_pipeline_arguments")

  lg <- lgr::get_logger("glm_pipeline/l3_setup")
  lg$set_threshold(gpa$lgr_threshold)

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
  session <- l3_df$session[1L]
  l1_model <- l3_df$l1_model[1L]
  l2_model <- l3_df$l2_model[1L]
  l3_model <- l3_df$l3_model[1L]
  l1_cope_name <- l1_contrast <- l3_df$l1_cope_name[1L] #use the double assign for synonyms used in glue() output expression
  l2_cope_name <- l2_contrast <- l3_df$l2_cope_name[1L]

  # tracking data frame for this model (column names should follow variable names)
  if (isTRUE(gpa$multi_run)) {
    feat_l3_df <- data.frame(l1_model, l1_cope_name, l2_model, l2_cope_name, l3_model)
  } else {
    feat_l3_df <- data.frame(l1_model, l1_cope_name, l3_model)
  }

  # we need to regenerate the l3 model for the inputs provided
  # l3_df should contain FEAT copes that have been vetted in setup_l3_models.R to exist and be complete
  # handle model respecification based on available data (e.g., if some subjects failed to run)
  mobj <- respecify_l3_model(gpa$l3_models$models[[l3_model]], new_data=l3_df)

  # now make sure that the l3_df aligns perfectly with the data being modeled (mostly dropping copes in l3_df that aren't in model)
  l3_df <- l3_df %>% dplyr::inner_join(mobj$metadata, by=c("id", "session"))

  # generate FSL EV syntax for these regressors
  ev_syntax <- fsl_generate_fsf_ev_syntax(inputs = l3_df$cope_file, dmat = mobj$model_matrix)

  # generate FSF contrast syntax for this setup
  contrast_syntax <- fsl_generate_fsf_contrast_syntax(mobj$contrasts)

  l3_fsf_syntax <- readLines(system.file("feat_lvl3_copefiles_template.fsf", package = "fmri.pipeline"))

  # Add EVs and contrasts into FSF
  l3_fsf_syntax <- c(l3_fsf_syntax, ev_syntax, contrast_syntax)

  # need to determine number of copes (contrasts) at level 1, which depends on the model being fit
  # FSL usually reads this from the .feat directories itself, but for batch processing, better to insert into the FSF ourselves
  # Need to put this just after the high pass filter cutoff line for Feat to digest it happily

  #n_l3_models <- length(gpa$l3_models$models)
  #n_l1_models <- length(gpa$l1_models$models)

  l3_outdir <- get_output_directory(
    l1_contrast=l1_cope_name, l1_model=l1_model,
    l2_contrast=l2_cope_name, l2_model=l2_model,
    l3_model=l3_model, gpa = gpa, what = "l3"
  )

  l3_feat_fsf <- file.path(
    l3_outdir,
    glue::glue(gpa$output_locations$feat_l3_fsf) # evaluate glue expression
  )

  l3_feat_dir <- sub(".fsf$", ".gfeat", l3_feat_fsf)

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

  # handle outlier deweighting flag (ifelse won't keep list names...)
  deweight <- if(isTRUE(mobj$fsl_outlier_deweighting)) list(robust_yn = 1) else list(robust_yn = 0)
  l3_fsf_syntax <- add_custom_feat_syntax(l3_fsf_syntax, deweight, lg)

  # skip re-creation of FSF and do not run below unless force==TRUE
  if (!file.exists(l3_feat_fsf) || isTRUE(gpa$glm_settings$fsl$force_l3_creation)) {
    lg$info("Writing L3 FSF syntax to: %s", l3_feat_fsf)
    cat(l3_fsf_syntax, file = l3_feat_fsf, sep = "\n")
  } else {
    lg$info("Skipping existing L3 FSF syntax: %s", l3_feat_fsf)
  }

  return(feat_l3_df)
}
