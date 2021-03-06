#This is a small helper function to validate the glm_model_arguments list structure.
#It adds a few details such as the output directory to make it less burdensome for to setup a pipeline
#N.B. gpa is a shorthand abbreviation for glm_model_arguments, to save typing

#' @param gpa A \code{glm_pipeline_arguments} object setup by \code{setup_glm_pipeline}
#' @importFrom string str_count fixed
finalize_pipeline_configuration <- function(gpa) {

  #new approach: use internal model names for creating output directories at subject level
  #default to <analysis_name>/<l1_model_name>
  #add suffix if using preconvolution approach
  gpa$l1_models$models <- lapply(gpa$l1_models$models, function(mm) {
    mm$outdir <- file.path(gpa$analysis_name, paste0(mm$name, ifelse(gpa$use_preconvolve, "_preconvolve", "")))
    return(mm)
  })

  if (!is.null(gpa$run_number_regex)) {
    if (stringr::str_count(gpa$run_number_regex, stringr::fixed("(")) != 1L) {
      stop("run_number_regex: ", gpa$run_number_regex,
        " must have exactly one opening parenthesis, denoting start of run number capture")
    }

    if (stringr::str_count(gpa$run_number_regex, stringr::fixed(")")) != 1L) {
      stop("run_number_regex: ", gpa$run_number_regex,
        " must have exactly one closing parenthesis, denoting end of run number capture")
    }
  }
  
  #setup l1 copes, cope names, and contrasts.
  gpa$n_l1_copes <- sapply(gpa$l1_models$models, function(mm) { nrow(mm$contrasts) }) #number of level 1 copes per model
  gpa$l1_cope_names <- lapply(gpa$l1_models$models, function(mm) { rownames(mm$contrasts) }) #names of level 1 copes for each model
  gpa$l1_working_directory <- file.path(gpa$working_directory, gpa$outdir) #temp folder for each analysis variant
  if (is.null(gpa$force_l1_creation)) { gpa$force_l1_creation <- FALSE } #whether to overwrite existing level 1 setup files (e.g., .fsf)

  # ---- PARALLELISM SETUP
  # pipeline_cores: number of cores used in push_pipeline when looping over l1 model variants
  if (is.null(gpa$parallel$pipeline_cores) || gpa$parallel$pipeline_cores == "default") {
    gpa$pipeline_cpus <- length(gpa$l1_models$models) #number of workers to setup at the pipeline level (i.e., over l1 model variants)
  }

  #l1_setup_cores defines how many cores to use when looping over subjects within a single l1 model setup
  if (is.null(gpa$parallel$l1_setup_cores) || gpal$parallel$l1_setup_cores == "default") {
    #default to serial execution within a single l1 model variant in setup_lvl1_models
    gpa$parallel$l1_setup_cores <- 1L
  } else {
    checkmate::assert_integerish(gpa$parallel$l1_setup_cores, lower=1)
  }
  
  if (is.null(gpa$l2_cpus)) { gpa$l2_cpus <- 20 } #number of cores to use in Feat LVL2 analyses (fixed effects combination of runs)

  #TODO: deprecate this -- should not be required when executing as an R package
  if (is.null(gpa$pipeline_home)) { gpa$pipeline_home <- "/proj/mnhallqlab/users/michael/fmri.pipeline" }
  if (is.null(gpa$group_output_directory) || gpa$group_output_directory == "default") {
    gpa$group_output_directory <- file.path(getwd(), "group_analyses", gpa$analysis_name)
  }
  if (is.null(gpa$center_l3_predictors)) { gpa$center_l3_predictors <- TRUE }
  if (is.null(gpa$bad_ids)) { gpa$bad_ids <- c() }
  if (is.null(gpa$scheduler)) { gpa$scheduler <- "slurm" } #HPC batch system
  
  if (is.null(gpa$zthresh)) { gpa$zthresh <- 3.09 }  #1-tailed p=.001 for z stat
  if (is.null(gpa$clustsize)) { gpa$clustsize <- 50 } #arbitrary reasonable lower bound on clusters
  if (is.null(gpa$glm_software)) { gpa$glm_software <- "fsl" } #default to FSL FEAT

  if (is.null(gpa$log_json)) { gpa$log_json <- TRUE } #whether to write JSON log files
  if (is.null(gpa$log_txt)) { gpa$log_txt <- TRUE } #whether to write text log files
  if (is.null(gpa$l1_setup_log)) { gpa$l1_setup_log <- paste0(names(gpa$l1_models$models), "_l1setup") %>% setNames(names(gpa$l1_models$models)) }
  if (is.null(gpa$l1_execution_log)) { gpa$l1_execution_log <- paste0(names(gpa$l1_models$models), "_l1execution") %>% setNames(names(gpa$l1_models$models)) }
  
  #remove bad ids before running anything further
  if (!is.null(gpa$bad_ids) && length(gpa$bad_ids) > 0L) {
    gpa$subject_data <- gpa$subject_data %>% filter(! (!!sym(gpa$vm["id"]) %in% gpa$bad_ids)) #remove bad ids
  }

  #build design matrix default arguments
  if (is.null(gpa$additional$bdm_args)) {
    gpa$additional$bdm_args <- list(baseline_coef_order=2, center_values=TRUE, plot=FALSE, convolve_wi_run=TRUE, output_directory="run_timing")
  }

  #default settings for feat l1
  if (is.null(gpa$additional$feat_l1_args)) {
    gpa$additional$feat_l1_args <- list(feat_l1_zthresh=1.96, feat_l1_pthresh=.05)
  }

  return(gpa)
}
