#This is a small helper function to validate the glm_model_arguments list structure.
#It adds a few details such as the output directory to make it less burdensome for to setup a pipeline
#N.B. gpa is a shorthand abbreviation for glm_model_arguments, to save typing

#' @param gpa A \code{glm_pipeline_arguments} object setup by \code{setup_glm_pipeline}
#' @importFrom string str_count fixed
finalize_pipeline_configuration <- function(gpa) {

  #gpa$outdir <- sapply(gpa$l1_model_variants, function(x) {
  #  paste0("l1-", paste(x, collapse="-"), #define output directory based on combination of signals requested
  #    ifelse(gpa$use_preconvolve, "-preconvolve", ""),
  #    gpa$model_suffix)
  #})

  #new approach: use internal model names
  gpa$outdir <- names(gpa$l1_models$models)
  if (isTRUE(gpa$use_preconvolve)) { gpa$outdir <- paste0(gpa$outdir, "-preconvolve") } #add suffix if using preconvolution approach

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
  #for now, always create a diagonal contrast matrix so that we start with 1 cope per EV
  gpa$n_l1_copes <- c() #number of level 1 copes per model
  gpa$l1_cope_names <- list() #names of level 1 copes per model

  gpa$n_l1_copes <- sapply(gpa$l1_models$models, function(mm) { nrow(mm$contrasts) }) #names of level 1 copes per model
  gpa$l1_cope_names <- lapply(gpa$l1_models$models, function(mm) { rownames(mm$contrasts) })

  gpa$l1_working_directory <- file.path(gpa$working_directory, gpa$outdir) #temp folder for each analysis variant

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
  if (is.null(gpa$pipeline_home)) { gpa$pipeline_home <- "/proj/mnhallqlab/clock_analysis/fmri/fsl_pipeline" }
  if (is.null(gpa$group_output_directory) || gpa$group_output_directory == "default") {
    gpa$group_output_directory <- file.path(getwd(), "group_analyses", gpa$analysis_name)
  }
  if (is.null(gpa$center_l3_predictors)) { gpa$center_l3_predictors <- TRUE }
  if (is.null(gpa$bad_ids)) { gpa$bad_ids <- c() }
  if (is.null(gpa$scheduler)) { gpa$scheduler <- "slurm" } #HPC batch system
  
  if (is.null(gpa$zthresh)) { gpa$zthresh <- 3.09 }  #1-tailed p=.001 for z stat
  if (is.null(gpa$clustsize)) { gpa$clustsize <- 50 } #arbitrary reasonable lower bound on clusters
  if (is.null(gpa$glm_software)) { gpa$glm_software <- "fsl" } #default to FSL FEAT

  if (is.null(gpa$l1_setup_logfile)) { gpa$l1_setup_logfile <- paste0(names(gpa$l1_models$models), "_l1setup.txt") }
  if (is.null(gpa$l1_execution_logfile)) { gpa$l1_setup_logfile <- paste0(names(gpa$l1_models$models), "_l1execution.txt") }
  
  #remove bad ids before running anything further
  if (!is.null(gpa$bad_ids) && length(gpa$bad_ids) > 0L) {
    gpa$subject_data <- gpa$subject_data %>% filter(! !!sym(gpa$vm["id"]) %in% gpa$bad_ids) #remove bad ids
  }
  
  return(gpa)
}
