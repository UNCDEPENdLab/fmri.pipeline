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
  final_l1_cmats <- list() #for holding contrast matrices after setup
  
  #populate model names for l1_model_variants
  if (is.null(names(gpa$l1_model_variants))) {
    names(gpa$l1_model_variants) <- sapply(gpa$l1_model_variants, function(x) { paste(x, collapse="-") })
  } else {
    names(gpa$l1_model_variants) <- sapply(1:length(gpa$l1_model_variants), function(i) {
      if (names(gpa$l1_model_variants)[i] == "") {
        return(paste(gpa$l1_model_variants[[i]], collapse="-")) #default name to collapse of individual EVs
      } else {
        return(names(gpa$l1_model_variants)[i]) #unchanged user nomenclature for model name
      }
    })
  }
  
  for (ii in 1:length(gpa$l1_model_variants)) {
    #generate a diagonal matrix of contrasts
    regressors <- gpa$l1_model_variants[[ii]]
    
    cmat <- diag(length(regressors))
    rownames(cmat) <- colnames(cmat) <- regressors

    mname <- names(gpa$l1_model_variants)[ii]

    #are there additional l1 contrasts?
    if (!is.null(gpa$l1_contrasts[[mname]])) {
      l1_contrasts <- gpa$l1_contrasts[[mname]]
      #l1_contrasts should be a list of named vectors
      #each element of the list is a new contrast to be added to the current l1 model
      #the names of the elements become the contrast names
      #the names of each vector refer to the non-zero elements of the contrast
      #the values of each vector are the contrast values to be set
      # for example:
      #   list(
      #    pe1h_gt_pe2h=c(pe_1h=1, pe_2h=-1),
      #    entropy1h_gt_entropy2h=c(v_entropy_1h=1, v_entropy_2h=-1)
      #   )
      # will generate two contrasts named pe1h_gt_pe2h and entropy1h_gt_entropy2h
      # the values of first contrast will be 1 for the EV pe_1h and -1 for pe_2h EV
      add_contrasts <- matrix(0, nrow=length(l1_contrasts), ncol=ncol(cmat), dimnames=list(names(l1_contrasts), colnames(cmat)))
      
      for(con in 1:length(l1_contrasts)) {
        add_contrasts[con,names(l1_contrasts[[con]])] <- l1_contrasts[[con]]
      }

      cmat <- rbind(cmat, add_contrasts)
    }

    final_l1_cmats[[mname]] <- cmat
    gpa$n_l1_copes[mname] <- nrow(cmat)
    l1_cope_names <- rownames(cmat)
    names(l1_cope_names) <- paste0("cope", 1:length(l1_cope_names)) #names attribute holds cope numbering scheme in FSL
    gpa$l1_cope_names[[mname]] <- l1_cope_names
  }

  #remove user-specified contrasts and supplant them with computed contrast matrices
  #use of final_l1_cmats temp variable ensures that the length of $l1_contrasts is the same as other elements like $l1_model_variants
  gpa$l1_contrasts <- final_l1_cmats
  
  gpa$workdir <- file.path(gpa$root_workdir, gpa$outdir) #temp folder for each analysis variant
  gpa$pipeline_cpus <- length(gpa$l1_model_variants) #number of workers to setup at the pipeline level (i.e., over run variants)

  #defaults

  #l1_setup_cores defines how many cores to use when looping over subjects within a single l1 model setup
  if (is.null(gpa$l1_setup_cores)) {
    #default to serial execution within a single l1 model variant in setup_lvl1_models
    gpa$l1_setup_cores <- 1L
  } else {
    checkmate::assert_integerish(gpa$l1_setup_cores, lower=1)
  }
  
  if (is.null(gpa$l2_cpus)) { gpa$l2_cpus <- 20 } #number of cores to use in Feat LVL2 analyses (fixed effects combination of runs)
  if (is.null(gpa$pipeline_home)) { gpa$pipeline_home <- "/proj/mnhallqlab/clock_analysis/fmri/fsl_pipeline" }
  if (is.null(gpa$group_output_dir)) { gpa$group_output_dir <- file.path(dirname(gpa$fmri_dir), "group_analyses", gpa$analysis_name) }
  if (is.null(gpa$center_l3_predictors)) { gpa$center_l3_predictors <- TRUE }
  if (is.null(gpa$bad_ids)) { gpa$bad_ids <- c() }
  if (is.null(gpa$scheduler)) { gpa$scheduler <- "slurm" } #HPC batch system
  
  if (is.null(gpa$zthresh)) { gpa$zthresh <- 3.09 }  #1-tailed p=.001 for z stat
  if (is.null(gpa$clustsize)) { gpa$clustsize <- 34 } #based on 3dClustSim using ACFs for first-level FEAT runs
  if (is.null(gpa$glm_software)) { gpa$glm_software <- "fsl" } #default to FSL FEAT
  
  #ensure that the user has specified some sort of clock event in the model
  for (v in gpa$l1_model_variants) {
    if (!any(c("clock", "clock_bs") %in% v)) {
      stop("No clock event is in the model: ", paste(v, collapse=","))
    }
    if (!any(c("feedback", "feedback_bs") %in% v)) {
      stop("No feedback event is in the model: ", paste(v, collapse=","))
    }
  }  

  #remove bad ids before running anything further
  if (!is.null(gpa$bad_ids) && length(gpa$bad_ids) > 0L) {
    gpa$subject_data <- gpa$subject_data %>% filter(!id %in% gpa$bad_ids) #remove bad ids
  }
  
  return(gpa)
}
