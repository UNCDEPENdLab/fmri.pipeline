#This is a small helper function to validate the glm_model_arguments list structure.
#It adds a few details such as the output directory to make it less burdensome for to setup a pipeline
#N.B. gma is a shorthand abbreviation for glm_model_arguments, to save typing
finalize_pipeline_configuration <- function(gma) {

  gma$outdir <- sapply(gma$l1_model_variants, function(x) {
    paste0("l1-", paste(x, collapse="-"), #define output directory based on combination of signals requested
      ifelse(gma$use_preconvolve, "-preconvolve", ""),
      gma$model_suffix)
  })

  #setup l1 copes, cope names, and contrasts.
  #for now, always create a diagonal contrast matrix so that we start with 1 cope per EV
  gma$n_l1_copes <- c() #number of level 1 copes per model
  gma$l1_cope_names <- list() #names of level 1 copes per model
  final_l1_cmats <- list() #for holding contrast matrices after setup
  
  #populate model names for l1_model_variants
  if (is.null(names(gma$l1_model_variants))) {
    names(gma$l1_model_variants) <- sapply(gma$l1_model_variants, function(x) { paste(x, collapse="-") })
  } else {
    names(gma$l1_model_variants) <- sapply(1:length(gma$l1_model_variants), function(i) {
      if (names(gma$l1_model_variants)[i] == "") {
        return(paste(gma$l1_model_variants[[i]], collapse="-")) #default name to collapse of individual EVs
      } else {
        return(names(gma$l1_model_variants)[i]) #unchanged user nomenclature for model name
      }
    })
  }
  
  for (ii in 1:length(gma$l1_model_variants)) {
    #generate a diagonal matrix of contrasts
    regressors <- gma$l1_model_variants[[ii]]
    
    cmat <- diag(length(regressors))
    rownames(cmat) <- colnames(cmat) <- regressors

    mname <- names(gma$l1_model_variants)[ii]

    #are there additional l1 contrasts?
    if (!is.null(gma$l1_contrasts[[mname]])) {
      l1_contrasts <- gma$l1_contrasts[[mname]]
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
    gma$n_l1_copes[mname] <- nrow(cmat)
    l1_cope_names <- rownames(cmat)
    names(l1_cope_names) <- paste0("cope", 1:length(l1_cope_names)) #names attribute holds cope numbering scheme in FSL
    gma$l1_cope_names[[mname]] <- l1_cope_names
  }

  #remove user-specified contrasts and supplant them with computed contrast matrices
  #use of final_l1_cmats temp variable ensures that the length of $l1_contrasts is the same as other elements like $l1_model_variants
  gma$l1_contrasts <- final_l1_cmats
  
  gma$workdir <- file.path(gma$root_workdir, gma$outdir) #temp folder for each analysis variant
  gma$pipeline_cpus <- length(gma$l1_model_variants) #number of workers to setup at the pipeline level (i.e., over run variants)
  if (is.null(gma$l2_cpus)) { gma$l2_cpus <- 20 } #number of cores to use in Feat LVL2 analyses (fixed effects combination of runs)
  if (is.null(gma$pipeline_home)) { gma$pipeline_home <- "/proj/mnhallqlab/clock_analysis/fmri/fsl_pipeline" }
  if (is.null(gma$group_output_dir)) { gma$group_output_dir <- file.path(dirname(gma$fmri_dir), "group_analyses", gma$analysis_name) }
  if (is.null(gma$center_l3_predictors)) { gma$center_l3_predictors <- TRUE }
  if (is.null(gma$bad_ids)) { gma$bad_ids <- c() }
  if (is.null(gma$scheduler)) { gma$scheduler <- "slurm" } #HPC batch system
  
    if (is.null(gma$zthresh)) { gma$zthresh <- 3.09 }  #1-tailed p=.001 for z stat
  if (is.null(gma$clustsize)) { gma$clustsize <- 34 } #based on 3dClustSim using ACFs for first-level FEAT runs
  if (is.null(gma$glm_software)) { gma$glm_software <- "fsl" } #default to FSL FEAT
  
  #ensure that the user has specified some sort of clock event in the model
  for (v in gma$l1_model_variants) {
    if (!any(c("clock", "clock_bs") %in% v)) {
      stop("No clock event is in the model: ", paste(v, collapse=","))
    }
    if (!any(c("feedback", "feedback_bs") %in% v)) {
      stop("No feedback event is in the model: ", paste(v, collapse=","))
    }
  }  

  #remove bad ids before running anything further
  if (!is.null(gma$bad_ids) && length(gma$bad_ids) > 0L) {
    gma$subject_data <- gma$subject_data %>% filter(!id %in% gma$bad_ids) #remove bad ids
  }
  
  return(gma)
}
