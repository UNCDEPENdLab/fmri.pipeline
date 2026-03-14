#' Harvest FSL L2 COPEs for 3dLMEr or other group-level backends
#'
#' @param gpa a glm_pipeline_arguments object
#' @param l3_model_names optional character vector of L3 model names to process
#'
#' @return a list of data.frames, one per L1/L2 contrast combination, for each L3 model
#' @keywords internal
harvest_fsl_copes <- function(gpa, l3_model_names = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  lg <- lgr::get_logger("glm_pipeline/harvest_fsl_copes")
  
  if (is.null(l3_model_names)) {
    l3_model_names <- names(gpa$l3_models$models)
  }
  
  # Get the compatible L2/L3 pairs
  pairs <- resolve_l2_l3_compatible_pairs(gpa, l3_model_names = l3_model_names)
  
  if (nrow(pairs) == 0L) {
    lg$warn("No compatible L2/L3 model pairs found for harvesting.")
    return(NULL)
  }
  
  # For each pair, we need to identify the L1 model and then the L1/L2 copes
  # The bridge between L1/L2/L3 cope numbers is usually defined during setup
  
  # We reuse the existing logic in get_fsl_l3_model_df and get_feat_l3_inputs
  # because they already handle the complex mapping of multi-level FSL hierarchies.
  
  all_harvested <- list()
  
  for (i in seq_len(nrow(pairs))) {
    this_l2 <- pairs$l2_model[i]
    this_l3 <- pairs$l3_model[i]
    
    # We need to know which L1 model(s) are associated with this L2 model
    # Usually this is defined in the L2 model object
    l1_model_names <- gpa$l2_models$models[[this_l2]]$l1_model_names
    if (is.null(l1_model_names)) {
      # Fallback: all L1 models in gpa
      l1_model_names <- names(gpa$l1_models$models)
    }
    
    model_df <- expand.grid(
      l1_model = l1_model_names,
      l2_model = this_l2,
      l3_model = this_l3,
      stringsAsFactors = FALSE
    )
    
    # Get the cope configuration (mapping of L1/L2/L3 cope numbers)
    # This involves identify_l1_l2_l3_copes which might be in specify_contrasts or setup_l3
    # For now, let's assume we can use the logic in get_fsl_l3_model_df
    
    # subj_df: only good subjects/sessions
    subj_df <- gpa$subject_data %>%
      dplyr::filter(exclude_subject == FALSE) %>%
      dplyr::select(id, session)
    
    l3_cope_config <- get_fsl_l3_model_df(gpa, model_df, subj_df)
    
    # This config has one row per id/session/l1_model/l2_model/l3_model/l1_cope/l2_cope
    # Now get the actual file paths
    feat_inputs <- get_feat_l3_inputs(gpa, l3_cope_config, lg = lg)
    
    # feat_inputs is a list of data.frames split by (l1_cope_name, l2_cope_name, l1_model, l2_model, l3_model, l3_session_partition)
    # Each data.frame has columns: id, session, cope_file, etc.
    
    # Rename cope_file to InputFile for 3dLMEr consistency
    for (nm in names(feat_inputs)) {
      feat_inputs[[nm]] <- feat_inputs[[nm]] %>%
        dplyr::rename(InputFile = cope_file)
    }
    
    all_harvested[[this_l3]] <- c(all_harvested[[this_l3]], feat_inputs)
  }
  
  return(all_harvested)
}
