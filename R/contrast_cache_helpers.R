extract_contrast_names <- function(contrasts, model_name = NULL, level_label = "model", allow_empty = TRUE) {
  checkmate::assert_string(level_label)
  checkmate::assert_flag(allow_empty)

  if (is.null(contrasts)) {
    if (isTRUE(allow_empty)) return(character(0))
    stop(
      sprintf("No %s contrasts available for model '%s'.", level_label, model_name),
      call. = FALSE
    )
  }

  if (checkmate::test_data_frame(contrasts)) {
    contrasts <- as.matrix(contrasts)
  }

  if (!is.matrix(contrasts)) {
    if (isTRUE(allow_empty)) return(character(0))
    stop(
      sprintf(
        "Expected %s contrasts for model '%s' to be a matrix or data.frame.",
        level_label, model_name
      ),
      call. = FALSE
    )
  }

  n_contrasts <- nrow(contrasts)
  if (is.null(n_contrasts) || n_contrasts == 0L) {
    if (isTRUE(allow_empty)) return(character(0))
    stop(
      sprintf("No %s contrasts available for model '%s'.", level_label, model_name),
      call. = FALSE
    )
  }

  contrast_names <- rownames(contrasts)
  if (is.null(contrast_names) || length(contrast_names) != n_contrasts) {
    contrast_names <- paste0("contrast_", seq_len(n_contrasts))
  }

  as.character(contrast_names)
}

refresh_l1_cope_names <- function(gpa, lg = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  if (is.null(lg)) lg <- lgr::get_logger()
  checkmate::assert_class(lg, "Logger", null.ok = TRUE)

  if (!checkmate::test_class(gpa$l1_models, "l1_model_set") ||
      is.null(gpa$l1_models$models) ||
      length(gpa$l1_models$models) == 0L) {
    gpa$l1_cope_names <- list()
    return(gpa)
  }

  existing_cache <- if (checkmate::test_list(gpa$l1_cope_names)) gpa$l1_cope_names else list()
  model_names <- names(gpa$l1_models$models)
  if (is.null(model_names)) {
    model_names <- as.character(seq_along(gpa$l1_models$models))
  }

  l1_cope_names <- setNames(vector("list", length(model_names)), model_names)
  for (mname in model_names) {
    this_model <- gpa$l1_models$models[[mname]]
    model_contrasts <- extract_contrast_names(
      this_model$contrasts,
      model_name = mname,
      level_label = "L1",
      allow_empty = TRUE
    )
    cache_contrasts <- existing_cache[[mname]]
    if (length(model_contrasts) == 0L && !is.null(cache_contrasts) && length(cache_contrasts) > 0L) {
      model_contrasts <- as.character(cache_contrasts)
    }
    l1_cope_names[[mname]] <- model_contrasts
  }

  gpa$l1_cope_names <- l1_cope_names
  gpa
}
