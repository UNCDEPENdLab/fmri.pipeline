longitudinal_l2_scopes <- function() {
  c("id_session", "id")
}

 longitudinal_l3_input_modes <- function() {
  c(
    "per_session",
    "pooled_sessions_subject_ev",
    "subject_rows",
    "3dlmer"
  )
}

normalize_l3_input_mode <- function(x, default = "per_session") {
  allowed <- longitudinal_l3_input_modes()
  checkmate::assert_string(default)
  checkmate::assert_subset(default, allowed)

  if (is.null(x)) return(default)
  if (!is.character(x)) x <- as.character(x)
  x <- trimws(x)
  x[is.na(x) | !nzchar(x)] <- default
  bad <- !x %in% allowed
  if (any(bad)) {
    stop(
      "Unknown l3_input_mode: ", paste(unique(x[bad]), collapse = ", "),
      ". Allowed: ", paste(allowed, collapse = ", "),
      call. = FALSE
    )
  }
  x
}

normalize_longitudinal_model_signatures <- function(gpa, lg = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  if (is.null(lg)) lg <- lgr::get_logger()
  checkmate::assert_class(lg, "Logger", null.ok = TRUE)

  l2_allowed <- longitudinal_l2_scopes()
  l3_allowed <- longitudinal_l3_input_modes()

  if (checkmate::test_class(gpa$l2_models, "hi_model_set")) {
    for (mname in names(gpa$l2_models$models)) {
      mm <- gpa$l2_models$models[[mname]]
      scope <- mm$l2_scope
      if (is.null(scope) || !is.character(scope) || length(scope) != 1L || !nzchar(scope)) {
        scope <- "id_session"
      }
      if (!scope %in% l2_allowed) {
        stop(
          "Unknown l2_scope for L2 model '", mname, "': ", scope,
          ". Allowed: ", paste(l2_allowed, collapse = ", "),
          call. = FALSE
        )
      }
      mm$l2_scope <- scope
      gpa$l2_models$models[[mname]] <- mm
    }
  }

  if (checkmate::test_class(gpa$l3_models, "hi_model_set")) {
    for (mname in names(gpa$l3_models$models)) {
      mm <- gpa$l3_models$models[[mname]]
      mode <- normalize_l3_input_mode(mm$l3_input_mode)
      if (!mode %in% l3_allowed) {
        stop(
          "Unknown l3_input_mode for L3 model '", mname, "': ", mode,
          ". Allowed: ", paste(l3_allowed, collapse = ", "),
          call. = FALSE
        )
      }
      mm$l3_input_mode <- mode
      gpa$l3_models$models[[mname]] <- mm
    }
  }

  gpa
}

l2_l3_signature_compatible <- function(l2_scope, l3_input_mode) {
  checkmate::assert_string(l2_scope)
  l3_input_mode <- normalize_l3_input_mode(l3_input_mode)

  if (l3_input_mode == "subject_rows") return(l2_scope == "id")
  if (l3_input_mode %in% c("per_session", "pooled_sessions_subject_ev", "3dlmer")) return(l2_scope == "id_session")
  FALSE
}

l2_l3_signature_incompatibility_reason <- function(l2_scope, l3_input_mode) {
  checkmate::assert_string(l2_scope)
  l3_input_mode <- normalize_l3_input_mode(l3_input_mode)

  if (l2_l3_signature_compatible(l2_scope, l3_input_mode)) return("")

  if (l3_input_mode == "subject_rows") {
    return("l3_input_mode='subject_rows' requires l2_scope='id'")
  }

  if (l3_input_mode %in% c("per_session", "pooled_sessions_subject_ev", "3dlmer")) {
    return(paste0("l3_input_mode='", l3_input_mode, "' requires l2_scope='id_session'"))
  }

  "Unknown L2/L3 signature mismatch"
}

enumerate_l2_l3_signature_pairs <- function(gpa, l2_model_names = NULL, l3_model_names = NULL, lg = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_character(l2_model_names, null.ok = TRUE)
  checkmate::assert_character(l3_model_names, null.ok = TRUE)
  if (is.null(lg)) lg <- lgr::get_logger()
  checkmate::assert_class(lg, "Logger", null.ok = TRUE)

  gpa <- normalize_longitudinal_model_signatures(gpa, lg)

  if (!checkmate::test_class(gpa$l2_models, "hi_model_set") ||
      !checkmate::test_class(gpa$l3_models, "hi_model_set")) {
    return(data.frame(
      pair_id = character(0),
      l2_model = character(0),
      l3_model = character(0),
      l2_scope = character(0),
      l3_input_mode = character(0),
      compatible = logical(0),
      reason = character(0),
      stringsAsFactors = FALSE
    ))
  }

  l2_all <- names(gpa$l2_models$models)
  l3_all <- names(gpa$l3_models$models)
  if (is.null(l2_model_names)) l2_model_names <- l2_all
  if (is.null(l3_model_names)) l3_model_names <- l3_all

  checkmate::assert_subset(l2_model_names, l2_all)
  checkmate::assert_subset(l3_model_names, l3_all)

  model_set <- expand.grid(l2_model = l2_model_names, l3_model = l3_model_names, stringsAsFactors = FALSE)
  model_set$l2_scope <- vapply(
    model_set$l2_model,
    function(mm) gpa$l2_models$models[[mm]]$l2_scope,
    character(1)
  )
  model_set$l3_input_mode <- vapply(
    model_set$l3_model,
    function(mm) gpa$l3_models$models[[mm]]$l3_input_mode,
    character(1)
  )
  model_set$compatible <- mapply(
    l2_l3_signature_compatible,
    l2_scope = model_set$l2_scope,
    l3_input_mode = model_set$l3_input_mode,
    SIMPLIFY = TRUE
  )
  model_set$reason <- mapply(
    l2_l3_signature_incompatibility_reason,
    l2_scope = model_set$l2_scope,
    l3_input_mode = model_set$l3_input_mode,
    SIMPLIFY = TRUE
  )
  model_set$pair_id <- paste(model_set$l2_model, model_set$l3_model, sep = "::")
  model_set <- model_set[, c(
    "pair_id", "l2_model", "l3_model", "l2_scope",
    "l3_input_mode", "compatible", "reason"
  ), drop = FALSE]

  model_set
}

format_l2_l3_incompatibilities <- function(pair_df, max_items = 8L) {
  checkmate::assert_data_frame(pair_df)
  checkmate::assert_subset(c("l2_model", "l3_model", "l2_scope", "l3_input_mode", "reason"), names(pair_df))
  checkmate::assert_count(max_items, positive = TRUE)

  if (nrow(pair_df) == 0L) return("")

  if ("compatible" %in% names(pair_df)) {
    pair_df <- pair_df[!pair_df$compatible, , drop = FALSE]
  }

  if (nrow(pair_df) == 0L) return("")

  n_show <- min(nrow(pair_df), max_items)
  detail <- sprintf(
    "%s/%s (l2_scope=%s, l3_input_mode=%s; %s)",
    pair_df$l2_model[seq_len(n_show)],
    pair_df$l3_model[seq_len(n_show)],
    pair_df$l2_scope[seq_len(n_show)],
    pair_df$l3_input_mode[seq_len(n_show)],
    pair_df$reason[seq_len(n_show)]
  )

  suffix <- if (nrow(pair_df) > n_show) {
    paste0("; ... +", nrow(pair_df) - n_show, " more incompatible pair(s)")
  } else {
    ""
  }

  paste0(paste(detail, collapse = "; "), suffix)
}

resolve_l2_l3_compatible_pairs <- function(gpa, l2_model_names = NULL, l3_model_names = NULL, lg = NULL) {
  pair_df <- enumerate_l2_l3_signature_pairs(
    gpa = gpa,
    l2_model_names = l2_model_names,
    l3_model_names = l3_model_names,
    lg = lg
  )
  if (is.null(lg)) lg <- lgr::get_logger()
  checkmate::assert_class(lg, "Logger", null.ok = TRUE)

  incompatible <- pair_df[!pair_df$compatible, , drop = FALSE]
  if (nrow(incompatible) > 0L) {
    lg$info(
      "Dropping %d incompatible L2+L3 pairs based on signatures.",
      nrow(incompatible)
    )
    lg$debug(
      "Dropped pairs: %s",
      paste(sprintf(
        "%s/%s (l2_scope=%s, l3_input_mode=%s; %s)",
        incompatible$l2_model, incompatible$l3_model,
        incompatible$l2_scope, incompatible$l3_input_mode, incompatible$reason
      ), collapse = "; ")
    )
  }

  model_set <- pair_df[pair_df$compatible, c("pair_id", "l2_model", "l3_model", "l2_scope", "l3_input_mode"), drop = FALSE]

  if (nrow(model_set) == 0L) {
    lg$warn("No compatible L2+L3 model pairs remain after applying signature rules.")
  }

  model_set
}
