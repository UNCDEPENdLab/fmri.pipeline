#' Internal function to add, edit, rename, delete contrasts
#'
#' @param mobj \code{l1_model_spec} or \code{hi_model_spec} object
#' @param signals a list of signals that provide information about regressor specification
#' @param spec_list a list of settings for this model based on an external YAML/JSON specification file
#' @keywords internal
#' @author Michael Hallquist
#' @importFrom magrittr %>%
#' @importFrom emmeans emmeans emtrends
#' @importFrom lgr get_logger
#' @importFrom checkmate assert_multi_class assert_integerish
specify_contrasts <- function(mobj = NULL, signals = NULL, spec_list = NULL) {
  checkmate::assert_multi_class(mobj, c("l1_model_spec", "l1_wi_spec", "hi_model_spec")) # verify that we have an object of known structure
  checkmate::assert_list(signals, null.ok = TRUE)

  checkmate::assert_integerish(mobj$level, lower = 1L, upper = 3L, len = 1L)
  lg <- lgr::get_logger(paste0("glm_pipeline/l", mobj$level, "_setup"))

  # ideally, we should expect a fitted model object to be passed in to support condition means and differences
  # default values for different contrast selections
  cond_means <- c() # factors for which we want model-predicted condition means
  pairwise_diffs <- c() # factors for which we want pairwise differences
  cell_means <- FALSE # whether to include cell means across all factors in design
  overall_response <- FALSE # whether to include predicted average response across all predictors
  simple_slopes <- list() # continuous variables over which to test slopes
  weights <- "cells"
  cat_vars <- c()
  cont_vars <- c()
  delete <- c() # contrasts to be deleted based on user specification or redundancy

  # subfunction to populate contrast_spec and create $contrast_list from spec
  # uses lexical scope of specify_contrasts
  populate_mobj_contrasts <- function(mobj) {
    mobj$contrast_spec$diagonal <- include_diagonal
    mobj$contrast_spec$cond_means <- cond_means
    mobj$contrast_spec$pairwise_diffs <- pairwise_diffs
    mobj$contrast_spec$cell_means <- cell_means
    mobj$contrast_spec$overall_response <- overall_response
    mobj$contrast_spec$simple_slopes <- simple_slopes
    mobj$contrast_spec$weights <- weights
    mobj$contrast_spec$cat_vars <- cat_vars
    mobj$contrast_spec$regressors <- mobj$regressors
    mobj$contrast_spec$delete <- delete
    mobj <- get_contrasts_from_spec(mobj)
  }

  # list used internally for caching generated contrasts
  if (is.null(mobj$contrast_list)) {
    mobj$contrast_list <- list()
  }

  prompt_contrasts <- FALSE
  if (!is.null(spec_list)) {
    if (!is.null(spec_list$diagonal)) {
      checkmate::assert_logical(spec_list$diagonal, len = 1L)
      include_diagonal <- spec_list$diagonal
    } else {
      lg$debug("In contrast specification from spec file, including diagonal contrasts by default")
      include_diagonal <- TRUE
    }

    if (!is.null(spec_list$cond_means)) {
      checkmate::assert_character(spec_list$cond_means)
      cond_means <- spec_list$cond_means
    }

    if (!is.null(spec_list$pairwise_diffs)) {
      checkmate::assert_character(spec_list$pairwise_diffs)
      pairwise_diffs <- spec_list$pairwise_diffs
    }

    if (!is.null(spec_list$cell_means)) {
      checkmate::assert_logical(spec_list$cell_means, len = 1L)
      cell_means <- spec_list$cell_means
    }

    if (!is.null(spec_list$overall_response)) {
      checkmate::assert_logical(spec_list$overall_response, len = 1L)
      overall_response <- spec_list$overall_response
    }

    if (!is.null(spec_list$simple_slopes)) {
      # checkmate::assert_character(spec_list$simple_slopes)
      simple_slopes <- spec_list$simple_slopes
    }

    if (!is.null(spec_list$weights)) {
      checkmate::assert_string(spec_list$weights)
      checkmate::assert_subset(spec_list$weights, c("equal", "cells", "proportional", "flat"))
      weights <- spec_list$weights
    } else {
      lg$debug("In contrast specficiation from spec file, setting emmeans weights to 'cells' by default")
      weights <- "cells"
    }

    if (!is.null(spec_list$delete)) {
      checkmate::assert_character(spec_list$delete)
      delete <- spec_list$delete
    }

    # populate categorical vars for within-subject contrast and higher-level models
    if (inherits(mobj, c("l1_wi_spec", "hi_model_spec"))) {
      term_labels <- attr(terms(mobj$lmfit), "term.labels") # names of regressors
      dclass <- attr(mobj$lmfit$terms, "dataClasses")[term_labels] # omit response variable
      cat_vars <- names(dclass[which(dclass %in% c("factor", "character"))])
    }

    # now handle population of within-subject contrasts (wi_models)

    if (!is.null(spec_list$wi_contrasts)) {
      checkmate::assert_list(spec_list$wi_contrasts)

      for (ff in seq_along(spec_list$wi_contrasts)) {
        sig <- names(spec_list$wi_contrasts)[ff]

        wi_obj <- list(
          level = 1L, lmfit = signals[[sig]]$wi_model, signal_name = sig,
          wi_factors = signals[[sig]]$wi_factors, wi_formula = signals[[sig]]$wi_formula
        )
        class(wi_obj) <- c("list", "l1_wi_spec")

        mobj$wi_models[[sig]] <- specify_contrasts(wi_obj, spec_list = spec_list$wi_contrasts[[ff]])
      }
    }

    mobj <- populate_mobj_contrasts(mobj)

  } else if (is.null(mobj$contrast_spec)) {
    mobj$contrast_spec <- list()
    prompt_contrasts <- TRUE
  } else {
    cat("\nCurrent model contrasts:\n\n")
    print(mobj$contrasts)
    cat("\n--------\n\n")

    res <- menu(c("No", "Yes"), title = "Do you want to modify model contrasts?")
    if (res == 2L) {
      prompt_contrasts <- TRUE
      mobj$contrast_list <- mobj$contrast_spec <- list() # clear prior contrasts that were calculated from the specification
    }
  }

  # if the user does not want to modify existing contrasts, return object unchanged
  if (isFALSE(prompt_contrasts)) {
    lg$debug("Returning model object without changing contrasts")
    return(mobj)
  }

  # handle diagonal contrasts and L1 within-factor modulation of signals
  if (inherits(mobj, "l1_wi_spec")) {
    include_diagonal <- FALSE # for within-subject factor, no need to re-ask about diagonal contrasts for a given factor
    cat("\n\n---------\n")
    cat(glue("The signal {mobj$signal_name} is modulated by the following within-subject factors: {c_string(mobj$wi_factors)}.\n", .trim=FALSE))
    cat(glue("  Please specify contrasts for this signal based on this model: {deparse(mobj$wi_formula, width.cutoff=500)}.\n\n", .trim=FALSE))
  } else if (mobj$level == 1L && inherits(mobj, "l1_model_spec")) {
    has_wi_contrasts <- FALSE # default overridden below if relevant
    # in case of level 1 object, look for any within-subject factors and prompt for relevant contrasts of each
    if (is.null(signals)) {
      lg$warn("An L1 model has been passed to specify_contrasts without a signals list. Cannot check for within-subject factors.")
    } else {
      wi_factors <- sapply(mobj$signals, function(x) !is.null(signals[[x]]$wi_factors))
      if (any(wi_factors)) {
        # setup models for each wi factor
        for (ss in names(wi_factors)[wi_factors]) {
          wi_obj <- list(
            level = 1L, lmfit = signals[[ss]]$wi_model, signal_name = ss,
            wi_factors = signals[[ss]]$wi_factors, wi_formula = signals[[ss]]$wi_formula
          )
          class(wi_obj) <- c("list", "l1_wi_spec")
          mobj$wi_models[[ss]] <- specify_contrasts(wi_obj)
        }
        has_wi_contrasts <- sapply(mobj$wi_models, function(m) {
          if (is.null(m$contrasts) || nrow(m$contrasts) == 0L) FALSE else TRUE
        })
      }
    }

    # add L1 diagonal contrasts automatically if there are no wi_factors or the user didn't ask for wi_contrasts
    if (!any(has_wi_contrasts)) {
      lg$info("Automatically adding diagonal contrasts")
      include_diagonal <- TRUE
    } else {
      # if we have wi contrasts, then prompt whether we also want overall diagonal contrasts
      include_diagonal <- menu(c("Yes", "No"), title = "Do you want to include diagonal contrasts for each regressor?")
      include_diagonal <- ifelse(include_diagonal == 1L, TRUE, FALSE)
    }
  } else {
    # at l2 and l3, ask about diagonal contrasts    
    include_diagonal <- menu(c("Yes", "No"), title = "Do you want to include diagonal contrasts for each regressor?")
    include_diagonal <- ifelse(include_diagonal == 1L, TRUE, FALSE)
  }

  format_prompt <- function(vars, signal_name=NULL) {
    if (is.null(signal_name)) {
      return(vars)
    } else {
      return(paste(vars, "modulation of", signal_name))
    }
  }

  # if we are working from a model specification in lm(), ask about factor means and pairwise comparisons
  if (is.null(mobj$lmfit)) {
    lg$debug("No $lmfit object present in mobj within specify_contrasts.")
  } else {
    # these apply to L2 and L3 model specification
    term_labels <- attr(terms(mobj$lmfit), "term.labels") # names of regressors
    int_terms <- attr(terms(mobj$lmfit), "order") # interactions in model
    dclass <- attr(mobj$lmfit$terms, "dataClasses")[term_labels] # omit response variable
    cat_vars <- names(dclass[which(dclass %in% c("factor", "character"))])
    cont_vars <- names(dclass[which(dclass %in% c("integer", "numeric"))])

    cat(
      "First, please decide on your weighting scheme for computing contrasts.",
      "For details, see the weights section of ?emmeans.",
      "In general, we recommend 'cell' weights in unbalanced designs to",
      "  approximate ordinary marginal means.\n",
      sep = "\n"
    )
    ii <- menu(c(
      "Equal (all cells weighted the same)",
      "Cell (frequencies of cells being averaged)",
      "Proportional (frequencies of factor combinations)",
      "Flat (frequencies of cells being averaged, ignoring empty cells)"
    ),
    title = paste0("How should contrasts be computed by emmeans?")
    )
    weights <- switch(EXPR = ii,
      `1` = "equal",
      `2` = "cells",
      `3` = "proportional",
      `4` = "flat"
    )
    if (is.null(weights)) weights <- "cells"

    # handle model-predicted means for each level of each factor (marginalizing over any other factors)
    for (vv in seq_along(cat_vars)) {
      title_str <- paste0(
        "Do you want to include model-predicted means for ",
        ifelse(inherits(mobj, "l1_wi_spec"), paste(cat_vars[vv], "modulation of", mobj$signal_name), cat_vars[vv]), "?"
      )

      ii <- menu(c("Yes", "No"), title = title_str)
      ii <- ifelse(ii == 1L, TRUE, FALSE)
      if (isTRUE(ii)) {
        cond_means <- c(cond_means, cat_vars[vv])
      }
    }

    # for categorical variables and their interactions, ask about pairwise differences, computed by emmeans
    emm_vars <- cat_vars
    if (any(int_terms > 1)) {
      cat_int <- sapply(term_labels[which(int_terms > 1)], function(tt) {
        term_split <- strsplit(tt, ":")[[1]] # split terms on colon
        all(term_split %in% cat_vars)
      })

      if (any(cat_int)) {
        emm_vars <- c(emm_vars, names(cat_int)[cat_int])
      }
    }
    
    # handle pairwise differences for all factors (and interactions, if specified)
    for (vv in seq_along(emm_vars)) {
      title_str <- paste0(
        "Do you want to include pairwise differences for ",
        ifelse(inherits(mobj, "l1_wi_spec"), paste(emm_vars[vv], "modulation of", mobj$signal_name), emm_vars[vv]), "?"
      )
      ii <- menu(c("Yes", "No"), title = title_str)
      ii <- ifelse(ii == 1L, TRUE, FALSE)
      if (isTRUE(ii)) {
        pairwise_diffs <- c(pairwise_diffs, emm_vars[vv])
      }
    }

    # handle model-predicted means for each *cell* of a factorial design
    if (length(cat_vars) > 1L) {
      title_str <- glue(
        ifelse(inherits(mobj, "l1_wi_spec"), "For {mobj$signal_name}, do ", "Do "),
        "you want to include model-predicted cell means across the factors: {c_string(cat_vars)}?"
      )
      ii <- menu(c("Yes", "No"), title = title_str)
      cell_means <- ifelse(ii == 1L, TRUE, FALSE)
    }

    title_str <- glue(
      "Do you want to include the overall average response",
      ifelse(inherits(mobj, "l1_wi_spec"), " across the levels of {c_string(cat_vars)}?", "?")
    )
    ii <- menu(c("Yes", "No"), title = title_str)
    overall_response <- ifelse(ii == 1L, TRUE, FALSE)

    # Handle model-predicted slopes across levels of the factors.
    # This should only fire if we have categorical x continuous interactions and
    # we wish to calculate simple slopes.
    if (length(cont_vars) > 0L && any(int_terms > 1)) {
      which_int <- int_terms > 1 # which terms are interactions (order > 1)
      term_split <- strsplit(term_labels[which_int], ":") # split terms on colon
      for (vv in cont_vars) {
        has_vv <- which(sapply(term_split, function(x) any(x == vv))) # does this term have the target var
        for (tt in seq_along(has_vv)) {
          this_term <- term_split[[has_vv[tt]]]
          ii <- menu(c("Yes", "No"), title = paste0(
            "Do you want to include model-predicted simple slopes across levels of the factors ",
            paste(this_term[!this_term == vv], collapse = ", "), "?"
          ))
          ii <- ifelse(ii == 1L, TRUE, FALSE)
          if (isTRUE(ii)) { # build list where top-level elements are continuous variables
            simple_slopes[[vv]] <- c(simple_slopes[[vv]], list(this_term[!this_term == vv]))
          }
        }
      }
    }
  }

  # populate contrast specification object
  mobj <- populate_mobj_contrasts(mobj)

  cmat_base <- matrix(numeric(0), nrow = 0, ncol = length(mobj$regressors), dimnames = list(NULL, mobj$regressors))
  cmat <- mobj$contrasts

  # N.B. For now, hand off full model matrix and contrasts to fMRI fitting function (e.g., feat) and
  # let it handle the rank degeneracy. We may need to come back here and drop bad contrasts and use
  # the reduced model matrix in the output. The challenge is that the contrasts calculated by emmeans
  # may be wrong for any column involving aliasing
  if (!is.null(mobj$aliased_terms)) {
    cat(
      "Aliased terms in model.",
      "We will drop aliased columns from the contrast matrix to match the model matrix.",
      "This may invalidate some contrast vectors. Check these manually, check these carefully!",
      "You will probably need to fix the contrasts by hand!",
      sep = "\n"
    )

    has_alias <- which(apply(cmat[, mobj$aliased_terms, drop=F], 1, function(x) {
      any(abs(x - 0) > 1e-5)
    }))

    if (length(has_alias) > 0L) {
      cat("These contrasts used at least one aliased column. Check them!\n")
      cat(strwrap(paste(names(has_alias), collapse = ", "), 70, exdent = 5), sep = "\n")
    }

    # Not used at present
    cmat_reduce <- cmat[, !colnames(cmat) %in% mobj$aliased_terms, drop=F]
    # write.csv(cmat_reduce, file = "cmat_test_reduce.csv")
  }

  # for a within-subject factor model, no need to request contrast editor since this should be handled at overall model level
  if (inherits(mobj, "l1_wi_spec")) {
    return(mobj)
  }

  # syntax of contrast
  cat("\n\n-----\nContrast editor\n\n")
  cat("When entering contrasts, you have two options:\n\n")
  cat("  1) enter a vector of numeric values, one for each regressor\n     Example: 1 0 -1 0 0\n\n")
  cat("  2) enter name=value pairs for relevant coefficients\n     Example: pe_1h = 1 pe_2h = -1\n\n")

  add_more <- 1
  while (add_more != 4) {
    # refresh model object with each iteration since some things may change
    mobj <- get_contrasts_from_spec(mobj)

    add_more <- menu(c("Add contrast", "Show contrasts", "Delete contrast", "Done with contrast setup"),
      title = "Contrast setup menu"
    )

    if (add_more == 1L) {
      cname <- readline("Enter contrast name: ")
      cat("Columns of contrast matrix (one per regressor)\n")
      cat(strwrap(paste(colnames(cmat_base), collapse = ", "), 70, exdent = 5), sep = "\n")
      centry <- readline("Enter contrast syntax: ")

      # Support column = value named syntax or numeric vector syntax
      if (grepl("=", centry)) {
        # split into pairs
        centry <- gsub("([^\\s]+)\\s*=\\s*([^\\s]+)", "\\1=\\2", centry, perl = TRUE)
        csplit <- do.call(rbind, strsplit(strsplit(centry, "\\s+")[[1L]], "=")) # first element is name, second is value
        cvec <- rep(0, ncol(cmat_base)) %>% setNames(colnames(cmat_base))

        for (ii in seq_len(nrow(csplit))) {
          if (!csplit[ii, 1] %in% colnames(cmat_base)) {
            lg$warn("Cannot find column: %s in contrast matrix. Ignoring contrast input.", csplit[ii, 1])
            cvec <- NULL # do not add
            break
          } else {
            cvec[csplit[ii, 1]] <- as.numeric(csplit[ii, 2])
          }
        }
      } else {
        cvec <- sapply(strsplit(centry, "(,\\s*|\\s+)")[[1L]], as.numeric)
        if (length(cvec) != ncol(cmat_base)) {
          lg$warn("Number of elements in entered contrast does not match number of columns in design. Not adding contrast.")
          cvec <- NULL # do not add
        } else {
          names(cvec) <- colnames(cmat_base)
        }
      }

      if (!is.null(cvec)) {
        mobj$contrast_list$custom <- rbind(mobj$contrast_list$custom, cvec)
        rownames(mobj$contrast_list$custom)[nrow(mobj$contrast_list$custom)] <- cname
      }
    } else if (add_more == 2L) {
      summarize_contrasts(mobj$contrasts)
    } else if (add_more == 3L) {
      whichdel <- menu(rownames(mobj$contrasts), title = "Which contrast should be deleted?")
      if (whichdel != 0) {
        mobj$contrast_spec$delete <- c(mobj$contrast_spec$delete, rownames(mobj$contrasts)[whichdel])
      }
    } else if (add_more == 4L) {
      cat("Exiting contrast setup\n")
    }
  }

  if (!is.null(mobj$aliased_terms)) mobj$contrasts_noalias <- cmat_reduce
  return(mobj)
}