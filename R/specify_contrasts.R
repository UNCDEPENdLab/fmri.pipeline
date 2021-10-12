#' Internal function to add, edit, rename, delete contrasts
#'
#' @param mobj \code{l1_model_spec} or \code{hi_model_spec} object
#' @param include_diagonal A logical indicating whether to include a diagonal matrix of contrasts
#'   in the model contrasts. Defaults to TRUE.
#' @keywords internal
#' @author Michael Hallquist
#' @importFrom magrittr %>%
#' @importFrom emmeans emmeans emtrends
specify_contrasts <- function(mobj = NULL) {
    checkmate::assert_multi_class(mobj, c("l1_model_spec", "hi_model_spec")) # verify that we have an object of known structure

    checkmate::assert_integerish(mobj$level, lower=1L, upper=3L, len=1L)
    lg <- lgr::get_logger(paste0("glm_pipeline/l", mobj$level, "_setup"))

    prompt_contrasts <- FALSE
    if (is.null(mobj$contrast_spec)) {
      mobj$contrast_spec <- list()
      prompt_contrasts <- TRUE
    } else {
      cat("\nCurrent model contrasts:\n\n")
      print(mobj$contrasts)
      cat("\n--------\n\n")

      res <- menu(c("No", "Yes"), title = "Do you want to modify model contrasts?")
      if (res == 2L) {
        prompt_contrasts <- TRUE
      }
    }

    # list used internally for caching generated contrasts
    if (is.null(mobj$contrast_list)) {
      mobj$contrast_list <- list()
    }

    # if the user does not want to modify existing contrasts, return object unchanged
    if (isFALSE(prompt_contrasts)) {
      lg$debug("Returning model object without changing contrasts")
      return(mobj)
    }

    # handle diagonal contrasts
    include_diagonal <- menu(c("Yes", "No"), title = "Do you want to include diagonal contrasts for each regressor?")
    include_diagonal <- ifelse(include_diagonal == 1L, TRUE, FALSE)

    # ideally, we should expect a fitted model object to be passed in to support condition means and differences
    cond_means <- c() # factors for which we want model-predicted condition means
    pairwise_diffs <- c() # factors for which we want pairwise differences
    cell_means <- FALSE # whether to include cell means across all factors in design
    overall_response <- FALSE # whether to include predicted average response across all predictors
    simple_slopes <- list() # continuous variables over which to test slopes
    weights <- "cells"
    cat_vars <- c()
    cont_vars <- c()

    # if we are working from a model specification in lm(), ask about factor means and pairwise comparisons
    if (is.null(mobj$lmfit)) {
      lg$debug("No $lmfit object present in mobj within specify_contrasts.")
    } else {
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
        "Flat (frequencies of cells being averaged, ignoring empty cells)"),
        title = paste0("How should contrasts be computed by emmeans?")
      )
      weights <- switch(EXPR = ii,
        `1` = "equal",
        `2` = "cells",
        `3` = "proportional",
        `4` = "flat"
      )
      if (is.null(weights)) weights <- "cells"

      # handle model-predicted means for each level of each factor
      for (vv in seq_along(cat_vars)) {
          ii <- menu(c("Yes", "No"), title = paste0("Do you want to include model-predicted means for ", cat_vars[vv], "?"))
          ii <- ifelse(ii == 1L, TRUE, FALSE)
          if (isTRUE(ii)) { cond_means <- c(cond_means, cat_vars[vv]) }

          ii <- menu(c("Yes", "No"), title = paste0("Do you want to include pairwise differences for ", cat_vars[vv], "?"))
          ii <- ifelse(ii == 1L, TRUE, FALSE)
          if (isTRUE(ii)) { pairwise_diffs <- c(pairwise_diffs, cat_vars[vv]) }
      }

      # handle model-predicted means for each cell of a factorial design
      if (length(cat_vars) > 1L) {
          ii <- menu(c("Yes", "No"), title = paste0(
              "Do you want to include model-predicted cell means across the factors ",
              paste(cat_vars, collapse = ", "), "?"
          ))
          cell_means <- ifelse(ii == 1L, TRUE, FALSE)
      }

      ii <- menu(c("Yes", "No"), title = "Do you want to include the overall average response?")
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
            this_term <- term_split[[ has_vv[tt] ]]
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

    #populate contrast specification object
    mobj$contrast_spec$diagonal <- include_diagonal
    mobj$contrast_spec$cond_means <- cond_means
    mobj$contrast_spec$pairwise_diffs <- pairwise_diffs
    mobj$contrast_spec$cell_means <- cell_means
    mobj$contrast_spec$overall_response <- overall_response
    mobj$contrast_spec$simple_slopes <- simple_slopes
    mobj$contrast_spec$weights <- weights
    mobj$contrast_spec$cat_vars <- cat_vars
    mobj$contrast_spec$regressors <- mobj$regressors

    cmat_base <- matrix(numeric(0), nrow = 0, ncol = length(mobj$regressors), dimnames = list(NULL, mobj$regressors))
    mobj <- get_contrasts_from_spec(mobj)
    cmat <- mobj$contrasts

    # N.B. For now, hand off full model matrix and contrasts to fMRI fitting function (e.g., feat) and
    # let it handled the rank degeneracy. We may need to come back here and drop bad contrasts and use
    # the reduced model matrix in the output. The challenge is that the contrasts calculated by emmeans
    # may be wrong for any column involving aliasing
    if (!is.null(mobj$aliased_terms)) {
      cat(
        "Aliased terms in model.",
        "We will drop aliased columns from the contrast matrix to match the model matrix.",
        "This may invalidate some contrast vectors. Check these manually, check these carefully!",
        "You will probably need to fix the contrasts by hand!",
        sep="\n"
      )
      has_alias <- which(apply(cmat[, mobj$aliased_terms], 1, function(x) {
        any(abs(x - 0) > 1e-5)
      }))

      if (length(has_alias) > 0L) {
        cat("These contrasts used at least one aliased column. Check them!\n")
        cat(strwrap(paste(names(has_alias), collapse = ", "), 70, exdent = 5), sep="\n")
      }

      # Not used at present
      cmat_reduce <- cmat[, !colnames(cmat) %in% mobj$aliased_terms]
      #write.csv(cmat_reduce, file = "cmat_test_reduce.csv")
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