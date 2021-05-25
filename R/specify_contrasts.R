#' Internal function to add, edit, rename, delete contrasts
#'
#' @param mobj \code{l1_model_spec} or \code{hi_model_spec} object
#' @param include_diagonal A logical indicating whether to include a diagonal matrix of contrasts
#'   in the model contrasts. Defaults to TRUE.
#' @internal
#' @author Michael Hallquist
#' @importFrom magrittr %>%
#' @importFrom emmeans emmeans emtrends
specify_contrasts <- function(mobj = NULL, include_diagonal = TRUE, include_factormeans = FALSE) {
    checkmate::assert_multi_class(mobj, c("l1_model_spec", "hi_model_spec")) # verify that we have an object of known structure
    checkmate::assert_logical(include_diagonal, len = 1)

    prompt_contrasts <- FALSE
    if (is.null(mobj$contrasts)) {
      # handle diagonal contrasts
      include_diagonal <- menu(c("Yes", "No"), title = "Do you want to include diagonal contrasts for each regressor?")
      include_diagonal <- ifelse(include_diagonal == 1L, TRUE, FALSE)
      prompt_contrasts <- TRUE

      cond_means <- c() # factors for which we want model-predicted condition means
      pairwise_diffs <- c() # factors for which we want pairwise differences
      cell_means <- FALSE # whether to include cell means across all factors in design

      # if we are working from a model specification in lm(), ask about factor means and pairwise comparisons
      if (!is.null(mobj$lmfit)) {
          term_labels <- attr(terms(mobj$lmfit), "term.labels") # names of regressors
          int_terms <- attr(terms(mobj$lmfit), "order") # interactions in model
          dclass <- attr(mobj$lmfit$terms, "dataClasses")[term_labels] # omit response variable
          cat_vars <- names(dclass[which(dclass %in% c("factor", "character"))])
          cont_vars <- names(dclass[which(dclass %in% c("integer", "numeric"))])

          cat(
            "First, please decide on your weighting scheme for computing contrasts.",
            "For details, see the weights section of ?emmeans.",
            "In general, we recommend cell weights in unbalanced designs to",
            "  approximate ordinary marginal means.\n",
            sep = "\n"
          )
          ii <- menu(c(
            "Equal (all cells weighted the same)",
            "Cells (frequencies of cells being averaged)",
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
          simple_slopes <- list()
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
    } else {
        cat("\nCurrent model contrasts:\n\n")
        print(mobj$contrasts)
        cat("\n--------\n\n")

        res <- menu(c("No", "Yes"), title = "Do you want to modify model contrasts?")
        if (res == 2L) {
            prompt_contrasts <- TRUE
        }
    }

    #if the user does not want to modify existing contrasts, return object unchanged
    if (isFALSE(prompt_contrasts)) { return(mobj) }

    #start with existing contrast matrix, if any
    cmat <- mobj$contrasts

    if (isTRUE(include_diagonal)) {
        diag_mat <- diag(length(mobj$model_regressors))
        rownames(diag_mat) <- paste0("EV_", mobj$model_regressors) # simple contrast naming for each individual regressor
        colnames(diag_mat) <- mobj$model_regressors # always have columns named by regressor

        cmat <- rbind(cmat, diag_mat)
        cmat <- cmat[!duplicated(cmat, MARGIN = 1), ] # don't add duplicate diagonal contrasts to matrix, if already present
    }

    if (is.null(cmat)) {
        cmat <- matrix(numeric(0), nrow = 0, ncol = length(mobj$model_regressors), dimnames = list(NULL, mobj$model_regressors))
    } else {
        # currently blow up if we change regressors (contrast matrix does not match new regressors)
        # TODO: build out zeros in added columns
        stopifnot(identical(colnames(cmat), mobj$model_regressors))
    }

    # add condition means, if requested
    if (length(cond_means) > 0L) {
      for (vv in cond_means) {
        ee <- emmeans(mobj$lmfit, as.formula(paste("~", vv)), weights=weights)
        edata <- summary(ee)
        econ <- ee@linfct
        enames <- as.character(edata[[vv]]) #names of contrasts

        # if any emmeans are not estimable, then this is likely due to aliasing. For now, drop
        which_na <- is.na(edata$emmean)
        if (any(which_na)) {
          econ <- econ[!which_na, , drop=FALSE]
          enames <- enames[!which_na]
        }

        # add contrast names to matrix
        rownames(econ) <- enames

        # add contrasts to matrix
        cmat <- rbind(cmat, econ)
        cmat <- cmat[!duplicated(cmat, MARGIN = 1), ] # don't add duplicate diagonal contrasts to matrix, if already present
      }
    }

    # add pairwise differences, if requested
    if (length(pairwise_diffs) > 0L) {
      for (vv in pairwise_diffs) {
        ee <- emmeans(mobj$lmfit, as.formula(paste("~", vv)), weights=weights)
        pp <- pairs(ee)
        edata <- summary(pp)
        econ <- pp@linfct
        enames <- make.names(sub("\\s+-\\s+", "_M_", edata$contrast, perl = TRUE)) # names of contrasts

        # if any emmeans are not estimable, then this is likely due to aliasing. For now, drop
        which_na <- is.na(edata$estimate)
        if (any(which_na)) {
          econ <- econ[!which_na, , drop = FALSE]
          enames <- enames[!which_na]
        }

        # add contrast names to matrix
        rownames(econ) <- enames

        cmat <- rbind(cmat, econ)
        cmat <- cmat[!duplicated(cmat, MARGIN = 1), ] # remove any duplicate contrasts from matrix

      }
    }

    # add cell means for all factors, if requested
    if (isTRUE(cell_means)) {
      # get model-predicted means for each factor
      ee <- emmeans(mobj$lmfit, as.formula(paste("~", paste(cat_vars, collapse = "*"))), weights=weights)
      #pp <- pairs(ee)
      edata <- summary(ee)
      econ <- ee@linfct
      enames <- apply(edata[, cat_vars, drop = FALSE], 1, paste, collapse = ".")

      # if any emmeans are not estimable, then this is likely due to aliasing. For now, drop
      which_na <- is.na(edata$emmean)
      if (any(which_na)) {
        econ <- econ[!which_na, , drop = FALSE]
        enames <- enames[!which_na]
      }

      rownames(econ) <- enames
      cmat <- rbind(cmat, econ)
    }

    #add overall response mean
    if (isTRUE(overall_response)) {
      ee <- emmeans(mobj$lmfit, ~1, weights=weights)
      econ <- ee@linfct
      rownames(econ) <- "overall"
      cmat <- rbind(cmat, econ)
    }

    #add per-cell slopes
    if (length(simple_slopes) > 0L) {
      # get model-predicted simple slopes
      for (vv in seq_along(simple_slopes)) {
        trend_var <- names(simple_slopes)[vv]
        for (comb in simple_slopes[[vv]]) {
          ee <- emtrends(mobj$lmfit,
            specs = as.formula(paste("~", paste(comb, collapse = "*"))),
            var = trend_var, weights = weights
          )
          edata <- summary(ee)
          econ <- ee@linfct
          enames <- paste(trend_var, apply(edata[, comb, drop = FALSE], 1, paste, collapse = "."), sep = ".")

          # if any emmeans are not estimable, then this is likely due to aliasing. For now, drop
          which_na <- is.na(edata[[paste(trend_var, "trend", sep=".")]])
          if (any(which_na)) {
            econ <- econ[!which_na, , drop = FALSE]
            enames <- enames[!which_na]
          }

          rownames(econ) <- enames
          cmat <- rbind(cmat, econ)
        }
      }
    }

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

    #write.csv(cmat, file="cmat_test.csv")

    # syntax of contrast
    cat("\n\n-----\nContrast editor\n\n")
    cat("When entering contrasts, you have two options:\n\n")
    cat("  1) enter a vector of numeric values, one for each regressor\n     Example: 1 0 -1 0 0\n\n")
    cat("  2) enter name=value pairs for relevant coefficients\n     Example: pe_1h = 1 pe_2h = -1\n\n")

    add_more <- 1
    while (add_more != 4) {
        add_more <- menu(c("Add contrast", "Show contrasts", "Delete contrast", "Done with contrast setup"),
            title = "Contrast setup menu"
        )

        if (add_more == 1L) {
            cname <- readline("Enter contrast name: ")
            cat("Columns of contrast matrix (one per regressor)\n")
            cat(strwrap(paste(colnames(cmat), collapse = ", "), 70, exdent = 5), sep="\n")
            centry <- readline("Enter contrast syntax: ")
            if (grepl("=", centry)) {
                # split into pairs
                centry <- gsub("([^\\s]+)\\s*=\\s*([^\\s]+)", "\\1=\\2", centry, perl = TRUE)
                csplit <- do.call(rbind, strsplit(strsplit(centry, "\\s+")[[1L]], "=")) # first element is name, second is value
                cvec <- rep(0, ncol(cmat)) %>% setNames(colnames(cmat))
                bad_con <- FALSE
                for (ii in seq_len(nrow(csplit))) {
                    if (!cvec[csplit[ii, 1]] %in% colnames(cmat)) {
                        warning("Cannot find column: ", cvec[csplit[ii, 1]], " in contrast matrix. Ignoring contrast input.")
                        bad_con <- TRUE
                    } else {
                        cvec[csplit[ii, 1]] <- as.numeric(csplit[ii, 2])
                    }
                }
                if (!bad_con) {
                    cmat <- rbind(cmat, cvec)
                    rownames(cmat)[nrow(cmat)] <- cname
                }
            } else {
                cvec <- sapply(strsplit(centry, "(,\\s*|\\s+)")[[1L]], as.numeric)
                if (length(cvec) != ncol(cmat)) {
                    warning("Number of elements in entered contrast does not match number of columns in design. Not adding contrast")
                } else {
                    cmat <- rbind(cmat, cvec)
                    rownames(cmat)[nrow(cmat)] <- cname
                }
            }
        } else if (add_more == 2L) {
            summarize_contrasts(cmat)
        } else if (add_more == 3L) {
            whichdel <- menu(rownames(cmat), title = "Which contrast should be deleted?")
            if (whichdel != 0) {
                cmat <- cmat[-whichdel, , drop = FALSE]
            }
        } else if (add_more == 4L) {
            cat("Exiting contrast setup\n")
        }
    }

    mobj$contrasts <- cmat
    mobj$contrasts_noalias <- cmat_reduce
    return(mobj)
}