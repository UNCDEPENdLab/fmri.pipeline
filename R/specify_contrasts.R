  #' Internal function to add, edit, rename, delete contrasts
  #'
  #' @param mobj \code{l1_model_spec} or \code{hi_model_spec} object
  #' @param include_diagonal A logical indicating whether to include a diagonal matrix of contrasts
  #'   in the model contrasts. Defaults to TRUE.
  #' @internal
  #' @author Michael Hallquist
  #' @importFrom magrittr %>%
  specify_contrasts <- function(mobj = NULL, include_diagonal = TRUE) {
      checkmate::assert_multi_class(mobj, c("l1_model_spec", "hi_model_spec")) # verify that we have an object of known structure
      checkmate::assert_logical(include_diagonal, len = 1)

      cmat <- mobj$contrasts

      if (isTRUE(include_diagonal)) {
          diag_mat <- diag(length(mobj$model_regressors))
          rownames(diag_mat) <- mobj$model_regressors # simple contrast naming for each single-regressor contrast
          colnames(diag_mat) <- mobj$model_regressors # always have columns named by regressor

          cmat <- rbind(cmat, diag_mat)
          cmat <- cmat[!duplicated(cmat, MARGIN = 1), ] # don't add duplicate diagonal contrasts to matrix, if alread present
      }

      if (is.null(cmat)) {
          cmat <- matrix(numeric(0), nrow = 0, ncol = length(mobj$model_regressors), dimnames = list(NULL, mobj$regressors))
      } else {
          # currently blow up if we change regressors (contrast matrix does not match new regressors)
          # TODO: build out zeros in added columns
          stopifnot(identical(colnames(cmat), mobj$model_regressors))
      }

      # syntax of contrast
      cat("\n\n-----\nContrast editor\n\n")
      cat("When entering contrasts, you have two options:\n\n")
      cat("  1) enter a vector of numeric values, one for each regressor\n     Example: 1 0 -1 0 0\n\n")
      cat("  2) enter name=value pairs for relevant coefficients\n     Example: pe_1h = 1 pe_2h = -1\n\n")

      summarize_contrasts <- function(cmat) {
          sapply(seq_len(nrow(cmat)), function(x) {
              con_name <- rownames(cmat)[x]
              cvec <- cmat[x, ]
              nzcols <- which(cvec != 0)
              cols <- colnames(cmat)[nzcols]
              cat("Contrast: ", con_name, "\n")
              for (ii in seq_along(cols)) {
                  cat(cols[ii], "=", cvec[nzcols[ii]], "\n")
              }
              cat("-----\n")
          })
      }

      add_more <- 1
      while (add_more != 4) {
          add_more <- menu(c("Add contrast", "Show contrasts", "Delete contrast", "Done with contrast setup"),
              title = "Contrast setup menu"
          )

          if (add_more == 1L) {
              cname <- readline("Enter contrast name: ")
              cat("Columns of contrast matrix (one per regressor)\n")
              cat(strwrap(paste(colnames(cmat), collapse = ", "), 70, exdent = 5), "\n")
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
      return(mobj)
  }