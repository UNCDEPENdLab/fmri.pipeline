#' Mixed by runs a set of mixed-effects models for each combination of a set of factors. Its primary use is to
#'   run the same model on different splits of the data.
#'
#' @param data A data.frame or data.table object containing stacked data for each combination of the \code{split_on}
#'   variables. The function will run separate mixed-effect models for each combination. Alternatively, a vector of
#'   filenames can be passed, which will be read in sequentially and fit (.rds, .csv, .dat, and .txt
#'   supported at present).
#' @param outcomes A character vector of outcome variables to be analyzed
#' @param rhs_model_formulae A lme4-format formula specifying the exact model to be run for each data split.
#' @param split_on A character vector of columns in \code{data} used to split the analyses into separate models.
#' @param external_df An optional data.frame/data.table containing external data that should be joined with \code{data}
#'   prior to model fitting. Useful if \code{data} contains external time series data and \code{external_df} is a
#'   dataset of behavioral variables that do not vary by neural sensor/region. Optionally, this can be a single
#'   filename to a .rds file if you want to pass in the filename, not the data itself.
#' @param external_merge_by A character vector specifying which columns of \code{external_df} and \code{data} should
#'   be used for creating a combined dataset.
#' @param padjust_by A character vector or list consisting of one or more variables over which an adjusted p-value
#'   should be calculated. This defaults to adjusting by each term (fixed effect) in the model, which will adjust for
#'   all tests for a given term across variables in \code{split_on}.
#' @param padjust_method The adjustment method (see \code{?p.adjust}) for adjusting p-values. Multiple values
#'   can be passed as a character vector, in which case multiple corrections will be added as distinct columns.
#' @param outcome_transform A vectorized function that will be applied to the outcome variable prior to running
#'   the model. For example, \code{outcome_transform=function(x) { as.vector(scale(x)) } }.
#' @param ncores The number of compute cores to be used in the computation. Defaults to 1.
#' @param cl An optional external cl object (created by a variant of makeCluster) used for computation. Can
#'   save the overhead of starting and stopping many workers in a loop context.
#' @param refit_on_nonconvergence The number of times a model should be refit if it does not converge. Final estimates
#'   from one iteration are used as starting values for the next.
#' @param tidy_args A list of arguments passed to tidy.merMod for creating the coefficient data.frame. By default,
#'   the function only returns the fixed effects and computes confidence intervals using the Wald method.
#' @param lmer_control An lmerControl object specifying any optimization settings to be passed to lmer()
#'
#' @return A data.table object containing all coefficients for each model estimated separately by \code{split_on}
#'
#' @importFrom lmerTest lmer
#' @importFrom lme4 lmerControl
#' @importFrom checkmate assert_data_frame assert_character assert_subset assert_formula
#' @importFrom foreach %dopar% foreach registerDoSEQ
#' @importFrom iterators iter
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom broom.mixed tidy
mixed_by <- function(data, outcomes = NULL, rhs_model_formulae = NULL, split_on = NULL,
                     external_df = NULL, external_merge_by = NULL,
                     padjust_by = "term", padjust_method = "BY", outcome_transform = NULL,
                     ncores = 1L, cl = NULL, refit_on_nonconvergence = 3,
                     tidy_args = list(effects = "fixed", conf.int = TRUE),
                     lmer_control = lmerControl(optimizer = "nloptwrap")) {
  require(data.table) # remove for package
  require(dplyr)
  require(lme4)
  require(lmerTest)
  require(foreach)
  require(doParallel)
  require(iterators)
  require(broom.mixed)

  ## VALIDATE INPUTS
  # support data.frame input for single dataset execution or a vector of files that are imported and fit sequentially
  if (checkmate::test_data_frame(data)) {
    single_df <- TRUE
  } else {
    if (!checkmate::test_character(data)) {
      stop("If data is not a data.frame, it should be a set of file names.")
    }
    checkmate::assert_file_exists(data)
    single_df <- FALSE
  }

  checkmate::assert_character(outcomes, null.ok = FALSE)
  checkmate::assert_character(split_on, null.ok = TRUE, unique = TRUE)
  if (!is.null(padjust_by)) {
    if (!is.list(padjust_by)) {
      padjust_by <- list(padjust_by)
    } # convert to list for consistency
    sapply(padjust_by, checkmate::assert_character, null.ok = TRUE)
  }
  checkmate::assert_string(padjust_method)
  checkmate::assert_integerish(ncores, lower = 1L)
  checkmate::assert_class(cl, "cluster", null.ok = TRUE)

  # turn off refitting if user specifies 'FALSE'
  if (is.logical(refit_on_nonconvergence) && isFALSE(refit_on_nonconvergence)) {
    refit_on_nonconvergence <- 0L
  }
  checkmate::assert_integerish(refit_on_nonconvergence, null.ok = FALSE)

  if (inherits(rhs_model_formulae, "formula")) {
    rhs_model_formulae <- list(rhs_model_formulae)
  } # wrap as list
  lapply(rhs_model_formulae, checkmate::assert_formula)
  if (is.null(names(rhs_model_formulae))) {
    nm <- paste0("model", seq_along(rhs_model_formulae))
    message("Using default model names of ", nm)
    names(rhs_model_formulae) <- nm
  }

  # handle external_df
  if (!is.null(external_df)) {
    if (checkmate::test_string(external_df)) {
      checkmate::assert_file_exists(external_df)
      external_df <- readRDS(external_df) # only supports .rds at the moment
    } else {
      checkmate::assert_data_frame(external_df)
      if (!is.data.table(external_df)) {
        setDT(external_df)
      }
      checkmate::assert_subset(external_merge_by, names(external_df))
      setkeyv(external_df, external_merge_by) # key external data by merge columns
    }
  }

  # worker subfunction to fit a given model to a data split
  model_worker <- function(data, model_formula, lmer_control, outcome_transform = NULL) {
    if (!is.null(outcome_transform)) { # apply transformation to outcome
      lhs <- all.vars(model_formula)[1]
      data[[lhs]] <- outcome_transform(data[[lhs]])
    }

    md <- lmerTest::lmer(model_formula, data, control = lmer_control)

    if (refit_on_nonconvergence > 0L) {
      rfc <- 0
      while (any(grepl("failed to converge", md@optinfo$conv$lme4$messages)) && rfc < refit_on_nonconvergence) {
        # print(md@optinfo$conv$lme4$conv)
        ss <- getME(md, c("theta", "fixef"))
        lmod <- lmer_control
        lmod$optimizer <- "bobyqa" # produces convergence more reliably
        md <- update(md, start = ss, control = lmod)
        rfc <- rfc + 1 # increment refit counter
      }
    }
    return(md)
  }

  # convert to data.table and nest
  if (isTRUE(single_df)) {
    if (!is.data.table(data)) {
      setDT(data)
    } # convert to data.table by reference to avoid RAM copy
    # bad idea -- will copy data and then original dataset is maintained in calling environment after rm() since it is not fully dereferenced
    # dt <- data.table(data)
    # rm(data) #avoid any lingering RAM demand
    df_set <- c("internal")
  } else {
    df_set <- data
    # even though we aren't splitting on this, adding it will propagate the filename to the output structure
    if (!".filename" %in% split_on) {
      split_on <- c(".filename", split_on)
    }
  }

  # setup parallel compute
  if (ncores > 1L && is.null(cl)) {
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    on.exit(try(stopCluster(cl)))
  } else if (!is.null(cl)) {
    message("Using external cl object for parallelism. mixed_by will not stop the cluster upon completion.")
    message("Use stopCluster() yourself when done with all computation.")
    registerDoParallel(cl, cores = ncores)
  } else {
    registerDoSEQ()
  }

  model_set <- expand.grid(outcome = outcomes, rhs = rhs_model_formulae)
  mresults <- vector(mode = "list", length(df_set)) # preallocate list

  # loop over each dataset to be fit
  for (i in seq_along(df_set)) {
    df_i <- df_set[i]

    # read each dataset if operating in multiple data scenario
    if (isFALSE(single_df)) {
      message("Reading file: ", df_i)
      if (grepl(".rds$", df_i, ignore.case = TRUE)) {
        data <- readRDS(df_i)
        if (!is.data.table(data)) {
          setDT(data)
        }
      } else if (grepl("(.csv|.csv.gz|.csv.bz2|.dat|.txt|.txt.gz|.txt.bz2)", df_i, ignore.case = TRUE, perl = TRUE)) {
        data <- fread(df_i, data.table = TRUE)
      } else {
        stop("Unable to sort out this data input: ", df_i)
      }
      data[, .filename := basename(df_i)] # add filename
    }

    # validate structure of data against models to be fit
    if (is.null(split_on)) {
      split_on <- "split" # dummy split to make code function consistently
      data[, split := factor(1)]
      has_split <- FALSE
    } else {
      has_split <- TRUE
    }

    # handle external_df, if requested
    if (!is.null(external_df)) {
      checkmate::assert_subset(external_merge_by, names(data))
      data <- merge(data, external_df, by = external_merge_by)
    }

    # verify that outcomes and split variables are present in data
    checkmate::assert_subset(outcomes, names(data))
    checkmate::assert_subset(split_on, names(data))

    # nest data.tables for each combination of split factors
    setkeyv(data, split_on)
    data <- data[, .(dt = list(.SD)), by = split_on]

    # loop over outcomes and rhs formulae within each chunk to maximize compute time by chunk (reduce worker overhead)
    mresults[[i]] <- foreach(
      dt_split = iter(data, by = "row"), .packages = c("lme4", "lmerTest", "data.table", "dplyr", "broom.mixed"),
      .noexport = "data", .export = "split_on", .inorder = FALSE, .combine = rbind
    ) %dopar% {
      split_results <- lapply(seq_len(nrow(model_set)), function(mm) {
        ff <- update.formula(model_set$rhs[[mm]], paste(model_set$outcome[[mm]], "~ .")) # replace LHS
        ret <- copy(dt_split)
        thism <- model_worker(ret$dt[[1]], ff, lmer_control, outcome_transform)
        ret[, dt := NULL] # drop original data.table for this split
        ret[, outcome := model_set$outcome[[mm]]]
        ret[, model_name := names(model_set$rhs)[mm]]
        ret[, rhs := as.character(model_set$rhs[mm])]
        # ret[, model := list(list(thism))] #not sure why a double list is needed, but a single list does not make a list-column
        ret[, coef_df := list(do.call(tidy, append(tidy_args, x = thism)))]
        return(ret)
      })

      split_results <- rbindlist(split_results)
      coef_df <- split_results[, coef_df[[1]], by = .(outcome, model_name, rhs)] # unnest coefficients
      coef_df <- cbind(dt_split[, ..split_on], coef_df) # add back metadata for this split

      coef_df # return
    }

    rm(data) # cleanup big datasets and force memory release before next iteration
    gc()
  }

  # need to put this above, but trying to avoid tmp objects
  # data[, filename:=df_i] #tag for later

  mresults <- rbindlist(mresults) # combine results from each df (in the multiple df case)
  setorderv(mresults, split_on) # since we allow out-of-order foreach, reorder coefs here.

  # compute adjusted p values
  if (!is.null(padjust_by)) {
    checkmate::assert_subset(padjust_method, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
    for (ff in padjust_by) {
      checkmate::assert_subset(ff, names(mresults))
      cname <- paste0("padj_", padjust_method, "_", paste(ff, collapse = "_"))
      mresults <- mresults[, (cname) := p.adjust(p.value, method = padjust_method), by = ff]
    }
  }

  # drop off dummy split if irrelevant
  if (isFALSE(has_split)) {
    mresults[, split := NULL]
  }

  return(mresults)
}
