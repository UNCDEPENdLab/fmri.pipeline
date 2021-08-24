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
#' @param calculate A character vector specifying what calculations should be returned by the function. 
#'   The options are: "parameter_estimates_reml", "parameter_estimates_ml", and "fit_statistics".
#' @param emmeans_spec A named list of emmeans calls to be run for each model to obtain model predictions.
#'   Any arguments that are valid for emmean can be passed through the list structure.
#' @param return_models A boolean indicating whether to return fitted model objects, which can be used for
#'   post hoc contrasts, visualization and statistics. Note that model objects can get very large, so be
#'   careful with this option since it could generate a massive data object.
#'   
#' @details In general, restricted maximum likelihood (REML) should be used for making inferences about parameter
#'   estimates, whereas ML should be used for model comparisons based on log-likelihood (e.g., AIC). 
#'   If "parameter_estimates_reml" are requested in \code{calculate}, then models will
#'   be fitted with REML. If "parameter_estimates_ml" and/or "fit_statistics" are requested, models will be
#'   fitted with ML. Note that if both REML and ML are requested, each model is fit twice since the estimators
#'   each have advantages and disadvantages noted above.
#'   
#' Example of emmeans_spec usage:
#'   mixed_by(data, emmeans_spec=list(
#'     em1=list(specs = ~ memory | noise_level, adjust = "sidak", weights = "cells"),
#'     em2=list(specs = ~ memory * noise_level, weights = "equal"),
#'     em3=list(specs = ~ memory)
#'   ))
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
#' @importFrom data.table fread setDT setkeyv
mixed_by <- function(data, outcomes = NULL, rhs_model_formulae = NULL, split_on = NULL,
                     external_df = NULL, external_merge_by = NULL,
                     padjust_by = "term", padjust_method = "BY", outcome_transform = NULL,
                     ncores = 1L, cl = NULL, refit_on_nonconvergence = 3,
                     tidy_args = list(effects = "fixed", conf.int = TRUE),
                     lmer_control = lmerControl(optimizer = "nloptwrap"),
                     calculate=c("parameter_estimates_reml", "parameter_estimates_ml", "fit_statistics"),
                     return_models = FALSE, emmeans_spec = NULL) {
  
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
  checkmate::assert_subset(calculate, c("parameter_estimates_reml", "parameter_estimates_ml", "fit_statistics"))
  checkmate::assert_logical(return_models, len = 1L)
  checkmate::assert_list(emmeans_spec, null.ok = TRUE)

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
  model_worker <- function(data, model_formula, lmer_control, outcome_transform = NULL, REML = TRUE) {
    if (!is.null(outcome_transform)) { # apply transformation to outcome
      lhs <- all.vars(model_formula)[1]
      data[[lhs]] <- outcome_transform(data[[lhs]])
    }

    md <- lmerTest::lmer(model_formula, data, control = lmer_control, REML = REML, na.action=na.exclude)

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
    # bad idea -- will copy data and then original dataset is maintained in calling
    # environment after rm() since it is not fully dereferenced
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

  # process emmeans specification setup
  if (!is.null(emmeans_spec)) {
    #has the outcome and model name to lookup whether to run a given emmeans specification
    if (is.null(names(emmeans_spec))) {
      names(emmeans_spec) <- paste("emm", seq_along(emmeans_spec), sep="_")
    } else {
      empty_names <- which(names(emmeans_spec) == "")
      names(emmeans_spec)[empty_names] <- paste("emm", empty_names, sep="_")
    }
    
    emm_metadata <- rbindlist(lapply(emmeans_spec, function(ee) { data.frame(ee[c("outcome", "model_name")])})) %>%
      mutate(emm_number = 1:n(), emm_label=names(emmeans_spec))
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
    message("Starting processing of data splits")
    mresults[[i]] <- foreach(
      dt_split = iter(data, by = "row"), .packages = c("lme4", "lmerTest", "data.table", "dplyr", "broom.mixed", "emmeans"),
      .noexport = "data", .inorder = FALSE
    ) %dopar% {
      split_results <- lapply(seq_len(nrow(model_set)), function(mm) {
        ff <- update.formula(model_set$rhs[[mm]], paste(model_set$outcome[[mm]], "~ .")) # replace LHS
        ret <- copy(dt_split)
        ret[, outcome := model_set$outcome[[mm]]]
        ret[, model_name := names(model_set$rhs)[mm]]
        ret[, rhs := as.character(model_set$rhs[mm])]
        
        if ("parameter_estimates_reml" %in% calculate) {
          thism <- model_worker(ret$dt[[1]], ff, lmer_control, outcome_transform)
          ret[, coef_df_reml := list(do.call(tidy, append(tidy_args, x = thism)))]
        }
        
        if (any(c("parameter_estimates_ml", "fit_statistics") %in% calculate)) {
          thism_ml <- model_worker(ret$dt[[1]], ff, lmer_control, outcome_transform, REML=FALSE) #refit with ML for AIC/BIC
          ret[, fit_df := list(glance(thism_ml))]
          ret[, coef_df_ml := list(do.call(tidy, append(tidy_args, x = thism_ml)))]
        }
        
        # no option to return model object at present
        # ret[, model := list(list(thism))] #not sure why a double list is needed, but a single list does not make a list-column
        
        # process emmeans
        if (!is.null(emmeans_spec)) {
          emm_torun <- emm_metadata %>% 
            filter(outcome == !!model_set$outcome[[mm]] & model_name == !!names(model_set$rhs)[mm])
          
          if (nrow(emm_torun) > 0L) {
            this_emmspec <- emmeans_spec[emm_torun %>% pull(emm_number)] #subset list to relevant elements
            emms <- lapply(seq_along(this_emmspec), function(emm_i) {
              tidy(do.call(emmeans, c(this_emmspec[[emm_i]], object=thism))) %>%
                dplyr::bind_cols(emm_torun[emm_i,])
            })
            names(emms) <- names(this_emmspec)
            
            #N.B. Need to double list wrap for data.table to keep list class for singleton list
            ret[, emm := list(list(emms))] 
          }
        }
        
        ret[, dt := NULL] # drop original data.table for this split from data
        
        return(ret)
      })
      
      split_results <- rbindlist(split_results, fill=TRUE)

      #need to unnest
      
      coef_df_reml <- NULL
      if ("parameter_estimates_reml" %in% calculate) {
        coef_df_reml <- split_results[, coef_df_reml[[1]], by = .(outcome, model_name, rhs)] # unnest coefficients
        coef_df_reml <- cbind(dt_split[, ..split_on], coef_df_reml) # add back metadata for this split
      }
      
      coef_df_ml <- NULL
      if ("parameter_estimates_ml" %in% calculate) {
        coef_df_ml <- split_results[, coef_df_ml[[1]], by = .(outcome, model_name, rhs)] # unnest coefficients
        coef_df_ml <- cbind(dt_split[, ..split_on], coef_df_ml) # add back metadata for this split
      }

      fit_df <- NULL
      if ("fit_statistics" %in% calculate) {
        fit_df <- split_results[, fit_df[[1]], by = .(outcome, model_name, rhs)] # unnest coefficients
        fit_df <- cbind(dt_split[, ..split_on], fit_df) # add back metadata for this split
      }
      
      emm_data <- NULL
      if (!is.null(emmeans_spec)) {
        emm_data <- subset(split_results, sapply(emm, function(x) !is.null(x)), 
                           select=c(split_on, "outcome", "model_name", "rhs", "emm"))
      }
      
      list(coef_df_reml=coef_df_reml, coef_df_ml=coef_df_ml, fit_df=fit_df, emm_data=emm_data) # return
    }
    
    rm(data) # cleanup big datasets and force memory release before next iteration
    gc()
  }

  # need to put this above, but trying to avoid tmp objects
  # data[, filename:=df_i] #tag for later
  
  #helper subfunction to pull out and rbind element from doubly nested structure above
  extract_df <- function(nested_list, element, split_on) {
    result_df <- rbindlist(lapply(nested_list, function(df_set) {
      rbindlist(lapply(df_set, "[[", element))
    })) # combine results from each df (in the multiple df case)
    setattr(result_df, "split_on", split_on) #tag split variables for secondary analysis
    setorderv(result_df, split_on) # since we allow out-of-order foreach, reorder coefs here.
    return(result_df)
  }
  
  emm_data <- NULL
  coef_results_reml <- NULL
  coef_results_ml <- NULL
  fit_results <- NULL
  if ("parameter_estimates_reml" %in% calculate) { coef_results_reml <- extract_df(mresults, "coef_df_reml", split_on) }
  if ("parameter_estimates_ml" %in% calculate) { coef_results_ml <- extract_df(mresults, "coef_df_ml", split_on) }
  if ("fit_statistics" %in% calculate) { fit_results <- extract_df(mresults, "fit_df", split_on) }
  
  emmeans_list <- NULL
  if (!is.null(emmeans_spec)) {
    emm_data <- extract_df(mresults, "emm_data", split_on)
    
    emmeans_list <- lapply(seq_along(emmeans_spec), function(aa) {
      em_sub <- subset(emm_data, outcome==emmeans_spec[[aa]]$outcome & model_name==emmeans_spec[[aa]]$model_name)
      nn <- names(em_sub)
      other_keys <- nn[nn != "emm"]
      em_name <- names(emmeans_spec)[aa]
      em_sub[, emm[[1]][[em_name]], by=other_keys]  
    })
    names(emmeans_list) <- names(emmeans_spec)
  }
  
  #helper subfunction to adjust p-values
  adjust_dt <- function(dt, padjust_by) {
    for (ff in padjust_by) {
      checkmate::assert_subset(ff, names(dt))
      cname <- paste0("padj_", padjust_method, "_", paste(ff, collapse = "_"))
      dt <- dt[, (cname) := p.adjust(p.value, method = padjust_method), by = ff]
      if (isTRUE(all.equal(dt[["p.value"]], dt[[cname]]))) {
        warning("p-value adjustment: ", paste(ff, collapse=", "), 
                " yields the same result as the uncorrected p-value. Setting ",
                cname, " to NA.")
        dt[[cname]] <- NA_real_
      }
    }
    return(dt)
  }
  
  # compute adjusted p values
  if (!is.null(padjust_by)) {
    checkmate::assert_subset(padjust_method, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
    if (!is.null(coef_results_reml)) coef_results_reml <- adjust_dt(coef_results_reml, padjust_by)
    if (!is.null(coef_results_ml)) coef_results_ml <- adjust_dt(coef_results_ml, padjust_by)
  }

  # drop off dummy split if irrelevant
  if (isFALSE(has_split)) {
    if (!is.null(coef_results_reml)) { coef_results_reml[, split := NULL] }
    if (!is.null(coef_results_ml)) { coef_results_ml[, split := NULL] }
    if (!is.null(fit_results)) { fit_results[, split := NULL] }
  }

  return(list(coef_df_reml=coef_results_reml, coef_df_ml=coef_results_ml, fit_df=fit_results, emmeans_list=emmeans_list))
}
