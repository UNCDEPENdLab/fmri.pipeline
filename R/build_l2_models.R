#' Interactive function to build an l1 model specification for setup_glm_pipeline
#'
#' @param data a data.frame containing trial-level data for one or more subjects
#' @param model_set optional existing model_set to be modified
#' @param variable_mapping a vector of mappings between columns in \code{data} and internal constructs
#' @param regressor_cols an optional character vector of columns in \code{data} that should be considered
#'   as possible regressors
#'
#' @return a \code{hi_model_set} object containing a list of higher-level regression models
#' @author Michael Hallquist
#' @importFrom checkmate assert_data_frame assert_class assert_subset
#' @importFrom magrittr %>%
#' @importFrom tidyselect all_of
#' @importFrom dplyr select
#' @export
#'
build_l2_models <- function(data, model_set = NULL,
                            variable_mapping = c(id = "id", run = "run", trial = "trial", run_trial = "trial", mr_dir = "mr_dir"),
                            regressor_cols = NULL) {

    # Maybe allow glm object to be passed in that would have data and variable_mapping.
    # I guess that would be like "add_l1_model"
    checkmate::assert_data_frame(data) # yeah, move toward allowing the broader model specification object here
    checkmate::assert_class(model_set, "hi_model_set", null.ok = TRUE)
    checkmate::assert_subset(regressor_cols, names(data)) # make sure all parametric regressor columns are in the data frame

    possible_cols <- names(data)
    possible_cols <- possible_cols[!possible_cols %in% variable_mapping[c("mr_dir")]]
    data <- data %>% dplyr::select(all_of(possible_cols))

    if (!is.null(regressor_cols)) { # select only columns requested by user
        data <- data %>% dplyr::select(regressor_cols)
    }

    if (is.null(model_set)) {
        ## initialize overall higher-level model set object
        model_set <- list(models = NULL)
        class(model_set) <- c("list", "hi_model_set")
    }

    # regressor manager
    get_regressors <- function(data, regressor_cols = NULL) {
        done_regressors <- FALSE
        while (isFALSE(done_regressors)) {
            cat("Current regressors for this model:\n\n  ", paste(regressor_cols, collapse = ", "), "\n\n")
            action <- menu(c("Add/modify regressors", "Delete regressors", "Done with regressor selection"),
                title = "Would you like to modify the model regressors?"
            )

            if (action == 1L) { # Add/modify
                regressor_cols <- select.list(names(data),
                    multiple = TRUE, preselect = regressor_cols,
                    title = "Choose all model regressors\n(Command/Control-click to select multiple)"
                )
            } else if (action == 2L) { # Delete
                if (length(regressor_cols) == 0L) {
                    message("No regressors yet. Please add at least one.")
                } else {
                    which_del <- menu(regressor_cols, title = "Which regressor would you like to remove?")
                    if (which_del > 0) {
                        proceed <- menu(c("Proceed", "Cancel"),
                            title = paste0("Are you sure you want to delete ", regressor_cols[which_del], "?")
                        )
                        if (proceed == 1) {
                            cat("  Deleting ", regressor_cols[which_del], "\n")
                            regressor_cols <- regressor_cols[-which_del]
                        } else {
                            cat("  Not deleting ", regressor_cols[which_del], "\n")
                        }
                    }
                }
            } else if (action == 3L) { # Done
                if (length(regressor_cols) == 0L) { 
                    message("No regressors yet. Please add at least one.")
                } else {
                    done_regressors <- TRUE
                    cat("The following columns were chosen as model regressors.\n\n")
                    cat("  ", paste(regressor_cols, collapse = ", "), "\n\n")
                }
            }
        }
        return(regressor_cols)
    }

    ############### BUILD MODELS

    summarize_models <- function(ml) {
        if (length(ml) == 0L) {
            return(invisible(NULL))
        }
        lapply(seq_along(ml), function(ii) {
            this <- ml[[ii]]
            cat("--------\nModel ", ii, "\n\n")
            cat("  Name:", this$name, "\n")
            cat("  Regressors:", paste(this$model_variables, collapse = ", "), "\n")
            if (ncol(this$contrasts) < 30) {
                cat("  Contrasts:\n\n")
                print(round(this$contrasts, 3))
                cat("\n--------\n\n")
            } else {
                cat("More than 30 regressor columns. Not printing contrasts\n")
            }
        })
    }

    create_new_model <- function(data, to_modify = NULL) {
        checkmate::assert_class(to_modify, "hi_model_spec", null.ok = TRUE)
        if (is.null(to_modify)) {
            mm <- list()
            class(mm) <- c("list", "hi_model_spec")
            modify <- FALSE
        } else {
            mm <- to_modify
            modify <- TRUE
        }

        ### ------ model name ------
        if (isTRUE(modify)) {
            cat("Current model name:", mm$name, "\n")
            res <- menu(c("No", "Yes"), title = "Change model name?")
            if (res == 2) {
                mm$name <- NULL
            } # clear out so that it is respecified
        }

        while (is.null(mm$name) || mm$name == "") {
            res <- trimws(readline("Enter the model name: "))
            if (res != "") {
                res <- make.names(res)
                if (res %in% names(model_list)) {
                    cat("\nModel name:", res, "already exists. Names must be unique.\n")
                    cat("Current models:", paste(names(model_list), collapse = ", "), "\n")
                } else {
                    mm$name <- res
                }
            }
        }

        if (isTRUE(modify)) {
            cat("Current model regressors:", paste(names(mm$model_variables), collapse = ", "), "\n")
            res <- menu(c("No", "Yes"), title = "Respecify model regressors (and contrasts)?")
            if (res == 2) { # clear out so that it is respecified
                cat("Okay, resetting model regressors and contrasts\n")
                mm$model_variables <- mm$contrasts <- NULL
            }
        }

        #let user decide LM formula approach versus walkthrough
        model_approach <- 0L
        while (model_approach == 0L) {
            model_approach <- menu(c("Specify formula", "Model builder"),
                title = "Do you want to specify the model formula or walk through the model builder?"
            )
        }
        if (model_approach == 1L) {
            # formula approach
            res <- NULL
            while (is.null(res)) {
                cat(c(
                    "\nSpecify the right-hand side of the model you wish to fit.",
                    "Use variable names in the dataset provided to this function.",
                    "Note that this syntax follows standard R model syntax. See ?lm for details.",
                    "Example: ~ emotion * wmload + run_number\n",
                    "Available column names: \n"
                ), sep = "\n")
                cat(paste(names(data), collapse = ", "), sep = "\n")
                res <- trimws(readline("Enter the model formula: "))
                # always trim any LHS specification
                res <- sub("^[^~]*", "", res, perl=TRUE)
                res <- tryCatch(as.formula(res), error = function(e) {
                    print(e)
                    cat("Problem converting your syntax to formula. Try again\n")
                    return(NULL)
                })
            }

            #ensure that all variables in formula are present in the data frame
            mm$model_variables <- all.vars(res)
            checkmate::assert_subset(mm$model_variables, names(data), empty.ok = FALSE)
        } else if (model_approach == 2L) {
            # walkthrough approach (only support additive model for now)
            mm$model_variables <- get_regressors(data, regressor_cols = mm$model_variables)
            res <- as.formula(paste("~", paste(mm$model_variables, collapse = " + ")))
        }

        # convert integer-like variables to integer
        for (vv in mm$model_variables) {
            #convert integer-like numbers to integers for possible conversion below
            if (is.numeric(data[[vv]]) && checkmate::test_integerish(data[[vv]])) {
                data[[vv]] <- as.integer(data[[vv]])
            }
            if (is.integer(data[[vv]]) && length(unique(data[[vv]]) < 20)) {
                cat(
                    "\n",
                    vv, " appears to be an integer with < 20 levels.\n",
                    "If you convert it to a factor, R will generate dummy code variables in the model.\n",
                    "If you leave it as a number, R will include the linear trend in the model.\n\n",
                    sep = ""
                )
                conv <- model_approach <- menu(c("Convert to factor", "Leave as number"),
                    title = paste("Do you want to convert", vv, "to a factor?")
                )
                if (conv == 1L) {
                    cat("Converting", vv, "into a factor.\n")
                    data[[vv]] <- as.factor(data[[vv]])
                }
            }
        }

        #need to build a model LM-style
        cat("Summary of variables included in model:\n")
        print(summary(data[, mm$model_variables]))

        modelmat <- model.matrix(res, data)
        data$dummy <- rnorm(nrow(data))
        ffit <- update.formula(res, "dummy ~ .") # add LHS
        mm$lmfit <- lm(ffit, data)
        mm$model_regressors <- colnames(modelmat) # actual regressors after expanding categorical variables
        mm$model_matrix <- modelmat

        # handle coefficient aliasing
        al <- alias(mm$lmfit)
        if (!is.null(al$Complete)) {
            cat("Problems with aliased (redundant) terms in model.\n\n")
            bad_terms <- rownames(al$Complete)
            cat(paste(bad_terms, collapse = ", "), "\n\n")

            # find unaliased (good) terms in the model
            good_terms <- colnames(modelmat)[!(colnames(modelmat) %in% bad_terms)]
            good_terms <- good_terms[!good_terms == "(Intercept)"]

            # build new model formula with only good terms
            newf <- as.formula(paste("dummy ~", paste(good_terms, collapse = " + ")))

            # also generate a model data.frame that expands dummy codes for terms, retaining only good variables
            mm$model_matrix_noalias <- modelmat[, grep(":", good_terms, fixed = TRUE, value = TRUE, invert = TRUE)]
            modeldf <- as.data.frame(mm$model_matrix_noalias)
            modeldf$dummy <- data$dummy # copy across dummy DV for fitting
            mm$lmfit_noalias <- lm(newf, modeldf)
            mm$aliased_terms <- bad_terms

            # N.B. emmeans needs to calculate contrasts on the original design to see the factor structure
            # So, we also need to drop out columns from the emmeans linfct
        }

        # walk through contrast generation for this model
        mm <- specify_contrasts(mm)

        return(mm)
    }

    model_list <- model_set$models
    add_more <- 1
    while (add_more != 4) {
        summarize_models(model_list)

        add_more <- menu(c("Add model", "Modify model", "Delete model", "Done with higher-level model setup"),
            title = "Higher-level model setup menu"
        )

        if (add_more == 1L) { # add
            mm <- create_new_model(data)
            if (mm$name %in% names(model_set)) {
                warning("An model with the same name exists: ", mm$name, ". Overwriting it.")
            }
            model_list[[mm$name]] <- mm # add to set
        } else if (add_more == 2L) { # modify
            if (is.null(model_list)) {
                message("No models available to modify. Add at least one model first")
            } else {
                res <- 0L
                while (res == 0L) {
                    res <- menu(names(model_list), title = "Which model do you want to modify?")
                }
                model_list[[res]] <- create_new_model(data, to_modify = model_list[[res]])
            }
        } else if (add_more == 3L) { # delete
            which_del <- menu(names(model_list), title = "Which model would you like to delete?")
            if (which_del > 0) {
                proceed <- menu(c("Proceed", "Cancel"),
                    title = paste0("Are you sure you want to delete ", names(model_list)[which_del], "?")
                )
                if (proceed == 1) {
                    cat("  Deleting ", names(model_list)[which_del], "\n")
                    model_list[[which_del]] <- NULL
                } else {
                    cat("  Not deleting ", names(model_list)[which_del], "\n")
                }
            }
        }
    }

    model_set$models <- model_list

    return(model_set)
}

#' @rdname build_l3_models
#' @export
build_l3_models <- build_l2_models