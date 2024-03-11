#' Helper function to output l1_models into '.yaml' object
#' @param gpa the \code{gpa} object
#' @importFrom yaml as.yaml
#' @return a yaml document of l1_models for an object
#' @export
l1_yaml <- function(gpa) {
  gpa2 <- gpa$l1_models
  ## Model Specifications
  end_yaml <- list(
  onsets = gpa2$onset,
  durations = gpa2$durations,
  isis = gpa2$isis,
  wi_factors = gpa2$wi_factors,
  values = gpa2$values
  )
  eobj <- list()
  for (ee in gpa2$events) {
    checkmate::assert_string(ee$name)
    checkmate::assert_string(ee$onset)
    if (!checkmate::test_string(ee$duration)) {
      checkmate::assert_number(ee$duration, lower=0)
    }
    checkmate::assert_string(ee$isi, null.ok = TRUE)
    name_e <- ee$name
    eobj$onset <- ee$onset
    eobj$duration <- ee$duration
    eobj$isi <- ee$isi
    end_yaml$events[[name_e]] <- eobj # this will overwrite existing specification
  }
  sobj <- list()
  for (ss in gpa2$signals) {
    name_s <- ss$name
    sobj$event <- ss$event
    sobj$trial_subset_expression <- ss$trial_subset_expression
    if (!is.null(ss$normalization)) {
      checkmate::assert_subset(ss$normalization, c("none", "evtmax_1", "durmax_1"))
      sobj$normalization <- ss$normalization
    }
    if (!is.null(ss$parametric_modulator)) {
      stopifnot(ss$parametric_modulator %in% gpa2$values)
      sobj$parametric_modulator <- ss$parametric_modulator
      sobj$value_type <- "parametric"
      sobj$value_fixed <- NULL
    } else if (!is.null(ss$value_fixed)) {
      checkmate::assert_number(ss$value_fixed)
      sobj$value_fixed <- ss$value_fixed
      if (abs(ss$value_fixed - 1) < 1e-5) {
        sobj$value_type <- "unit"
      } else {
        sobj$value_type <- "number"
      }
    }
    end_yaml$signals[[name_s]] <- sobj
  }
  mobj <- list()
    for (mm in gpa2$models) {
    name_m <- mm$name
    mobj$diagonal <- mm$contrast_spec$diagonal
    #diag(mm$contrasts) <- ifelse(mobj$diagonal==TRUE, 0, diag(mm$contrasts))
    mobj$cell_means <- mm$contrast_spec$cell_means
    mobj$overall_response <- mm$contrast_spec$overall_response
    mobj$weights <- mm$contrast_spec$weights
    mobj$signals <- mm$signals
    cobj <- list()
    for (cc in seq_along(colnames(mm$contrasts))) {
      name_c <- NULL
      name_c <- colnames(mm$contrasts)[cc]
      cobj$row <- rownames(mm$contrasts)[which(mm$contrasts[cc,]!=0)]
      cobj$contrast <- mm$contrasts[cc,which(mm$contrasts[cc,]!=0)]
      mobj$contrasts[[name_c]] <- cobj
    }
    end_yaml$l1_models[[name_m]] <- mobj
    mobj$contrasts <- NULL
  }
  endyaml3 <- as.yaml(end_yaml)
  yaml_choice <- menu(c(
          "Console Output",
          "File Output",
          "Both",
          "Exit"), title= "How would you like to receive the YAML file?")
          if (yaml_choice == 1) {
            return(cat(endyaml3, "\nGoodbye.\n"))
          } else if (yaml_choice == 2) {
            var <- readline(prompt = "Enter Output File Name: ")
            var = as.character(var)
            var2 <- paste(c(var, "yaml"), collapse=".")
            writeLines(endyaml3, var2)
            return(cat("\nFile should be seen as", var2, "\nGoodbye.\n"))
          } else if (yaml_choice == 3) {
            var <- readline(prompt = "Enter Output File Name: ")
            var = as.character(var)
            var2 <- paste(c(var, "yaml"), collapse=".")
            writeLines(endyaml3, var2)
            return(cat(endyaml3, "\nFile should be seen as", var2, "\nGoodbye.\n"))
          } else (
            return(cat("\nGoodbye.\n"))
          )
}

l2_yaml <- function(gpa) {
  gpa2 <- gpa$l2_models
  ## Model Specifications
  end_yaml <- list()
  mobj <- list()
    for (mm in gpa2$models) {
    name_m <- mm$name
    mobj$diagonal <- mm$contrast_spec$diagonal
    #diag(mm$contrasts) <- ifelse(mobj$diagonal==TRUE, 0, diag(mm$contrasts))
    mobj$cell_means <- mm$contrast_spec$cell_means
    mobj$overall_response <- mm$contrast_spec$overall_response
    mobj$weights <- mm$contrast_spec$weights
    mobj$regressors <- mm$regressors
    mobj$signals <- mm$signals
    cobj <- list()
    for (cc in seq_along(colnames(mm$contrasts))) {
      name_c <- NULL
      name_c <- colnames(mm$contrasts)[cc]
      cobj$row <- rownames(mm$contrasts)[which(mm$contrasts[cc,]!=0)]
      cobj$contrast <- mm$contrasts[cc,which(mm$contrasts[cc,]!=0)]
      mobj$contrasts[[name_c]] <- cobj
    }
    end_yaml$l2_models[[name_m]] <- mobj
    mobj$contrasts <- NULL
  }
  endyaml3 <- as.yaml(end_yaml)
  yaml_choice <- menu(c(
          "Console Output",
          "File Output",
          "Both",
          "Exit"), title= "How would you like to receive the YAML file?")
          if (yaml_choice == 1) {
            return(cat(endyaml3, "\nGoodbye.\n"))
          } else if (yaml_choice == 2) {
            var <- readline(prompt = "Enter Output File Name: ")
            var = as.character(var)
            var2 <- paste(c(var, "yaml"), collapse=".")
            writeLines(endyaml3, var2)
            return(cat("\nFile should be seen as", var2, "\nGoodbye.\n"))
          } else if (yaml_choice == 3) {
            var <- readline(prompt = "Enter Output File Name: ")
            var = as.character(var)
            var2 <- paste(c(var, "yaml"), collapse=".")
            writeLines(endyaml3, var2)
            return(cat(endyaml3, "\nFile should be seen as", var2, "\nGoodbye.\n"))
          } else (
            return(cat("\nGoodbye.\n"))
          )
}

