#' Fit a linear mixed model with Julia MixedModels
#'
#' \code{jlmer()} is an \code{lme4::lmer()}-style wrapper around Julia's
#' \code{MixedModels.fit()}. It fits the model in Julia, then converts the result
#' back to an R \code{lmerMod} object through JellyMe4 so downstream R tooling
#' such as \code{broom.mixed}, \code{emmeans}, and \code{lme4} methods can be
#' used at the call site.
#'
#' @param formula A two-sided mixed model formula.
#' @param data A data frame containing variables used in \code{formula}.
#' @param REML Logical scalar. If \code{TRUE}, estimate with restricted maximum
#'   likelihood; otherwise use maximum likelihood.
#' @param JULIA_HOME Optional path to the Julia installation. Passed to
#'   \code{JuliaCall::julia_setup()}.
#' @param na.action Optional missing-data action for formula variables. Supports
#'   \code{stats::na.omit}, \code{stats::na.exclude}, \code{stats::na.fail}, or
#'   \code{NULL}. Omitted or excluded rows are removed before sending data to
#'   Julia.
#' @param lmer_test Logical scalar. If \code{TRUE}, convert the returned
#'   \code{lmerMod} object to \code{lmerTest}'s \code{lmerModLmerTest} class.
#'   This matches the default \code{mixed_by(engine = "lme4")} call site, which
#'   uses \code{lmerTest::lmer()}.
#' @param lmer_test_tol Numeric tolerance passed to
#'   \code{lmerTest::as_lmerModLmerTest()} when \code{lmer_test = TRUE}.
#' @param julia_setup_args Named list of additional arguments passed to
#'   \code{JuliaCall::julia_setup()}.
#' @param julia_packages Character vector of Julia packages to load. Defaults to
#'   \code{c("MixedModels", "RCall", "JellyMe4")}.
#' @param fit_args Named list of additional keyword arguments passed to Julia's
#'   \code{fit(MixedModel, ...)} call.
#' @param keep_julia_model Logical scalar. If \code{TRUE}, the temporary Julia
#'   model variable is left in Julia's \code{Main} module and its name is stored
#'   in \code{attr(result, "jlmer")}. By default all temporary Julia variables
#'   are cleared after conversion to R.
#' @param verbose Logical scalar. If \code{TRUE}, print the Julia fit call.
#' @param ... Additional named keyword arguments passed to Julia's
#'   \code{fit(MixedModel, ...)} call. These are combined with \code{fit_args}.
#'
#' @return An \code{lmerMod} object converted by JellyMe4, with a \code{"jlmer"}
#'   attribute containing backend metadata. By default, the object is further
#'   converted to \code{lmerTest}'s \code{lmerModLmerTest} class.
#'
#' @details
#' This function requires Julia packages \code{MixedModels}, \code{RCall}, and
#' \code{JellyMe4} in the active Julia environment. It does not install Julia
#' packages automatically.
#'
#' @export
jlmer <- function(formula,
                  data,
                  REML = TRUE,
                  JULIA_HOME = NULL,
                  na.action = stats::na.omit,
                  lmer_test = TRUE,
                  lmer_test_tol = 1e-8,
                  julia_setup_args = list(),
                  julia_packages = c("MixedModels", "RCall", "JellyMe4"),
                  fit_args = list(),
                  keep_julia_model = FALSE,
                  verbose = FALSE,
                  ...) {
  if (!requireNamespace("JuliaCall", quietly = TRUE)) {
    stop("Package \"JuliaCall\" must be installed to use jlmer().", call. = FALSE)
  }

  if (is.character(formula)) {
    formula <- stats::as.formula(formula)
  }
  checkmate::assert_formula(formula)
  if (length(formula) != 3L) {
    stop("formula must be a two-sided mixed model formula.", call. = FALSE)
  }
  checkmate::assert_data_frame(data)
  data <- as.data.frame(data)
  checkmate::assert_flag(REML)
  checkmate::assert_string(JULIA_HOME, null.ok = TRUE)
  if (!is.null(na.action) && !is.function(na.action)) {
    stop("na.action must be a function such as stats::na.omit, stats::na.exclude, stats::na.fail, or NULL.", call. = FALSE)
  }
  checkmate::assert_flag(lmer_test)
  checkmate::assert_number(lmer_test_tol, lower = 0, finite = TRUE)
  checkmate::assert_list(julia_setup_args, names = "unique")
  checkmate::assert_character(julia_packages, min.len = 1L, any.missing = FALSE)
  checkmate::assert_list(fit_args, names = "unique")
  checkmate::assert_flag(keep_julia_model)
  checkmate::assert_flag(verbose)

  dot_args <- list(...)
  if (length(dot_args) > 0L) {
    if (is.null(names(dot_args)) || any(names(dot_args) == "")) {
      stop("Additional arguments in ... must be named Julia fit() keyword arguments.", call. = FALSE)
    }
    fit_args <- c(fit_args, dot_args)
    if (anyDuplicated(names(fit_args))) {
      stop("Duplicate Julia fit() keyword arguments were supplied.", call. = FALSE)
    }
  }

  if (length(fit_args) > 0L) {
    valid_julia_names <- grepl("^[A-Za-z_][A-Za-z0-9_]*$", names(fit_args))
    if (!all(valid_julia_names)) {
      stop(
        "All Julia fit() keyword arguments must be valid Julia identifiers. Invalid: ",
        paste(names(fit_args)[!valid_julia_names], collapse = ", "),
        call. = FALSE
      )
    }
  }

  model_vars <- all.vars(formula)
  missing_vars <- setdiff(model_vars, names(data))
  if (length(missing_vars) > 0L) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }

  omitted_rows <- integer()
  if (!is.null(na.action)) {
    model_data <- data[, model_vars, drop = FALSE]
    incomplete_rows <- which(!stats::complete.cases(model_data))
    if (length(incomplete_rows) > 0L) {
      if (identical(na.action, stats::na.fail)) {
        stop("Missing values in model variables and na.action = stats::na.fail.", call. = FALSE)
      }
      data <- data[-incomplete_rows, , drop = FALSE]
      omitted_rows <- incomplete_rows
    }
  }

  julia_setup_call <- c(list(JULIA_HOME = JULIA_HOME), julia_setup_args)
  julia_setup_call <- julia_setup_call[!vapply(julia_setup_call, is.null, logical(1L))]
  tryCatch(
    do.call(JuliaCall::julia_setup, julia_setup_call),
    error = function(e) {
      stop("Unable to initialize Julia with JuliaCall::julia_setup(): ", conditionMessage(e), call. = FALSE)
    }
  )

  julia_using_call <- sprintf("using %s", paste(julia_packages, collapse = ", "))
  tryCatch(
    JuliaCall::julia_command(julia_using_call),
    error = function(e) {
      stop(
        "Unable to load required Julia package(s): ",
        paste(julia_packages, collapse = ", "),
        ". Install them in the active Julia environment before using jlmer(). Original error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  prefix <- paste0("fmri_pipeline_jlmer_", as.integer(stats::runif(1L, 1e8, 1e9 - 1L)))
  data_var <- paste0(prefix, "_data")
  formula_var <- paste0(prefix, "_formula")
  model_var <- paste0(prefix, "_model")
  fit_arg_vars <- character()

  cleanup_vars <- function() {
    vars <- c(data_var, formula_var, fit_arg_vars)
    if (!isTRUE(keep_julia_model)) {
      vars <- c(vars, model_var)
    }
    cleanup_call <- paste0(vars, " = nothing", collapse = "; ")
    try(JuliaCall::julia_command(paste0(cleanup_call, "; GC.gc()")), silent = TRUE)
  }
  on.exit(cleanup_vars(), add = TRUE)

  JuliaCall::julia_assign(data_var, data)
  JuliaCall::julia_assign(formula_var, formula)

  fit_arg_expr <- character()
  if (length(fit_args) > 0L) {
    fit_arg_expr <- vapply(names(fit_args), function(arg_name) {
      arg_var <- paste0(prefix, "_arg_", arg_name)
      fit_arg_vars <<- c(fit_arg_vars, arg_var)
      JuliaCall::julia_assign(arg_var, fit_args[[arg_name]])
      paste0(arg_name, "=", arg_var)
    }, character(1L))
  }

  reml_arg <- paste0("REML=", if (isTRUE(REML)) "true" else "false")
  fit_call <- sprintf(
    "%s = fit(MixedModel, %s, %s, %s);",
    model_var,
    formula_var,
    data_var,
    paste(c(reml_arg, fit_arg_expr), collapse = ", ")
  )
  if (isTRUE(verbose)) {
    message("Julia jlmer call: ", fit_call)
  }

  tryCatch(
    JuliaCall::julia_command(fit_call),
    error = function(e) {
      stop("Julia MixedModels fit failed: ", conditionMessage(e), call. = FALSE)
    }
  )
  model <- tryCatch(
    JuliaCall::julia_eval(
      sprintf("robject(:lmerMod, (%s, %s));", model_var, data_var),
      need_return = "R"
    ),
    error = function(e) {
      stop("JellyMe4 conversion to an R lmerMod failed: ", conditionMessage(e), call. = FALSE)
    }
  )

  if (isTRUE(lmer_test)) {
    if (!requireNamespace("lmerTest", quietly = TRUE)) {
      stop("Package \"lmerTest\" must be installed when lmer_test = TRUE.", call. = FALSE)
    }
    model <- tryCatch(
      lmerTest::as_lmerModLmerTest(model, tol = lmer_test_tol),
      error = function(e) {
        stop("lmerTest conversion failed: ", conditionMessage(e), call. = FALSE)
      }
    )
  }

  attr(model, "jlmer") <- list(
    engine = "jlmer",
    REML = REML,
    lmer_test = lmer_test,
    lmer_test_tol = if (isTRUE(lmer_test)) lmer_test_tol else NA_real_,
    julia_model_var = if (isTRUE(keep_julia_model)) model_var else NA_character_,
    omitted_rows = omitted_rows,
    fit_args = fit_args
  )

  model
}
