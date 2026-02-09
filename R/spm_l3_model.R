#' Setup a third-level SPM model using lower-level con images
#'
#' @param l3_df a data.frame containing one row per subject with fields including
#'   id, session, l1_model, l3_model, l1_cope_name, and con_file
#' @param gpa a gpa (glm_pipeline_arguments) object containing model specification
#' @param execute_spm whether to execute SPM setup/estimation immediately. Default: FALSE
#' @param model_type optional override for SPM L3 model type. Options:
#'   "flexible_factorial", "mreg", "one_sample"
#'
#' @importFrom checkmate assert_data_frame assert_subset assert_string assert_class assert_logical
#' @importFrom lgr get_logger
#' @importFrom dplyr left_join
#' @author Michael Hallquist
#' @export
#'
spm_l3_model <- function(l3_df = NULL, gpa, execute_spm = FALSE, model_type = NULL) {
  checkmate::assert_data_frame(l3_df)
  checkmate::assert_subset(c("id", "session", "l1_model", "l3_model", "l1_cope_name", "con_file"), names(l3_df))
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_logical(execute_spm, len = 1L)
  checkmate::assert_string(model_type, null.ok = TRUE)

  lg <- lgr::get_logger("glm_pipeline/l3_setup")
  lg$set_threshold(gpa$lgr_threshold)

  if (length(unique(l3_df$session)) > 1L) {
    msg <- "spm_l3_model is designed for execution on a single session data.frame"
    lg$error(msg)
    stop(msg)
  }

  if (length(unique(l3_df$l1_model)) > 1L) {
    msg <- "spm_l3_model is designed for execution on a single l1 model"
    lg$error(msg)
    stop(msg)
  }

  if (length(unique(l3_df$l3_model)) > 1L) {
    msg <- "spm_l3_model is designed for execution on a single l3 model"
    lg$error(msg)
    stop(msg)
  }

  # metadata
  session <- l3_df$session[1L]
  l1_model <- l3_df$l1_model[1L]
  l3_model <- l3_df$l3_model[1L]
  l1_cope_name <- l1_contrast <- l3_df$l1_cope_name[1L]

  spm_defaults <- list(
    l3_model_type = "flexible_factorial",
    l3_use_one_sample_when_intercept_only = TRUE,
    spm_execute_l3_setup = FALSE,
    spm_execute_l3_glm = FALSE,
    spm_execute_l3_contrasts = FALSE,
    spm_execute_setup = FALSE,
    spm_execute_glm = FALSE,
    spm_execute_contrasts = FALSE,
    estimation_method = "Classical",
    write_residuals = FALSE,
    print_spm_run_instructions = FALSE,
    spm_path = "/gpfs/group/mnh5174/default/lab_resources/spm12",
    matlab_cmd = "matlab",
    matlab_args = "-batch",
    force_l3_creation = FALSE
  )

  spm_settings <- populate_defaults(gpa$glm_settings$spm, spm_defaults)

  if (is.null(model_type)) model_type <- spm_settings$l3_model_type
  if (!model_type %in% c("flexible_factorial", "mreg", "one_sample")) {
    lg$warn("Unknown SPM L3 model_type: %s. Defaulting to multiple regression.", model_type)
    model_type <- "mreg"
  }

  if (isTRUE(execute_spm)) {
    spm_settings$spm_execute_l3_setup <- TRUE
    spm_settings$spm_execute_l3_glm <- TRUE
    spm_settings$spm_execute_l3_contrasts <- TRUE
  }

  # tracking data frame for this model (column names should follow variable names)
  spm_l3_df <- data.frame(
    l1_model = l1_model, l1_cope_name = l1_cope_name,
    l3_model = l3_model
  )

  # respecify model based on available subjects
  mobj <- respecify_l3_model(gpa$l3_models$models[[l3_model]], new_data = l3_df)

  # align l3_df to model metadata order
  l3_df <- dplyr::left_join(mobj$metadata, l3_df, by = c("id", "session"))

  missing_con <- is.na(l3_df$con_file) | !file.exists(l3_df$con_file)
  if (any(missing_con)) {
    lg$warn("Missing con files for %d subjects in spm_l3_model. Dropping.", sum(missing_con))
    l3_df <- l3_df[!missing_con, , drop = FALSE]
    mobj <- respecify_l3_model(gpa$l3_models$models[[l3_model]], new_data = l3_df)
    l3_df <- dplyr::left_join(mobj$metadata, l3_df, by = c("id", "session"))
  }

  con_files <- l3_df$con_file
  if (length(con_files) == 0L) {
    lg$warn("No con files available for SPM L3 model %s. Skipping.", l3_model)
    return(spm_l3_df)
  }

  # output directory for this l3 analysis
  spm_l3_output_dir <- get_output_directory(
    l1_contrast = l1_contrast, l1_model = l1_model,
    l3_model = l3_model, gpa = gpa, what = "l3", glm_software = "spm"
  )

  if (!dir.exists(spm_l3_output_dir)) {
    lg$debug("Creating SPM L3 output directory: %s", spm_l3_output_dir)
    dir.create(spm_l3_output_dir, recursive = TRUE)
  }

  spm_status <- get_spm_status(spm_l3_output_dir, lg = lg)
  force_l3_creation <- isTRUE(spm_settings$force_l3_creation) || isTRUE(spm_status$spm_failed)
  if (isTRUE(spm_status$spm_mat_exists) && isFALSE(force_l3_creation)) {
    lg$info("SPM.mat exists in %s. Skipping SPM L3 setup (force_l3_creation = FALSE).", spm_l3_output_dir)
  } else {
    # determine whether to use one-sample or multiple regression
    regressors <- colnames(mobj$model_matrix)
    intercept_only <- length(regressors) == 1L && regressors[1L] == "(Intercept)"
    if (model_type == "one_sample" ||
      (isTRUE(spm_settings$l3_use_one_sample_when_intercept_only) && intercept_only && model_type != "mreg")) {
      model_type <- "one_sample"
    } else if (model_type == "flexible_factorial") {
      # currently handled via multiple regression to preserve exact design matrix and contrasts
      lg$info("SPM L3 model_type 'flexible_factorial' is currently implemented via multiple regression.")
      model_type <- "mreg"
    }

    spm_preamble <- c(
      ifelse(is.null(spm_settings$spm_path), "", paste0("addpath('", spm_settings$spm_path, "');")),
      "spm('defaults', 'fmri');",
      "spm_jobman('initcfg');",
      ""
    )
    build_matlab_call <- function(script_path) {
      cmd_str <- paste0(
        "try; run('", script_path, "'); ",
        "catch ME; disp(getReport(ME,'extended')); exit(1); end; exit(0);"
      )
      paste(spm_settings$matlab_cmd, spm_settings$matlab_args, shQuote(cmd_str))
    }
    build_shell_call <- function(script_path) {
      compute_env <- get_compute_environment(gpa, c("spm"))
      cmd <- build_matlab_call(script_path)
      if (!is.null(compute_env) && length(compute_env) > 0L) {
        return(paste(c(compute_env, cmd), collapse = " && "))
      }
      cmd
    }

    baseobj <- paste0("matlabbatch{1}.spm.stats.factorial_design")
    m_string <- c(
      spm_preamble,
      "matlabbatch = []; %initialize empty structure",
      paste0(baseobj, ".dir = {'", spm_l3_output_dir, "'};")
    )

    scan_block <- c("{", paste0("'", con_files, ",1'"), "};")
    scan_block_first <- scan_block[1L]
    scan_block_rest <- scan_block[-1L]

    if (model_type == "one_sample") {
      m_string <- c(
        m_string,
        paste0(baseobj, ".des.t1.scans = ", scan_block_first),
        scan_block_rest
      )
    } else {
      # multiple regression: use model_matrix columns as covariates to preserve exact design
      m_string <- c(
        m_string,
        paste0(baseobj, ".des.mreg.scans = ", scan_block_first),
        scan_block_rest
      )
      for (cc in seq_along(regressors)) {
        cname <- make.names(regressors[cc])
        cvec <- mobj$model_matrix[, cc]
        cvals <- paste(format(cvec, scientific = FALSE, digits = 6), collapse = " ")
        m_string <- c(
          m_string,
          paste0(baseobj, ".des.mreg.mcov(", cc, ").c = [ ", cvals, " ];"),
          paste0(baseobj, ".des.mreg.mcov(", cc, ").cname = '", cname, "';"),
          paste0(baseobj, ".des.mreg.mcov(", cc, ").iCC = 5;")
        )
      }
      m_string <- c(m_string, paste0(baseobj, ".des.mreg.incint = 0;"))
    }

    # common design settings
    m_string <- c(
      m_string,
      paste0(baseobj, ".cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});"),
      paste0(baseobj, ".multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});"),
      paste0(baseobj, ".masking.tm.tm_none = 1;"),
      paste0(baseobj, ".masking.im = 1;"),
      paste0(baseobj, ".masking.em = {''};"),
      paste0(baseobj, ".globalc.g_omit = 1;"),
      paste0(baseobj, ".globalm.gmsca.gmsca_no = 1;"),
      paste0(baseobj, ".globalm.glonorm = 1;")
    )

    exec_string <- c(
      spm_preamble,
      ifelse(file.exists(file.path(spm_l3_output_dir, "SPM.mat")), paste0("delete('", file.path(spm_l3_output_dir, "SPM.mat"), "');"), ""),
      paste0("run('", file.path(spm_l3_output_dir, "glm_l3_design.m"), "');"),
      "spm_jobman('run',matlabbatch);"
    )

    cat(m_string, file = file.path(spm_l3_output_dir, "glm_l3_design.m"), sep = "\n")
    cat(exec_string, file = file.path(spm_l3_output_dir, "setup_l3_design.m"), sep = "\n")

    if (isTRUE(spm_settings$spm_execute_l3_setup) || isTRUE(spm_settings$spm_execute_setup)) {
      system(build_shell_call(file.path(spm_l3_output_dir, "setup_l3_design.m")))
    } else {
      cmd <- build_shell_call(file.path(spm_l3_output_dir, "setup_l3_design.m"))
      if (isTRUE(spm_settings$print_spm_run_instructions)) {
        message("To run SPM L3 design setup, execute: ", cmd)
      } else {
        lg$debug("To run SPM L3 design setup, execute: %s", cmd)
      }
    }

    # estimation
    method_label <- normalize_spm_estimation_method(spm_settings$estimation_method)
    write_residuals <- normalize_spm_write_residuals(spm_settings$write_residuals, method_label)

    baseobj_est <- paste0("matlabbatch{1}.spm.stats.fmri_est")
    m_est <- c(
      spm_preamble,
      "matlabbatch = []; %initialize empty structure",
      "% ESTIMATE MODEL",
      paste0(baseobj_est, ".spmmat = { [ '", spm_l3_output_dir, "' filesep 'SPM.mat']};"),
      paste0(baseobj_est, ".write_residuals = ", ifelse(isTRUE(write_residuals), 1, 0), ";"),
      paste0(baseobj_est, ".method.", method_label, " = 1;"),
      "spm_jobman('run',matlabbatch);"
    )
    cat(m_est, file = file.path(spm_l3_output_dir, "run_l3_glm.m"), sep = "\n")

    if (isTRUE(spm_settings$spm_execute_l3_glm) || isTRUE(spm_settings$spm_execute_glm)) {
      system(build_shell_call(file.path(spm_l3_output_dir, "run_l3_glm.m")))
    } else {
      cmd <- build_shell_call(file.path(spm_l3_output_dir, "run_l3_glm.m"))
      if (isTRUE(spm_settings$print_spm_run_instructions)) {
        message("To estimate SPM L3 GLM, execute: ", cmd)
      } else {
        lg$debug("To estimate SPM L3 GLM, execute: %s", cmd)
      }
    }

    # contrasts
    contrast_names <- rownames(mobj$contrasts)
    if (!is.null(contrast_names) && length(contrast_names) > 0L) {
      baseobj_con <- paste0("matlabbatch{1}.spm.stats.con")
      m_con <- c(
        spm_preamble,
        "matlabbatch = []; %initialize empty structure",
        paste0(baseobj_con, ".spmmat = { [ '", spm_l3_output_dir, "' filesep 'SPM.mat']};")
      )

      for (cc in seq_along(contrast_names)) {
        weights <- paste(format(mobj$contrasts[cc, ], scientific = FALSE, digits = 6), collapse = " ")
        m_con <- c(
          m_con,
          paste0(baseobj_con, ".consess{", cc, "}.tcon.name = '", contrast_names[cc], "';"),
          paste0(baseobj_con, ".consess{", cc, "}.tcon.weights = [ ", weights, " ];"),
          paste0(baseobj_con, ".consess{", cc, "}.tcon.sessrep = 'none';")
        )
      }

      m_con <- c(
        m_con,
        paste0(baseobj_con, ".delete = 0;"),
        "spm_jobman('run',matlabbatch);"
      )
      cat(m_con, file = file.path(spm_l3_output_dir, "estimate_l3_contrasts.m"), sep = "\n")

      if (isTRUE(spm_settings$spm_execute_l3_contrasts) || isTRUE(spm_settings$spm_execute_contrasts)) {
        system(build_shell_call(file.path(spm_l3_output_dir, "estimate_l3_contrasts.m")))
      } else {
        cmd <- build_shell_call(file.path(spm_l3_output_dir, "estimate_l3_contrasts.m"))
        if (isTRUE(spm_settings$print_spm_run_instructions)) {
          message("To estimate SPM L3 contrasts, execute: ", cmd)
        } else {
          lg$debug("To estimate SPM L3 contrasts, execute: %s", cmd)
        }
      }
    }
  }

  spm_status <- get_spm_status(spm_l3_output_dir, lg = lg)
  spm_l3_df$spm_dir <- spm_l3_output_dir
  spm_l3_df <- cbind(spm_l3_df, spm_status)
  spm_l3_df$to_run <- !spm_l3_df$spm_complete

  return(spm_l3_df)
}
