#' Internal function to build the data table for AFNI 3dLMEr
#'
#' @param subject_data a data.frame containing subject-level covariates
#' @param input_files a data.frame with columns 'id', 'session', and 'InputFile'
#' @param model_variables a character vector of variables that must be included in the table
#'
#' @return a data.frame formatted for 3dLMEr -dataTable
#' @keywords internal
build_3dlmer_datatable <- function(subject_data, input_files, model_variables) {
  checkmate::assert_data_frame(subject_data)
  checkmate::assert_data_frame(input_files)
  checkmate::assert_character(model_variables)
  checkmate::assert_subset(c("id", "session", "InputFile"), names(input_files))

  subject_data <- subject_data %>%
    dplyr::select(dplyr::all_of(c("id", "session", model_variables))) %>%
    dplyr::distinct()

  input_files <- input_files %>%
    dplyr::select(id, session, InputFile) %>%
    dplyr::distinct()

  # Merge input files with subject data to get covariates
  # subject_data should contain 'id' and 'session' for longitudinal
  dt <- input_files %>%
    dplyr::inner_join(subject_data, by = c("id", "session"))

  # Keep only necessary columns: Subj, model_variables, InputFile
  # 3dLMEr requires the column name 'Subj' for the subject ID
  dt <- dt %>%
    dplyr::rename(Subj = id) %>%
    dplyr::select(dplyr::all_of(c("Subj", "session", model_variables, "InputFile")))

  return(dt)
}

build_3dlmer_intersection_mask <- function(mask_files, output_file, lg = NULL) {
  checkmate::assert_character(mask_files, min.len = 1L)
  checkmate::assert_file_exists(mask_files)
  checkmate::assert_string(output_file)
  checkmate::assert_class(lg, "Logger", null.ok = TRUE)

  mask_files <- unique(mask_files)
  if (is.null(lg)) lg <- lgr::get_logger("fmri.pipeline/afni_3dlmer/mask")
  if (!dir.exists(dirname(output_file))) dir.create(dirname(output_file), recursive = TRUE)

  mask_imgs <- lapply(mask_files, RNifti::readNifti)
  mask_dims <- vapply(mask_imgs, function(img) paste(dim(img), collapse = "x"), character(1))
  if (length(unique(mask_dims)) != 1L) {
    stop(
      sprintf(
        "Cannot build AFNI 3dLMEr intersection mask from mismatched FSL masks: %s",
        paste(unique(mask_dims), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  mask_4d <- do.call(abind::abind, c(mask_imgs, list(along = 4)))
  mask_sum <- rowSums(mask_4d, dims = 3)
  intersect_img <- 1L * (mask_sum == dim(mask_4d)[4L])

  header <- RNifti::niftiHeader(mask_files[1L])
  intersect_nifti <- RNifti::asNifti(intersect_img, header)
  RNifti::writeNifti(intersect_nifti, file = output_file)

  lg$info(
    "Computed AFNI 3dLMEr intersection mask from %d FSL L2 masks: %s",
    length(mask_files), output_file
  )

  output_file
}

resolve_3dlmer_mask <- function(target_dir, harvested_inputs, explicit_mask = NULL, requires_l2 = FALSE, lg = NULL) {
  checkmate::assert_string(target_dir)
  checkmate::assert_data_frame(harvested_inputs)
  checkmate::assert_string(explicit_mask, null.ok = TRUE)
  checkmate::assert_flag(requires_l2)
  checkmate::assert_class(lg, "Logger", null.ok = TRUE)

  if (!is.null(explicit_mask)) return(explicit_mask)
  if (!isTRUE(requires_l2) || !"feat_dir" %in% names(harvested_inputs)) return(NULL)

  mask_files <- unique(file.path(as.character(harvested_inputs$feat_dir), "mask.nii.gz"))
  if (length(mask_files) == 0L || anyNA(mask_files) || any(!file.exists(mask_files))) {
    if (!is.null(lg)) {
      lg$warn("Could not derive AFNI 3dLMEr intersection mask because one or more FSL L2 masks were missing.")
    }
    return(NULL)
  }

  build_3dlmer_intersection_mask(
    mask_files = mask_files,
    output_file = file.path(target_dir, "intersection_mask.nii.gz"),
    lg = lg
  )
}

#' Internal function to construct the 3dLMEr command string
#'
#' @param prefix output prefix for 3dLMEr
#' @param model_formula fixed effects formula string
#' @param qVars character vector of quantitative variables
#' @param glt_codes list of named GLT code strings
#' @param data_table_file path to the data table text file
#' @param mask path to the brain mask file
#' @param njobs number of parallel jobs for 3dLMEr
#' @param ss_type sum of squares type (default 3)
#'
#' @return a character string containing the 3dLMEr command
#' @keywords internal
# Collapse whitespace because AFNI counts spaces inside -model formulas.
sanitize_3dlmer_model_formula <- function(model_formula) {
  checkmate::assert_string(model_formula)
  gsub("[[:space:]]+", "", model_formula)
}

build_3dlmer_command <- function(prefix, model_formula, qVars = NULL, glt_codes = NULL,
                                 data_table_file, mask = NULL, njobs = 1, ss_type = 3) {
  checkmate::assert_string(prefix)
  checkmate::assert_string(model_formula)
  checkmate::assert_character(qVars, null.ok = TRUE)
  checkmate::assert_list(glt_codes, null.ok = TRUE)
  checkmate::assert_string(data_table_file)
  checkmate::assert_string(mask, null.ok = TRUE)
  checkmate::assert_integerish(njobs, lower = 1)
  checkmate::assert_integerish(ss_type, lower = 1, upper = 3)

  model_formula <- sanitize_3dlmer_model_formula(model_formula)

  cmd <- glue::glue("3dLMEr -prefix {prefix} -model '{model_formula}'")

  if (!is.null(qVars) && length(qVars) > 0) {
    cmd <- paste(cmd, glue::glue("-qVars '{paste(qVars, collapse=\",\")}'"))
  }

  if (!is.null(mask)) {
    cmd <- paste(cmd, glue::glue("-mask {mask}"))
  }

  if (njobs > 1) {
    cmd <- paste(cmd, glue::glue("-jobs {njobs}"))
  }

  cmd <- paste(cmd, glue::glue("-SS_type {ss_type}"))

  if (!is.null(glt_codes) && length(glt_codes) > 0) {
    for (i in seq_along(glt_codes)) {
      glt_name <- names(glt_codes)[i]
      glt_code <- glt_codes[[i]]
      # glt_code should already be in AFNI format: "Factor : 1*Level1 -1*Level2"
      cmd <- paste(cmd, glue::glue("-gltCode {glt_name} '{glt_code}'"))
    }
  }

  cmd <- paste(cmd, glue::glue("-dataTable @{data_table_file}"))

  return(cmd)
}

#' Helper to write 3dLMEr script and data table
#'
#' @param output_dir directory where files should be written
#' @param dt the data table data.frame
#' @param cmd the 3dLMEr command string
#' @param script_name name of the shell script to write
#'
#' @return list with paths to script and data table
#' @keywords internal
write_3dlmer_files <- function(output_dir, dt, cmd, script_name = "run_3dlmer.sh") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  dt_file <- file.path(output_dir, "dataTable.txt")
  # write.table with proper formatting for AFNI (space separated, header)
  write.table(dt, file = dt_file, row.names = FALSE, quote = FALSE, sep = "\t")

  script_file <- file.path(output_dir, script_name)
  # Basic script structure
  script_content <- c(
    "#!/bin/bash",
    "",
    "# Load AFNI if needed (placeholder, setup_l3_run will handle compute env)",
    "script_dir=\"$(cd \"$(dirname \"$0\")\" && pwd)\"",
    "cd \"$script_dir\"",
    "",
    cmd
  )
  writeLines(script_content, script_file)
  Sys.chmod(script_file, "0755")

  return(list(script = script_file, data_table = dt_file))
}

#' Status checker for 3dLMEr
#'
#' @param output_file expected output file from 3dLMEr
#' @param lg optional logger
#' @param prefix optional prefix for returned column names
#' @return data.frame indicating output existence and completion status
#' @keywords internal
get_3dlmer_status <- function(output_file, lg = NULL, prefix = NULL) {
  checkmate::assert_string(output_file)
  checkmate::assert_class(lg, "Logger", null.ok = TRUE)
  checkmate::assert_string(prefix, null.ok = TRUE)
  if (is.null(lg)) lg <- lgr::get_logger()

  resolved_output <- output_file
  candidates <- unique(c(
    output_file,
    paste0(output_file, ".HEAD"),
    paste0(output_file, ".BRIK"),
    paste0(output_file, "+tlrc.HEAD"),
    paste0(output_file, "+orig.HEAD")
  ))

  if (grepl("\\+tlrc$|\\+orig$", output_file)) {
    candidates <- unique(c(
      output_file,
      paste0(output_file, ".HEAD"),
      paste0(output_file, ".BRIK")
    ))
  } else if (grepl("\\.(nii|nii\\.gz|HEAD|BRIK)$", output_file, perl = TRUE)) {
    sans_ext <- file_sans_ext(output_file)
    if (!is.na(sans_ext)) {
      candidates <- unique(c(
        output_file,
        paste0(sans_ext, "+tlrc.HEAD"),
        paste0(sans_ext, "+orig.HEAD")
      ))
    }
  }

  existing <- candidates[file.exists(candidates)]
  output_exists <- length(existing) > 0L
  if (output_exists) {
    resolved_output <- existing[1L]
  } else {
    lg$debug("No AFNI 3dLMEr output found for '%s'. Checked: %s", output_file, paste(candidates, collapse = ", "))
  }

  out <- data.frame(
    output_file = output_file,
    output_file_resolved = resolved_output,
    output_file_exists = output_exists,
    feat_complete = output_exists,
    feat_failed = if (output_exists) FALSE else NA,
    stringsAsFactors = FALSE
  )
  if (!is.null(prefix)) names(out) <- paste0(prefix, names(out))
  out
}

#' Translator from emmeans-style contrast spec to 3dLMEr -gltCode
#'
#' @param mobj high-level model specification object (hi_model_spec)
#' @param data the data.frame used for the model (to check factor levels)
#'
#' @return a list of 3dLMEr gltCode strings
#' @keywords internal
emmeans_to_3dlmer_glt <- function(mobj, data) {
  checkmate::assert_class(mobj, "hi_model_spec")
  checkmate::assert_data_frame(data)
  lg <- lgr::get_logger("fmri.pipeline/afni_3dlmer/emmeans_to_3dlmer_glt")

  if (is.null(mobj$lmfit)) {
    lg$warn("No lmfit in mobj. Cannot generate emmeans-to-3dlmer GLTs.")
    return(NULL)
  }

  spec <- mobj$contrast_spec
  if (is.null(spec)) return(NULL)

  glt_list <- list()

  # 1. Diagonal contrasts
  if (isTRUE(spec$diagonal)) {
    # 3dLMEr can handle simple regressors by name if they are in the dataTable
    # But for factors, it's better to use GLTs.
    # Actually, for 3dLMEr, 'diagonal' usually means we want a test for each beta.
    # For now, let's focus on the emmeans-based ones as they are the "hard work".
  }

  # 2. EMMeans-based contrasts (cond_means, pairwise_diffs, etc.)
  # We'll reproduce the emmeans calls from get_contrasts_from_spec
  # but instead of taking @linfct, we'll build strings.

  process_emm <- function(vv, type = "means") {
    # vv is a variable or interaction, e.g., "Group" or "Group:Session"
    ee <- emmeans::emmeans(mobj$lmfit, as.formula(paste("~", vv)), weights = spec$weights)
    edata <- summary(ee)
    
    # Identify factor columns
    fac_cols <- strsplit(vv, ":")[[1]]
    
    if (type == "means") {
      # One GLT per row of edata
      for (i in seq_len(nrow(edata))) {
        # Name: e.g., Group.patient or Group.Session.patient.post
        con_name <- paste(vv, paste(unname(unlist(edata[i, fac_cols])), collapse="."), sep=".")
        
        # Code: "Group : 1*patient" or "Group : 1*patient Session : 1*post"
        code_parts <- c()
        for (f in fac_cols) {
          code_parts <- c(code_parts, paste0(f, " : 1*", edata[i, f]))
        }
        glt_list[[con_name]] <<- paste(code_parts, collapse = " ")
      }
    } else if (type == "pairwise") {
      # Pairwise differences
      pw <- pairs(ee)
      pw_data <- summary(pw)
      pw_linfct <- pw@linfct # This is against the emmeans grid rows
      
      grid_data <- summary(ee)
      
      for (i in seq_len(nrow(pw_data))) {
        con_name <- as.character(pw_data$contrast[i])
        # Clean up name for AFNI (no spaces, etc.)
        con_name <- gsub("[^[:alnum:]_]", ".", con_name)
        con_name <- paste0(vv, ".pw.", con_name)
        
        # Build code from non-zero weights in pw_linfct[i,]
        weights <- pw_linfct[i, ]
        non_zero <- which(abs(weights) > 1e-8)
        
        code_parts <- c()
        # 3dLMEr GLT syntax for multi-factor interaction contrasts can be tricky.
        # If it's a single factor, it's "Factor : 1*Level1 -1*Level2"
        # If it's an interaction, it's harder.
        
        if (length(fac_cols) == 1) {
          # Simple case
          f <- fac_cols[1]
          term_parts <- c()
          for (j in non_zero) {
            term_parts <- c(term_parts, paste0(weights[j], "*", grid_data[j, f]))
          }
          glt_list[[con_name]] <<- paste0(f, " : ", paste(term_parts, collapse = " "))
        } else {
          # Interaction pairwise: "Factor1 : 1*F1L1 -1*F1L2 Factor2 : 1*F2L1" etc.
          # Actually, emmeans pairwise on interaction is usually (F1L1 F2L1) - (F1L2 F2L1) etc.
          # We can represent any linear combination of grid rows:
          # "Factor1 : 1*F1L1 Factor2 : 1*F2L1 ++ Factor1 : -1*F1L2 Factor2 : 1*F2L1"
          # Wait, 3dLMEr supports weights on combinations:
          # "Factor1 : 1*F1L1 Factor2 : 1*F2L1"
          # Let's use the '++' syntax or just combine if they share factors.
          
          # More robust way for 3dLMEr:
          # Each non-zero row in the grid gets its own "Factor1 : W*L1 Factor2 : 1*L2" snippet
          full_code <- c()
          for (j in non_zero) {
            snippet <- c()
            for (f in fac_cols) {
              # The weight only needs to be on one factor's level for that combination
              val <- if (f == fac_cols[1]) weights[j] else 1
              snippet <- c(snippet, paste0(f, " : ", val, "*", grid_data[j, f]))
            }
            full_code <- c(full_code, paste(snippet, collapse = " "))
          }
          glt_list[[con_name]] <<- paste(full_code, collapse = " ")
        }
      }
    }
  }

  if (length(spec$cond_means) > 0L) {
    for (vv in spec$cond_means) process_emm(vv, "means")
  }
  
  if (length(spec$pairwise_diffs) > 0L) {
    for (vv in spec$pairwise_diffs) process_emm(vv, "pairwise")
  }

  # 3. Custom contrasts from contrast_list$custom
  if (!is.null(mobj$contrast_list$custom)) {
    # Custom contrasts are against regressors (betas).
    # 3dLMEr supports this too: "-gltCode mycon 'reg1 : 1 reg2 : -1'"
    cmat <- mobj$contrast_list$custom
    for (i in seq_len(nrow(cmat))) {
      con_name <- rownames(cmat)[i]
      weights <- cmat[i, ]
      non_zero <- which(abs(weights) > 1e-8)
      
      code_parts <- c()
      for (j in non_zero) {
        # 3dLMEr syntax for quantitative/beta contrasts: "Regname : Weight"
        code_parts <- c(code_parts, paste0(names(weights)[j], " : ", weights[j]))
      }
      glt_list[[con_name]] <- paste(code_parts, collapse = " ")
    }
  }

  return(glt_list)
}
