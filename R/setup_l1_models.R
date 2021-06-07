#' Function for setting up level 1 model specifications and corresponding files
#'
#' @param gpa A glm_model_arguments function setup by the \code{setup_glm_pipeline} function
#' @param to_setup A character vector of model names within \code{gpa$l1_models} whose l1
#'    inputs should be generated. If omitted, all models within \code{gpa$l1_models} will be generated.
#'
#' @details The \code{to_setup} argument allows the creation of multiple l1 models to be parallelized at
#'   a superordinate level or, even better, to be spread across independent jobs on a cluster. This function
#'   already provides the option to parallelize over subjects for a single model if \code{gpa$l1_setup_cpus}
#'   is greater than 1.
#'
#' @importFrom checkmate assert_class assert_character
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach registerDoSEQ
#' @importFrom parallel makeCluster stopCluster
#' @importFrom RNifti niftiHeader
#' @importFrom dplyr filter pull
#' @importFrom rlang sym
#' @importFrom magrittr %>%
setup_l1_models <- function(gpa, to_setup=NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_character(to_setup, null.ok=TRUE)
  checkmate::assert_data_frame(gpa$subject_data)

  #if no model subset is requested, output all models
  if (is.null(to_setup)) { to_setup <- names(gpa$l1_models$models) }

  lg <- lgr::get_logger("glm_pipeline/l1_setup")
  if (isTRUE(gpa$log_txt)) { lg$add_appender(lgr::AppenderFile$new("setup_l1_models.txt"), name="txt") }

  #setup parallel worker pool, if requested
  if (gpa$parallel$l1_setup_cores > 1L) {
    lg$info("Initializing l1 setup cluster with %d cores", gpa$parallel$l1_setup_cores)
    cl <- parallel::makeCluster(gpa$parallel$l1_setup_cores)
    doParallel::registerDoParallel(cl)
    on.exit(try(parallel::stopCluster(cl))) #cleanup pool upon exit of this function
  } else {
    lg$info("Initializing l1 setup with serial execution")
    foreach::registerDoSEQ() #formally register a sequential 'pool' so that dopar is okay
  }

  #FOR TESTING
  gpa$subject_data <- gpa$subject_data[1:3,]

  # loop over each subject, identify relevant fMRI data, and setup level 1 analysis files
  all_subj_l1_list <- foreach(subj_df = iter(gpa$subject_data, by="row"), .inorder=FALSE, .packages=c("dependlab", "dplyr"),
    .export=c("lg", "gpa", "truncateRuns", "fsl_l1_model", "spm_l1_model", "runFSLCommand")) %dopar% {
      subj_df <- subj_df # avoid complaints about visible global binding in R CMD check
      subj_id <- subj_df$id
      subj_session <- subj_df$session

      #find the run data for analysis
      rdata <- gpa$run_data %>% dplyr::filter(id == !!subj_id & session == !!subj_session & run_nifti_present == TRUE)
      run_nifti <- get_mr_abspath(rdata, "run_nifti")
      mr_run_nums <- rdata %>% dplyr::pull("run_number")

      subj_mr_dir <- subj_df$mr_dir

      # process files
      if (length(run_nifti) == 0L) {
        lg$warn("Unable to find any preprocessed fMRI files in dir: %s", subj_mr_dir)
        return(NULL)
      } else if (any(!file.exists(run_nifti))) {
        lg$warn("Could not find some of the expected preprocessed fMRI files. These will be dropped.")
        lg$warn("Missing: %s", run_nifti[!file.exists(run_nifti)])
        run_nifti <- run_nifti[file.exists(run_nifti)]
      } else {
        lg$debug(paste("MR files to analyze:", run_nifti)) #log files that were found
      }

      ## read number of volumes from NIfTI header
      # RNifti is unexpectedly slow compared to oro.nifti
      mr_dims <- lapply(run_nifti, function(x) { oro.nifti::readNIfTI(x, read_data = FALSE)@dim_ })
      run_lengths <- sapply(mr_dims, "[[", 5) # number of volumes is 4th dimension
      nvoxels <- sapply(mr_dims, function(x) { prod(x[2:5]) })

      #we also need xyz to get number of voxels
      #run_lengths <- unname(sapply(run_nifti, function(x) { oro.nifti::readNIfTI(x, read_data=FALSE)@dim_[5L] }))
      lg$debug("Run lengths of run_nifti: %s", paste(run_lengths, collapse=", "))

      ## create truncated run files to end analysis 12s after last ITI (or big head movement)
      ## also handle removal of N volumes from the beginning of each run due to steady state magnetization

      #mrdf <- truncate_runs(b, run_nifti, mr_run_nums, run_lengths, drop_volumes=drop_volumes)
      mrdf <- data.frame(
        id = subj_id, session=subj_session,
        run_nifti = run_nifti, run_number = mr_run_nums,
        last_volume = run_lengths, drop_volumes = gpa$drop_volumes
      )

      run_nifti <- mrdf$run_nifti
      run_lengths <- mrdf$last_volume

      # determine whether to include each run
      mrdf$exclude_run <- sapply(seq_len(nrow(mrdf)), function(rr) {
        ll <- as.list(mrdf[rr, , drop = FALSE]) # rrth row of mrdf
        ll[["gpa"]] <- gpa
        ll[["run_nifti"]] <- NULL
        ex <- do.call(get_l1_confounds, ll)
        return(ex$exclude_run)
      })

      # Tracking list containing data.frames for each software, where we expect one row per run-level model (FSL)
      # or subject-level model (AFNI). The structure varies because FSL estimates per-run GLMs, while AFNI concatenates.
      l1_file_setup <- list(fsl = list(), spm = list(), afni = list(), metadata = mrdf)

      #loop over models to output
      for (ii in seq_along(to_setup)) {
        this_model <- to_setup[ii]
        if (isTRUE(gpa$log_json)) { lg$add_appender(lgr::AppenderJson$new(paste0(gpa$l1_setup_log[ii], ".json")), name="json") }
        if (isTRUE(gpa$log_txt)) { lg$add_appender(lgr::AppenderFile$new(paste0(gpa$l1_setup_log[ii], ".txt")), name="txt") }

        # setup design matrix for any given software package
        m_events <- data.table::rbindlist(
          lapply(gpa$l1_models$events, function(this_event) {
            this_event %>% dplyr::filter(id == !!subj_id & session == !!subj_session)
          })
        )

        m_signals <- lapply(gpa$l1_models$signals[gpa$l1_models$models[[this_model]]$signals], function(this_signal) {
          # filter down to this id if the signal is a data.frame
          if (inherits(this_signal$value, "data.frame")) {
            this_signal$value <- this_signal$value %>% dplyr::filter(id == !!subj_id & session == !!subj_session)
          }
          return(this_signal)
        })

        subj_out <- file.path(subj_mr_dir, gpa$l1_models$models[[this_model]]$outdir)

        if (!dir.exists(subj_out)) {
          lg$info("Creating subject output directory: %s", subj_out)
          dir.create(subj_out, showWarnings=FALSE, recursive=TRUE)
        }

        t_out <- gpa$glm_software
        if (isTRUE(gpa$use_preconvolve)) { t_out <- c("convolved", t_out) } #compute preconvolved regressors
        bdm_args <- gpa$additional$bdm_args
        bdm_args$events <- m_events
        bdm_args$signals <- m_signals
        bdm_args$tr <- gpa$tr
        bdm_args$write_timing_files <- t_out
        bdm_args$drop_volumes <- gpa$drop_volumes
        bdm_args$run_volumes <- run_lengths
        bdm_args$run_4d_files <- run_nifti
        bdm_args$runs_to_output <- mr_run_nums
        bdm_args$output_directory <- file.path(subj_out, "timing_files")
        d_obj <- tryCatch(do.call(build_design_matrix, bdm_args), error=function(e) {
          lg$error("Failed build_design_matrix for id: %s, session: %s, model: %s", subj_id, subj_session, this_model)
          lg$error("Error message: %s", as.character(e))
          return(NULL)
        })

        if (is.null(d_obj)) { next } #skip to next iteration on error

        save(d_obj, bdm_args, mrdf, mr_run_nums, subj_mr_dir, run_nifti, run_lengths, subj_id, subj_session, this_model,
          file=file.path(subj_out, paste0(gpa$l1_models$models[[this_model]]$name, "_bdm_setup.RData")))

        if ("fsl" %in% gpa$glm_software) {
          #Setup FSL run-level models for each combination of signals
          #Returns a data.frame of feat l1 inputs and the fsf file
          feat_l1_df <- tryCatch({fsl_l1_model(id=subj_id, session=subj_session, d_obj=d_obj,
          gpa = gpa, this_model, nvoxels = nvoxels, lg=lg)},
            error=function(e) {
              lg$error("Problem with fsl_l1_model. Model: %s, Subject: %s, Session: %s", this_model, subj_id, subj_session)
              lg$error("Error message: %s", as.character(e))
              return(NULL)
            })

          if (!is.null(feat_l1_df)) {
            #add to tracking data.frame (simple, inefficient append)
            l1_file_setup$fsl <- rbind(l1_file_setup$fsl, feat_l1_df)
          }
        }

        #TODO: finalize SPM approach
        if ("spm" %in% gpa$glm_software) {
          #Setup spm run-level models for each combination of signals
          spm_files <- tryCatch(spm_l1_model(d_obj, gpa, this_model, run_nifti),
            error=function(e) {
              lg$error("Problem running spm_l1_model. Model: %s, Subject: %s, Session: %s", this_model, subj_id, subj_session)
              lg$error("Error message: %s", as.character(e))
              return(NULL)
            })

        }
      }

      lg$info("Completed processing of subject: %s", subj_id)
      return(l1_file_setup)
    }

  # N.B. rbindlist converts the bound elements into a single data.table object
  all_subj_l1_combined <- list(
    fsl=rbindlist(lapply(all_subj_l1_list, "[[", "fsl")),
    metadata=rbindlist(lapply(all_subj_l1_list, "[[", "metadata"))
    #spm = rbindlist(lapply(all_subj_l1_list, "[[", "spm"))
  )

  class(all_subj_l1_combined) <- c("list", "l1_setup")

  #append l1 setup to gpa
  gpa$l1_model_setup <- all_subj_l1_combined

  return(gpa)

}
