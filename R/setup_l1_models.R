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
    on.exit(try(stopCluster(cl))) #cleanup pool upon exit of this function
  } else {
    lg$info("Initializing l1 setup with serial execution")
    foreach::registerDoSEQ() #formally register a sequential 'pool' so that dopar is okay
  }

  idvec <- gpa$subject_data[[ gpa$vm["id"] ]]

  # loop over each subject, identify relevant fMRI data, and setup level 1 analysis files
  ll <- foreach(subid = iter(idvec), .inorder=FALSE, .packages=c("dependlab", "dplyr"),
    .export=c("lg", "gpa", "truncateRuns", "fsl_l1_model", "spm_l1_model", "runFSLCommand")) %dopar% {
      subid <- subid # avoid complaints about visible global binding in R CMD check
      subj_mr_dir <- gpa$subject_data %>%
        dplyr::filter(!!sym(gpa$vm["id"]) == !!subid) %>%
        pull(!!gpa$vm["mr_dir"])

      #use specific run NIfTIs included in run_data, rather finding these by regex
      if (gpa$vm["run_nifti"] %in% names(gpa$run_data)) {
        lg$info("Using run_data to identify NIfTI files for analysis")
        subj_runs <- gpa$run_data %>% filter(!!sym(gpa$vm["id"]) == !!subid)
        if (nrow(subj_runs) != 1L) {
          print(subj_runs);
          log$fatal("Found an unexpected number of runs (%d) for subject %s", nrow(subj_runs), subid)
          stop("Unexpected number of runs in subj_runs")
        }

        mrfiles <- subj_runs[[ gpa$vm["run_nifti"] ]]
        mr_run_nums <- subj_runs[[ gpa$vm["run_number"] ]]

        mr_found <- file.exists(mrfiles)
        if (any(mr_found != TRUE)) {
          lg$debug(
            "Cannot find some files: %s. Prepend mr_dir and try again.",
            paste(mrfiles[!mr_found], collapse = ", ")
          )
          mrdirs <- subj_runs[[ gpa$vm["mr_dir"] ]]

          #try prepending the mr_dir path if run_nifti is relative
          mrfiles[!mr_found] <- file.path(mrdirs[!mr_found], mrfiles[!mr_found])
        }

        mr_found <- file.exists(mrfiles)
        if (any(mr_found != TRUE)) {
          lg$warn(
            "Could not find the following run files: %s. Dropping from analysis",
            paste(mrfiles[!mr_found], collapse = ", ")
          )
        }

        #pass forward variables for analysis
        mrfiles <- mrfiles[mr_found]
        mr_run_nums <- mr_run_nums[mr_found]
      } else {
        lg$info("Using regex-based find approach to identify run NIfTIs")

        #find run nifti files based on directory and regular expression settings

        if (length(subj_mr_dir) != 1L) {
          lg$warn("Unable to find fMRI directory record in subject_data for subid: %s", subid)
          return(NULL)
        }

        if (!dir.exists(file.path(subj_mr_dir))) {
          lg$warn("Unable to find subject data directory: %s for subid: %s", subj_mr_dir, subid)
          return(NULL)
        }

        ## Find processed fMRI run-level data for this subject
        # mrfiles <- list.files(subj_mr_dir, pattern=gpa$fmri_file_regex, full.names=TRUE, recursive=TRUE)
        ## cat(paste0("command: find ", subj_mr_dir, " -iname '", expectfile, "' -ipath '*", expectdir, "*' -type f\n"))

        # -ipath '*", expectdir, "*' -type f | sort -n"), intern=TRUE)
        if (!is.null(gpa$fmri_path_regex)) {
          addon <- paste0(" -ipath '*/", gpa$fmri_path_regex, "/*'")
        } else {
          addon <- ""
        }
        find_string <- paste0(
          "find ", subj_mr_dir, " -regextype posix-egrep -iregex '.*",
          gpa$fmri_file_regex, "'", addon, " -type f | sort -n"
        )
        lg$debug("mrfiles find syntax: %s", find_string)
        mrfiles <- system(find_string, intern=TRUE)

        #extract run number from file name
        mr_run_nums <- as.integer(sub(paste0(gpa$run_number_regex), "\\1", mrfiles, perl=TRUE))

        ## NB. If we reorder the mrfiles, then the run numbers diverge unless we sort(mr_run_nums). 
        ## Remove for now for testing
        ## mrfiles <- mrfiles[order(mr_run_nums)] #make absolutely sure that runs are ordered ascending
      }

      #process files
      if (length(mrfiles) == 0L) {
        lg$warn("Unable to find any preprocessed fMRI files in dir: %s", subj_mr_dir)
        return(NULL)
      } else {
        lg$debug(paste("MR files to analyze:", mrfiles[mr_found])) #log files that were found
      }

      ## read number of volumes from NIfTI header
      # RNifti is unexpectedly slow compared to oro.nifti
      mr_dims <- lapply(mrfiles, function(x) { oro.nifti::readNIfTI(x, read_data = FALSE)@dim_ })
      run_lengths <- sapply(mr_dims, "[[", 5) # number of volumes is 4th dimension
      nvoxels <- sapply(mr_dims, function(x) { prod(x[2:5]) })

      #we also need xyz to get number of voxels
      #run_lengths <- unname(sapply(mrfiles, function(x) { oro.nifti::readNIfTI(x, read_data=FALSE)@dim_[5L] }))
      lg$debug("Run lengths of mrfiles: %s", paste(run_lengths, collapse=", "))

      ## create truncated run files to end analysis 12s after last ITI (or big head movement)
      ## also handle removal of N volumes from the beginning of each run due to steady state magnetization

      #mrdf <- truncate_runs(b, mrfiles, mr_run_nums, run_lengths, drop_volumes=drop_volumes)
      mrdf <- data.frame(
        mrfile_to_analyze = mrfiles, run = mr_run_nums,
        last_vol_analysis = run_lengths, drop_volumes = gpa$drop_volumes
      )

      mrfiles <- mrdf$mrfile_to_analyze
      run_lengths <- mrdf$last_vol_analysis

      # Tracking list containing data.frames for each software, where we expect one row per run-level model (FSL)
      # or subject-level model (AFNI). The structure varies because FSL estimates per-run GLMs, while AFNI concatenates.
      l1_file_setup <- list(fsl=list(), spm=list(), afni=list())

      #loop over models to output
      for (ii in seq_along(to_setup)) {
        this_model <- to_setup[ii]
        if (isTRUE(gpa$log_json)) { lg$add_appender(lgr::AppenderJson$new(paste0(gpa$l1_setup_log[ii], ".json")), name="json") }
        if (isTRUE(gpa$log_txt)) { lg$add_appender(lgr::AppenderFile$new(paste0(gpa$l1_setup_log[ii], ".txt")), name="txt") }

        # setup design matrix for any given software package
        m_events <- data.table::rbindlist(
          lapply(gpa$l1_models$events, function(this_event) {
            this_event %>% filter(!!sym(gpa$vm["id"]) == !!subid)
          })
        )

        m_signals <- lapply(gpa$l1_models$signals[gpa$l1_models$models[[this_model]]$model_signals], function(this_signal) {
          # filter down to this id if the signal is a data.frame
          if (inherits(this_signal$value, "data.frame")) {
            this_signal$value <- this_signal$value %>% filter(!!sym(gpa$vm["id"]) == !!subid)
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
        bdm_args$run_4d_files <- mrfiles
        bdm_args$runs_to_output <- mr_run_nums
        bdm_args$output_directory <- file.path(subj_out, "timing_files")
        d_obj <- tryCatch(do.call(build_design_matrix, bdm_args), error=function(e) {
          lg$error("Failed build_design_matrix for subject: %s, model: %s", subid, this_model)
          lg$error("Error message: %s", as.character(e))
          return(NULL)
        })

        if (is.null(d_obj)) { next } #skip to next iteration on error

        save(d_obj, bdm_args, mrdf, mr_run_nums, subj_mr_dir, mrfiles, run_lengths, subid, this_model,
          file=file.path(subj_out, paste0(gpa$l1_models$models[[this_model]]$name, "_bdm_setup.RData")))

        if ("fsl" %in% gpa$glm_software) {
          #Setup FSL run-level models for each combination of signals
          #Returns a data.frame of feat l1 inputs and the fsf file
          feat_files <- tryCatch(fsl_l1_model(id=subid, d_obj, gpa, this_model, nvoxels = nvoxels, lg=lg),
            error=function(e) {
              lg$error("Problem running fsl_l1_model. Model: %s, Subject: %s", this_model, subid)
              lg$error("Error message: %s", as.character(e))
              return(NULL)
            })

          if (!is.null(feat_files)) {
            #add to tracking data.frame (simple append)
            l1_file_setup$fsl <- c(l1_file_setup$fsl, list(feat_files))
          }
        }

        #TODO: finalize SPM approach
        if ("spm" %in% gpa$glm_software) {
          #Setup spm run-level models for each combination of signals
          spm_files <- tryCatch(spm_l1_model(d_obj, gpa, this_model, mr_files),
            error=function(e) {
              lg$error("Problem running spm_l1_model. Model: %s, Subject: %s", this_model, subid)
              lg$error("Error message: %s", as.character(e))
              return(NULL)
            })

        }
      }

      browser()
      lg$info("Completed processing of subject: %s", subid)
      return(l1_file_setup)
    }

  return(ll)
}
