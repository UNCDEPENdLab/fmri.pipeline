#cutting down on redundancy across glm setup scripts

#wrapper for running an AFNI command safely within R
#if AFNI does not have its environment setup properly, commands may not work
runAFNICommand <- function(args, afnidir=NULL, stdout=NULL, stderr=NULL, ...) {
  #look for AFNIDIR in system environment if not passed in
  if (is.null(afnidir)) {
    env <- system("env", intern=TRUE)
    if (length(afnidir <- grep("^AFNIDIR=", env, value=TRUE)) > 0L) {
      afnidir <- sub("^AFNIDIR=", "", afnidir)
    } else {
      warning("AFNIDIR not found in environment. Defaulting to ", paste0(normalizePath("~/"), "/afni"))
      afnidir <- paste0(normalizePath("~/"), "/afni")
    }
  }
  
  Sys.setenv(AFNIDIR=afnidir) #export to R environment
  afnisetup <- paste0("AFNIDIR=", afnidir, "; PATH=${AFNIDIR}:${PATH}; DYLD_FALLBACK_LIBRARY_PATH=${AFNIDIR}; ${AFNIDIR}/")
  afnicmd  <- paste0(afnisetup, args)
  if (!is.null(stdout)) { afnicmd <- paste(afnicmd, ">", stdout) }
  if (!is.null(stderr)) { afnicmd <- paste(afnicmd, "2>", stderr) }
  cat("AFNI command: ", afnicmd, "\n")
  retcode <- system(afnicmd, ...)
  return(retcode)
}


#wrapper for running an fsl command safely within R
#if FSL does not have its configuration setup properly, commands such as feat don't work, or hang strangely
runFSLCommand <- function(args, fsldir=NULL, stdout=NULL, stderr=NULL) {
  #look for FSLDIR in system environment if not passed in
  if (is.null(fsldir)) {
    #check for FSLDIR in sourced .bashrc
    bashrc_fsldir <- character(0)
    if (file.exists("~/.profile")) {
      bashrc_fsldir <- system("source ~/.profile && echo $FSLDIR", intern=TRUE)
    }
    
    #check for FSLDIR in current environment
    env <- system("env", intern=TRUE)
    if (length(fsldir <- grep("^FSLDIR=", env, value=TRUE)) > 0L) {
      fsldir <- sub("^FSLDIR=", "", fsldir)
    } else if (!identical(bashrc_fsldir, character(0))) {
      fsldir <- bashrc_fsldir      
    } else {
      warning("FSLDIR not found in environment. Defaulting to /usr/local/fsl.")
      fsldir <- "/usr/local/fsl"
    }
  }
  
  #Sys.setenv(LD_LIBRARY_PATH="/gpfs/group/mnh5174/default/sw/openblas/lib")
  Sys.setenv(FSLDIR=fsldir) #export to R environment
  fslsetup=paste0("FSLDIR=", fsldir, "; PATH=${FSLDIR}/bin:${PATH}; . ${FSLDIR}/etc/fslconf/fsl.sh; ${FSLDIR}/bin/")
  fslcmd=paste0(fslsetup, args)
  if (!is.null(stdout)) { fslcmd=paste(fslcmd, ">", stdout) }
  if (!is.null(stderr)) { fslcmd=paste(fslcmd, "2>", stderr) }
  cat("FSL command: ", fslcmd, "\n")
  retcode <- system(fslcmd)
  return(retcode)
}


#' internal function for getting a file extension that may include a compressed ending
#' this is adapted from neurobase::file_imgext, but extended for other file types
#' @keywords internal
file_ext <- function(file, withdot = TRUE) {
  file <- tolower(file)
  matches <- grepl("^.*\\.(csv|dat|hdr|img|brik|head|nii|txt|tsv)(\\.gz|\\.bz2|\\.zip|\\.xz)*$", file)
  ext <- rep(NA, length=length(file)) # return NA for inputs that can't be parsed
  ext[matches] <- sub("^(.*)\\.(csv|dat|hdr|img|brik|head|nii|txt|tsv)(\\.gz|\\.bz2|\\.zip|\\.xz)*$", "\\2\\3", file[matches])
  if (isTRUE(withdot)) ext[matches] <- paste0(".", ext[matches])
  return(ext)
}

#' helper function for generating motion regressors from raw 6-parameter motion coregistration
#'
#' @param motion_params_file file containing 6 motion parameters
#' @param col.names names of columns in \code{motion_params_file}, in order from left to right
#' @param drop_volumes number of volumes to drop from beginning of motion params
#' @param last_volume final volume to include from motion params (if truncated at end). If \code{NULL},
#'   then the end of the time series is not truncated.
#' @param rot_units The units of the rotation parameters. Default is "rad" for radians
#' @param tra_units The units of the translation parameter. Only support "mm" right now millimeters
#'
#' @importFrom data.table fread
#' @importFrom checkmate assert_file_exists assert_subset assert_integerish
#' @keywords internal
generate_motion_regressors <- function(motion_params_file = "motion.par",
                                       col.names = c("rx", "ry", "rz", "tx", "ty", "tz"),
                                       demean = FALSE, drop_volumes = 0L, last_volume=NULL,
                                       rot_units="rad", tra_units="mm",
                                       na.strings="NA", lg=NULL) {

  checkmate::assert_class(lg, "Logger")
  checkmate::assert_file_exists(motion_params_file)
  if (is.null(col.names)) {
    # explicit defaults in case of null
    lg$debug("In generate_motion_regressors, defaulting column order rx, ry, rz, tx, ty, tz")
    col.names <- c("rx", "ry", "rz", "tx", "ty", "tz")
  }
  checkmate::assert_subset(col.names, c("rx", "ry", "rz", "tx", "ty", "tz"))
  checkmate::assert_logical(demean, len = 1L)
  checkmate::assert_integerish(drop_volumes, lower = 0)
  checkmate::assert_subset(rot_units, c("rad", "deg"))
  checkmate::assert_subset(tra_units, c("mm"))
  checkmate::assert_character(na.strings)

  # read raw motion parameters file
  mot <- data.table::fread(motion_params_file, na.strings = na.strings)

  if (ncol(mot) != 6L) {
    lg$error("In generate_motion_regressors, motion file %s has %d columns instead of 6!", motion_params_file, ncol(mot))
    return(invisible(NULL))
  } else {
    data.table::setnames(mot, col.names)
  }

  if (is.null(last_volume)) { last_volume <- nrow(mot) }
  checkmate::assert_integerish(last_volume, upper=nrow(mot))

  # subset motion regressors to match drops and truncations in fmri data
  mot <- mot[(1 + drop_volumes):last_volume, ]

  mot_deriv <- mot[, lapply(.SD, function(x) { c(0, diff(x)) })]
  data.table::setnames(mot_deriv, paste0("d", names(mot))) #add delta to names
  mot <- cbind(mot, mot_deriv)

  mot_quad <- mot[, lapply(.SD, function(x) { x^2 })]
  data.table::setnames(mot_quad, paste0("q", names(mot))) #add delta to names
  mot <- cbind(mot, mot_quad)

  # https://wiki.cam.ac.uk/bmuwiki/FMRI
  if (rot_units=="rad") {
    FD <- apply(mot_deriv[, c("drx", "dry", "drz")], 1, function(x) 50 * sum(abs(x))) +
      apply(mot_deriv[, c("dtx", "dty", "dtz")], 1, function(x) sum(abs(x)))
  } else if (rot_units=="deg") {
    FD <- apply(mot_deriv[, c("drx", "dry", "drz")], 1, function(x) 50 * (pi / 180) * sum(abs(x))) +
      apply(mot_deriv[, c("dtx", "dty", "dtz")], 1, function(x) sum(abs(x)))
  } else {
    stop("unknown")
  }
  mot <- cbind(mot, FD = FD)

  # Demean all columns, if requested
  # Note that this is not used internally at present to avoid surprising results when testing motion thresholds (e.g., FD)
  if (isTRUE(demean)) {
    mot <- mot[, lapply(.SD, function(x) { x - mean(x, na.rm=TRUE) }) ]
  }

  ##just PCA motion on the current run
  ##mregressors <- pca_motion(run_nifti[r], runlengths[r], motion_parfile="motion.par", numpcs=3, drop_volumes=drop_volumes)$motion_pcs_concat

  return(mot)
}

#' compute spike regressors from a data.frame or matrix of motion parameters and a set of
#'   expressions that are evaluated against this parameter matrix.
#' 
#' @param mot a data.frame or matrix with volumes on rows and named motion parameters on columns
#' @param spike_volume a character vector of expressions to evaluate against \code{mot}. Resulting
#'   columns will be prefixed with the names of each expression. If the expressions are unnamed,
#'   prefixes will be expr1_, expr2_, etc.
#'
#' @return a matrix of spike regressors (volumes on rows, spike regressors on columns)
#' @importFrom dplyr lead lag
#' @keywords internal
compute_spike_regressors <- function(mot = NULL, spike_volume = NULL, lg=NULL) {
  if (is.null(mot)) return(NULL)
  if (is.null(spike_volume)) return(NULL)

  checkmate::assert_character(spike_volume, null.ok = TRUE)
  checkmate::assert_class(lg, "Logger")

  spikes <- do.call(cbind, lapply(seq_along(spike_volume), function(ii) {
    has_bounds <- grepl(";", spike_volume[ii], fixed = TRUE)
    if (isTRUE(has_bounds)) {
      esplit <- strsplit(spike_volume[ii], "\\s*;\\s*", perl = TRUE)[[1L]]
      stopifnot(length(esplit) == 2L)
      spike_bounds <- as.integer(eval(parse(text = paste0("c(", esplit[1], ")"))))
      expr <- esplit[2L]
    } else {
      spike_bounds <- 0L # only volume where expr is TRUE
      expr <- spike_volume[ii]
    }

    spike_vec <- tryCatch(with(mot, eval(parse(text = expr))),
      error = function(e) {
        lg$error(
          "Problem evaluating spike regressors for file %s, expr: %s",
          motion_params_file, expr
        )
        return(NULL)
      }
    )

    which_spike <- which(spike_vec)
    if (length(which_spike) == 0L) return(NULL) #don't attempt to generate spikes if no TRUE values

    # will evaluate to NULL if which_spike is empty
    spike_df <- do.call(cbind, lapply(which_spike, function(xx) {
      vec <- rep(0, nrow(mot))
      vec[xx] <- 1
      return(vec)
    }))

    colnames(spike_df) <- paste0("spike_", seq_len(ncol(spike_df)))

    if (!is.null(spike_df) && !identical(spike_bounds, 0L)) {
      # apply shifts
      shifts <- spike_bounds[spike_bounds != 0L]
      res <- do.call(cbind, lapply(shifts, function(ss) {
        shift_mat <- apply(spike_df, 2, function(col) {
          if (ss < 0) {
            dplyr::lead(col, abs(ss), default = 0)
          } else {
            dplyr::lag(col, ss, default = 0)
          }
        })
        colnames(shift_mat) <- paste0("spike_", ifelse(ss < 0, "m", "p"), abs(ss), "_", 1:ncol(shift_mat))
        return(shift_mat)
      }))

      spike_df <- cbind(spike_df, res)
    }

    if (is.null(names(spike_volume)[ii]) || names(spike_volume)[ii] == "") {
      colnames(spike_df) <- paste0("expr", ii, "_", colnames(spike_df))
    } else {
      colnames(spike_df) <- paste0(names(spike_volume)[ii], "_", colnames(spike_df))
    }
    return(spike_df)

  }))

  # we don't want to generate any duplicate columns (linear dependency problems) 
  # due to overlapping expressions or ranges of evaluation
  spikes <- spikes[, !duplicated(spikes, MARGIN=2)]

  return(spikes)
}


#' small helper function for compressing motion parameters using PCA
#'
#' @param motion_df volumes x motion parameters data frame for PCA compression
#' @param num_pcs number of principal components to extract
#' @param zscore whether to standardize motion parameters prior to PCA (this is a good idea)
#' @param verbose whether to print the variance explained by PCs and the multiple correlations of
#'   among motion parameters
#'
#' @return A matrix of PCA-compressed motion regressors
#' @importFrom psych smc
#' @keywords internal
pca_motion <- function(motion_df, num_pcs=3L, zscore=TRUE, verbose=FALSE) {
  checkmate::assert_data_frame(motion_df)
  checkmate::assert_integerish(num_pcs, lower=1, upper=50)
  checkmate::assert_logical(zscore)
  checkmate::assert_logical(verbose)

  #compute the PCA decomposition of motion parameters and their derivatives
  pc <- prcomp(motion_df, retx=TRUE, scale.=zscore)
  cumvar <- cumsum(pc$sdev^2/sum(pc$sdev^2))

  if (isTRUE(verbose)) message("first", num_pcs, "motion principal components account for: ", round(cumvar[num_pcs], 3))
  mregressors <- pc$x[, 1:num_pcs] #cf Churchill et al. 2012 PLoS ONE
  attr(mregressors, "variance.explained") <- cumvar[num_pcs]

  if (isTRUE(verbose)) {
    cat("Multiple correlation of motion parameters:\n\n")
    print(round(sqrt(psych::smc(motion_df)), 2))
  }

  return(as.data.frame(mregressors))
}

generateRunMask <- function(mr_files, outdir=getwd(), outfile="runmask") {
  if (file.exists(file.path(outdir, paste0(outfile, ".nii.gz")))) { return(invisible(NULL)) }
  ##generate mask of mr_files where temporal min is > 0 for all runs
  for (f in seq_along(mr_files)) {
    runFSLCommand(paste0("fslmaths ", mr_files[f], " -Tmin -bin ", outdir, "/tmin", f))#, fsldir="/usr/local/ni_tools/fsl")
  }

  ##sum mins together over runs and threshold at number of runs
  runFSLCommand(paste0("fslmaths ", paste(paste0(outdir, "/tmin", seq_along(mr_files)), collapse=" -add "), " ", outdir, "/tminsum"))#, fsldir="/usr/local/ni_tools/fsl")
  runFSLCommand(paste0("fslmaths ", outdir, "/tminsum -thr ", length(mr_files), " -bin ", outdir, "/", outfile))#, fsldir="/usr/local/ni_tools/fsl")
  runFSLCommand(paste0("imrm ", outdir, "/tmin*"))#, fsldir="/usr/local/ni_tools/fsl") #cleanup 
}

visualizeDesignMatrix <- function(d, outfile=NULL, runboundaries=NULL, events=NULL, includeBaseline=TRUE) {
  require(ggplot2)
  require(reshape2)

  if (!includeBaseline) {
    d <- d[,!grepl("run[0-9]+base", colnames(d))]
  }

  print(round(cor(d), 3))
  d <- as.data.frame(d)
  d$volume <- seq_len(nrow(d))
  d.m <- melt(d, id.vars="volume")
  g <- ggplot(d.m, aes(x = volume, y = value)) +
    geom_line(size = 1.2) +
    theme_bw(base_size = 15) +
    facet_grid(variable ~ ., scales = "free_y")
  
  colors <- c("black", "blue", "red", "orange") #just a hack for color scheme right now

  if (!is.null(runboundaries)) {
    g <- g + geom_vline(xintercept=runboundaries, color=colors[1L])
  }

  if (!is.null(events)) {
    for (i in seq_along(events)) {
      g <- g + geom_vline(xintercept=events[[i]], color=colors[i+1])
    }
  }

  if (!is.null(outfile)) {
    ggsave(filename=outfile, plot=g, width=21, height=9)
  }
  return(invisible(g))
}

#compute the mean of each cluster in roimask
#within a 4d array: x, y, z, subbrik
#
#roimask is an x,y,z array with integer mask values
#ni4d is a concatenated set of cope values, either across subjects
# or, in the case of beta series, within a subject (over trials)
get_cluster_means <- function(roimask, ni4d) {
  #ni4d should be a 4s array in which the fourth dimension is either subject (single beta) or a beta series (single subject)
  maskvals <- sort(unique(as.vector(roimask)))
  maskvals <- maskvals[!maskvals == 0]

  roimats <- sapply(maskvals, function(v) {
    mi <- which(roimask==v, arr.ind=TRUE)
    nsubbriks <- dim(ni4d)[4]
    nvox <- nrow(mi)
    mi4d <- cbind(pracma::repmat(mi, nsubbriks, 1), rep(1:nsubbriks, each=nvox))

    mat <- matrix(ni4d[mi4d], nrow=nvox, ncol=nsubbriks) #need to manually reshape into matrix from vector
    #for each subject, compute huber m-estimator of location/center Winsorizing at 2SD across voxels (similar to voxel mean)
    #clusavg <- apply(mat, 2, function(x) { MASS::huber(x, k=2)$mu })
    clusavg <- apply(mat, 2, mean)

    return(clusavg)
  })
}

#function to extract mean beta series
get_beta_series <- function(inputs, roimask, n_bs=50) {
  #inputs <- inputs[1:5] #speed up testing

  beta_res <- foreach(i=iter(1:length(inputs)), .packages=c("reshape2", "oro.nifti", "dplyr", "abind"), .export="get_cluster_means") %dopar% {
  #beta_res <- lapply(1:length(inputs), function(i) {
    run_dirs <- list.files(path=inputs[i], pattern="FEAT_LVL1_run\\d+\\.feat", recursive=FALSE, full.names=TRUE)
    run_betas <- list()
    
    for (r in seq_along(run_dirs)) {
      runnum <- as.numeric(sub(".*/FEAT_LVL1_run(\\d+)\\.feat", "\\1", run_dirs[r], perl=TRUE))
      copes <- file.path(run_dirs[r], "stats", paste0("cope", 1:n_bs, ".nii.gz"))

      if (!all(file.exists(copes))) {
        cat("File listing in dir:", run_dirs[r], list.files(run_dirs[r], recursive=FALSE), sep="\n", append=TRUE)
        cat("Cannot find all ", length(copes), " beta series copes in: ",
          run_dirs[r], "\n",
          file = "beta_series_errors.txt", append = TRUE
        )
        next #don't error for now, just omit
      }

      #in testing, it is faster to use readNIfTI to read each cope file individually (7.5s) than to use fslmerge + readNIfTI on the 4d (13s)
      #concat_file <- tempfile()
      #system.time(runFSLCommand(paste("fslmerge -t", concat_file, paste(copes, collapse=" "))))
      #system.time(cout <- readNIfTI(concat_file, reorient=FALSE)@.Data)

      cout <- do.call(abind, list(along=4, lapply(copes, function(x) { readNIfTI(x, reorient=FALSE)@.Data })))
      beta_series_cluster_means <- get_cluster_means(roimask, cout)
      beta_melt <- reshape2::melt(beta_series_cluster_means, value.name="bs_value", varnames=c("trial", "cluster_number")) %>%
        mutate(feat_input_id=i, run=runnum)

      run_betas[[r]] <- beta_melt
    }

    return(do.call(rbind, run_betas))

  }

  return(do.call(rbind, beta_res))
}

#' Internal helper function for printing a summary of a contrast matrix
#' @param cmat a contrast matrix where row names are the contrast names and column
#'   names are the column names of the corresponding design matrix.
#' @keywords internal
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

#' helper function to rename columns of input data.frame to internal nomenclature
#'   based on the variable mapping (vm) vector
#' 
#' @param df a data.frame containing columns to be renamed to internal standards
#' @param vm a named vector of columns in \code{df} that identify internal constructs
#'   such as id, session, and run.
#'
#' @return a modified version of \code{df} with column names modified to use internal names
#' @importFrom data.table setnames
#' @keywords internal
names_to_internal <- function(df, vm) {
  # look for naming collisions
  nm <- names(df)
  nv <- names(vm)
  conflict <- intersect(nm[nm %in% names(vm)], nv[nv != vm])

  if (length(conflict) > 0L) {
    new_names <- paste0("..", conflict, "..")
    cat(
      "To avoid naming conflicts, renaming these columns:", paste(conflict, collapse = ", "),
      "to:", paste0(new_names, collapse = ", "), "\n"
    )

    nm[which(nm %in% conflict)] <- new_names
  }

  # convert variable names to internal constructs via the variable mapping
  vm_present <- vm[which(vm %in% nm)]
  vm_pos <- sapply(vm_present, function(x) { which(nm == x)[1] })
  nm[vm_pos] <- names(vm_present)

  return(nm)

}

#' small helper function to pull absolute paths to a given column in run_data or trial_data
#'
#' @param mr_df a data.frame from a gpa object, following standard variable mapping nomenclature
#' @param col a character string denoting the column in \code{mr_df} to be used for looking up
#'   absolute paths
#'
#' @details Note that if a given value of the requested column is an absolute path, it will
#'   not be combined with $mr_dir to generate the combined path. Thus, the \code{col} in
#'   \code{mr_df} can contain a mixture of relative an absolute paths. The relative paths will
#'   be combined with 
#' @return Absolute paths to all files in the specified \code{col}
#' @importFrom R.utils getAbsolutePath
#' @keywords internal
get_mr_abspath <- function(mr_df, col="run_nifti") {
  checkmate::assert_data_frame(mr_df)
  checkmate::assert_string(col)
  checkmate::assert_subset(col, names(mr_df))

  sapply(seq_len(nrow(mr_df)), function(ii) {
    if (!"mr_dir" %in% names(mr_df)) {
      # this should be rare (or ideally, not happen at all)
      # but if there is no mr_dir in the mr_df, we can only use the column itself
      wd <- NULL # defaults to getwd()
    } else {
      this_dir <- mr_df$mr_dir[ii]
      if (is.na(this_dir)) {
        wd <- getwd() # defaults to getwd()
      } else {
        wd <- this_dir
      }
    }

    R.utils::getAbsolutePath(mr_df[[col]][ii], workDirectory = wd, expandTilde = TRUE)
  })

}



#' small helper function to populate subject exclusions in $run_data and $subject_data
#'
#' @param gpa a \code{glm_pipeline_arguments} object that already has the \code{$exclude_run}
#'   column populated in \code{$l1_model_setup}.
#' @return a modified copy of \code{gpa} where $exclude_subject has been added to $run_data
#'   and $subject_data. This function also adds $n_good_runs to $subject_data which is useful
#'   if we want to enforce a lower bound on the number of runs used to define subject exclusion.
#' @keywords internal
#' @importFrom dplyr filter count left_join
#' @importFrom checkmate assert_class assert_data_table
#' @importFrom lgr get_logger
#' @author Michael Hallquist
calculate_subject_exclusions <- function(gpa) {
  lg <- lgr::get_logger("glm_pipeline/setup_glm_pipeline")
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_data_frame(gpa$run_data)
  checkmate::assert_subset(c("id", "session", "run_number", "exclude_run"), names(gpa$run_data)) # verify that exclude_run is populated

  if ("exclude_subject" %in% names(gpa$subject_data) && "exclude_subject" %in% names(gpa$run_data)) {
    lg$debug("exclude_subject already calculated and populated in $subject_data and $run_data")
    return(gpa)
  }

  # if no subject exclusion string is provided, then keep all subjects
  if (is.null(gpa$confound_settings$exclude_subject)) {
    gpa$run_data$exclude_subject <- FALSE
    gpa$subject_data$exclude_subject <- FALSE
    return(gpa)
  }

  n_good_runs_df <- gpa$run_data %>%
    dplyr::group_by(id, session) %>%
    dplyr::summarize(n_good_runs = sum(exclude_run == FALSE)) %>% # sum of non-excluded runs
    dplyr::ungroup() %>%
    dplyr::select(id, session, n_good_runs)

  gpa$subject_data <- gpa$subject_data %>%
    dplyr::left_join(n_good_runs_df, by = c("id", "session"))

  gpa$run_data <- gpa$run_data %>%
    dplyr::left_join(n_good_runs_df, by = c("id", "session"))

  # evaluate subject exclusion in the context of the broader $subject_data data.frame
  gpa$subject_data$exclude_subject <- sapply(seq_len(nrow(gpa$subject_data)), function(ss) {
    s_info <- gpa$subject_data[ss, , drop = FALSE]
    tryCatch(with(s_info, eval(parse(text = gpa$confound_settings$exclude_subject))),
      error = function(e) {
        lg$error(
          "Problem evaluating subject exclusion for subject: %s, session: %s, expr: %s",
          s_info$id, s_info$session,
          gpa$confound_settings$exclude_subject
        )
        lg$error("Defaulting to retaining this subject.")
        return(FALSE)
      }
    )
  })

  # propagate subject exclusion down to $run_data
  gpa$run_data <- gpa$run_data %>%
    dplyr::left_join(
      gpa$subject_data %>% dplyr::select(id, session, exclude_subject),
      by = c("id", "session")
    )

  return(gpa)
}

#' Internal function to add last_onset and last_offset columns to gpa$run_data based on events in l1_models
#' @param gpa a \code{glm_pipeline_arguments} object containing valid $run_data and $l1_models objects
#' @details The last_onset and last_offset columns are calculated for each run based on the timing of all
#'   events in the $l1_models$events list. These are then used to facilitate run truncation if the user
#'   requests truncation after a final onset or offset.
#' @return a modified copy of gpa with $run_data populated with last_offset and last_onset columns (times in seconds)
#' @keywords internal
#' @importFrom dplyr summarize group_by left_join
populate_last_events <- function(gpa) {
  # get all events as a long data.frame
  m_events <- data.table::rbindlist(lapply(gpa$l1_models$events, function(this_event) this_event$data))

  last_events <- m_events %>%
    group_by(id, session, run_number) %>%
    dplyr::summarize(last_onset = max(onset, na.rm = TRUE), last_offset = max(onset + duration, na.rm = TRUE), .groups="drop")

  gpa$run_data <- gpa$run_data %>%
    dplyr::left_join(last_events, by=c("id", "session", "run_number"))

  return(gpa)
}

# slurm_job_array <- function(job_name = "slurm_array") {
# }

#' helper function to convert the $sched_args field to a vector
#'  of directives that can be included dynamically in the header of
#'  PBS or SBATCH scripts.
#'
#' @param gpa a \code{glm_pipeline_arguments} object containing the
#'   $parallel$sched_args field
#'
#' @return a character vector where the $sched_args elements are converted to
#'   corresponding # SBATCH or # PBS directives for inclusion in scripts
#' @keywords internal
sched_args_to_header <- function(gpa) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")

  #if there are no scheduler arguments, return nothing
  if (is.null(gpa$parallel$sched_args)) return(NULL)

  is_slurm <- gpa$scheduler == "slurm"
  directives <- sapply(gpa$parallel$sched_args, function(x) {
    ifelse(isTRUE(is_slurm), paste("#SBATCH", x), paste("#PBS", x))
  }, USE.NAMES = FALSE)
  return(directives)
}

#' helper function to generate a contrast matrix from an lm() object
#'   and a set of user-specified contrasts using emmeans
#'
#' @param mobj a model object created by build_l<X>_models
#' @param lmfit an optional lm() object used for emmeans calculations. If provided, this object
#'   will be used instead of mobj$lmfit (the parent lm on the overall dataset).
#'
#' @return a modified copy of the model object \code{mobj} with $contrast_list and $contrasts
#'   fully populated 
#' @keywords internal
#' @importFrom emmeans emmeans emtrends
get_contrasts_from_spec <- function(mobj, lmfit=NULL) {
  checkmate::assert_multi_class(mobj, c("l1_model_spec", "l1_wi_spec", "hi_model_spec")) # verify that we have an object of known structure
  if (is.null(lmfit)) {
    lmfit <- mobj$lmfit # use parent lmfit
  } else {
    # override existing contrasts to force lmfit-specific respecification
    mobj$contrast_list <- list()
    mobj$contrasts <- NULL
  }
  spec <- mobj$contrast_spec
  contrast_list <- mobj$contrast_list

  checkmate::assert_class(lmfit, "lm", null.ok=TRUE)
  checkmate::assert_list(spec)
  checkmate::assert_list(contrast_list)

  if (is.null(lmfit)) {
    c_colnames <- spec$regressors
  } else {
    c_colnames <- names(coef(lmfit)) # prefer model-specific regressors (if they vary by subject)
  }

  ### add diagonal contrasts
  c_diagonal <- contrast_list$diagonal
  if (isTRUE(spec$diagonal) && is.null(c_diagonal)) {

    diag_mat <- diag(length(c_colnames))
    rownames(diag_mat) <- paste0("EV_", c_colnames) # simple contrast naming for each individual regressor
    colnames(diag_mat) <- c_colnames # always have columns named by regressor

    c_diagonal <- diag_mat
  }

  ### add condition means, if requested
  c_cond_means <- contrast_list$cond_means
  if (length(spec$cond_means) > 0L && is.null(c_cond_means)) {
    for (vv in spec$cond_means) {
      ee <- emmeans(lmfit, as.formula(paste("~", vv)), weights = spec$weights)
      edata <- summary(ee)
      econ <- ee@linfct
      enames <- as.character(edata[[vv]]) # names of contrasts

      # if any emmeans are not estimable, then this is likely due to aliasing. For now, drop
      which_na <- is.na(edata$emmean)
      if (any(which_na)) {
        econ <- econ[!which_na, , drop = FALSE]
        enames <- enames[!which_na]
      }

      # add contrast names to matrix
      rownames(econ) <- paste(vv, enames, sep=".")

      # add contrasts to matrix
      c_cond_means <- rbind(c_cond_means, econ)
    }
    colnames(c_cond_means) <- c_colnames
  }

  ### add cell means for all factors, if requested
  c_cell_means <- contrast_list$cell_means
  if (isTRUE(spec$cell_means) && is.null(c_cell_means)) {
    # get model-predicted means for each factor
    ee <- emmeans(lmfit, as.formula(paste("~", paste(spec$cat_vars, collapse = "*"))), weights = spec$weights)
    # pp <- pairs(ee)
    edata <- summary(ee)

    econ <- ee@linfct
    enames <- apply(edata[, spec$cat_vars, drop = FALSE], 1, function(x) { paste(names(x), x, sep=".", collapse = "_") })

    # if any emmeans are not estimable, then this is likely due to aliasing. For now, drop
    which_na <- is.na(edata$emmean)
    if (any(which_na)) {
      econ <- econ[!which_na, , drop = FALSE]
      enames <- enames[!which_na]
    }

    rownames(econ) <- enames
    c_cell_means <- rbind(c_cell_means, econ)
    colnames(c_cell_means) <- c_colnames
  }

  ### add overall response mean
  c_overall <- contrast_list$overall
  if (isTRUE(spec$overall_response) && is.null(c_overall)) {
    ee <- emmeans(lmfit, ~1, weights = spec$weights)
    econ <- ee@linfct
    rownames(econ) <- "overall"
    colnames(econ) <- c_colnames
    c_overall <- econ
  }

  ### add pairwise differences, if requested
  c_pairwise_diffs <- contrast_list$pairwise_diffs
  if (length(spec$pairwise_diffs) > 0L && is.null(c_pairwise_diffs)) {
    for (vv in spec$pairwise_diffs) {
      ee <- emmeans(lmfit, as.formula(paste("~", vv)), weights = spec$weights)
      pp <- pairs(ee)
      edata <- summary(pp)
      econ <- pp@linfct
      enames <- paste(vv, make.names(sub("\\s+-\\s+", "_M_", edata$contrast, perl = TRUE)), sep=".") # names of contrasts

      # if any emmeans are not estimable, then this is likely due to aliasing. For now, drop
      which_na <- is.na(edata$estimate)
      if (any(which_na)) {
        econ <- econ[!which_na, , drop = FALSE]
        enames <- enames[!which_na]
      }

      # add contrast names to matrix
      rownames(econ) <- enames

      # add pairwise differences for this factor to the pairwise matrix
      c_pairwise_diffs <- rbind(c_pairwise_diffs, econ)
    }
    colnames(c_pairwise_diffs) <- c_colnames
  }

  ### add per-cell model-predicted simple slopes
  c_simple_slopes <- contrast_list$simple_slopes
  if (length(spec$simple_slopes) > 0L && is.null(c_simple_slopes)) {
    for (vv in seq_along(spec$simple_slopes)) {
      trend_var <- names(spec$simple_slopes)[vv]
      for (comb in spec$simple_slopes[[vv]]) {
        ee <- emtrends(lmfit,
          specs = as.formula(paste("~", paste(comb, collapse = "*"))),
          var = trend_var, weights = spec$weights
        )
        edata <- summary(ee)
        econ <- ee@linfct
        enames <- paste(trend_var, apply(edata[, comb, drop = FALSE], 1, paste, collapse = "."), sep = ".")

        # if any emmeans are not estimable, then this is likely due to aliasing. For now, drop
        which_na <- is.na(edata[[paste(trend_var, "trend", sep=".")]])
        if (any(which_na)) {
          econ <- econ[!which_na, , drop = FALSE]
          enames <- enames[!which_na]
        }

        rownames(econ) <- enames
        c_simple_slopes <- rbind(c_simple_slopes, econ)
      }
    }
    colnames(c_simple_slopes) <- c_colnames
  }

  c_custom <- contrast_list$custom

  # walk through all within-subject factor contrasts and return as a single matrix
  c_within <- get_wi_contrast_matrix(mobj, c_colnames)

  #combine each element
  cmat_full <- rbind(
    c_diagonal,
    c_cond_means,
    c_cell_means,
    c_overall,
    c_pairwise_diffs,
    c_simple_slopes,
    c_within,
    c_custom
  )

  if (!is.null(spec$delete)) {
    which_del <- match(spec$delete, rownames(cmat_full))
    if (!is.na(which_del[1L])) {
      cmat <- cmat_full[-1*which_del, , drop=FALSE]
    }
  } else {
    cmat <- cmat_full
  }

  dupes <- duplicated(cmat, MARGIN = 1)
  cmat <- cmat[!dupes, , drop=FALSE] # drop duplicated contrasts

  # parentheses in contrast names generate problems in running through FEAT (chokes on submission)
  # so far, this comes up only with intercept terms -- may need to do a more general substitution later
  if (!is.null(cmat)) { rownames(cmat) <- gsub("(Intercept)", "Intercept", rownames(cmat), fixed = T) }

  # populate any updates to the contrast_list object based on new calculations
  contrast_list$diagonal <- c_diagonal
  contrast_list$cond_means <- c_cond_means
  contrast_list$cell_means <- c_cell_means
  contrast_list$overall <- c_overall
  contrast_list$pairwise_diffs <- c_pairwise_diffs
  contrast_list$simple_slopes <- c_simple_slopes
  contrast_list$custom <- c_custom

  # populate contrasts back into model object (spec should not be updated)
  mobj$contrast_list <- contrast_list
  mobj$contrasts <- cmat
  return(mobj)
}

get_wi_contrast_matrix <- function(mobj, c_colnames) {
  # this function is only intended to crawl over within-subject contrasts of an l1 model
  if (!checkmate::test_class(mobj, "l1_model_spec") || is.null(mobj$wi_models)) {
    return(invisible(NULL)) # quietly return NULL if we don't have an l1 model
  }

  # get contrast matrix for each wi_model
  clist <- lapply(mobj$wi_models, "[[", "contrasts")

  cnames <- do.call(c, lapply(clist, rownames))
  cmat <- matrix(0, nrow = length(cnames), ncol = length(c_colnames), dimnames = list(cnames, c_colnames))
  
  # fill in relevant rows and columns in combined matrix for each wi model contrast matrix
  for (cc in clist) {
    cmat[rownames(cc), colnames(cc)] <- cc
  }

  return(cmat)
}

# get_run_lmfits <- function() {

# }

#' Helper function to run the requested GLM model for each subject+session separately
#'
#' @param mobj an \code{l1_model_spec} or \code{hi_model_spec} object containing the GLM model to run
#' @param data The run-level data frame containing data for all ids and sessions. This will be split into individual chunks
#'
#' @return a modified copy of \code{mobj} where the $by_subject field has been added
#'
#' @details The function adds the $by_subject field, which contains the design matrices and contrasts
#'   for each subject and session in \code{data} based on the available data for that session. For example, if
#'   a subject is missing a few runs (or these are dropped from analysis), then some contrasts may change or drop out of the model.
#'
#' The $by_subject field is a keyed data.table object containing list elements for the cope_list (mapping cope numbers to contrast names),
#'   the contrasts, and the design matrix for each session.
#'
#' @keywords internal
#' @importFrom checkmate assert_data_frame assert_multi_class assert_subset
#' @importFrom data.table data.table
respecify_l2_models_by_subject <- function(mobj, data) {
  checkmate::assert_multi_class(mobj, c("l1_model_spec", "hi_model_spec")) # verify that we have an object of known structure
  checkmate::assert_data_frame(data)
  checkmate::assert_subset(c("id", "session"), names(data))

  # create a nested data.table object for fitting each id + session separately
  data <- data.table(data, key = c("id", "session"))
  data$dummy <- rnorm(nrow(data))
  dsplit <- data[, .(dt = list(.SD)), by = c("id", "session")]

  # use model formula from parent object
  model_formula <- terms(mobj$lmfit)

  cope_list <- list()
  contrast_list <- list()
  model_matrix_list <- list()
  n_l2_copes <- rep(NA_integer_, nrow(dsplit))
  for (vv in seq_len(nrow(dsplit))) {
    lmfit <- lm(model_formula, data = dsplit[[vv, "dt"]])

    mm <- get_contrasts_from_spec(mobj, lmfit)
    cope_df <- data.frame(
      id = dsplit$id[vv], session = dsplit$session[vv],
      l2_cope_number = seq_len(nrow(mm$contrasts)), l2_cope_name = rownames(mm$contrasts)
    )

    cope_list[[vv]] <- cope_df
    contrast_list[[vv]] <- mm$contrasts
    model_matrix_list[[vv]] <- model.matrix(lmfit)
    n_l2_copes[vv] <- ncol(mm$contrasts)
  }

  dsplit[, cope_list := cope_list]
  dsplit[, contrasts := contrast_list]
  dsplit[, model_matrix := model_matrix_list]
  dsplit[, n_l2_copes := n_l2_copes]
  dsplit[, dt := NULL] # no longer need original split data

  mobj$by_subject <- dsplit
  return(mobj)
}

#' Helper function to re-calculate the l3 design matrix for a model based on available data
#'
#' @param mobj a \code{hi_model_spec} object containing the L3 GLM model to run
#' @param data The run-level data frame containing data for all ids and sessions. This will be split into individual chunks
#'
#' @return a modified copy of \code{mobj} where the $by_subject field has been added
#'
#' @details The function adds the $by_subject field, which contains the design matrices and contrasts
#'   for each subject and session in \code{data} based on the available data for that session. For example, if
#'   a subject is missing a few runs (or these are dropped from analysis), then some contrasts may change or drop out
#'   of the model.
#'
#' The $by_subject field is a keyed data.table object containing list elements for the cope_list (mapping cope numbers
#'   to contrast names), the contrasts, and the design matrix for each session.
#'
#' @keywords internal
#' @importFrom checkmate assert_data_frame assert_multi_class assert_subset
#' @importFrom data.table data.table

mobj_refit_lm <- function(mobj, new_data) {

}

#' internal helper function to setup a linear model for a given l1, l2, or l3 model
#' 
#' @param mobj a model object to be populated or modified
#' @param model_formula a formula of the model to be fit
#' @param data a data.frame containing all columns used in model fitting
#' @param id_cols a character vector of column names in \code{data} that identify the observations
#'   and can be used for merging the model against related datasets
#' @return a model object containing the fitted model
#' @keywords internal
mobj_fit_lm <- function(mobj=NULL, model_formula=NULL, data, id_cols=NULL, lg=NULL) {
  if (is.null(lg)) { lg <- lgr::get_logger() }
  # verify that we have an object of known structure
  checkmate::assert_multi_class(mobj, c("l1_model_spec", "hi_model_spec"), null.ok=TRUE)
  if (is.null(mobj)) {
    mobj <- list()
    class(mobj) <- c("list", "hi_model_spec")
  } else {
    # if mobj is passed in, but other elements are NULL, we are refitting a model to a new/updated dataset
    if (is.null(model_formula)) {
      model_formula <- mobj$model_formula
    }

    if (is.null(id_cols)) {
      id_cols <- names(mobj$metadata)
    }
  }

  checkmate::assert_formula(model_formula)
  checkmate::assert_data_frame(data)
  checkmate::assert_subset(id_cols, names(data))

  # verify that there is a dummy DV for lm fitting
  if (!"dummy" %in% names(data)) {
    data$dummy <- rnorm(nrow(data))
  }

  model_formula <- update.formula(model_formula, "dummy ~ .") # add LHS

  # use model formula from parent object
  # model_formula <- terms(mobj$lmfit)
  model_vars <- all.vars(model_formula)

  checkmate::assert_subset(model_vars, names(data)) # verify that all predictors are present

  # look for missingness on any predictor variable
  # https://stackoverflow.com/questions/64287986/create-variable-that-captures-if-there-are-missing-fields-in-4-string-variables
  data <- data %>%
    dplyr::select(!!model_vars, !!id_cols) %>% # just keep model-relevant variables
    mutate(any_miss = rowSums(is.na(select(., any_of(!!model_vars)))) > 0)

  miss_data <- data %>%
    dplyr::filter(any_miss == TRUE) %>%
    dplyr::select(-any_miss)

  # retain non-missing data
  data <- data %>%
    dplyr::filter(any_miss == FALSE) %>%
    dplyr::select(-any_miss)

  if (nrow(miss_data) > 0L) {
    lg$warn("Model data contain missing values for one or more covariates.")
    lg$warn("These observations will be dropped from model outputs!")
    lg$warn("%s", capture.output(print(miss_data)))
  }

  # fit model and populate model information
  mobj$lmfit <- lm(model_formula, data)
  mobj$model_formula <- model_formula
  mobj$model_matrix <- model.matrix(mobj$lmfit)
  mobj$regressors <- colnames(mobj$model_matrix) # actual regressors after expanding categorical variables
  mobj$metadata <- data %>% dplyr::select(!!id_cols) # for merging model matrix against identifying columns
  mobj$model_data <- data %>% dplyr::select(!!model_vars) #keep track of data used for fitting model for refitting in case of missing data

  # handle coefficient aliasing
  al <- alias(mobj$lmfit)
  if (!is.null(al$Complete)) {
    cat("Problems with aliased (redundant) terms in model.\n\n")
    bad_terms <- rownames(al$Complete)
    cat(paste(bad_terms, collapse = ", "), "\n\n")

    # find unaliased (good) terms in the model
    good_terms <- colnames(modelmat)[!(colnames(modelmat) %in% bad_terms)]
    good_terms <- good_terms[!good_terms == "(Intercept)"]

    if (length(good_terms) == 0L) {
      warning(
        "No unaliased (good) terms in this model, suggesting that all covariates are constant or dependent. ",
        "Reverting to intercept-only model."
      )
      good_terms <- "1"
    }

    # build new model formula with only good terms
    newf <- as.formula(paste("dummy ~", paste(good_terms, collapse = " + ")))

    # also generate a model data.frame that expands dummy codes for terms, retaining only good variables
    # mobj$model_matrix_noalias <- modelmat[, grep(":", good_terms, fixed = TRUE, value = TRUE, invert = TRUE)]
    # modeldf <- as.data.frame(mobj$model_matrix_noalias)
    # modeldf$dummy <- data$dummy # copy across dummy DV for fitting
    mobj$lmfit_noalias <- lm(newf, data)
    mobj$aliased_terms <- bad_terms

    # N.B. emmeans needs to calculate contrasts on the original design to see the factor structure
    # So, we also need to drop out columns from the emmeans linfct
  }

  return(mobj)

}


#OLD version
# respecify_l3_model <- function(mobj, data) {
#   checkmate::assert_multi_class(mobj, c("l1_model_spec", "hi_model_spec")) # verify that we have an object of known structure
#   checkmate::assert_data_frame(data)
#   checkmate::assert_subset(c("id", "session"), names(data))

#   # use model formula from parent object
#   model_formula <- terms(mobj$lmfit)

#   mobj$lmfit <- lm(model_formula, data = data)
#   mobj$regressors <- colnames(mobj$model_matrix)
#   mobj$model_matrix <- model.matrix(mobj$lmfit)
#   mobj$metadata <- data %>% dplyr::select(id, session)
#   mobj$model_data <- data

#   mobj <- get_contrasts_from_spec(mobj, mobj$lmfit)

#   return(mobj)
# }


# new version, but depends on mobj being specified under new approach
respecify_l3_model <- function(mobj, new_data) {
  checkmate::assert_class(mobj, "hi_model_spec") # verify that we have an object of known structure
  checkmate::assert_data_frame(new_data)
  checkmate::assert_subset(c("id", "session"), names(new_data))

  # use id and session in new_data to subset the full data from initial model specification
  full_data <- cbind(mobj$metadata, mobj$model_data)
  new_data <- new_data %>% dplyr::select(id, session) #just metadata of interest
  new_data <- new_data %>%
    dplyr::left_join(full_data, by=c("id", "session"))

  mobj <- mobj_fit_lm(mobj=mobj, data=new_data)

  mobj <- get_contrasts_from_spec(mobj, mobj$lmfit)

  return(mobj)
}

#' Helper function to obtain the number of volumes in a 4D nifti file
#' 
#' @details This function prefers to use fslval instead of an internal R library
#'   because both oro.nifti and RNifti are rather slow to obtain a single value from
#'   the NIfTI header
#' 
#' @param nifti a 4D nifti file
#' @return the number of volumes in \code{nifti}
#' @keywords internal
lookup_run_volumes <- function(nifti) {
  if (!file.exists(nifti)) return(NA_integer_)

  # fslval is much faster than any internal R command. Use it, if possible
  # TODO: Make this more robust, more like runFSLCommand with path expectations
  has_fslval <- system2("which", "fslval", stdout = FALSE, stderr = FALSE)
  has_fslval <- ifelse(has_fslval == 0L, TRUE, FALSE)
  if (isTRUE(has_fslval)) {
    nvol <- as.integer(system2("fslval", args = c(nifti, "dim4"), stdout = TRUE))
  } else {
    nvol <- oro.nifti::readNIfTI(nifti, read_data = FALSE)@dim_[5L]
  }

  # system.time(run_volumes <- oro.nifti::readNIfTI(run_nifti[nn], read_data = FALSE)@dim_[5L])
  # system.time(run_volumes <- RNifti::readNifti(run_nifti[nn], internal = TRUE))

  return(nvol)
}


#' Helper function to obtain the number of voxels in a 4D file for populating FEAT FSF
#'
#' @details This function prefers to use fslval instead of an internal R library
#'   because both oro.nifti and RNifti are rather slow to obtain a single value from
#'   the NIfTI header
#'
#' @param nifti a 4D nifti file
#' @return the number of voxels in \code{nifti} as calculated by x * y * z * t
#' @keywords internal
lookup_nvoxels <- function(nifti) {
  if (!file.exists(nifti)) return(NA_integer_)

  # fslval is much faster than any internal R command. Use it, if possible
  # TODO: Make this more robust, more like runFSLCommand with path expectations
  # TODO: should probably have a single global lookup for fslval in gpa (finalize step)
  has_fslval <- system2("which", "fslval", stdout = FALSE, stderr = FALSE)
  has_fslval <- ifelse(has_fslval == 0L, TRUE, FALSE)
  if (isTRUE(has_fslval)) {
    ndim <- as.integer(system2("fslval", args = c(nifti, "dim0"), stdout = TRUE))
    nvoxels <- prod(sapply(1:ndim, function(xx) {
      as.integer(system2("fslval", args = c(nifti, paste0("dim", xx)), stdout = TRUE))
    }))
  } else {
    nihead <- oro.nifti::readNIfTI(nifti, read_data = FALSE)
    ndim <- nihead@dim_[1L]
    nvoxels <- prod(nihead@dim_[2:(ndim + 1)])
  }

  return(nvoxels)
}

#' small helper function to parse duration syntax of days-hours:minutes:seconds
#'   into lubridate duration object
#'
#' @param str string containing a duration that may include a days specification
#' @importFrom lubridate hms
#' @keywords internal
dhms <- function(str) {
  checkmate::assert_string(str)
  if (grepl("^\\d+-", str, perl = TRUE)) {
    split_hyphen <- strsplit(str, "-", fixed = TRUE)[[1]]
    days <- as.numeric(split_hyphen[1])
    period <- lubridate::hms(split_hyphen[2:length(split_hyphen)])
    period@day <- days
  } else {
    period <- lubridate::hms(str)
  }
  return(period)
}

#' convert a number of hours to a days, hours, minutes, seconds format
#' 
#' @importFrom lubridate day hour minute second seconds_to_period dhours
#' @keywords internal
hours_to_dhms <- function(hours, frac=FALSE) {
  checkmate::assert_number(hours, lower = 0)
  dur <- lubridate::dhours(hours)
  period <- seconds_to_period(dur)

  if (isTRUE(frac)) {
    str <- sprintf("%02d:%02d:%.03f", hour(period), minute(period), second(period))
  } else {
    str <- sprintf("%02d:%02d:%02d", hour(period), minute(period), round(second(period)))
  }

  if (day(period) > 0) {
    str <- paste0(sprintf("%d-", day(period)), str)
  }

  return(str)
}

#' helper function to refresh l3 model status and save gpa object from batch pipeline back to its cache
#' 
#' @param gpa a glm_pipeline_arguments object
#' @return a refreshed version of the gpa object
#' @importFrom lgr get_logger
#' @export
cleanup_glm_pipeline <- function(gpa) {
  lg <- lgr::get_logger("glm_pipeline/cleanup_glm_pipeline")
  gpa <- refresh_feat_status(gpa, level = 1L, lg = lg)
  gpa <- refresh_feat_status(gpa, level = 2L, lg = lg)
  gpa <- refresh_feat_status(gpa, level = 3L, lg = lg)
  res <- tryCatch(saveRDS(gpa, file = gpa$output_locations$object_cache), error = function(e) {
    lg$error("Could not save gpa object to file: %s", gpa$output_locations$object_cache)
    return(NULL)
  })

  # replace copied L1 niftis with symlinks to the corresponding file
  if (isTRUE(gpa$fsl$replace_l1_nifti_symlink)) {
    feat_files <- file.path(
      gpa$l1_model_setup$fsl$feat_dir,
      "filtered_func_data.nii.gz"
    )
    feat_files <- feat_files[file.exists(feat_files)]
    input_files <- gpa$l1_model_setup$fsl$run_nifti[file.exists(feat_files)]
    for (ff in seq_along(feat_files)) {
      unlink(feat_files[ff])
      file.symlink(input_files[ff], feat_files[ff])
    }
  }

  return(gpa)
}

#' helper function to copy any missing fields in target
#' @param target a named list to be populated by defaults if fields are missing
#' @param defaults a named list containing default values
populate_defaults <- function(target = NULL, defaults) {
  if (is.null(target)) target <- list()
  if (is.data.frame(target) && nrow(target) == 1L && is.data.frame(defaults) && nrow(defaults) == 1L) {
    return_df <- TRUE
    target <- as.list(target)
    defaults <- as.list(defaults)
  } else {
    return_df <- FALSE
  }

  miss_fields <- setdiff(names(defaults), names(target))
  if (length(miss_fields) > 0L) {
    for (mm in miss_fields) {
      target[[mm]] <- defaults[[mm]]
    }
  }

  if (isTRUE(return_df)) target <- as.data.frame(target)
  return(target)
}

# little helper function to create named list from objects
named_list <- function(...) {
  vnames <- as.character(match.call())[-1]
  return(setNames(list(...), vnames))
}

#' Helper function to create named data.frame from a set of objects
#' 
#' @details This function helps with the problem of having several vectors
#'   in the workspace that you want to combine with data.fram
# example:
# 
# named_data_frame <- function(...) {
#   vnames <- as.character(match.call())[-1]
#   return(setNames(data.frame(...), vnames))
# }


enforce_glms_complete <- function(gpa, level=1L, lg=NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_integerish(level, lower=1L, upper=3L)
  checkmate::assert_class(lg, "Logger")

  obj_name <- glue("l{level}_model_setup")
  expect_class <- glue("l{level}_setup")
  obj <- gpa[[obj_name]]

  if (is.null(obj) || !inherits(obj, expect_class)) {
    msg <- sprintf("No l%d_model_setup found in the glm pipeline object.", level)
    lg$error(msg)
    stop(msg)
  } else {
    if ("fsl" %in% gpa$glm_software) {
      nmiss <- sum(obj$fsl$feat_complete == FALSE)
      nruns <- nrow(obj$fsl)
      if (nmiss == nruns) {
        msg <- sprintf("All feat runs in %s$fsl are incomplete.", obj_name)
        lg$error(msg)
        stop(msg)
      } else if (nmiss > 0) {
        lg$warn(
          "There are %d missing runs in %s$fsl. Using complete %d runs.",
          nmiss, obj_name, nruns - nmiss
        )
      }
    }
  }

  return(invisible(NULL))
}


#' helper function to ask user to choose models at a given level for further processing
#' @param gpa a \code{glm_pipeline_arguments} with models specified at a given level
#' @param model_names a user-specified string of models to process/use
#' @param level the level of GLM analysis to be specified (1, 2, or 3)
#' @return a character vector of user-specified model names at this level
#' @keywords internal
choose_glm_models <- function(gpa, model_names, level, lg=NULL) {
  checkmate::assert_integerish(level, min = 1, max = 3)
  all_m_names <- names(gpa[[paste0("l", level, "_models")]]$models)
  checkmate::assert_subset(model_names, c("prompt", "all", "none", all_m_names))
  if (is.null(lg)) lg <- lgr::get_logger("")
  checkmate::assert_class(lg, "Logger")

  if (is.null(model_names)) {
    chosen_models <- NULL # happens when user de-selects all models at one level
  } else if (model_names[1L] == "all") {
    chosen_models <- all_m_names
  } else if (model_names[1L] == "none") {
    chosen_models <- NULL
  } else if (model_names[1L] == "prompt") {
    cat("\n")
    chosen_models <- select.list(all_m_names,
      multiple = TRUE,
      title = paste("Choose all level", level, "models to include: ")
    )
    if (identical(chosen_models, character(0))) {
      lg$info(paste("No level", level, "models were selected."))
      chosen_models <- NULL
    }
  } else {
    chosen_models <- model_names # user-specified set
  }
  return(chosen_models)
}