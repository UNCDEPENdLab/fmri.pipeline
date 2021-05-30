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
  fslsetup=paste0("FSLDIR=", fsldir, "; PATH=${FSLDIR}/bin:${PATH}; . ${FSLDIR}/etc/fslconf/fsl.sh; LD_LIBRARY_PATH=/gpfs/group/mnh5174/default/sw/openblas/lib ${FSLDIR}/bin/")
  fslcmd=paste0(fslsetup, args)
  if (!is.null(stdout)) { fslcmd=paste(fslcmd, ">", stdout) }
  if (!is.null(stderr)) { fslcmd=paste(fslcmd, "2>", stderr) }
  cat("FSL command: ", fslcmd, "\n")
  retcode <- system(fslcmd)
  return(retcode)
}

gen_emo_interaction_regressors <- function(examplefile, regressors, emotions=c("fear","scram","happy"), timingdir, mrrunnums, runlengths, dropVolumes) {
    ##make between-session regressors based on emotion condition.
    ##fmriDir <- "/Volumes/Serena/MMClock/MR_Raw"
    fmriDir <- "/Volumes/Serena/MMClock/MR_Proc"
    ##fmriDir <- "/Volumes/Serena/SPECC/MR_Proc"
    fitDir <- file.path(getMainDir(), "clock_analysis", "fmri", "fmri_fits")
    ##/Volumes/Serena/MMClock/MR_Raw/10997_20140308/MBclock_recon/clock1/nfswudktm_clock1_5_trunc282.nii.gz
    subid <- factor(sub(paste0(fmriDir, "/([0-9]{5})_\\d+/(?:mni_5mm_wavelet|native_nosmooth)/.*$"), "\\1", examplefile, perl=TRUE))
    ##subid <- as.integer(sub(paste0(fmriDir, "/([0-9]{3})[A-z]{2}_.*/mni_5mm_wavelet/.*$"), "\\1", examplefile, perl=TRUE)) #for SPECC
    
    loc <- local({load(file.path(fitDir, paste0(as.character(subid), "_fitinfo.RData"))); environment()})$f #time-clock fit object (load as local var)
    emocon <- data.frame(emotion=loc$run_condition[mrrunnums], contingency=loc$rew_function[mrrunnums]) #vector of emotion and contingency

    ##generate interactions for ev, rpe_neg, and rpe_pos with emotion
    csum <- cumsum(runlengths - dropVolumes)
    
    for (reg in regressors) {
        for (emo in emotions) {
            vec <- rep(0, max(csum))
            rmatch <- which(emocon$emotion == emo)
            for (r in rmatch) {
                f <- read.table(file.path(timingdir, paste0("run", r, "_", reg, ".1D")))$V1
                start <- ifelse(r > 1, csum[r-1] + 1, 1)
                end <- csum[r] #if (r < nrow(emocon)) csum[r] else 
                vec[start:end] <- f
            }
            write.table(vec, file=file.path(timingdir, paste0(reg, "_", emo, "_concat.1D")), row.names=FALSE, col.names=FALSE)
        }
    }
}


truncate_runs <- function(s, mrfiles, mrrunnums, niftivols, drop_volumes=0) {
  ##Identify the last valid volume acquired in a given run.
  ##Subjects often exhibit head movement after run ends (MATLAB closes), but scan hasn't stopped
  ##This occurs because the MB raw transfer of the prior run is occurring, but does not finish before the current run
  ##Thus, truncate mr files to be 12 seconds after final feedback presentation, which is how the paradigm timing files are setup
  ##note that all of this would need to be reworked if TR were not 1.0 (i.e., 1 second = 1 volume)

  require(dplyr)

  mrdf <- do.call(rbind, lapply(seq_along(mrfiles), function(r) {
    #iti_durations <- s$runs[[ mrrunnums[r] ]]$orig_data_frame$iti_ideal #for clockfit objects
    #last_iti <- s$runs[[ mrrunnums[r] ]]$iti_onset[length(s$runs[[ mrrunnums[r] ]]$iti_onset)]

    iti_durations <- s %>% dplyr::filter(run == mrrunnums[r]) %>% dplyr::pull(iti_ideal) #for outputs from parse_sceptic_outputs (2018+)
    last_iti <- s %>% dplyr::filter(run == mrrunnums[r]) %>% dplyr::pull(iti_onset) %>% tail(n=1)

    last_vol_behavior <- floor(last_iti + iti_durations[length(iti_durations)]) #use floor to select last vol in the iti window
    first_vol <- drop_volumes #first volume to use for analysis 

    if (last_vol_behavior < niftivols[r]) {
      ##more vols were acquired than presented in paradigm. Thus, truncation may be needed
      ##check framewise displacement and truncate earlier than 12 second ITI if a big movement occurred
      fd <- read.table(file.path(dirname(mrfiles[r]), "motion_info", "fd.txt"))$V1
      badfd <- do.call(c, sapply(seq_along(fd), function(x) { if (x >= last_iti && fd[x] > 0.9) x else NULL })) #flag volumes after last_iti with high FD
      if (length(badfd) == 0L) {
        ##no frames flagged in last volumes
        last_vol_analysis <- last_vol_behavior
      } else {
        ##use either the last volume of the task or the volume before the earliest bad movement 
        last_vol_analysis <- min(last_vol_behavior, (min(badfd) - 1))
      }

      #length of truncated file
      trunc_length <- last_vol_analysis - first_vol

      #generate filename for truncated volume
      if (first_vol > 0) {
        truncfile <- sub("(^.*/[a-z]+_clock[0-9](?:_5)*)\\.nii\\.gz$",
          paste0("\\1_drop", drop_volumes, "_trunc", last_vol_analysis, ".nii.gz"), mrfiles[r],
          perl = TRUE
        )
      } else {
        truncfile <- sub("(^.*/[a-z]+_clock[0-9](?:_5)*)\\.nii\\.gz$",
          paste0("\\1_trunc", last_vol_analysis, ".nii.gz"), mrfiles[r],
          perl = TRUE
        )
      }
      
      if (!file.exists(truncfile)) { runFSLCommand(paste("fslroi", mrfiles[r], truncfile, first_vol, trunc_length)) } #create truncated volume
      mrfile_to_analyze <- truncfile
    } else {
      last_vol_analysis <- niftivols[r] 
      if (drop_volumes > 0) {
        trunc_length <- niftivols[r] - drop_volumes
        truncfile <- sub("(^.*/[a-z]+_clock[0-9](?:_5)*)\\.nii\\.gz$",
          paste0("\\1_drop", drop_volumes, ".nii.gz"), mrfiles[r],
          perl = TRUE
        )
        if (!file.exists(truncfile)) { runFSLCommand(paste("fslroi", mrfiles[r], truncfile, first_vol, trunc_length)) } #create truncated volume
        mrfile_to_analyze <- truncfile
      } else {
        mrfile_to_analyze <- mrfiles[r] #just use original file  
      }

    }
    #cat(paste0(paste(mrfiles[r], niftivols[r], floor(last_iti), trunc_length, sep="\t"), "\n"), file="trunclog", append=TRUE)
    return(data.frame(last_vol_analysis, mrfile_to_analyze, stringsAsFactors=FALSE))
  }))

  mrdf

}

#' helper function for generating motion regressors from raw 6-parameter motion coregistration
#'
#' @param motion_params file containing 6 motion parameters
#' @param col.names names of columns in \code{motion_params}, in order from left to right
#' @param regressors names of regressors to generate and return to caller
#' @param drop_volumes number of volumes to drop from beginning of motion params
#' @param last_volume final volume to include from motion params (if truncated at end). If \code{NULL},
#'   then the end of the time series is not truncated.
#' @param rot_units The units of the rotation parameters. Default is "rad" for radians
#' @param tra_units The units of the translation parameter. Default is "mm" for millimeters
#'
#' @importFrom data.table fread
#' @importFrom checkmate assert_file_exists assert_subset assert_integerish
#' @keywords internal
generate_motion_regressors <- function(motion_params = "motion.par",
                                       col.names = c("rx", "ry", "rz", "tx", "ty", "tz"),
                                       regressors = c("rx", "ry", "rz", "tx", "ty", "tz"),
                                       demean = TRUE, drop_volumes = 0L, last_volume=NULL,
                                       rot_units="rad", tra_units="mm") {

  checkmate::assert_file_exists(motion_params)
  if (is.null(col.names)) { col.names <- c("rx", "ry", "rz", "tx", "ty", "tz") } #explicit defaults in case of null
  checkmate::assert_subset(col.names, c("rx", "ry", "rz", "tx", "ty", "tz"))
  checkmate::assert_logical(demean)
  checkmate::assert_integerish(drop_volumes, lower=0)

  # Regressors that are not motion-related can be passed in from outside (since we have a combined 
  # l1_confound_regressors argument). These could include CSF, WM, or whatever. Subset down to just
  # the values that can be calculated from motion params alone so that the logic below of parameter naming holds up.
  regressors <- intersect(regressors,
    c("fd", "rx", "ry", "rz", "tx", "ty", "tz",
    "drx", "dry", "drz", "dtx", "dty", "dtz",
    "qdrx", "qdry", "qdrz", "qdtx", "qdty", "qdtz")
  )

  #if none of the confound regressors is a motion parameter, then quietly return NULL
  if (length(regressors) == 0L) { return(invisible(NULL)) }

  mot <- data.table::fread(motion_params, col.names=col.names)

  if (is.null(last_volume)) { last_volume <- nrow(mot) }
  checkmate::assert_integerish(last_volume, upper=nrow(mot))

  mot <- mot[(1+drop_volumes):last_volume, ]

  if ("fd" %in% regressors || any(derivcols <- grepl("^q?d{1}.*", regressors, perl=TRUE))) {
    motderiv <- mot[, lapply(.SD, function(x) { c(0, diff(x)) })]
    setnames(motderiv, paste0("d", names(mot))) #add delta to names
    mot <- cbind(mot, motderiv)
  }

  #quadratics always computed after derivative calculation
  if (any(quadcols <- grepl("^q{1}.*", regressors, perl=TRUE))) {
    motquad <- mot[, lapply(.SD, function(x) { x^2 })]
    data.table::setnames(motquad, paste0("q", names(mot))) #add delta to names
    mot <- cbind(mot, motquad)
  }

  if ("fd" %in% regressors) {
    #need to adapt in case of degrees (this is based on radians)
    #https://wiki.cam.ac.uk/bmuwiki/FMRI
    if (rot_units=="rad") {
      fd <- apply(mot.deriv[,c("drx", "dry", "drz")], 1, function(x) sum(2*pi*50*(abs(x)/360))) +
        apply(mot.deriv[,c("dtx", "dty", "dtz")], 1, function(x) sum(abs(x)))
    } else {
      stop("not done yet")
    }
    mot <- cbind(mot, fd=fd)
  }

  mot <- mot[, ..regressors] #keep regressors of interest
  if (isTRUE(demean)) {
    mot <- mot[, lapply(.SD, function(x) { x - mean(x, na.rm=TRUE) }) ] #demean all columns
  }

  ##just PCA motion on the current run
  ##mregressors <- pca_motion(mr_files[r], runlengths[r], motion_parfile="motion.par", numpcs=3, drop_volumes=drop_volumes)$motion_pcs_concat

  #need to adapt this implementation -- create spike regressors based on a motion threshold
  ## if (spikeregressors) { #incorporate spike regressors if requested (not used in conventional AROMA)
  ##   censorfile <- file.path(dirname(mr_files[rr]), "motion_info", "fd_0.9.mat")
  ##   if (file.exists(censorfile) && file.info(censorfile)$size > 0) {
  ##     censor <- read.table(censorfile, header=FALSE)
  ##     censor <- censor[(1+drop_volumes):runlengths[rr],,drop=FALSE] #need no drop here in case there is just a single volume to censor
  ##     #if the spikes fall outside of the rows selected above, we will obtain an all-zero column. remove these
  ##     censor <- censor[,sapply(censor, sum) > 0,drop=FALSE]
  ##     if (ncol(censor) == 0L) { censor <- NULL } #no volumes to censor within valid timepoints
  ##     mregressors <- censor
  ##   }
  ## }

  ##add CSF and WM regressors (with their derivatives)
  ## nuisancefile <- file.path(dirname(mr_files[rr]), "nuisance_regressors.txt")
  ## if (file.exists(nuisancefile)) {
  ##   nuisance <- read.table(nuisancefile, header=FALSE)
  ##   nuisance <- nuisance[(1+drop_volumes):runlengths[rr],,drop=FALSE]
  ##   nuisance <- as.data.frame(lapply(nuisance, function(col) { col - mean(col) })) #demean
  ##   #cat("about to cbind with nuisance\n")
  ##   #print(str(mregressors))
  ##   #print(str(nuisance))
  ##   if (!is.null(mregressors)) { mregressors <- cbind(mregressors, nuisance) #note that in R 3.3.0, cbind with NULL or c() is no problem...
  ##   } else { mregressors <- nuisance }
  ## }
  
  return(mot)
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
#' 
pca_motion <- function(motion_df, num_pcs=3L, zscore=TRUE, verbose=FALSE) {
  checkmate::assert_data_frame(motion_df)
  checkmate::assert_integerish(num_pcs, lower=1, upper=50)
  checkmate::assert_logical(zscore)
  checkmate::assert_logical(verbose)

  #compute the PCA decomposition of motion parameters and their derivatives
  pc <- prcomp(motion_df, retx=TRUE, scale.=zscore)
  cumvar <- cumsum(pc$sdev^2/sum(pc$sdev^2))

  if (isTRUE(verbose)) message("first", num_pcs, "motion principal components account for: ", round(cumvar[num_pcs], 3))
  mregressors <- pc$x[,1:num_pcs] #cf Churchill et al. 2012 PLoS ONE
  attr(mregressors, "variance.explained") <- cumvar[num_pcs]

  if (isTRUE(verbose)) {
    cat("Multiple correlation of motion parameters:\n\n")
    print(round(sqrt(psych::smc(motion_df)), 2))
  }

  return(as.data.frame(mregressors))
}


#' helper function to generate confounds txt file for inclusion as additional regressors
#' @param id The subject id
#' @param session The session number
#' @param run_number The run number
#' @param gpa a \code{glm_pipeline_arguments} object containing pipeline specification
get_confound_txt <- function(id = NULL, session = NULL, run_number = NULL, gpa, drop_volumes=0L, last_volume=NULL, demean=TRUE) {
  if (checkmate::test_null(id)) { stop("get_confound_txt requires a specific id for lookup") }
  checkmate::assert_integerish(session, null.ok=TRUE)
  if (is.null(session)) session <- 1
  checkmate::assert_integerish(run_number, lower=1, null.ok = FALSE)
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_integerish(drop_volumes, lower = 0)
  checkmate::assert_integerish(last_volume, null.ok = TRUE)

  lg <- lgr::get_logger("glm_pipeline/l1_setup")

  rinfo <- gpa$run_data %>% dplyr::filter(id == !!id & session == !!session & run_number == !!run_number)
  if (nrow(rinfo) > 1L) {
    print(rinfo)
    lg$error("Multiple matches for a single run in get_confound_txt.")
    return(NULL)
  } else if (nrow(rinfo) == 0L) {
    lg$error("Unable to locate a record in gpa$run_data for id %s, session %s, run_number %s.", id, session, run_number)
    return(NULL)
  }

  # Park these in an analysis-level subfolder, not a particular model, since they are re-used across models.
  # Note: normalizePath will fail to evaluate properly if directory does not exist
  analysis_outdir <- get_l1_directory(id=id, session=session, gpa=gpa, create_if_missing = TRUE)

  expect_file <- file.path(analysis_outdir, paste0("run", run_number, "_l1_confounds.txt"))
  if (file.exists(expect_file)) {
    lg$debug("Returning extant file: %s in get_confound_txt", expect_file)
    return(expect_file)
  }

  confound_df <- NULL
  if (isTRUE(rinfo$confound_file_present[1])) {
    cfile <- get_mr_abspath(rinfo, "confound_file")[1]
    lg$debug("Reading confound file: %s", cfile)
    confound_df <- tryCatch(data.table::fread(cfile), error = function(e) {
      lg$error("Failed to read confound file: %s with error %s", cfile, as.character(e))
      return(NULL)
    })

    if (!is.null(gpa$confound_settings$confound_columns) && !is.null(confound_df)) {
      if (length(gpa$confound_settings$confound_columns) != ncol(confound_df)) {
        lg$warn(
          "Mismatch in number of columns in confound file: %s relative to $confound_settings$confound_columns",
          rinfo$confound_file[1]
        )
      }
      data.table::setnames(confound_df, gpa$confound_settings$confound_columns)
    }

    confound_df <- confound_df[(1 + drop_volumes):last_volume, ]
  }

  motion_df <- NULL
  if (isTRUE(rinfo$motion_params_present[1])) {
    mfile <- get_mr_abspath(rinfo, "motion_params")[1]
    lg$debug("Reading motion file: %s", mfile)
    motion_df <- tryCatch({
      generate_motion_regressors(
        mfile,
        col.names = gpa$confound_settings$motion_params_columns,
        regressors = gpa$confound_settings$l1_confound_regressors,
        drop_volumes = drop_volumes, last_volume = last_volume
      )}, error = function(e) {
      lg$error("Failed to read motion file: %s with error %s", rinfo$motion_params[1], as.character(e))
      return(NULL)
    })

    if (!is.null(gpa$confound_settings$motion_params_columns) && !is.null(motion_df)) {
      if (length(gpa$confound_settings$motion_params_columns) != ncol(motion_df)) {
        lg$warn(
          "Mismatch in number of columns in confound file: %s relative to $confound_settings$confound_columns",
          mfile
        )
      }
      data.table::setnames(motion_df, gpa$confound_settings$motion_params_columns)
    }

  }

  if (is.null(motion_df) && is.null(confound_df)) {
    lg$info("Neither confounds nor motion parameters are available for %s", mr_file)
    return(NULL)
  } else if (is.null(motion_df)) {
    confounds <- confound_df
  } else if (is.null(confound_df)) {
    confounds <- motion_df
  } else {
    if (nrow(motion_df) != nrow(confound_df)) {
      lg$error("Number of rows in motion_df is: %d and in confound_df is %d. Cannot combine", nrow(motion_df), nrow(confound_df))
      return(NULL)
    } else {
      # prefer confound_df as authoritative in cases where motion_df has overlapping columns
      overlap_names <- intersect(names(motion_df), names(confound_df))
      if (length(overlap_names) > 0L) {
        lg$info(
          "Motion parameters have overlapping columns with confounds file: %s. Preferring confounds to motion params",
          rinfo$confound_file[1]
        )
        lg$info("Overlap: %s", overlap_names)
        data.table::setnames(motion_df, old = overlap_names, new = paste0(overlap_names, ".mot"))
      }

      confounds <- dplyr::bind_cols(confound_df, motion_df)
    }
  }

  if (!all(gpa$confound_settings$l1_confound_regressors %in% names(confounds))) {
    lg$warn("Missing confound columns for subject: %s, session: %s", rinfo$id[1], rinfo$session[1])
    lg$warn("Column: %s", setdiff(gpa$confound_settings$l1_confound_regressors, names(confounds)))

    confounds <- confounds[, intersect(gpa$confound_settings$l1_confound_regressors, names(confounds))]
  }

  if (isTRUE(demean)) {
    lg$debug("Demeaning columns of confounds matrix")
    confounds <- as.data.frame(apply(confounds, 2, function(x) {
      x - mean(x, na.rm = TRUE)
    }))
  }

  lg$debug("Writing l1 confounds to file: %s", expect_file)
  write.table(confounds, file=expect_file, row.names=FALSE, col.names=FALSE)

  return(expect_file)
}

generateRunMask <- function(mrfiles, outdir=getwd(), outfile="runmask") {
  if (file.exists(file.path(outdir, paste0(outfile, ".nii.gz")))) { return(invisible(NULL)) }
  ##generate mask of mrfiles where temporal min is > 0 for all runs
  for (f in seq_along(mrfiles)) {
    runFSLCommand(paste0("fslmaths ", mrfiles[f], " -Tmin -bin ", outdir, "/tmin", f))#, fsldir="/usr/local/ni_tools/fsl")
  }

  ##sum mins together over runs and threshold at number of runs
  runFSLCommand(paste0("fslmaths ", paste(paste0(outdir, "/tmin", seq_along(mrfiles)), collapse=" -add "), " ", outdir, "/tminsum"))#, fsldir="/usr/local/ni_tools/fsl")
  runFSLCommand(paste0("fslmaths ", outdir, "/tminsum -thr ", length(mrfiles), " -bin ", outdir, "/", outfile))#, fsldir="/usr/local/ni_tools/fsl")
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
        wd <- NULL # defaults to getwd()
      } else {
        wd <- this_dir
      }
    }

    R.utils::getAbsolutePath(mr_df[[col]][ii], workDirectory = wd, expandTilde = TRUE)
  })

}

#TODO: This function returns irregular results... the desired results for different
# cases is not entirely clear.
# consider whether we want a 'what' type argument, like subject directory, analysis directory,
# run directory, model-specific analysis directory, etc.

#' small helper function to return the location of an l1 directory based on
#'   id, session, and run number
#' 
#' @param id The id of a participant
#' @param session The session number to lookup
#' @param run_number The run number to lookup
#' @param gpa A \code{glm_pipeline_arguments} object
#' @param glm_software which software is being used for the analysis (since directories may vary)
#' @param create_if_missing whether to create the directory if it does not exist
get_l1_directory <- function(id = NULL, session = NULL, run_number = NULL, model_name = NULL,
                             gpa, glm_software = "fsl", create_if_missing = FALSE) {
  if (checkmate::test_null(id)) {
    stop("get_l1_directory requires a specific id for lookup")
  }
  checkmate::assert_integerish(session, null.ok = TRUE)
  if (is.null(session)) session <- 1
  checkmate::assert_integerish(run_number, lower = 1, null.ok = TRUE)
  checkmate::assert_string(model_name, null.ok=TRUE)
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  stopifnot("run_data" %in% names(gpa))
  checkmate::assert_data_frame(gpa$run_data)
  checkmate::assert_logical(create_if_missing)
  if (!is.null(model_name)) checkmate::assert_subset(model_name, names(gpa$l1_models$models))

  lg <- lgr::get_logger("glm_pipeline/l1_setup")
  if (is.null(run_number)) {
    rinfo <- gpa$run_data %>% dplyr::filter(id == !!id & session == !!session)
  } else {
    rinfo <- gpa$run_data %>% dplyr::filter(id == !!id & session == !!session & run_number == !!run_number)
  }

  # if (nrow(rinfo) > 1L) {
  #   print(rinfo)
  #   lg$error("Multiple matches for a single run in get_l1_directory.")
  #   return(NULL)
  if (nrow(rinfo) == 0L) {
    lg$error("Unable to locate a record in gpa$run_data for id %s, session %s, run_number %s.", id, session, run_number)
    return(NULL)
  }

  l1_dir <- NULL
  if (gpa$output_settings$l1_directory == "local") {
    if ("mr_dir" %in% names(rinfo)) {
      lg$debug("Using mr_dir l1 directory lookup: %s", rinfo$mr_dir[1L])
      if (is.null(model_name)) {
        lg$debug("Lookup from analysis_name: %s", gpa$analysis_name)
        l1_dir <- file.path(normalizePath(file.path(rinfo$mr_dir[1L], "..")), gpa$analysis_name)
      } else {
        lg$debug("Lookup from model outdir: %s", gpa$l1_models$models[[model_name]]$outdir)
        l1_dir <- file.path(normalizePath(file.path(rinfo$mr_dir[1L], "..")), gpa$l1_models$models[[model_name]]$outdir)
      }
    } else {
      # look in parent folder of relevant run nifti and place l1 model there
      rn <- get_mr_abspath(rinfo[1, , drop = F], "run_nifti")
      lg$debug("Using run_nifti l1 directory lookup: %s", rn)

      if (is.null(model_name)) {
        lg$debug("Lookup from analysis_name: %s", gpa$analysis_name)
        l1_dir <- file.path(
          normalizePath(file.path(dirname(rn), "..")),
          gpa$analysis_name
        )
      } else {
        lg$debug("Lookup from model outdir: %s", gpa$l1_models$models[[model_name]]$outdir)
        l1_dir <- file.path(
          normalizePath(file.path(dirname(rn), "..")),
          gpa$l1_models$models[[model_name]]$outdir
        )
      }
    }
  } else {
    stop("not yet implemented in l1_get_directory")
  }

  if (isTRUE(create_if_missing) && !dir.exists(l1_dir)) {
    lg$debug("Create l1 directory: %s", l1_dir)
    dir.create(l1_dir, recursive = TRUE)
  }

  return(l1_dir)
}