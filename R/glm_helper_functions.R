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
      run_nifti <- truncfile
    } else {
      last_vol_analysis <- niftivols[r]
      if (drop_volumes > 0) {
        trunc_length <- niftivols[r] - drop_volumes
        truncfile <- sub("(^.*/[a-z]+_clock[0-9](?:_5)*)\\.nii\\.gz$",
          paste0("\\1_drop", drop_volumes, ".nii.gz"), mrfiles[r],
          perl = TRUE
        )
        if (!file.exists(truncfile)) { # create truncated volume
          runFSLCommand(paste("fslroi", mrfiles[r], truncfile, first_vol, trunc_length))
        } 
        run_nifti <- truncfile
      } else {
        run_nifti <- mrfiles[r] #just use original file  
      }

    }
    #cat(paste0(paste(mrfiles[r], niftivols[r], floor(last_iti), trunc_length, sep="\t"), "\n"), file="trunclog", append=TRUE)
    return(data.frame(last_vol_analysis, run_nifti, stringsAsFactors=FALSE))
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
                                       demean = FALSE, drop_volumes = 0L, last_volume=NULL,
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
    c("FD", "rx", "ry", "rz", "tx", "ty", "tz",
    "drx", "dry", "drz", "dtx", "dty", "dtz",
    "qdrx", "qdry", "qdrz", "qdtx", "qdty", "qdtz")
  )

  #if none of the confound regressors is a motion parameter, then quietly return NULL
  if (length(regressors) == 0L) { return(invisible(NULL)) }

  mot <- data.table::fread(motion_params, col.names=col.names)

  if (is.null(last_volume)) { last_volume <- nrow(mot) }
  checkmate::assert_integerish(last_volume, upper=nrow(mot))

  mot <- mot[(1 + drop_volumes):last_volume, ]

  if ("FD" %in% regressors || any(derivcols <- grepl("^q?d{1}.*", regressors, perl=TRUE))) {
    mot_deriv <- mot[, lapply(.SD, function(x) { c(0, diff(x)) })]
    setnames(mot_deriv, paste0("d", names(mot))) #add delta to names
    mot <- cbind(mot, mot_deriv)
  }

  #quadratics always computed after derivative calculation
  if (any(quadcols <- grepl("^q{1}.*", regressors, perl=TRUE))) {
    mot_quad <- mot[, lapply(.SD, function(x) { x^2 })]
    data.table::setnames(mot_quad, paste0("q", names(mot))) #add delta to names
    mot <- cbind(mot, mot_quad)
  }

  if ("FD" %in% regressors) {
    #need to adapt in case of degrees (this is based on radians)
    #https://wiki.cam.ac.uk/bmuwiki/FMRI
    if (rot_units=="rad") {
      FD <- apply(mot_deriv[, c("drx", "dry", "drz")], 1, function(x) 50 * sum(abs(x))) +
        apply(mot_deriv[, c("dtx", "dty", "dtz")], 1, function(x) sum(abs(x)))
    } else if (rot_units=="deg") {
      FD <- apply(mot_deriv[, c("drx", "dry", "drz")], 1, function(x) 50 * (pi / 180) * sum(abs(x))) +
        apply(mot_deriv[, c("dtx", "dty", "dtz")], 1, function(x) sum(abs(x)))
    } else {
      stop("not done yet")
    }
    mot <- cbind(mot, FD=FD)
  }

  mot <- mot[, ..regressors] #keep regressors of interest
  if (isTRUE(demean)) {
    mot <- mot[, lapply(.SD, function(x) { x - mean(x, na.rm=TRUE) }) ] #demean all columns
  }

  ##just PCA motion on the current run
  ##mregressors <- pca_motion(run_nifti[r], runlengths[r], motion_parfile="motion.par", numpcs=3, drop_volumes=drop_volumes)$motion_pcs_concat

  #need to adapt this implementation -- create spike regressors based on a motion threshold
  ## if (spikeregressors) { #incorporate spike regressors if requested (not used in conventional AROMA)
  ##   censorfile <- file.path(dirname(run_nifti[rr]), "motion_info", "fd_0.9.mat")
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
  ## nuisancefile <- file.path(dirname(run_nifti[rr]), "nuisance_regressors.txt")
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
  lg <- lgr::get_logger("glm_pipeline/l2_setup")
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_class(gpa$l1_model_setup, "l1_setup")
  checkmate::assert_data_table(gpa$l1_model_setup$metadata)
  if (!"exclude_run" %in% names(gpa$l1_model_setup$metadata)) {
    msg <- "calculate_subject_exclusions depends on exclude_run being populated by setup_l1_models."
    lg$error(msg)
    stop(msg)
  }

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

  n_good_runs_df <- gpa$l1_model_setup$metadata %>%
    dplyr::filter(exclude_run == FALSE) %>%
    dplyr::count(id, session, name = "n_good_runs")

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

  # propagate subject exclusion down to $run_data for simplicity
  gpa$run_data <- gpa$run_data %>%
    dplyr::left_join(
      gpa$subject_data %>% dplyr::select(id, session, exclude_subject),
      by = c("id", "session")
    )

  # TODO: cleanup redundancy, or at least be more thoughtful about this
  # For now, also populate exclude_subject into gpa$l1_model_setup$metadata to avoid
  # additional joins in l2 model setup
  gpa$l1_model_setup$metadata <- gpa$l1_model_setup$metadata %>%
    dplyr::left_join(
      gpa$subject_data %>% dplyr::select(id, session, exclude_subject),
      by = c("id", "session")
    )

  return(gpa)
}

slurm_job_array <- function(job_name = "slurm_array") {

}

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

generate_subject_l2_models <- function(gpa) {

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
#' @importFrom emmeans emmeans emtrends pairs
get_contrasts_from_spec <- function(mobj, lmfit=NULL) {
  checkmate::assert_multi_class(mobj, c("l1_model_spec", "hi_model_spec")) # verify that we have an object of known structure
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

  ### add diagonal contrasts
  c_diagonal <- contrast_list$diagonal
  if (isTRUE(spec$diagonal) && is.null(c_diagonal)) {
    cnames <- spec$regressors
    diag_mat <- diag(length(cnames))
    rownames(diag_mat) <- paste0("EV_", cnames) # simple contrast naming for each individual regressor
    colnames(diag_mat) <- cnames # always have columns named by regressor

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
      rownames(econ) <- enames

      # add contrasts to matrix
      c_cond_means <- rbind(c_cond_means, econ)
    }
  }

  ### add cell means for all factors, if requested
  c_cell_means <- contrast_list$cell_means
  if (isTRUE(spec$cell_means) && is.null(c_cell_means)) {
    # get model-predicted means for each factor
    ee <- emmeans(lmfit, as.formula(paste("~", paste(spec$cat_vars, collapse = "*"))), weights = spec$weights)
    # pp <- pairs(ee)
    edata <- summary(ee)
    econ <- ee@linfct
    enames <- apply(edata[, spec$cat_vars, drop = FALSE], 1, paste, collapse = ".")

    # if any emmeans are not estimable, then this is likely due to aliasing. For now, drop
    which_na <- is.na(edata$emmean)
    if (any(which_na)) {
      econ <- econ[!which_na, , drop = FALSE]
      enames <- enames[!which_na]
    }

    rownames(econ) <- enames
    c_cell_means <- rbind(c_cell_means, econ)
  }

  ### add overall response mean
  c_overall <- contrast_list$overall
  if (isTRUE(spec$overall_response) && is.null(c_overall)) {
    ee <- emmeans(lmfit, ~1, weights = spec$weights)
    econ <- ee@linfct
    rownames(econ) <- "overall"
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
      enames <- make.names(sub("\\s+-\\s+", "_M_", edata$contrast, perl = TRUE)) # names of contrasts

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
  }

  c_custom <- contrast_list$custom

  #combine each element
  cmat_full <- rbind(
    c_diagonal,
    c_cond_means,
    c_cell_means,
    c_overall,
    c_pairwise_diffs,
    c_simple_slopes,
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

# get_run_lmfits <- function() {

# }


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
  for (vv in seq_len(nrow(dsplit))) {
    lmfit <- lm(model_formula, data = dsplit[[vv, "dt"]])

    mm <- get_contrasts_from_spec(mobj, lmfit)
    cope_df <- data.frame(
      id = dsplit$id[1L], session = dsplit$session[1L],
      l2_cope = seq_len(nrow(mm$contrasts)), l2_cope_names = rownames(mm$contrasts)
    )

    cope_list[[vv]] <- cope_df
    contrast_list[[vv]] <- mm$contrasts
    model_matrix_list[[vv]] <- model.matrix(lmfit)
  }

  dsplit[, cope_list := cope_list]
  dsplit[, contrasts := contrast_list]
  dsplit[, model_matrix := model_matrix_list]
  dsplit[, dt := NULL] # no longer need original split data

  mobj$by_subject <- dsplit
  return(mobj)
}