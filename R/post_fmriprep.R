## simple script to handle post-fmriprep processing
require(glue)
require(oro.nifti)
#nfsa

# matrix must be time x units/regions
mat_to_nii <- function(mat, ni_out="mat") {
  require(oro.nifti)
  if (is.data.frame(mat)) { mat <- as.amtrix(mat) }

  # this always puts regressors along the x dimension; y and z are singletons
  ydim <- zdim <- 1 # size of y and z dimensions
  xsz <- ysz <- zsz <- 1 # voxel size in x y z
  tr <- 1
  xorigin <- yorigin <- zorigin <- 0

  system(glue("fslcreatehd {nrow(mat)} {ydim} {zdim} {ncol(mat)} {xsz} {ysz} {zsz} {tr} {xorigin} {yorigin} {zorigin} 64 {ni_out}"))

  ## read empty NIfTI into R
  nif <- readNIfTI(ni_out, reorient = FALSE)
  nif <- drop_img_dim(nif) # need to cleanup dim_ attribute to avoid writeNIfTI failure

  # populate nifti
  nif@.Data <- array(mat, dim = c(nrow(mat), 1, 1, ncol(mat))) # add singleton dimensions for y and z

  # write NIfTI with regressors back to file
  writeNIfTI(nif, filename = ni_out) # this returns the filename to the caller
}

nii_to_mat <- function(ni_in) {
  checkmate::assert_file_exists(ni_in)

  nii <- readNIfTI(ni_in, reorient = FALSE, rescale_data = FALSE)
  mat <- nii[ , 1, 1, ] # x and z
  return(mat)
}

runFSLCommand <- function(args, fsldir=NULL, echo=TRUE, run=TRUE, stdout=NULL, stderr=NULL, log_file=NULL, intern=FALSE) {
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
  fslsetup <- paste0("FSLDIR=", fsldir, "; PATH=${FSLDIR}/bin:${PATH}; . ${FSLDIR}/etc/fslconf/fsl.sh; ${FSLDIR}/bin/")
  fslcmd <- paste0(fslsetup, args)
  if (!is.null(stdout)) { fslcmd <- paste(fslcmd, ">", stdout) }
  if (!is.null(stderr)) { fslcmd <- paste(fslcmd, "2>", stderr) }
  #cat("FSL command: ", fslcmd, "\n")
  if (!is.null(log_file)) { cat(args, file=log_file, append=TRUE, sep="\n") }
  if (isTRUE(echo)) { cat(args, "\n") }
  if (isTRUE(run)) {
    retcode <- system(fslcmd, intern=intern)
    return(retcode)
  } else {
    return(0) # no run retcode
  }
}

out_file_exists <- function(in_file, prefix, overwrite=TRUE) {
  # helper subfunction to enforce hyphen after initial postprocessing prefix
  p <- function(in_file, prefix) {
    has_prefix <- grepl("^\\w+-sub-.*", in_file, perl = TRUE)
    if (isTRUE(has_prefix)) {
      return(prefix)
    } else {
      return(paste0(prefix, "-")) # need to append hyphen
    }
  }

  # handle extant file
  out_file <- glue("{p(in_file, prefix)}{in_file}")
  skip <- FALSE
  if (checkmate::test_file_exists(out_file)) {
    if (isFALSE(overwrite)) {
      message(glue("Processed image already exists: {out_file}. Skipping this step."))
      skip <- TRUE
    } else {
      message(glue("Overwriting image: {out_file}."))
    }
  }
  return(list(out_file=out_file, skip=skip))
}

temporal_filter <- function(in_file, prefix="f", low_pass_hz=0, high_pass_hz=1/120, tr=NULL, overwrite=FALSE, log_file=NULL) {
  #checkmate::assert_file_exists(in_file)
  checkmate::assert_string(prefix)
  checkmate::assert_number(low_pass_hz)
  checkmate::assert_number(high_pass_hz)
  checkmate::assert_number(tr, lower = 0.01, upper = 30)
  stopifnot(low_pass_hz < high_pass_hz)

  # handle extant file
  res <- out_file_exists(in_file, prefix, overwrite)
  if (isTRUE(res$skip)) {
    return(res$out_file) # skip out
  } else  {
    out_file <- res$out_file
  }

  # bptf specifies its filter cutoffs in terms of volumes, not frequencies
  fwhm_to_sigma <- sqrt(8 * log(2)) # Details here: https://www.mail-archive.com/hcp-users@humanconnectome.org/msg01393.html

  if (is.infinite(high_pass_hz)) {
    #message("Low-pass filtering")
    hp_volumes <- -1 # do not apply high-pass
  } else {
    hp_volumes <- 1 / (high_pass_hz * fwhm_to_sigma * tr)
  }

  if (is.infinite(low_pass_hz) || low_pass_hz==0) {
    #message("High-pass filtering")
    lp_volumes <- -1 # do not apply low-pass
  } else {
    lp_volumes <- 1 / (low_pass_hz * fwhm_to_sigma * tr)
  }

  temp_tmean <- tempfile()
  runFSLCommand(glue("fslmaths {in_file} -Tmean {temp_tmean}"), log_file=log_file)
  runFSLCommand(glue("fslmaths {in_file} -bptf {hp_volumes} {lp_volumes} -add {out_file} "), log_file=log_file)
  return(out_file)
}

apply_aroma <- function(in_file, brain_mask=NULL, prefix="a", mixing_file, noise_file, overwrite=FALSE, log_file=NULL) {
  #checkmate::assert_file_exists(in_file)
  checkmate::assert_string(prefix)
  checkmate::assert_file_exists(mixing_file)
  checkmate::assert_file_exists(noise_file)

  # handle extant file
  res <- out_file_exists(in_file, prefix, overwrite)
  if (isTRUE(res$skip)) {
    return(res$out_file) # skip out
  } else {
    out_file <- res$out_file
  }

  #just read in the comma-separated noise ICs
  noise_ics <- readLines(noise_file, warn=FALSE)

  cmd <- glue("fsl_regfilt -i {in_file} -o {out_file} -d {mixing_file} -f {noise_ics}")
  if (!is.null(brain_mask) && checkmate::test_file_exists(brain_mask)) {
    cmd <- glue("{cmd} -m {brain_mask}")
  }

  runFSLCommand(cmd, log_file=log_file)
  return(out_file)
}

spatial_smooth <- function(in_file, prefix="s", fwhm_mm=6, brain_mask=NULL, overwrite=FALSE, log_file=NULL) {
  #checkmate::assert_file_exists(in_file)

  # handle extant file
  res <- out_file_exists(in_file, prefix, overwrite)
  if (isTRUE(res$skip)) {
    return(res$out_file) # skip out
  } else {
    out_file <- res$out_file
  }

  fwhm_to_sigma <- sqrt(8 * log(2)) # Details here: https://www.mail-archive.com/hcp-users@humanconnectome.org/msg01393.html
  sigma <- fwhm_mm / fwhm_to_sigma

  p2_intensity <- get_image_quantile(in_file, brain_mask, 2, log_file=log_file)
  median_intensity <- get_image_quantile(in_file, brain_mask, 50, log_file = log_file)
  susan_thresh <- (median_intensity - p2_intensity) * .75  # also see featlib.tcl

  # compute mean functional image used in susan
  temp_tmean <- tempfile()
  runFSLCommand(glue("fslmaths {in_file} -Tmean {temp_tmean}"), log_file=log_file) # save tmean to temporary file
  runFSLCommand(glue("susan {in_file} {sigma} 3 1 1 {temp_tmean} {susan_thresh} {out_file}"), log_file=log_file)
  if (checkmate::test_file_exists(temp_tmean))  { unlink(temp_tmean) } # cleanup
  return(out_file)

}

get_image_quantile <- function(in_file, brain_mask=NULL, quantile=50, log_file=NULL) {
  #checkmate::assert_file_exists(in_file)
  checkmate::assert_number(quantile, lower=0, upper=100)
  if (is.null(brain_mask)) {
     # median of non-zero voxels
     quantile_value <- as.numeric(runFSLCommand(glue("fslstats {in_file} -P {quantile}"), intern = TRUE, log_file = log_file))
  } else {
    checkmate::assert_file_exists(brain_mask)
    # median of all voxels in mask
    quantile_value <- as.numeric(runFSLCommand(glue("fslstats {in_file} -k {brain_mask} -p {quantile}"), intern = TRUE, log_file = log_file))
  }
  return(quantile_value)
}

intensity_normalize <- function(in_file, prefix="n", brain_mask=NULL, global_median=10000, overwrite=FALSE, log_file=NULL) {
  #checkmate::assert_file_exists(in_file)
  checkmate::assert_string(prefix)
  checkmate::assert_number(global_median)

  # handle extant file
  res <- out_file_exists(in_file, prefix, overwrite)
  if (isTRUE(res$skip)) {
    return(res$out_file) # skip out
  } else {
    out_file <- res$out_file
  }

  median_intensity <- get_image_quantile(in_file, brain_mask, 50, log_file=log_file)
  rescaling_factor <- global_median / median_intensity

  runFSLCommand(glue("fslmaths {in_file} -mul {rescaling_factor} {out_file} -odt float"), log_file=log_file)
  return(out_file)
}

confound_regression <- function(in_file, to_regress=NULL, prefix="r", brain_mask=NULL, overwrite=FALSE, log_file=NULL) {
  #checkmate::assert_file_exists(in_file)
  checkmate::assert_file_exists(to_regress)
  checkmate::assert_string(prefix)
  checkmate::assert_number(global_median)

  # handle extant file
  res <- out_file_exists(in_file, prefix, overwrite)
  if (isTRUE(res$skip)) {
    return(res$out_file) # skip out
  } else {
    out_file <- res$out_file
  }

  runFSLCommand(glue("fsl_glm -i {in_file} -d {to_regress} -m {brain_mask} --out_res={out_file}"), log_file=log_file)
  return(out_file)
}

get_fmriprep_outputs <- function(in_file) {
  first_char <- sub("(sub-\\d+_task-[^_]+_run-\\d+).*", "\\1", in_file, perl=TRUE)
  bold <- Sys.glob(glue("{first_char}*preproc_bold*nii*"))
  brain_mask <- Sys.glob(glue("{first_char}*_desc-brain_mask*nii*"))
  confounds <- glue("{first_char}_desc-confounds_regressors.tsv")
  melodic_mix <- glue("{first_char}_desc-MELODIC_mixing.tsv")
  noise_ics <- glue("{first_char}_AROMAnoiseICs.csv")
  ret_list <- list(bold = bold, brain_mask = brain_mask, confounds = confounds, melodic_mix = melodic_mix, noise_ics = noise_ics)
  ret_list <- lapply(ret_list, function(x) {
    ifelse(checkmate::test_file_exists(x), x, NULL)
  }) # NULL out missing files
  return(ret_list)
}

process_subject <- function(in_file, config_file="post_fmriprep.yaml") {
  checkmate::assert_file_exists(in_file)
  checkmate::assert_file_exists(config_file)
  #checkmate::assert_list(processing_sequence)
  proc_files <- get_fmriprep_outputs(in_file)

  cfg <- yaml::read_yaml(config_file)

  # currently insist on yaml
  # if (checkmate::test_file_exists(config_file)) {
  #   cfg <- yaml::read_yaml(config_file)
  # } else {
  #   processing_sequence <- c("spatial_smooth", "apply_aroma", "temporal_filter", "intensity_normalize")
  # }

  cur_file <- proc_files$bold

  if ("confound_regression" %in% cfg$processing_sequence || isTRUE(cfg$confound_calculate$compute)) {
    confounds <- data.table::fread(proc_files$confounds, na.strings = c("n/a", "NA"))
    confound_cols <- union(cfg$confound_regression$columns, cfg$confound_calculate$columns)
    confounds <- subset(confounds, select = confound_cols)
    confound_nii <- mat_to_nii(confounds, ni_out = tempfile(pattern = "confounds"))
    if ("apply_aroma" %in% cfg$processing_sequence) {
      confound_nii <- apply_aroma(confound_nii,
        brain_mask = proc_files$brain_mask, mixing_file = proc_files$melodic_mix,
        noise_file = proc_files$noise_ics, overwrite=cfg$overwrite, log_file=cfg$log_file
      )
    }
    if ("temporal_filter" %in% cfg$processing_sequence) {
      confound_nii <- temporal_filter(confound_nii,
        tr = cfg$tr, low_pass_hz = cfg$temporal_filter$low_pass_hz,
        high_pass_hz = cfg$temporal_filter$high_pass_hz, overwrite=cfg$overwrite, log_file=cfg$log_file
      )
    }

    if (isTRUE(cfg$confound_calculate$compute)) {
      filtered_confounds <- nii_to_mat(confound_nii)
      data.table::fwrite(subset(filtered_confounds, select = cfg$confound_calculate$columns), file = cfg$confound_calculate$output_file)
    }

    if ("confound_regression" %in% cfg$processing_sequence) {
      to_regress <- tempfile(pattern = "filtered_confounds")
      data.table::fwrite(subset(filtered_confounds, select = cfg$confound_regression$columns), file = to_regress)
    }

  }

  for (step in cfg$processing_sequence) {
    if (step == "spatial_smooth") {
      cur_file <- spatial_smooth(cur_file,
        brain_mask = proc_files$brain_mask, prefix = cfg$spatial_smooth$prefix,
        fwhm_mm = cfg$spatial_smooth$fwhm_mm, overwrite = cfg$overwrite, log_file = cfg$log_file
      )
    } else if (step == "apply_aroma") {
      cur_file <- apply_aroma(cur_file,
        brain_mask = proc_files$brain_mask, mixing_file = proc_files$melodic_mix,
        noise_file = proc_files$noise_ics,
        overwrite=cfg$overwrite, log_file=cfg$log_file
      )
    } else if (step == "temporal_filter") {
      cur_file <- temporal_filter(cur_file,
        tr = cfg$tr, low_pass_hz = cfg$temporal_filter$low_pass_hz,
        high_pass_hz = cfg$temporal_filter$high_pass_hz,
        overwrite=cfg$overwrite, log_file=cfg$log_file
      )
    } else if (step == "intensity_normalize") {
      cur_file <- intensity_normalize(cur_file,
        brain_mask = proc_files$brain_mask,
        global_median = cfg$intensity_normalize$global_median,
        overwrite=cfg$overwrite, log_file=cfg$log_file
      )
    } else if (step == "confound_regression") {
      cur_file <- confound_regression(cur_file,
        brain_mask = proc_files$brain_mask,
        to_regress = to_regress,
        overwrite=cfg$overwrite, log_file=cfg$log_file
      )
    }
  }

}

sdir <- "/proj/mnhallqlab/studies/bsocial/clpipe/data_fmriprep/fmriprep/sub-221256/func"
setwd(sdir)
process_subject("sub-221256_task-clock_run-2_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz",
  config_file = "/proj/mnhallqlab/users/michael/fmri.pipeline/R/post_fmriprep.yaml"
)

