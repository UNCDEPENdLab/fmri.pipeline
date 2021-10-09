## simple script to handle post-fmriprep processing

#nfsa

# matrix must be time x units/regions
mat_to_nii <- function(mat, ni_out="mat") {
  require(oro.nifti)

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
  writeNIfTI(nif, filename = ni_out)
}



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
  fslsetup <- paste0("FSLDIR=", fsldir, "; PATH=${FSLDIR}/bin:${PATH}; . ${FSLDIR}/etc/fslconf/fsl.sh; ${FSLDIR}/bin/")
  fslcmd <- paste0(fslsetup, args)
  if (!is.null(stdout)) { fslcmd=paste(fslcmd, ">", stdout) }
  if (!is.null(stderr)) { fslcmd=paste(fslcmd, "2>", stderr) }
  #cat("FSL command: ", fslcmd, "\n")
  cat(args, "\n")
  #retcode <- system(fslcmd)
  #return(retcode)
  return(NULL)
}

out_file_exists <- function(in_file, prefix, overwrite=TRUE) {
  # handle extant file
  out_file <- glue("{prefix}{in_file}")
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

temporal_filter <- function(in_file, prefix="f", low_pass_hz=0, high_pass_hz=1/120, tr=NULL, overwrite=FALSE) {
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
  runFSLCommand(glue("fslmaths {in_file} -Tmean {temp_tmean}"))
  runFSLCommand(glue("fslmaths {in_file} -bptf {hp_volumes} {lp_volumes} -add {out_file} "))
  return(out_file)
}

apply_aroma <- function(in_file, brain_mask=NULL, prefix="a", mixing_file, noise_file, overwrite=FALSE) {
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

  runFSLCommand(cmd)
  return(out_file)
}

spatial_smooth <- function(in_file, prefix="s", fwhm_mm=6, brain_mask=NULL) {
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

  p2_intensity <- get_image_quantile(in_file, brain_mask, 2)
  median_intensity <- get_image_quantile(in_file, brain_mask, 50)
  susan_thresh <- (median_intensity - p2_intensity) * .75  # also see featlib.tcl

  # compute mean functional image used in susan
  temp_tmean <- tempfile()
  runFSLCommand(glue("fslmaths {in_file} -Tmean {temp_tmean}")) # save tmean to temporary file
  runFSLCommand(glue("susan {in_file} {sigma} 3 1 1 {temp_tmean} {susan_thresh} {out_file}"))
  if (checkmate::test_file_exists(temp_tmean))  { unlink(temp_tmean) } # cleanup
  return(out_file)

}

get_image_quantile <- function(in_file, brain_mask=NULL, quantile=50) {
  #checkmate::assert_file_exists(in_file)
  checkmate::assert_number(quantile, lower=0, upper=100)
  if (is.null(brain_mask)) {
    quantile_value <- runFSLCommand(glue("fslstats {in_file} -P {quantile}")) # median of non-zero voxels
  } else {
    checkmate::assert_file_exists(brain_mask)
    quantile_value <- runFSLCommand(glue("fslstats {in_file} -k {brain_mask} -p {quantile}")) # median of all voxels in mask
  }
  return(quantile_value)
}

intensity_normalize <- function(in_file, prefix="n", brain_mask=NULL, global_median=10000, overwrite=FALSE) {
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

  median_intensity <- get_image_quantile(in_file, brain_mask, 50)
  rescaling_factor <- global_median / median_intensity
  
  runFSLCommand(glue("fslmaths {in_file} -mul {rescaling_factor} {prefix}{in_file} -odt float"))
  return(out_file)
}

get_fmriprep_outputs <- function(in_file) {
  prefix <- sub("(sub-\\d+_task-[^_]+_run-\\d+).*", "\\1", in_file, perl=TRUE)
  bold <- Sys.glob(glue("{prefix}*preproc_bold*nii*"))
  mask <- Sys.glob(glue("{prefix}*_desc-brain_mask*nii*"))
  confounds <- glue("{prefix}_desc-confounds_regressors.tsv")
  melodic_mix <- glue("{prefix}_desc-MELODIC_mixing.tsv")
  noise_ics <- glue("{prefix}_AROMAnoiseICs.csv")
  ret_list <- list(bold = bold, mask = mask, confounds = confounds, melodic_mix = melodic_mix, noise_ics = noise_ics)
  ret_list <- lapply(ret_list, function(x) {
    if (!checkmate::test_file_exists(x)) {
      return(NULL)
    } else {
      return(x)
    }
  })
  return(ret_list)
}

process_subject <- function(in_file, processing_steps=list()) {
  checkmate::assert_file_exists(in_file)
  checkmate::assert_list(processing_steps)
  proc_files <- get_fmriprep_outputs(in_file)

  res <- spatial_smooth(proc_files$bold, brain_mask=proc_files$mask)
  res <- apply_aroma(res, brain_mask=proc_files$mask, mixing_file=proc_files$melodic_mix, noise_file=proc_files$noise_ics)
  res <- temporal_filter(res, tr=1)
  res <- intensity_normalize(res, brain_mask=proc_files$mask)

  # for (step in processing_steps) {
  #   if (step$name == "spatial_smooth") {

  #   } else if (step) {

  #   }

  # }
}

process_subject("sub-221256_task-clock_run-2_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz")

