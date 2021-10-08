## simple script to handle post-fmriprep processing

#nfsa

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

temporal_filter <- function(in_file, prefix="f", low_hz=1/120, high_hz=Inf tr=NULL) {
  checkmate::assert_file_exists(in_file)
  checkmate::assert_string(prefix)
  checkmate::assert_number(low_hz)
  checkmate::assert_number(high_hz)
  checkmate::assert_number(tr, lower=0.01, upper=30)

  if (is.infinite(high_hz)) {
    message("No high-pass filtering")
  }

  # bptf specifies its filter cutoffs in terms of volumes, not frequencies
  fwhm_to_sigma <- sqrt(8 * log(2)) # Details here: https://www.mail-archive.com/hcp-users@humanconnectome.org/msg01393.html
  lp_volumes <- 1 / (low_hz * fwhm_to_sigma * tr)
  

}

apply_aroma <- function(in_file, out_file, in_confounds, out_confounds) {

}

spatial_smooth <- function(in_file, prefix="s", fwhm_mm=6, brain_mask=NULL) {
  checkmate::assert_file_exists(in_file)

  # handle extant file
  out_file <- glue("{prefix}{in_file}")
  if (checkmate::test_file_exists(out_file)) {
    if (isFALSE(overwrite)) {
      message(glue("Spatially smoothed image already exists: {out_file}. Skipping this step.")
      return(out_file)
    } else {
      message(glue("Overwriting spatially smoothed image: {out_file}.")
    }
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
  checkmate::assert_file_exists(in_file)
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
  checkmate::assert_file_exists(in_file)
  checkmate::assert_string(prefix)
  checkmate::assert_number(global_median)

  # handle extant file
  out_file <- glue("{prefix}{in_file}")
  if (checkmate::test_file_exists(out_file)) {
    if (isFALSE(overwrite)) {
      message(glue("Intensity normalized image already exists: {out_file}. Skipping this step.")
      return(out_file)
    } else {
      message(glue("Overwriting intensity-normalized image: {out_file}.")
    }
  }

  median_intensity <- get_image_quantile(in_file, brain_mask, 50)
  rescaling_factor <- global_median / median_intensity
  
  runFSLCommand(glue("fslmaths {in_file} -mul {rescaling_factor} {prefix}{in_file} -odt float"))
  return(out_file)
}

process_subject <- function(in_file, processing_steps) {
  checkmate::assert_file_exists(in_file)
  checkmate::assert_list(processing_steps)

  for (step in processing_steps) {
    if (step$name == "spatial_smooth") {

    } else if (step)

  }
}