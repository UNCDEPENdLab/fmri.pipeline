#!/usr/bin/env Rscript
#
# This is a modification of the original pTFCE shell script: https://spisakt.github.io/pTFCE/
# The goal is to apply pTFCE to a single zstat image based on relevant inputs

print_help <- function() {
  cat(paste("This script runs pTFCE (probabilistic threshold-free cluster enhancement) on a single",
    "z-statistic image, as produced by FSL FEAT, for example.",
    "Options:",
    "  --zstat <z_img>: A NIfTI file containing a z-statistic map.",
    "  --mask <mask_img>: A required NIfTI file containing the brain mask for the --zstat input",
    "  --help: print the help",
    "  --residuals <residuals_img>: Estimate smoothness from an image containing GLM residuals.",
    "  --dof <integer>: The degrees of freedom for the corresponding GLM analysis.",
    "  --fsl_smoothest <smoothness_file>: A smoothness file produced by fsl smoothest, ideally from res4D.nii.gz",
    "  --twosided: If specified, then both the positive and negative z-statistics are corrected (DEFAULT).",
    "  --onesided: If specified, only the positive z-statistics are FWE-corrected and negative voxels are dropped.",
    "  --fwep <.05>: The p-values for which the map is corrected for multiple comparisons. Multiple values may be specified.",
    "  --write_thresh_imgs: If specified, then each p-value supplied in --fwep will be applied to the pTFCE-enhanced image, then saved",
    "  --verbose: If specified, the verbose option of pTFCE will be turned on, which prints out progress and diagnostics.",
    "\n\n",
    sep = "\n"
  ))
}

#read in command line arguments
args <- commandArgs(trailingOnly = FALSE)

scriptpath <- dirname(sub("--file=", "", grep("--file=", args, fixed=TRUE, value=TRUE), fixed=TRUE))
argpos <- grep("--args", args, fixed=TRUE)
if (length(argpos) > 0L) {
   args <- args[(argpos+1):length(args)]
} else {
  args <- c()
}

if (length(args) == 0L) {
  message("ptfce_zstat expects a single z-statistc image from a GLM analysis using --zstat <z_img> and a corresponding --mask <mask_img>.\n")
  print_help()
  quit(save = "no", 1, FALSE)
}

if (!require("pacman")) {
  install.packages("pacman")
  library(pacman)
}

pacman::p_load(oro.nifti, checkmate, pTFCE, dplyr)

z_img <- NA_character_
mask_img <- NA_character_
residuals_img <- NA_character_
dof <- NA_integer_
fsl_smoothest <- NA_character_
two_sided <- TRUE
fwe_p <- .05
write_thresh_imgs <- FALSE # whether to write hard-thresholded images at each FWE p-value
verbose <- FALSE

argpos <- 1
while (argpos <= length(args)) {
  #print(args[argpos])
  if (args[argpos] == "--zstat") {
    z_img <- args[argpos + 1] # name of z-stat image
    argpos <- argpos + 2
  } else if (args[argpos] == "--mask") {
    mask_img <- args[argpos + 1] # name of mask image
    argpos <- argpos + 2
  } else if (args[argpos] == "--residuals") {
    residuals_img <- args[argpos + 1] # name of residuals image
    checkmate::assert_file_exists(residuals_img)
    argpos <- argpos + 2
  } else if (args[argpos] == "--fsl_smoothest") {
    fsl_smoothest <- args[argpos + 1]
    checkmate::assert_file_exists(fsl_smoothest)
    argpos <- argpos + 2
  } else if (args[argpos] == "--help") {
    print_help()
    quit(save = "no", 0, FALSE)
  } else if (args[argpos] == "--dof") {
    dof_arg <- args[argpos + 1]
    if (checkmate::test_integerish(dof_arg)) {
      dof <- as.integer(dof_arg)
    } else if (checkmate::test_file_exists(dof_arg)) {
      dof <- as.integer(readLines(dof_arg))
    }
    argpos <- argpos + 2
  } else if (args[argpos] == "--twosided") {
    two_sided <- TRUE
    argpos <- argpos + 1
  } else if (args[argpos] == "--onesided") {
    two_sided <- FALSE
    argpos <- argpos + 1
  } else if (args[argpos] == "--fwep") {
    named_args <- grep("^--", args)
    if (any(named_args > argpos)) {
      next_arg <- named_args[named_args > argpos][1L] # first named argument after --fwep
      if (next_arg == argpos + 1) {
        stop("No valid --fwep values provided")
      }
      last_el <- next_arg - 1
    } else {
      last_el <- length(args)
    }
    fwep_inp <- as.numeric(args[(argpos + 1):last_el])
    checkmate::assert_numeric(fwep_inp, lower = 1e-10, upper = .9999, any.missing = FALSE)
    argpos <- argpos + 1 + length(fwep_inp)
    fwe_p <- fwep_inp
  } else if (args[argpos] == "--write_thresh_imgs") {
    write_thresh_imgs <- TRUE
    argpos <- argpos + 1
  } else if (args[argpos] == "--verbose") {
    verbose <- TRUE
    argpos <- argpos + 1
  } else {
    stop("Not sure what to do with argument: ", args[argpos])
  }
}

ptfce_worker <- function(args_list) {
  do.call(pTFCE::ptfce, args_list)
}

# for testing
# z_img <- "/proj/mnhallqlab/users/michael/mmclock_pe/mmclock_nov2021/feat_l3/L1m-abspe/L2m-l2_l2c-overall/L3m-int_only/FEAT_l1c-EV_abspe.gfeat/cope1.feat/stats/zstat1.nii.gz"
# mask_img <- "/proj/mnhallqlab/users/michael/mmclock_pe/mmclock_nov2021/feat_l3/L1m-abspe/L2m-l2_l2c-overall/L3m-int_only/FEAT_l1c-EV_abspe.gfeat/cope1.feat/mask.nii.gz"
# residuals_img <- "/proj/mnhallqlab/users/michael/mmclock_pe/mmclock_nov2021/feat_l3/L1m-abspe/L2m-l2_l2c-overall/L3m-int_only/FEAT_l1c-EV_abspe.gfeat/cope1.feat/stats/res4d.nii.gz"
# dof_arg <- "/proj/mnhallqlab/users/michael/mmclock_pe/mmclock_nov2021/feat_l3/L1m-abspe/L2m-l2_l2c-overall/L3m-int_only/FEAT_l1c-EV_abspe.gfeat/cope1.feat/stats/dof"
# fsl_smoothest <- "/proj/mnhallqlab/users/michael/mmclock_pe/mmclock_nov2021/feat_l3/L1m-abspe/L2m-l2_l2c-overall/L3m-int_only/FEAT_l1c-EV_abspe.gfeat/cope1.feat/stats/smoothness"

# both z and mask are required inputs
checkmate::assert_file_exists(z_img)
checkmate::assert_file_exists(mask_img)

z_dir <- normalizePath(dirname(z_img))
ext <- sub(".*(\\.nii(\\.gz)?)$", "\\1", z_img, perl = TRUE)
base <- sub(ext, "", basename(z_img), fixed = TRUE)

if (!is.na(residuals_img)) {
  if (is.na(dof)) {
    stop("If --residuals image is specified, then --dof must be provided, too!")
  }
  cat("Using residuals image to calculate smoothness\n")
  call_list <- list(residual=residuals_img, dof=dof, verbose = verbose)
} else if (!is.na(fsl_smoothest)) {
  cat("Using FSL smoothest file:", fsl_smoothest, "\n")
  smooth_data <- read.table(fsl_smoothest, nrow = 3) # DLH, VOLUME, RESELS
  V <- smooth_data[2, 2] # VOLUME
  Rd <- smooth_data[1, 2] * smooth_data[2, 2] # DLH * VOLUME
  resels <- smooth_data[3, 2] # RESELS

  call_list <- list(V = V, Rd = Rd, resels = resels, verbose = verbose)
} else {
  cat("Estimating smoothness internally from z-stat image\n")
  call_list <- list(verbose = verbose)
}

# run ptfce
cat("pTFCE running: ", as.character(Sys.time()), "\n")

Z <- oro.nifti::readNIfTI(z_img, reorient = FALSE)
mask <- oro.nifti::readNIfTI(mask_img, reorient = FALSE)

call_list[["mask"]] <- mask

# always use the positive z-stats for FWE (both one- and two-sided)
#z_pos <- Z # for now, we use the whole image (based on Smith and Nichols 2009)
#z_pos@.Data[z_pos@.Data < 0] <- 0 # zero negative values
call_list[["img"]] <- Z #z_pos
ptfce_pos <- ptfce_worker(call_list)

# The pTFCE correction produces massive negative values for the other tail of the distribution, including outside of
# the mask. I think this is just a coding oversight on their end since the positive values are the only ones of interest.
# Still, it makes for weird images, and if we want to combine two-sided images, we can't simply add them.
# Thus, zero out anything below zero in both cases for clarity.
ptfce_pos$Z@.Data[ptfce_pos$Z@.Data < 0] <- 0

# There is also a problem where some images (including mine) produce NaNs in the TFCE-enhanced image. This is documented here:
# https://github.com/spisakt/pTFCE/issues/8. But it is not entirely resolved. For now, rather than leave NAs in the image, set
# them to zero and produce a warning.
nmiss <- sum(is.na(ptfce_pos$Z))
if (nmiss > 0L) {
  warning("NAs produced in ", nmiss, " voxels for the positive z-statistic TFCE. These will be set to zero.")
  ptfce_pos$Z@.Data[is.na(ptfce_pos$Z@.Data)] <- 0
}

if (isTRUE(two_sided)) {
  z_neg <- -1*Z # the literal interpretation of Smith and Nichols is to use the whole image
  #z_neg@.Data[z_neg@.Data > 0] <- 0 # zero positive values
  #z_neg@.Data[z_neg@.Data < 0] <- z_neg@.Data[z_neg@.Data < 0] * -1 # invert negative stats
  call_list[["img"]] <- z_neg
  ptfce_neg <- ptfce_worker(call_list)
  ptfce_neg$Z@.Data[ptfce_neg$Z@.Data < 0] <- 0

  nmiss <- sum(is.na(ptfce_neg$Z))
  if (nmiss > 0L) {
    warning("NAs produced in ", nmiss, " voxels for the negative z-statistic TFCE. These will be set to zero.")
    ptfce_neg$Z@.Data[is.na(ptfce_neg$Z@.Data)] <- 0
  }

  ptfce_neg$Z <- -1 * ptfce_neg$Z # make negative z-stats negative again

  # create combined pTFCE-corrected image with positive and negative results
  # ptfce_obj@.Data[ptfce_neg$Z > 0] <- -1 * ptfce_neg$Z[ptfce_neg$Z > 0] # fill in negative z-stats
  ptfce_obj <- ptfce_pos$Z + ptfce_neg$Z
} else {
  ptfce_obj <- ptfce_pos$Z # only positive values
}

p_list <- list()

for (pv in seq_along(fwe_p)) {
  this_p <- ifelse(isTRUE(two_sided), fwe_p[pv] / 2, fwe_p[pv])
  z_thresh <- pTFCE::fwe.p2z(ptfce_pos$number_of_resels, FWEP = this_p)
  if (isTRUE(write_thresh_imgs)) {
    zi <- ptfce_obj
    zi[abs(zi) < z_thresh] <- 0 # zero out below-threshold values
    writeNIfTI(zi, file.path(z_dir, paste0(base, "_ptfce_fwep_", round(fwe_p[pv], 3))))
  }
  p_list[[pv]] <- data.frame(
    z_ptfce = round(z_thresh, 4), p_value = fwe_p[pv], two_sided = two_sided,
    number_of_resels = round(ptfce_pos$number_of_resels, 4)
  )
}

p_df <- dplyr::bind_rows(p_list)

# write overall pTFCE image here
writeNIfTI(ptfce_obj, file.path(z_dir, paste0(base, "_ptfce")))
write.csv(p_df, file.path(z_dir, paste0(base, "_ptfce_zthresh.csv")), row.names=FALSE)

cat("pTFCE completed: ", as.character(Sys.time()), "\n")

# if we want to build out the fsleyes call to look at the image.
#echo "fsleyes  $FSLDIR/data/standard/MNI152_T1_1mm_brain.nii.gz -cm greyscale pTFCE_Z_$BASE  -dr `cat thres_z_$BASE.txt` `fslstats pTFCE_Z_$BASE -p 100` -cm red-yellow &"