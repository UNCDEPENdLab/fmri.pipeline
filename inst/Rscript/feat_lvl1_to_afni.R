#!/usr/bin/env Rscript

#script to pull together FSL FEAT level 1 runs into an AFNI BRIK+HEAD format for review
print_help <- function() {
  cat(paste("This script converts a first-level feat directory to an AFNI-compatible stats file with",
    "the afni graphical viewer. This avoids the challenge of unlabeled and unintuitive individual",
    "outputs from feat. Here are the basic options:",
    "  --feat_dir <dir.feat>: the location of the feat output.",
    "  --help: print the help",
    "  --no_varcope: do not include the varcope (variance) estimates in the output.",
    "  --no_auxstats: do not output a more detailed file with auxiliary statistics, including zfstats, sigmasquareds, and pes.",
    "  --stat_outfile <file prefix>: The filename prefix to be used for stats outputs.",
    "  --aux_outfile <file prefix>: The filname prefix for auxiliary statistics.",
    "\n\n",
    sep="\n"))
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

# contains run_afni_command
if (!requireNamespace("fmri.pipeline", quietly = TRUE)) stop("Must install fmri.pipeline R package for this script to work")

if (length(args) == 0L) {
  message("feat_lvl1_to_afni.R expects a single .feat directory from a level 1 analysis --feat_dir <directory>.\n")
  print_help()
  quit(save="no", 1, FALSE)
}

outfilename <- "feat_stats"
auxfilename <- "feat_aux"
output_varcope <- TRUE
output_auxstats <- TRUE #thresh zstat f-stats etc.

argpos <- 1
while (argpos <= length(args)) {
  if (args[argpos] == "--feat_dir") {
    featdir <- args[argpos + 1] #name of preprocessed fMRI data
    stopifnot(file.exists(featdir))
    argpos <- argpos + 2
  } else if (args[argpos] == "--help") {
    print_help()
    quit(save="no", 0, FALSE)
  } else if (args[argpos] == "--no_varcope") {
    output_varcope <- FALSE
    argpos <- argpos + 1
  } else if (args[argpos] == "--no_auxstats") {
    output_auxstats <- FALSE
    argpos <- argpos + 1
  } else if (args[argpos] == "--stat_outfile") { #name of output file for main stats file
    outfilename <- args[argpos + 1]
    argpos <- argpos + 2
  } else if (args[argpos] == "--aux_outfile") {
    auxfilename <- args[argpos + 1]
    argpos <- argpos + 2
  } else {
    stop("Not sure what to do with argument: ", args[argpos])
  }
}

setwd(featdir)

#inside the stats directory we will have pes, copes, varcopes, and zstats
zfiles <- list.files("stats", pattern="zstat[0-9]+\\.nii.*", full.names=TRUE)
statnums <- as.numeric(sub(".*zstat(\\d+)\\.nii.*", "\\1", zfiles, perl=TRUE))
zfiles <- zfiles[order(statnums)] #need to sort stat files numerically for labeling to work appropriately
nstats <- length(zfiles)

##lookup names of parameter estimates

tcatcall <- paste("3dTcat -overwrite -prefix", outfilename)
zbriks <- c()
for (s in 1:nstats) {
  copefile <- sub("zstat", "cope", zfiles[s], fixed=TRUE)
  if (output_varcope) {
    varcopefile <- sub("zstat", "varcope", zfiles[s], fixed=TRUE)
    zbriks <- c(zbriks, s*3 - 2) #triplets with zstat in second position. subtract 1 because AFNI uses 0-based indexing and 1 to put zstat as second of triplet
    tcatcall <- paste(tcatcall, copefile, zfiles[s], varcopefile)
  } else {
    zbriks <- c(zbriks, s*2 - 1) #couplets with zstat in second position. subtract 1 because AFNI uses 0-based indexing
    tcatcall <- paste(tcatcall, copefile, zfiles[s])
  }
}

#concatenate stat images
fmri.pipeline::run_afni_command(tcatcall)

#design.con contains names of contrasts
dcon <- readLines("design.con")
connames <- sub("/ContrastName\\d+\\s+([\\w_.]+).*", "\\1", grep("/ContrastName", dcon, value=TRUE), perl=TRUE)
connames <- gsub("\\s", "_", connames, perl=TRUE) #replace spaces with underscores to make labels accurate in AFNI

#now rework AFNI header to correct the labels and add z-stat info
if (output_varcope) {
  briknames <- paste(rep(connames, each=3), c("coef", "z", "var"), sep="_", collapse=" ")
} else {
  briknames <- paste(rep(connames, each=2), c("coef", "z"), sep="_", collapse=" ")
}

#add this to zstat images
refitcall <- paste0("3drefit -fbuc ", paste("-substatpar", zbriks, "fizt", collapse=" "), " -relabel_all_str '", briknames, "' ", outfilename, "+tlrc")

fmri.pipeline::run_afni_command(refitcall)

if (output_auxstats) {

  ##read auxiliary files (PEs + error)
  zbriks <- c()
  pefiles <- list.files("stats", pattern="^pe.*\\.nii.*", full.names=TRUE)
  penums <- as.numeric(sub(".*pe(\\d+)\\.nii.*", "\\1", pefiles, perl=TRUE))
  pefiles <- pefiles[order(penums)]
  findex <- 0 #0-based indexing
  tcatcall <- paste("3dTcat -overwrite -prefix", auxfilename)
  for (file in seq_along(pefiles)) {
    tcatcall <- paste(tcatcall, pefiles[file])
    findex <- findex + 1
  }

  ##add in other files
  ##thresh_zstat files
  threshzfiles <- list.files(pattern="^thresh_zstat.*\\.nii.*", full.names=TRUE)
  if (length(threshzfiles) > 0L) {
    threshznums <- as.numeric(sub(".*thresh_zstat(\\d+)\\.nii.*", "\\1", threshzfiles, perl=TRUE))
    threshzfiles <- threshzfiles[order(threshznums)]
    for (file in seq_along(threshzfiles)) {
      tcatcall <- paste(tcatcall, threshzfiles[file])
      zbriks <- c(zbriks, findex)
      findex <- findex + 1
    }
  }

  zfstatfiles <- list.files("stats", pattern="zfstat.*\\.nii.*", full.names=TRUE)
  if (length(zfstatfiles) > 0L) {
    for (file in seq_along(zfstatfiles)) {
      tcatcall <- paste(tcatcall, zfstatfiles[file])
      zbriks <- c(zbriks, findex)
      findex <- findex + 1
    }
  }

  #add residuals
  if (file.exists("stats/sigmasquared.nii.gz")) {
    tcatcall <- paste(tcatcall, "stats/sigmasquareds.nii.gz")
    findex <- findex + 1
  }

  fmri.pipeline::run_afni_command(tcatcall)

  briknames <- paste(c(paste0("pe", sort(penums)), paste0("thresh_zstat", sort(threshznums)), paste0("zfstat", seq_along(zfstatfiles)), "sigmasquareds"), collapse=" ")

  #add this to zstat images
  refitcall <- paste0("3drefit -fbuc ", paste("-substatpar", zbriks, "fizt", collapse=" "), " -relabel_all_str '", briknames, "' ", auxfilename, "+tlrc")
  fmri.pipeline::run_afni_command(refitcall)
}
