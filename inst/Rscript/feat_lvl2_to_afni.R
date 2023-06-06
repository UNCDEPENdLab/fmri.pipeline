#!/usr/bin/env Rscript

##script to pull together FSL FEAT level 2 runs into an AFNI BRIK+HEAD format for review
##Note: this is largely a simple wrapper around the first-level conversion script
##For L2 analyses, the individual cope*.feat directories have the L1 structure and can be digested accordingly

print_help <- function() {
  cat("feat_lvl2_to_afni.R is a script that converts a .gfeat directory from FSL FEAT to an AFNI-compatible BRIK and HEAD file.",
    "The script labels all of the sub-briks of the output using the contrast names embedded in FSL's output structure.",
    "Altogether, the script is intended to make viewing various FSL outputs much easier and self-describing than",
    "the default structure, which places many NIfTI images in various directories with names like 'cope1.nii.gz'.",
    "",
    "Required inputs are:",
    "  --gfeat_dir <directory name>: The .gfeat directory to be parsed into an AFNI file.",
    "",
    "Optional inputs are:",
    "  --stat_outfile: The filename prefix (no suffix) for the AFNI stats files output by this script. Default: gfeat_stats",
    "  --no_varcope: Do not include the varcope (variance) estimate for each cope in the output.",
    "  --no_subjstats: Do not copy the between-subjects cope and varcope images from the .gfeat as individual NIfTIs.",
    "",
    sep="\n")
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
if (!requireNamespace("fmri.pipeline", quietly=TRUE)) stop("Must install fmri.pipeline R package for this script to work")

if (length(args) == 0L) {
  message("feat_lvl2_to_afni.R expects a single .gfeat directory from a level 2 analysis --gfeat_dir <directory>.\n")
  print_help()
  quit(save="no", 1, FALSE)
}

outfilename <- "gfeat_stats"
auxfilename <- "gfeat_aux"

argpos <- 1
output_varcope <- TRUE #whether to include group varcope (variance estimate) in the statistics output
output_subjstats <- TRUE

while (argpos <= length(args)) {
  if (args[argpos] == "--gfeat_dir") {
    gfeatdir <- args[argpos + 1] #name of preprocessed fMRI data
    argpos <- argpos + 2
  } else if (args[argpos] == "--help") {
    print_help()
    quit(save="no", 0, FALSE)
  } else if (args[argpos] == "--aux_outfile") {
    auxfilename <- args[argpos + 1] #not used at present
    argpos <- argpos + 2
  } else if (args[argpos] == "--no_varcope") {
    output_varcope <- FALSE
    argpos <- argpos + 1
  } else if (args[argpos] == "--no_subjstats") {
    output_subjstats <- FALSE
    argpos <- argpos + 1
  } else if (args[argpos] == "--stat_outfile") { #name of output file for main stats file
    outfilename <- args[argpos + 1]
    argpos <- argpos + 2
  } else {
    stop("Not sure what to do with argument: ", args[argpos])
  }
}

if (!file.exists(gfeatdir)) {
  stop("Unable to locate --gfeatdir:", gfeatdir, ". Did you pass in the absolute path?")
} else {
  gfeatdir <- normalizePath(gfeatdir) #convert to absolute so that commands below work as expected
}

#setwd(gfeatdir)

##find cope directories
copedirs <- grep("/cope[0-9]+\\.feat", list.dirs(path=gfeatdir, full.names=TRUE, recursive=FALSE), value=TRUE, perl=TRUE)

#sort copes in numeric order to match user expectation
copenums <- as.numeric(sub(".*/cope(\\d+).feat", "\\1", copedirs, perl=TRUE))
copedirs <- copedirs[order(copenums)]

cat("Located the following cope directories:\n")
print(copedirs)
cat("\n\n")

copeafni <- c()
l1_script <- system.file("Rscript/feat_lvl1_to_afni.R", package = "fmri.pipeline")
for (d in seq_along(copedirs)) {
  ##run the L1 -> AFNI conversion for each separate cope
  cat("Processing: ", copedirs[d], "\n\n")
  vc_suffix <- ifelse(output_varcope, "", "--no_varcope")
  system(paste(l1_script, "--no_auxstats --feat_dir", copedirs[d], vc_suffix))
  copename <- readLines(file.path(copedirs[d], "design.lev")) #contains the L2 effect name (e.g., clock_onset)
  copename <- gsub("\\s", "_", copename, perl=TRUE) #convert spaces to underscores for accurate sub-brik labels

  afniout <- file.path(copedirs[d], "feat_stats+tlrc")
  briklabels <- fmri.pipeline::run_afni_command(paste("3dinfo -label", afniout), intern=TRUE)
  briklabels <- paste(copename, strsplit(briklabels, "|", fixed=TRUE)[[1]], sep="_", collapse=" ")

  ##need to add prefix for each cope to make the stats unique
  retcode <- fmri.pipeline::run_afni_command(paste0("3drefit -relabel_all_str '", briklabels, "' ", afniout))

  ##for now, eliminate the aux file (now handled by --no_auxstats above)
  #system(paste("rm", file.path(copedirs[d], "feat_aux+tlrc*")))

  if (output_subjstats) {
    ##filtered_func_data contains the cope from the lower level. Copy to output directory and rename
    retcode <- fmri.pipeline::run_afni_command(paste("3dcopy -overwrite", file.path(copedirs[d], "filtered_func_data.nii.gz"), paste0(outfilename, "_", copename, "_cope.nii.gz")))
    retcode <- fmri.pipeline::run_afni_command(paste("3dcopy -overwrite", file.path(copedirs[d], "var_filtered_func_data.nii.gz"), paste0(outfilename, "_", copename, "_varcope.nii.gz")))
    #retcode <- fmri.pipeline::run_afni_command(paste("3dcopy -overwrite", file.path(copedirs[d], "tdof_filtered_func_data.nii.gz"), paste0(outfilename, "_", copename, "_tdof.nii.gz"))) #have no use for this at the moment
  }

  copeafni <- c(copeafni, afniout)
}

#glue together the stats files
retcode <- fmri.pipeline::run_afni_command(paste("3dTcat -overwrite -prefix", outfilename, paste(copeafni, collapse=" ")))

if (file.exists(paste0(outfilename, "+tlrc.BRIK"))) {
  system(paste0("gzip ", outfilename, "+tlrc.BRIK"))
}

#cleanup ingredients of individual cope aggregation
system(paste("rm", paste0(copeafni, "*", collapse=" ")))
