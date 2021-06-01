#this version of the FSL LVL1 feat estimation creates multiple qsub scripts in a temporary directory
##where each script has a number of feat calls that are forked, then the script waits for completion
##This circumvents the ICS limit on multiple nodes in R using a SNOW cluster.
##The primary calculation is how many files there are to run relative to processors and files per processor (chunking)

## to_run <- Sys.getenv("fsl_pipeline_file")

## run_model_index <- as.numeric(Sys.getenv("run_model_index")) #which variant to execute
## if (nchar(to_run) == 0L) { stop("Cannot locate environment variable fsl_pipeline_file") }
## if (!file.exists(to_run)) { stop("Cannot locate configuration file", to_run) }
## if (is.na(run_model_index)) { stop("Couldn't identify usable run_model_index variable.") }

## load(to_run)
## wait_for <- Sys.getenv("WAIT_FOR")

## #call function below
## run_feat_lvl1_sepqsub(feat_model_arguments, run_model_index, rerun=FALSE, wait_for=wait_for)

run_feat_lvl1_sepqsub <- function(gpa, model_name=NULL, rerun=FALSE, wait_for="") {

  #This function is now called within run_fsl_pipeline, rather than being run in its own qsub
  #to_run <- Sys.getenv("fsl_pipeline_file")
  #run_model_index <- as.numeric(Sys.getenv("run_model_index")) #which variant to execute
  #if (nchar(to_run) == 0L) { stop("Cannot locate environment variable fsl_pipeline_file") }
  #if (!file.exists(to_run)) { stop("Cannot locate configuration file: ", to_run) }
  #if (is.na(run_model_index)) { stop("Couldn't identify usable run_model_index variable.") }

  #load(to_run)

  subject_data <- gpa$subject_data
  expectdir <- gpa$expectdir
  model_match <- gpa$outdir[run_model_index]
  l1_working_directory <- gpa$l1_working_directory[run_model_index]

  cpusperjob <- 8 #number of cpus per qsub
  runsperproc <- 3 #number of feat calls per processor

  #look in the subfolder for each subject for fsf files
  fsf_files <- do.call(c, lapply(subject_data$mr_dir, function(s) {
    system(paste0(
      "find ", file.path(s, expectdir),
      " -mindepth 2 -iname \"FEAT_LVL1_*.fsf\" -ipath \"*/", model_match, "/*\" -type f"
    ), intern = TRUE)
  }))

  #figure out which fsf files have already been run
  dir_expect <- gsub("\\.fsf$", ".feat", fsf_files, perl=TRUE)
  torun <- c()

  for (f in seq_len(fsf_files)) {
    if (file.exists(dir_expect[f])) {
      if (rerun) {
        cmd <- paste0("rm -rf \"", dir_expect[f], "\"")
        cat("Removing old directory: ", cmd, "\n")
        system(cmd)
        torun <- c(torun, fsf_files[f]) #add to queue
      } else {
        cat("Skipping existing directory: ", dir_expect[f], "\n")
      }
    } else {
      torun <- c(torun, fsf_files[f])
    }
  }

  if (length(torun) == 0) {
    cat("No LVL1 .fsf files to execute.\n\n")
    return(NULL)
  }

  cat("About to run the following fsf files in parallel:\n\n")
  cat(torun, sep="\n")

  preamble <- c(
    "#PBS -A mnh5174_c_g_sc_default",
    paste0("#PBS -l nodes=1:ppn=", cpusperjob),
    paste0("#PBS -l pmem=8gb"),
    ifelse(wait_for != "", paste0("#PBS -W depend=afterok:", wait_for), ""), #allow job dependency on upstream setup
    "#PBS -l walltime=14:00:00",
    "#PBS -j oe",
    "#PBS -m n",
    "#PBS -W group_list=mnh5174_collab",
    "",
    "",
    "source /gpfs/group/mnh5174/default/lab_resources/ni_path.bash",
    "module unload fsl", #make sure that the ni_path version of FSL is unloaded
    "#module load \"openblas/0.2.20\" >/dev/null 2>&1",
    "module load \"fsl/6.0.1\" >/dev/null 2>&1",
    "module load gsl/2.5", #for dependlab R package to work (some new dependency)
    "",
    "cd $PBS_O_WORKDIR"
  )


  if (!file.exists(l1_working_directory)) { dir.create(l1_working_directory, recursive=TRUE) }

  njobs <- ceiling(length(torun)/(cpusperjob*runsperproc))

  #use length.out on rep to ensure that the vectors align even if chunks are uneven wrt files to run
  df <- data.frame(fsf=torun, job=rep(1:njobs, each=cpusperjob*runsperproc, length.out=length(torun)), stringsAsFactors=FALSE)
  df <- df[order(df$job), ]

  joblist <- rep(NA_character_, njobs)
  for (j in 1:njobs) {
    outfile <- paste0(l1_working_directory, "/qsub_featsep_", j, "_", basename(tempfile()), ".pbs")
    cat(preamble, file=outfile, sep="\n")
    thisrun <- with(df, fsf[job==j])
    cat(paste("feat", thisrun, "&"), file=outfile, sep="\n", append=TRUE)
    cat("wait\n\n", file=outfile, append=TRUE)
    cat(paste0("bash /gpfs/group/mnh5174/default/clock_analysis/fmri/gen_feat_reg_dir.bash ", unique(dirname(thisrun))), sep="\n", file=outfile, append=TRUE)
    joblist[j] <- qsub_file(outfile) #system(paste0("qsub ", outfile))
  }

  #write the list of separate feat qsub jobs that are now queued (so that LVL2 can wait on these)
  #should also return this to the caller as a function?
  writeLines(joblist, con=file.path(l1_working_directory, "sepqsub_lvl1_jobs.txt"))

  return(joblist)
}
