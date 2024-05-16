spm_sceptic_model <- function(subj_data, sceptic_signals, mrfiles, runlengths, mrrunnums, execute=FALSE, force=FALSE,
                              drop_volumes=0, outdir=NULL, usepreconvolve=FALSE, spikeregressors=FALSE, model_suffix="", ...) {

    
  if (is.null(outdir)) {
    outdir=paste0("sceptic-", paste(sceptic_signals, collapse="-")) #define output directory based on combination of signals requested
    if (usepreconvolve) { outdir=paste(outdir, "preconvolve", sep="-") }

    outdir <- paste0(outdir, model_suffix) #add any model suffix, if requested
  }

  #note: normalizePath will fail to evaluate properly if directory does not exist
  spm_run_output_dir <- file.path(normalizePath(file.path(dirname(mrfiles[1L]), "..")), outdir)

  if (dir.exists(file.path(spm_run_output_dir, "SPM.mat")) && force==FALSE) { message(file.path(spm_run_output_dir, "SPM.mat"), " exists. Skipping."); return(0) }
  cat("spm_run_output_dir create: ", spm_run_output_dir, "\n")
  dir.create(spm_run_output_dir, showWarnings=FALSE) #one directory up from a given clock run
  timingdir <- file.path(spm_run_output_dir, "run_timing_sceptic")

  #create the events and signals structures for the build_design_matrix call
  #this assumes that we have received a data.frame with the structure from the _trial_statistics.csv.gz generated
  #thus, we have $clock_onset and $feedback_onset, and $iti_onset available
  events <- subj_data %>% dplyr::select(id, run, trial, clock_onset, feedback_onset, iti_onset, rt_csv) %>%
    dplyr::mutate(clock_duration=rt_csv/1000, feedback_duration=iti_onset - feedback_onset) %>%
    dplyr::select(-iti_onset, -rt_csv) %>% tidyr::gather(key="key", value="value", -id, -run, -trial) %>%
    tidyr::separate(col = key, into = c("event", "onset_duration")) %>%
    tidyr::spread(key=onset_duration, value=value) %>% dplyr::select(event, run, trial, onset, duration)

  signals <- populate_sceptic_signals(sceptic_signals, subj_data)

  #NB. The tr argument should be passed in as part of ...
  d <- build_design_matrix(events=events, signals=signals, baseline_coef_order=2, write_timing_files=NULL, #since generate_spm_mat uses the d object itself
    center_values=TRUE, plot=FALSE, convolve_wi_run=TRUE, drop_volumes=drop_volumes,
    run_volumes=mrfiles, runs_to_output=mrrunnums, ...)

  all_nuisance <- list()
  for (r in 1:length(mrfiles)) {
    stopifnot(file.exists(file.path(dirname(mrfiles[r]), "motion.par"))) #can't find motion parameters
    
    runnum <- sub("^.*/clock(\\d+)$", "\\1", dirname(mrfiles[r]), perl=TRUE)
    mregressors <- NULL #start with NULL
    
    ##add CSF and WM regressors (with their derivatives)
    nuisancefile <- file.path(dirname(mrfiles[r]), "nuisance_regressors.txt")
    if (file.exists(nuisancefile)) {
      nuisance <- read.table(nuisancefile, header=FALSE)
      nuisance <- nuisance[(1+drop_volumes):runlengths[r],,drop=FALSE]
      nuisance <- as.data.frame(lapply(nuisance, function(col) { col - mean(col) })) #demean
      if (!is.null(mregressors)) { mregressors <- cbind(mregressors, nuisance) #note that in R 3.3.0, cbind with NULL or c() is no problem...
      } else { mregressors <- nuisance }
    }

    all_nuisance[[r]] <- mregressors
  }

  #put everything into one big data.frame for a concatenated run
  nuisance_concat <- do.call(rbind, all_nuisance)

  save(d, subj_data, events, signals, timingdir, runlengths, mrrunnums, all_nuisance, nuisance_concat, file=file.path(spm_run_output_dir, "designmatrix.RData"))

  spm_syntax <- generate_spm_mat(d, output_dir=spm_run_output_dir, nifti_tmpdir=spm_run_output_dir, concatenate_runs=TRUE,
    condition_contrasts=TRUE, effects_of_interest_F=TRUE, generate_qsub=TRUE, execute_qsub=execute)

}
