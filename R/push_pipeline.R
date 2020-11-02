#push each model variant through the full pipeline
#eventually, when we have l3 covariates, should probably allow l1 x l3 combinations

#LOOP OVER MODEL VARIANTS IN PARALLEL

push_pipeline <- function(fsl_model_arguments, ncpus=1) { #should this read from the fsl list, which also has ncpus??
  #this helper script walks through each variant of the level1 model and calls scripts to run the full analysis pipeline
  #this amounts to iterating over each of the l1_model_variants
  stopifnot(length(fsl_model_arguments$outdir) == length(fsl_model_arguments$l1_model_variants))
  stopifnot(is.numeric(ncpus) && ncpus >= 1)
  
  require(parallel)
  require(doParallel)
  require(dependlab) #has qsub_file and wait_for_job

  source(file.path(fsl_model_arguments$pipeline_home, "functions", "glm_helper_functions.R"))
  source(file.path(fsl_model_arguments$pipeline_home, "functions", "run_feat_lvl1_sepqsub.R")) #executes FSF files in parallel batches

  setwd(fsl_model_arguments$pipeline_home) #to make qsub_file calls below happy with local paths
  
  #setup parallel worker pool, if requested
  if (ncpus > 1) {
    cl <- makeCluster(ncpus)
    registerDoParallel(cl)    
    on.exit(try(stopCluster(cl))) #cleanup pool upon exit of this function
  } else {
    registerDoSEQ() #formally register a sequential 'pool' so that dopar is okay
  }

  #iterate over lvl1/run variants
  nothing <- foreach(run_model_index=iter(1:length(fsl_model_arguments$outdir)), .packages="dependlab", .export=c("run_feat_lvl1_sepqsub")) %dopar% {

    #setup FSF files for all runs in one qsub

    #execute_fsl_lvl1_pipeline.R is just a thin qsub wrapper around model_clock_fmri_lvl1 to aid in parallelism
    setup_fsf_jobid <- qsub_file(script="qsub_fmri_r_script.bash",
      pbs_args=c("-l nodes=1:ppn=18", "-l walltime=10:00:00"),
      env_variables=c(R_SCRIPT="execute_fsl_lvl1_pipeline.R",
        run_model_index=run_model_index,
        fsl_pipeline_file=file.path(fsl_model_arguments$pipeline_home, "configuration_files", paste0(fsl_model_arguments$analysis_name, ".RData")))
    )

    # Okay, here's the qsub dilemma: We need for model_clock_fmri_lvl1 to complete before we try to run
    # run_feat_lvl1_sepqsub. We could submit this as a job that depends on setup_fsf_jobid above
    # But! We need all of the subsidiary qsub jobs that are spawned by run_feat_lvl1_sepqsub to be available/defined
    # when we call execute_fsl_lvl2_pipeline.R. So, if we queue a qsub job to run_feat_lvl1_sepqsub that depends
    # on setup_fsf_jobid, the code will continue to the lvl2 pipeline before the multiple lvl1 pipeline jobids
    # are even defined. Thus, my provisional, if slightly hacky solution is to wait here in the code until the
    # model_clock_fmri_lvl1 job completes so that we can execute run_feat_lvl1_sepqsub directly as a function
    # and recover the separate jobids that need to be passed to the level 2 script.

    wait_for_job(setup_fsf_jobid) #pause R script here until model_clock_fmri_lvl1 above completes

    #then, call run_feat_lvl1_sepqsub.R to execute all .fsf files using separate qsub instances
    
    #do not pass a wait_for here, as it will lead to an 'invalid job dependency' error during qsub (i.e., the parent job is already complete!).
    #this occurs because we now pause in the execution while the model_clock_fmri_lvl1 finishes in the previous step
    sep_lvl1_jobs <- run_feat_lvl1_sepqsub(fsl_model_arguments, run_model_index, rerun=FALSE) #, wait_for=setup_fsf_jobid)

    #bring the LVL1 sep qsub out to here.
    #Create batched qsubs for all FSL LVL1 analyses after FSF setup completes
    ##Sys.setenv(TARGET=fsl_model_arguments$fmri_dir, MODEL_MATCH=fsl_model_arguments$outdir) #, WAIT_FOR=setup_fsf_jobid)
    
    ## lvl1_setup_jobid <- qsub_file(script="qsub_fmri_r_script.bash",
    ##   pbs_args=c("-l nodes=1:ppn=20", "-l walltime=10:00:00")
    ##   env_variables=c("R_SCRIPT=run_feat_lvl1_sepqsub.R",
    ##     paste0("WAIT_FOR=", setup_fsf_jobid),
    ##     paste0("run_model_index=", run_model_index),
    ##     paste0("fsl_pipeline_file=configuration_files/", fsl_model_arguments$analysis_name, ".RData"))
    ## )
    
    #NB. run_feat_lvl1_sepqsub.R will create a file called sepqsub_lvl1_jobs.txt in the relevant temporary directory in the scratch folder
    
    # Run lvl2 analyses for the current model. The execute_fsl_lvl2_pipeline.R script will look for lvl1 jobs that are still running.
    # In terms of dependencies, because the qsub above fires immediately, the lvl2 needs to wait on:
    #    a) the job id from execute_fsl_lvl1_pipeline.bash
    #    b) all subsidiary job ids from each batch created by run_feat_lvl1_sepqsub.R
    #
    # Notably, the subsidiary job ids for the separate qsubs should be available/spawned as part of the parent job (execute_fsl_lvl1_pipeline.bash),
    # so waiting on both sets for LVL 2 should be a safe approach.

    #the problem is that run_feat_lvl1_sepqsub.R will be queued, but the subidiary jobs have not been created yet. Thus, the check below for separate qsub jobs will
    #probably be invalid.
    
    #this is now handled as a function return above
    #check in on the separate l1 jobs to wait on
    ## if (file.exists(l1_jobfile <- file.path(fsl_model_arguments$workdir[m], "sepqsub_lvl1_jobs.txt"))) {
    ##   sep_lvl1_jobs <- readLines(l1_jobfile)
    ## } else {
    ##   sep_lvl1_jobs <- NULL
    ## }

    #Level 2 analyses must wait on the completion of all parallel Level 1 (run-level) jobs
    if (!is.null(sep_lvl1_jobs)) {
      wait_string <- paste0("-W depend:afterok:", paste(sep_lvl1_jobs, collapse=":")) #setup_fsf_jobid, 
    } else { wait_string <- NULL }
    
    l2_execution_jobid <- qsub_file(script="qsub_fmri_r_script.bash",
      pbs_args=c("-l nodes=1:ppn=20", "-l walltime=10:00:00", wait_string),
      env_variables=c(
        R_SCRIPT="execute_fsl_lvl2_pipeline.R",
        run_model_index=run_model_index,
        fsl_pipeline_file=file.path(fsl_model_arguments$pipeline_home, "configuration_files", paste0(fsl_model_arguments$analysis_name, ".RData")))
    )

    if (!is.null(l2_execution_jobid)) {
      wait_string <- paste0("-W depend:afterok:", l2_execution_jobid)
    } else { wait_string <- NULL }

    #only need 1 CPU for l3 because it mostly sets up the .fsf files, then qsubs them.
    l3_execution_jobid <- qsub_file(script="qsub_fmri_r_script.bash",
      pbs_args=c("-l nodes=1:ppn=1", "-l walltime=2:00:00", wait_string),
      env_variables=c(
        R_SCRIPT="execute_fsl_lvl3_pipeline.R",
        run_model_index=run_model_index,
        fsl_pipeline_file=file.path(fsl_model_arguments$pipeline_home, "configuration_files", paste0(fsl_model_arguments$analysis_name, ".RData")))
    )


    return(NA_character_) #explicit loop return
  }

}
