## mm <- model.matrix(designmat)
## write.table(mm, file="fsl_LVL3_emo_design.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
## cat(feat_l2_inputs_df$feat_dir, sep="\n", file="feat_dirlist.txt")
## print(dimnames(mm)[2]) #column names (for FSL)

## library(lsmeans)

##generate per-subject second-level FE analyses to get contrasts of interest for group analysis
#saved file does not include relevel for ref of scram (hence copy from above -- redundant)

run_feat_lvl2 <- function(feat_l2_inputs_df, run=TRUE, force=FALSE, ncpus=8) {
  allFeatRuns <- list()
  require(plyr)
  require(parallel)
  require(dependlab)
  
  if (is.null(feat_l2_inputs_df$n_l1_copes)) { stop("run_feat_lvl2 requires n_l1_copes column to be present in feat_l2_inputs_df") }    
  stopifnot(length(unique(feat_l2_inputs_df$n_l1_copes)) == 1)
  ncopes <- feat_l2_inputs_df$n_l1_copes[1] #number of lvl1 copes to combine for this model
  
  d_ply(feat_l2_inputs_df, .(subid, model), function(subdf) {
    subdf <- subdf[order(subdf$run_num),] #verify that we have ascending runs
    subdf$y_dummy <- rnorm(nrow(subdf))
    subdf$run_num <- subdf$run_num - min(subdf$run_num) #have intercept represent activation on earliest modeled run (usually 1).
    n_runs <- nrow(subdf) #used for FSF
    lm_dummy <- lm(y_dummy ~ emotion + run_num, subdf)
    mm <- model.matrix(lm_dummy)

    library(emmeans)
    v <- emmeans(lm_dummy, list(pairwise ~ emotion))
    condmeans <- v[[1]]@linfct #condition means
    rownames(condmeans) <- summary(v[[1]])$emotion
    contrasts <- v[[2]]@linfct #emo diffs (pairwise)
    rownames(contrasts) <- sub(" - ", "_gt_", summary(v[[2]])$contrast, fixed=TRUE)

    #get coefficients for grand mean contrast
    #gm_base <- emmeans(lm_dummy, ~ 1) #equally weights emotions even in case of unbalanced design (which applies here)

    #lsm <- emmeans(lm_dummy, "emotion")
    #contrast(lsm, method="eff")@linfct #cell versus gm

    #unbalanced design, so need to set weights based on relative frequency
    props <- prop.table(table(subdf$emotion))

    #at present, we can't easily work around subjects that are missing an emotion altogether (11343 only in the MMClock sample)
    #this would require an alternative approach to propagating contrasts to LVL3 since some contrasts would be undefined in LVL2.
    if (any(props < .01)) {
      warning("Missing at least one emotion altogether in combining the LVL1 runs: ", subdf$subid[1])
      return(NULL)
    }

    #1 for intercept/reference, proportions for other cells, add mean run value from emmeans computation
    gm_coef <- c(`(Intercept)`=1, props[2:length(props)], run_num=condmeans[1,"run_num"])
    names(gm_coef)[2:3] <- colnames(contrasts)[2:3] #should be safe because contrasts in lm() and table() ordering follow factor levels

    #run effect contrast
    run_coef <- as.vector(emtrends(lm_dummy, ~1, v="run_num")@linfct)
    
    #generate overall contrast matrix
    cmat <- round(rbind(condmeans, overall=gm_coef, contrasts, run=run_coef), 5) #the emtrends value tends to be .99999 something; round for clarity

    #generate FSF contrast syntax for this setup
    contrast_syntax <- generate_fsf_contrast_syntax(cmat)
    
    #generate and run lvl2 for this subject
    #fsfTemplate <- readLines(file.path(getMainDir(), "clock_analysis", "fmri", "fsf_templates", "feat_lvl2_clock_template.fsf"))

    #template with run trend modeled
    fsfTemplate <- readLines(file.path(getMainDir(), "clock_analysis", "fmri", "fsf_templates", "feat_lvl2_clock_template_runtrend.fsf"))
    #depending on lower-level model (e.g., TC versus value, will have different number of copes to compute
    
    # Template nomenclature:
    # .NINPUTS.   : number of .feat directories to be combined. This populates the fields: npts and multiple
    # .OUTPUTDIR. : the feat output location
    
    thisTemplate <- fsfTemplate #copy general template for adaptation in this subject

    #need to determine number of copes (contrasts) at level 1, which depends on the model being fit
    #FSL usually reads this from the .feat directories itself, but for batch processing, better to insert into the FSF ourselves
    #Need to put this just after the high pass filter cutoff line for Feat to digest it happily

    thisTemplate <- c(thisTemplate,
      "# Number of lower-level copes feeding into higher-level analysis",
      paste0("set fmri(ncopeinputs) ", ncopes))

    #tell FSL to analyze all lower-level copes in LVL2
    for (n in 1:ncopes) {
      thisTemplate <- c(thisTemplate,
        paste0("# Use lower-level cope ", n, " for higher-level analysis"),
        paste0("set fmri(copeinput.", n, ") 1"), "")      
    }

    ##thisTemplate <- gsub(".OUTPUTDIR.", file.path(dirname(subdf$feat_dir[1L]), "FEAT_LVL2"), thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".OUTPUTDIR.", file.path(dirname(subdf$feat_dir[1L]), "FEAT_LVL2_runtrend"), thisTemplate, fixed=TRUE)

    thisTemplate <- gsub(".NINPUTS.", n_runs, thisTemplate, fixed=TRUE)

    #add .feat directories to analyze
    for (i in 1:nrow(subdf)) { 
      thisTemplate <- gsub(paste0(".INPUT", i, "."), subdf$feat_dir[i], thisTemplate, fixed=TRUE)
    }

    #write the values for the EVs (run conditions)
    #EV1 is intercept (scrambled is reference)
    #EV2 is emofear
    #EV3 is emohappy
    #EV4 is linear run trend
    for (i in seq_len(nrow(mm))) {
      for (j in 1:ncol(mm)) {
        thisTemplate <- gsub(paste0(".I", i, "EV", j, "."), mm[i,j], thisTemplate, fixed=TRUE)
      }
    }

    #handle cleanup for less than 8 runs
    runs_to_dump <- (1:8)[1:8 > nrow(subdf)]
    if (length(runs_to_dump) > 0) {
      #The EV lines for each of the irrelevant runs should be omitted
      badfields <- c(paste0("evg", runs_to_dump, ".", rep(1:4, each=length(runs_to_dump))),
        paste0(".INPUT", runs_to_dump, "."), #unused .feat inputs
        paste0("groupmem", runs_to_dump)) #group membership for unused input slots

      badlines <- grep(sub(".", "\\.", paste0("(", paste(badfields, collapse="|"), ")"), fixed=TRUE), thisTemplate)
      if (length(badlines) > 0) {
        badlines <- c(badlines - 1, badlines) #also drop the preceding comment lines
        thisTemplate <- thisTemplate[-1*badlines]
      }
    }

    #add contrast information
    thisTemplate <- c(thisTemplate, contrast_syntax)

    featOutDir <- file.path(dirname(subdf$feat_dir[1L]), "FEAT_LVL2_runtrend.gfeat")
    featFile <- file.path(dirname(subdf$feat_dir[1L]), "FEAT_LVL2_runtrend.fsf")
    if (file.exists(featOutDir) && force==FALSE) { return(NULL) } #skip re-creation of FSF and do not run below unless force==TRUE 
    cat(thisTemplate, file=featFile, sep="\n")      
    
    allFeatRuns[[featFile]] <<- featFile
  })

  print(allFeatRuns)
  
  if (run == TRUE) {
    cl_fork <- makeForkCluster(nnodes=ncpus)
    runfeat <- function(fsf) {
      runname <- basename(fsf)
      runFSLCommand(paste("feat", fsf), stdout=file.path(dirname(fsf), paste0("feat_stdout_", runname)), stderr=file.path(dirname(fsf), paste0("feat_stderr_", runname)))
      system(paste0("feat_lvl2_to_afni.R --gfeat_dir ", sub(".fsf", ".gfeat", fsf, fixed=TRUE), " --no_subjstats --no_varcope --stat_outfile ", sub(".fsf", "_gfeat_stats", fsf, fixed=TRUE))) #aggregate FEAT statistics into a single file
    }
    clusterApply(cl_fork, allFeatRuns, runfeat)
    stopCluster(cl_fork)
  } 
  
}
