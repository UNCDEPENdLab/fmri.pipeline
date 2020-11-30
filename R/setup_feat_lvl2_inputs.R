#This function generates the inputs for an FSL level 2 analysis, where multiple runs for a subject are combined using
#fixed effects estimation. The latter is currently achieved using run_lvl2...

# The steps for this script are:
# 1) Identify .feat directories for the specified configuration
# 2) Examine head motion statistics and exclude runs that exceed specified thresholds
# 3) Load emotion conditions from first-level design matrices
# 4) Save all valid inputs into a single data.frame .RData object to be digested/run

setup_feat_lvl2_inputs <- function(fsl_model_arguments, run_model_index) {
  require(plyr)
 
  #define run-level model folder name for analysis
  expectdir <- fsl_model_arguments$expectdir #the subfolder for each subject such as mni_5mm_aroma
  odir <- fsl_model_arguments$outdir[run_model_index] #the expect subfolder name for outputs of this model
  n_l1_copes <- fsl_model_arguments$n_l1_copes[run_model_index]
  subject_data <- fsl_model_arguments$subject_data
  id_col <- fsl_model_arguments$id_col
  
  #setup inputs for all LVL2 analyses for the current model
  feat_runs <- do.call(rbind, lapply(1:nrow(subject_data), function(s) {
    rr <- system(paste0("find ", file.path(subject_data$mr_dir[s], expectdir), " -mindepth 2 -iname \"FEAT_LVL1_run*.feat\" -ipath \"*/", odir, "/*\" -type d"), intern=TRUE)
    if (length(rr) < 1L) { return(NULL) } # subjects with no results will fall out in the rbind
    dd <- data.frame(subid=subject_data[[id_col]][s], mr_dir=subject_data$mr_dir[s], feat_dir=rr, run_num=as.integer(sub(".*FEAT_LVL1_run(\\d+)\\.feat.*", "\\1", rr, perl=TRUE)), stringsAsFactors=FALSE)
    return(dd)
  }))
  
  cat("All runs identified:\n")
  print(feat_runs$feat_dir)

  #flag runs with more than 15% volumes with FD 0.9mm or greater
  #find fd.txt files corresponding to each FEAT run

  feat_runs$fd_file <- file.path(feat_runs$mr_dir, expectdir, paste0("clock", feat_runs$run_num), "motion_info", "fd.txt")

  #identify the length of runs used in the analysis, accounting for both initial dropped volumes and truncation at the end
  feat_runs$drop_volumes <- fsl_model_arguments$drop_volumes #replicate initial drops for each .feat directory (used to index FD files)

  feat_runs$trunc_lengths <- unname(sapply(feat_runs$feat_dir, function(feat_dir) {
    design_fsf <- readLines(file.path(feat_dir, "design.fsf"))
    nvol <- as.numeric(sub("set fmri(npts)", "", grep("set fmri(npts)", design_fsf, fixed=TRUE, value=TRUE)[1], fixed=TRUE))
    nvol <- nvol + fsl_model_arguments$drop_volumes  #in FD indexing, this is the end point of the vector
  }))

  if (any(is.na(feat_runs$trunc_lengths))) {
    cat("Cannot read volumes from L1 design.fsf\nFiles affected:\n")
    print(feat_runs %>% filter(is.na(trunc_lengths)))
  }
  
  cat("Excluding runs exceeding 10% frames with FD >= 0.9mm OR any movement > 5mm\n")
  #read each FD file, index it based on modeled fMRI data, then return FD statistics
  motexclude <- plyr::ldply(1:nrow(feat_runs), function(i) {
    fd <- read.table(feat_runs$fd_file[i], header=FALSE)$V1
    fd <- fd[(feat_runs$drop_volumes[i]+1):feat_runs$trunc_lengths[i]] #only include volumes within run
    propSpikes_0p9 <- sum(as.integer(fd > 0.9))/length(fd)
    spikeExclude <- if (propSpikes_0p9 > .10) 1 else 0
    maxFD <- max(fd)
    meanFD <- mean(fd)
    maxMotExclude <- if (maxFD > 5) 1 else 0
    data.frame(fd_file=feat_runs$fd_file[i], propSpikes_0p9, spikeExclude, meanFD, maxFD, maxMotExclude, stringsAsFactors=FALSE)
  })

  feat_runs <- inner_join(feat_runs, motexclude, by="fd_file")

  #drop subject altogether if fewer than 4 analyzeable runs
  feat_runs <- plyr::ddply(feat_runs, .(subid), function(subdf) {
    if (nrow(subset(subdf, maxMotExclude == 0 & spikeExclude == 0)) < 4) {
      subdf$lt4runs <- 1
    } else {
      subdf$lt4runs <- 0
    }
    subdf
  })

  feat_runs$anyExclude <- with(feat_runs, as.integer(spikeExclude | maxMotExclude | lt4runs))

  #2018: leaving this as comment for now, but not enforcing since it violates the algorithmic approach
  #10637 has pretty bad head movement in runs 5-8... in runs 7 and 8, it falls just below 10% FD > 0.9mm, so exclude subject altogether
  #feat_runs[which(feat_runs$subid == "10637"),"anyExclude"] <- 1

  nrow(feat_runs[which(feat_runs$anyExclude == 1),])
  badruns <- droplevels(subset(feat_runs, anyExclude==1))

  #table(badruns$subid)
  #feat_runs[which(feat_runs$subid == "10711"),] #runs 6, 7, 8 are bad
  #feat_runs[which(feat_runs$subid == "11324"),] #a lot of movement in runs 3 and 4, but otherwise very still...
  #feat_runs[which(feat_runs$subid == "11336"),] #run 1 has a 14.5 mm FD, run 4 has an 8.5mm movement, but otherwise still

  #generate data frame of runs to analyze

  feat_l2_inputs_df <- droplevels(feat_runs[which(feat_runs$anyExclude==0),]) #only retain good runs

  #figure out emotion and rew contingency for all runs
  run_conditions <- do.call(rbind, lapply(1:nrow(feat_l2_inputs_df), function(i) {
    designmat <- file.path(dirname(feat_l2_inputs_df[i,"feat_dir"]), "designmatrix.RData")
    loc <- local({load(designmat); environment()})$subj_data #load subj_data data.frame from designmatrix.RData
    head(subset(loc, run==feat_l2_inputs_df$run_num[i], select=c(emotion, rewFunc)), n=1) #just get emotion and contingency as a single-row data.frame
  }))

  names(run_conditions) <- c("emotion", "contingency") #rename 'rewFunc' -> 'contingency'
  feat_l2_inputs_df <- cbind(feat_l2_inputs_df, run_conditions)
  feat_l2_inputs_df$emotion <- relevel(feat_l2_inputs_df$emotion, ref="scram")
  feat_l2_inputs_df$model <- odir
  feat_l2_inputs_df$n_l1_copes <- n_l1_copes #number of level 1 copes to propagate to analyze/combine at L2 (E.g., clock, feedback, and pe)  
  feat_l2_inputs_df <- droplevels(arrange(feat_l2_inputs_df, subid, run_num)) #arrange by subid, run
  #feat_l2_inputs_df[sample(1:nrow(feat_l2_inputs_df), 30), ] #verify match across columns

  feat_l2_inputs_df <- feat_l2_inputs_df %>% select(subid, run_num, contingency, emotion, model, feat_dir, everything())

  save(feat_l2_inputs_df, feat_runs, file=file.path(fsl_model_arguments$pipeline_home, "configuration_files", paste0(paste(fsl_model_arguments$analysis_name, odir, "lvl2_inputs", sep="_"), ".RData")))

  return(feat_l2_inputs_df) #return run-level information for passing onto run_feat_lvl2
}
