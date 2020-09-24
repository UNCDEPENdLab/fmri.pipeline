#the goal of this script is to run an entire fmri analysis for SCEPTIC data, including levels 1-3 in FSL
library(dependlab)
library(foreach)
library(parallel)
library(doParallel)
library(plyr)
library(tidyverse)

scripts_dir <- "/proj/mnhallqlab/projects/clock_analysis/fmri/fsl_pipeline"
setwd(scripts_dir)

source(file.path(scripts_dir, "functions", "push_pipeline.R"))
source(file.path(scripts_dir, "functions", "get_mmy3_trial_df.R"))
source(file.path(scripts_dir, "functions", "finalize_pipeline_configuration.R"))

trial_df <- get_mmy3_trial_df(model="selective", groupfixed=TRUE)

subject_df <- read.table("/proj/mnhallqlab/projects/clock_analysis/fmri/data/mmy3_demographics.tsv", header=TRUE) %>%
  rename(id=lunaid, Age=age, Female=female, ScanDate=scandate) %>%
  mutate(mr_dir = paste0("/proj/mnhallqlab/studies/MMClock/MR_Proc/", id, "_", format((as.Date(ScanDate, format="%Y-%m-%d")), "%Y%m%d")), #convert to Date, then reformat YYYYMMDD
    I_Age = -1000*1/Age,
    I_Age_c = I_Age - mean(I_Age, na.rm=TRUE),
    Age_c = Age - mean(Age, na.rm=TRUE),
    Q_Age = Age_c^2,
    Q_Age_c = Q_Age - mean(Q_Age, na.rm=TRUE)
  )

#Jun2017: further ICAs on the MMClock data suggest a short steady-state problem. Drop 2 volumes for good measure.

#from 2017:
##results from Mean SCEPTIC regressor correlation.pdf indicate that regressors for vchosen, ventropy_decay_matlab, dauc, and pemax are
##reasonably uncorrelated. The worst is dauc with vchosen (mean r = -0.31), which makes sense that as learning progresses, chosen values
##are higher and there is less residue to decay. These 4 regressors are also of greatest theoretical interest

#Setup the global configuration for the full FSL pipeline
fsl_model_arguments <- list(
  #analysis_name="MMClock_aroma_preconvolve_fse",
  analysis_name="MMClock_aroma_preconvolve_fse_groupfixed",
  #analysis_name="MMClock_fixed_aroma_preconvolved",
  #analysis_name="MMClock_aroma_preconvolve_fse_groupfixed_unsmoothed",
  trial_statistics = trial_df,
  subject_covariates = subject_df,
  id_col = "id",
  fmri_dir = "/proj/mnhallqlab/studies/MMClock/MR_Proc",
  expectdir = "mni_5mm_aroma", #subfolder name for processed data
  expectfile = "nfaswuktm_clock[0-9]_5.nii.gz", #expected file name for processed clock data
  #expectdir = "mni_nosmooth_aroma",
  #expectfile = "nfawuktm_clock[0-9].nii.gz", #expected file name for processed clock data
  usepreconvolve=TRUE,
  ncpus=20,
  drop_volumes=2, #to handle steady state concerns
  tr=1.0, #seconds
  spikeregressors=FALSE, #don't include spike regressors in nuisance variables since we are using AROMA
  sceptic_run_variants=list(
      c("clock", "feedback", "u_chosen_quantile"),
      c("clock", "feedback", "v_chosen"), #individual regressors
      c("clock", "feedback", "v_entropy"), #clock-aligned
      c("clock", "feedback", "rt_vmax_change"),
      c("clock", "feedback", "pe_max")
#    c("clock", "feedback", "u_chosen"),
#    c("clock", "feedback", "u_chosen_sqrt"),
#    c("clock", "feedback", "u_chosen_change"),
#    c("clock", "feedback", "u_chosen_z"),
#    c("clock", "feedback", "run_trial", "u_chosen"),
#    c("clock", "feedback", "run_trial", "u_chosen_z"),
#    c("clock", "feedback", "run_trial", "u_chosen_change")    
#    c("clock", "feedback_bs")
#    c("clock_bs", "feedback")
#    c("clock", "feedback", "v_chosen", "v_entropy", "d_auc", "pe_max"), #all signals with entropy of weights
#    c("clock", "feedback", "v_chosen", "v_entropy_func", "d_auc", "pe_max"), #all signals with entropy of evaluated function
#    c("clock", "feedback", "v_entropy_feedback"), #feedback-aligned
#    c("clock", "feedback", "v_entropy_func"),
#    c("clock", "feedback", "d_auc"), #feedback-aligned
#    c("clock", "feedback", "d_auc_clock"), #clock-aligned
#    c("clock", "feedback", "v_entropy_no5"),
#    c("clock", "feedback", "v_auc"),
#    c("clock", "feedback", "d_auc_sqrt"),
#    c("clock", "feedback", "rt_swing"),
#    c("clock", "feedback", "rt_swing_sqrt"),
#    c("clock", "feedback", "v_max"),
#    c("clock", "feedback", "mean_kld"),
#    c("clock", "feedback", "intrinsic_discrepancy"),
#    c("clock", "feedback", "mean_kld_feedback"),
#    c("clock", "feedback", "intrinsic_discrepancy_feedback"),

#    c("clock", "feedback", "rt_vmax_change_dir"),
#    c("clock", "feedback", "v_trial_fixed"),
#    c("clock", "feedback", "pe_trial_fixed")    
#    c("clock", "feedback", "v_entropy_change"),
#    c("clock", "feedback", "v_entropy_change_pos"),
#    c("clock", "feedback", "v_entropy_change_neg")
#    c("clock", "feedback", "rew_om"),
#    c("clock", "feedback", "pe_max", "rew_om")
#    m1=c("clock", "feedback", "pe_1h", "pe_2h"), #use model names to cross-reference in add'l l1_contrasts
#    m2=c("clock", "feedback", "v_entropy_1h", "v_entropy_2h"),
#    c("clock", "feedback", "v_entropy", "pe_max"), #simultaneous model
#    c("clock", "feedback", "pe_trial_fixed_p05"), #fixed lr v .05
#    c("clock", "feedback", "pe_trial_fixed_p10"), #fixed lr v .1
#    c("clock", "feedback", "pe_trial_fixed_p15"), #fixed lr v .15
#    c("clock", "feedback", "pe_trial_fixed_p20") #fixed lr v .2
  ),
  group_model_variants=list(
    c("Intercept"),
    c("Intercept", "Age")
#    c("Intercept", "Age", "Female"),
#    c("Intercept", "I_Age"),
#    c("Intercept", "I_Age", "Female")
  ),
  l1_contrasts=list( #these are always in addition to the diagonal matrix of contrasts for each regressor
    m1=list(
      pe1h_gt_pe2h=c(pe_1h=1, pe_2h=-1),
      peavg=c(pe_1h=0.5, pe_2h=0.5) #closest to earlier whole-run analysis
    ),
    m2=list(
      entropy1h_gt_entropy2h=c(v_entropy_1h=1, v_entropy_2h=-1),
      entropyavg=c(v_entropy_1h=0.5, v_entropy_2h=0.5)
    )
  ),
  execute_feat=FALSE, #passed through to fsl_sceptic_model to create fsf, but not run the model
  #model_suffix="_fse", #factorized, selective, equal generalization width
  model_suffix="_fse_groupfixed", #factorized, selective, equal generalization width
  #model_suffix="_fixed_groupfixed", #fixed lr v model at group mean params
  root_workdir="/pine/scr/m/n/mnhallq/run_fsl_pipeline_qsub_tmp",
  n_cluster_beta_cpus=8, #should be number of l2 contrasts, or lower
  badids = c(11335, #low IQ, ADHD Hx, loss of consciousness
    11332, #should be excluded, but scan was terminated early due to repeated movement
    11282, #RTs at the floor for essentially all runs. Not appropriate
    11246, #huge movement and RTs at floor
    #10637, #large and many movements in later runs (need to revisit to confirm) ### OCT2018: 6 of 8 runs pass our algorithmic thresholds for motion
    10662  #I think there are reconstruction problems here -- need to revisit
  )
)

#validate and populate any other pipeline details before execution
fsl_model_arguments <- finalize_pipeline_configuration(fsl_model_arguments)

save(fsl_model_arguments, file=paste0("configuration_files/", fsl_model_arguments$analysis_name, ".RData"))

#this pushes the full analysis pipeline in parallel, where parallelism is across sceptic_run_variants
push_pipeline(fsl_model_arguments, ncpus=fsl_model_arguments$pipeline_cpus)


## Other pipelines go here
