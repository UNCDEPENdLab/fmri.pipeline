## Leverage the existing pipeline to run SPM GLMs followed by DCMs
#the goal of this script is to run an entire fmri analysis for SCEPTIC data, including levels 1-3 in FSL

library(dependlab)
library(foreach)
library(parallel)
library(doParallel)
library(plyr)
library(tidyverse)

scripts_dir <- "/gpfs/group/mnh5174/default/clock_analysis/fmri/fsl_pipeline"
setwd(scripts_dir)

source(file.path(scripts_dir, "functions", "push_pipeline.R"))
source(file.path(scripts_dir, "functions", "finalize_pipeline_configuration.R"))

#Jun2017: further ICAs on the MMClock data suggest a short steady-state problem. Drop 2 volumes for good measure.

###
# SCEPTIC MMClock Y3

#trial_df <- read.csv("/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs/mmclock_fmri_decay_factorize_mfx_trial_statistics.csv.gz")
#trial_df <- read.csv("/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs/mmclock_fmri_decay_mfx_trial_statistics.csv.gz")


#factorized, selective maintenance, equal basis-generalization width
#trial_df <- read.csv("/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs/mmclock_fmri_decay_factorize_selective_psequate_mfx_trial_statistics.csv.gz") %>%
trial_df <- read.csv("/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs/mmclock_fmri_decay_factorize_selective_psequate_fixedparams_ffx_trial_statistics.csv.gz") %>%
  mutate(
    run_trial=case_when(
    trial >= 1 & trial <= 50 ~ trial,
    trial >= 51 & trial <= 100 ~ trial - 50L, #dplyr/rlang has gotten awfully picky about data types!!
    trial >= 101 & trial <= 150 ~ trial - 100L,
    trial >= 151 & trial <= 200 ~ trial - 150L,
    trial >= 201 & trial <= 250 ~ trial - 200L,
    trial >= 251 & trial <= 300 ~ trial - 250L,
    trial >= 301 & trial <= 350 ~ trial - 300L,
    trial >= 351 & trial <= 400 ~ trial - 350L,
    TRUE ~ NA_integer_),
    v_entropy_no5=if_else(run_trial <= 5, NA_real_, v_entropy),
    d_auc_sqrt=if_else(d_auc > 0, NA_real_, sqrt(-1*d_auc)), #only compute the sqrt of d_auc for negative (i.e., reasonable) observations
    v_entropy_sqrt=sqrt(v_entropy),
    rew_om=if_else(score_vba > 0, 1, 0)
  ) %>% #for win/loss maps
  group_by(id, run) %>%
  dplyr::mutate(   #compute rt_swing within run and subject
    rt_vmax_lag = dplyr::lag(rt_vmax, 1, order_by=trial),
    rt_vmax_change = abs(rt_vmax - rt_vmax_lag),
    v_entropy_lag = dplyr::lag(v_entropy, 1, order_by=trial),
    v_entropy_change = v_entropy - v_entropy_lag, #change in entropy
    v_entropy_change_pos = v_entropy_change*(v_entropy_change > 0),
    v_entropy_change_neg = abs(v_entropy_change*(v_entropy_change < 0)),
    rt_swing = abs( c(NA, diff(rt_csv)))/1000,
    rt_swing_sqrt=sqrt(rt_swing)) %>%
  ungroup()


subject_df <- read.table("/gpfs/group/mnh5174/default/clock_analysis/fmri/data/mmy3_demographics.tsv", header=TRUE) %>%
  rename(id=lunaid, Age=age, Female=female, ScanDate=scandate) %>%
  mutate(mr_dir = paste0("/gpfs/group/mnh5174/default/MMClock/MR_Proc/", id, "_", format((as.Date(ScanDate, format="%Y-%m-%d")), "%Y%m%d")), #convert to Date, then reformat YYYYMMDD
    I_Age = -1000*1/Age,
    I_Age_c = I_Age - mean(I_Age, na.rm=TRUE),
    Age_c = Age - mean(Age, na.rm=TRUE),
    Q_Age = Age_c^2,
    Q_Age_c = Q_Age - mean(Q_Age, na.rm=TRUE)
  )

#from 2017:
##results from Mean SCEPTIC regressor correlation.pdf indicate that regressors for vchosen, ventropy_decay_matlab, dauc, and pemax are
##reasonably uncorrelated. The worst is dauc with vchosen (mean r = -0.31), which makes sense that as learning progresses, chosen values
##are higher and there is less residue to decay. These 4 regressors are also of greatest theoretical interest

#Setup the global configuration for the full FSL pipeline
fsl_model_arguments <- list(
  analysis_name="MMClock_aroma_preconvolve_fse_groupfixed",
  trial_statistics = trial_df,
  subject_covariates = subject_df,
  id_col = "id",
  fmri_dir = "/gpfs/group/mnh5174/default/MMClock/MR_Proc",
  expectdir = "mni_5mm_aroma", #subfolder name for processed data
  expectfile = "nfaswuktm_clock[0-9]_5.nii.gz", #expected file name for processed clock data
  usepreconvolve=TRUE,
  ncpus=20,
  drop_volumes=2, #to handle steady state concerns
  tr=1.0, #seconds
  spikeregressors=FALSE, #don't include spike regressors in nuisance variables since we are using AROMA
  glm_software="spm", #will run SPM GLMs
  sceptic_run_variants=list(
    c("clock", "feedback", "v_entropy", "pe_max") #core DCM model
  ),
  group_model_variants=list(
    c("Intercept"),
    c("Intercept", "Age")
  ),    
  execute_feat=FALSE, #passed through to fsl_sceptic_model to create fsf, but not run the model
  model_suffix="_fse_groupfixed", #factorized, selective, equal generalization width
  root_workdir="/gpfs/scratch/mnh5174/run_fsl_pipeline_qsub_tmp",
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
#push_pipeline(fsl_model_arguments, ncpus=fsl_model_arguments$pipeline_cpus)


#execute_fsl_lvl1_pipeline.R is just a thin qsub wrapper around model_clock_fmri_lvl1 to aid in parallelism
setup_fsf_jobid <- qsub_file(script="qsub_fmri_r_script.bash",
  pbs_args=c("-l nodes=1:ppn=20", "-l walltime=10:00:00"),
  env_variables=c(R_SCRIPT="execute_fsl_lvl1_pipeline.R",
    run_model_index=1,
    fsl_pipeline_file=file.path(fsl_model_arguments$pipeline_home, "configuration_files", paste0(fsl_model_arguments$analysis_name, ".RData")))
)

## Other pipelines go here


# masks <- c("/gpfs/group/mnh5174/default/MMClock/MR_Proc/11347_20141205/mni_5mm_aroma/sceptic-clock-feedback-d_auc-preconvolve_fse_groupfixed/spm_out/masks/l_v1_2.3mm.nii.gz", "/gpfs/group/mnh5174/default/MMClock/MR_Proc/11347_20141205/mni_5mm_aroma/sceptic-clock-feedback-d_auc-preconvolve_fse_groupfixed/spm_out/masks/r_v1_2.3mm.nii.gz",
#    "/gpfs/group/mnh5174/default/MMClock/MR_Proc/11347_20141205/mni_5mm_aroma/sceptic-clock-feedback-d_auc-preconvolve_fse_groupfixed/spm_out/masks/r_hippo.nii.gz", "/gpfs/group/mnh5174/default/MMClock/MR_Proc/11347_20141205/mni_5mm_aroma/sceptic-clock-feedback-d_auc-preconvolve_fse_groupfixed/spm_out/masks/l_motor_2.3mm.nii"
#  )

# aa <- spm_extract_anatomical_rois("/gpfs/group/mnh5174/default/MMClock/MR_Proc/11347_20141205/mni_5mm_aroma/sceptic-clock-feedback-d_auc-preconvolve_fse_groupfixed/spm_out", masks, threshold=0.2, threshdesc='none', session=1, extent=0, adjust_F_index=1, contrast_index=2, ncores=1, spm_path="/gpfs/group/mnh5174/default/lab_resources/spm12", matlab_path="/opt/aci/sw/matlab/R2017b/bin")
