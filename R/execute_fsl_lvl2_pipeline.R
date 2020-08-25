## This script sets up the .fsf files to combine runs within subject using an
## FSL Feat Level 2 analysis -- that is, fixed effects combinations of runs.

#load("/gpfs/group/mnh5174/default/clock_analysis/fmri/fsl_pipeline/configuration_files/SPECC_aroma_preconvolve_fse_groupfixed.RData")

#load the master configuration file
to_run <- Sys.getenv("fsl_pipeline_file")

run_model_index <- as.numeric(Sys.getenv("run_model_index")) #which variant to execute
if (nchar(to_run) == 0L) { stop("Cannot locate environment variable fsl_pipeline_file") }
if (!file.exists(to_run)) { stop("Cannot locate configuration file", to_run) }
if (is.na(run_model_index)) { stop("Couldn't identify usable run_model_index variable.") }

load(to_run)

library(tidyverse)

source(file.path(fsl_model_arguments$pipeline_home, "functions", "run_feat_lvl2.R"))
source(file.path(fsl_model_arguments$pipeline_home, "functions", "setup_feat_lvl2_inputs.R"))
source(file.path(fsl_model_arguments$pipeline_home, "functions", "glm_helper_functions.R"))

#first, identify the inputs for LVL2 analysis (doesn't really benefit from parallel execution)
#this creates a data.frame with first-level runs that should be combined.
feat_l2_inputs_df <- setup_feat_lvl2_inputs(fsl_model_arguments, run_model_index)

#generate .fsf files for LVL2 (subject-level) analysis, then run the analyses using feat
run_feat_lvl2(feat_l2_inputs_df, run=TRUE, force=FALSE, ncpus=fsl_model_arguments$l2_cpus)
