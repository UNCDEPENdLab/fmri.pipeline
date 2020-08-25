library(dependlab)

#load("/gpfs/group/mnh5174/default/clock_analysis/fmri/fsl_pipeline/configuration_files/MMClock_aroma_preconvolve_fse.RData") #current arguments
load("/gpfs/group/mnh5174/default/clock_analysis/fmri/fsl_pipeline/configuration_files/MMClock_aroma_preconvolve_fse_groupfixed.RData") #current arguments
#load("/gpfs/group/mnh5174/default/clock_analysis/fmri/fsl_pipeline/configuration_files/MMClock_fixed_aroma_preconvolved.RData") #fixed lr v

for (run_model_index in 1:length(fsl_model_arguments$outdir)) {
  for (extract_z in c(0,1)) {
    qsub_file(script="qsub_fmri_r_script.bash",
      pbs_args=c("-l nodes=1:ppn=8", "-l walltime=3:00:00", "-l pmem=8gb"),
      env_variables=c(R_SCRIPT="extract_sceptic_map_betas_atlas.R",
        run_model_index=run_model_index,
        extract_z=extract_z, #betas versus z-statistics
        fsl_pipeline_file=file.path(fsl_model_arguments$pipeline_home, "configuration_files", paste0(fsl_model_arguments$analysis_name, ".RData")))
    )

  }


}
