# This script sets up the .fsf files to run a group analysis
## FSL Feat Level 3 analysis -- that is, mixed effects combination of subjects

#used for testing: group fixed entropy
#Sys.setenv(fsl_pipeline_file="/gpfs/group/mnh5174/default/clock_analysis/fmri/fsl_pipeline/configuration_files/MMClock_aroma_preconvolve_fse_groupfixed.RData")
#Sys.setenv(run_model_index=2)

#load the master configuration file
to_run <- Sys.getenv("fsl_pipeline_file")

run_model_index <- as.numeric(Sys.getenv("run_model_index")) #which variant to execute
if (nchar(to_run) == 0L) { stop("Cannot locate environment variable fsl_pipeline_file") }
if (!file.exists(to_run)) { stop("Cannot locate configuration file", to_run) }
if (is.na(run_model_index)) { stop("Couldn't identify usable run_model_index variable.") }

load(to_run)

library(tidyverse)
library(dependlab)

#verify that mr_dir is present as expected
subinfo <- fsl_model_arguments$subject_covariates
id_col <- fsl_model_arguments$id_col
feat_run_outdir <- fsl_model_arguments$outdir[run_model_index] #the name of the subfolder for the current run-level model
feat_lvl3_outdir <- file.path(fsl_model_arguments$group_output_dir, feat_run_outdir) #output directory for this run-level model
n_l1_copes <- fsl_model_arguments$n_l1_copes[run_model_index] #number of l1 copes determines number of FEAT LVL3 analyses to run (1 per LVL1 cope)
l1_cope_names <- fsl_model_arguments$l1_cope_names[[run_model_index]] #names of l1 copes (used for folder naming)

subinfo$dir_found <- file.exists(subinfo$mr_dir)

rerun <- FALSE #TODO: move to fsl_model_arguments list

cat("The following subjects were in the covariate file, but not the processed MRI data\n")
print(subset(subinfo, dir_found==FALSE))

dir.create(feat_lvl3_outdir, showWarnings=FALSE, recursive=TRUE)
setwd(feat_lvl3_outdir)

#cope structure for preconvolve models
#1 = clock_onset
#2 = feedback_onset
#3 = regressor of interest (in single-param models)

feat_lvl2_dirname <- "FEAT_LVL2_runtrend.gfeat" #should populate this to the structure at some point
models <- fsl_model_arguments$group_model_variants #different covariate models for the current run-level model (run_model_index)

##rework using subinfo structure as the authoritative guide (rather than repeated searches)
copedf <- c()
for (s in 1:nrow(subinfo)) {
  for (cope in 1:n_l1_copes) {
    expectdir <- file.path(subinfo[s,"mr_dir"], fsl_model_arguments$expectdir, feat_run_outdir, feat_lvl2_dirname, paste0("cope", cope, ".feat"))
    if (dir.exists(expectdir)) {
      copedf <- rbind(copedf, data.frame(id=subinfo[s,id_col], model=feat_run_outdir, cope=cope, fsldir=expectdir))
    } else {
      message("could not find expected directory: ", expectdir)
    }
  }
}

names(copedf)[1] <- id_col #for matching

mdf <- merge(subinfo, copedf, by=id_col, all.y=TRUE)

#remove bad ids
mdf <- mdf %>% filter(!id %in% fsl_model_arguments$badids)
mdf <- arrange(mdf, id, model, cope) #should really use !!id_col here?

##fsl constructs models by cope
bycope <- lapply(split(mdf, mdf$cope), droplevels)

#loop over group-level models, setup the design matrix and spawn a FSL Level 3 job
l3template <- readLines(file.path(getMainDir(), "clock_analysis", "fmri", "fsf_templates", "feat_lvl3_sceptic_template.fsf"))

#loop over copes and group models, setting up .fsf files for each combination
for (cope in 1:length(bycope)) {
  if (is.null(bycope[[cope]]$Intercept)) { bycope[[cope]]$Intercept <- 1 } #add the column of ones

  copename <- l1_cope_names[cope]
  
  #cope-level subfolder
  #currently organized by run_model_name/l1_cope/l3_model (promotes comparisons of alternative l3 models of a given run-level effect)
  #could reorganize as run_model_name/l3_model/l1_cope (promotes comparisons of maps within an l3 model)
  model_output_dir <- file.path(feat_lvl3_outdir, copename) # paste0("cope", cope))
  dir.create(model_output_dir, showWarnings=FALSE)
  
  for (this_model in models) {

    model_df <- bycope[[cope]]
    fsf_syntax <- l3template #copy shared ingredients
    fsf_syntax <- gsub(".OUTPUTDIR.", file.path(model_output_dir, paste0(copename, "-", paste(this_model, collapse="-"))), fsf_syntax, fixed=TRUE)
    if (!"Intercept" %in% this_model) { this_model <- c("Intercept", this_model) } #at present, force an intercept column
    
    if (fsl_model_arguments$center_l3_predictors) {
      for (p in this_model) {
        if (p != "Intercept" && is.numeric(model_df[[p]])) {
          model_df[[p]] <- model_df[[p]] - mean(model_df[[p]], na.rm=TRUE)
        }
      }
    }

    model_df$dummy_ <- rnorm(nrow(model_df))
    mform <- as.formula(paste("dummy_ ~ -1 + ", paste(this_model, collapse=" + ")))
    fit_lm <- lm(mform, model_df)
    dmat <- model.matrix(fit_lm) #eventually allow interactions and so on??
    model_df$dummy_ <- NULL #clean up
    
    #add design matrix
    fsf_syntax <- c(fsf_syntax, generate_fsf_ev_syntax(inputs=model_df$fsldir, dmat=dmat))

    #generate diagonal contrast matrix, one per EV
    cmat <- diag(length(this_model))
    rownames(cmat) <- this_model #just name the contrasts after the EVs themselves

    fsf_syntax <- c(fsf_syntax, generate_fsf_contrast_syntax(cmat))

    #write the FSF to file
    out_fsf <- file.path(model_output_dir, paste0(copename, "-", paste(this_model, collapse="-"), ".fsf"))

    if (!file.exists(out_fsf) || rerun) { writeLines(fsf_syntax, con=out_fsf) }

    if (!file.exists(dmat_file <- sub(".fsf", "_design.txt", out_fsf, fixed=TRUE)) || rerun) {
      #write the design matrix to file for matching with extracted betas later
      model_df$feat_input_id <- 1:nrow(model_df) #for matching with extracted betas
      model_df <- model_df %>% select(-dir_found, -mr_dir) %>% select(id, feat_input_id, model, cope, fsldir, everything())
      write.table(model_df, file=dmat_file, row.names=FALSE)
    }
    
    if (!file.exists(sub(".fsf", ".gfeat", out_fsf, fixed=TRUE)) || rerun) {
      #run the L3 analysis in parallel using qsub
      qsub_file(script=file.path(fsl_model_arguments$pipeline_home, "qsub_feat_lvl3.bash"), env_variables=c(torun=out_fsf))
    }

  }

}
