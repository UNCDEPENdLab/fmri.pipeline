
library(tidyverse)
library(readxl)
library(doParallel)
library(foreach)
library(oro.nifti)


mm <- "/Users/hallquist/Data_Analysis/clock_analysis/fmri/pfc_entropy/original_masks/Schaefer2018_200Parcels_7Networks_order_fonov_1mm_ants.nii.gz"
# img <- oro.nifti::readNIfTI(mm, reorient = FALSE)
# sort(unique(as.vector(img)))


load("~/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/medusa_results_summary_dataframes/clock_decode_output.Rdata")

orig <- readxl::read_excel("/Users/hallquist/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/MNH Dan Labels.xlsx") %>%
  rowwise() %>%
  mutate(label=paste(Visuomotor_Gradient, Stream_Gradient, sub("([LR])_(.*)", "\\2_\\1", MNHLabel, perl=TRUE), sep="_")) %>%
  select(roinum, label)

ddf <- ddf %>% left_join(orig) %>% filter(term != "(Intercept)") #not interesting

setdiff(unique(ddf$label), orig$label)

rdfiles <- list.files("~/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/medusa_results_summary_dataframes", 
                      pattern=".*\\.Rdata", full.names = TRUE)
for (ff in rdfiles) {
  rd <- tools::file_path_sans_ext(basename(ff))
  odir <- file.path("/Users/hallquist/OneDrive/collected_letters/papers/sceptic_fmri/dan/medusa_backproject", rd)
  dir.create(odir, showWarnings = FALSE)
  load(ff)
  if (grepl("predict", rd)) {
    to_analyze <- rdf
  } else {
    to_analyze <- ddf
  }
  
  to_analyze <- to_analyze %>% left_join(orig) %>% 
    filter(term != "(Intercept)" & substr(term,1,3) != "sd_") #not interesting

  backproject_medusa(to_analyze, brain_mask = mm, plot_cols=c("estimate", "statistic", "p.value", "p_fdr"),
                     parcel_col="roinum", ncpus=5,
                     time_col="t", effect_col="term", output_dir=odir)
}


