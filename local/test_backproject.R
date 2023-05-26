
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




## temp code for ggbrain:brain


mm <- "/Users/hallquist/Data_Analysis/clock_analysis/fmri/pfc_entropy/original_masks/Schaefer2018_200Parcels_7Networks_order_fonov_1mm_ants.nii.gz"
# img <- oro.nifti::readNIfTI(mm, reorient = FALSE)
# sort(unique(as.vector(img)))


load("~/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/medusa_results_summary_dataframes/clock_decode_output.Rdata")

library(tidyverse)
library(readxl)
library(doParallel)
library(foreach)
library(oro.nifti)
anat <- readNIfTI("~/OneDrive/collected_letters/papers/sceptic_fmri/dan/medusa_backproject/clock_decode_output/mni_icbm152_t1_tal_nlin_asym_09c_brain.nii.gz", reorient=FALSE)
nzpos <- which(anat != 0, arr.ind=TRUE)
apply(anat, c(1,2,3), min)
medusa4d <- readNIfTI("~/OneDrive/collected_letters/papers/sceptic_fmri/dan/medusa_backproject/clock_decode_output/estimate_scale.rt_vmax_change..nii.gz", reorient=FALSE)

overlay(anat, medusa4d[,,,1])

library(microbenchmark)

microbenchmark(
  baser={
    png("test.png", res=300, width=8, height=8, units="in")
    image(anat[,,70])
    dev.off()
  },
  ggplotr={
    vv <- reshape2::melt(anat)
    png("ggplot.png", res=300, width=8, height=8, units="in")
    ggplot(vv, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + scale_fill_gradient(low="#000000", high="#FFFFFF") + 
      coord_fixed() + theme_void() +
      geom_tile(data=overlay_df)
    dev.off()
  }
)

source("/Users/hallquist/Data_Analysis/r_packages/fmri.pipeline/R/ggbrain.R")

setwd("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas/L1m-entropy_wiz")
pdf("entropy_wb_ptfce.pdf", width=11, height=5)
#png("example_plot.png", width=11, height=5, units="in", res=600)
# ggbrain(underlay="template_brain.nii", overlay = "zstat6_ptfce_fwep_0.05_1mm.nii.gz", axial_slices = list(xyz=c(58)), 
#         remove_null_space = TRUE, pos_thresh = 5, neg_thresh = -5, background_color = "black", text_color = "white") # already thresholded by ptfce

g <- ggbrain(underlay="template_brain.nii", overlay = "zstat6_ptfce_fwep_0.05_1mm.nii.gz",
        slices = data.frame(coord = c("z = 58", "y = -7", "y = 50%"), coord_labels = TRUE),
        remove_null_space = TRUE, pos_thresh = 5.11, neg_thresh = -5.11, 
        background_color = "black", text_color = "white",
        panel_borders = FALSE, symmetric_legend = TRUE, base_size = 18,
        positive_colorscale = scale_fill_viridis_c(), underlay_contrast = "medium"
)

#slices = data.frame(coord = c("z = 58", "x = -7.45", "y = 10")),         
#positive_colorscale = scale_fill_viridis_c(), negative_colorscale = scale_fill_distiller(palette="Blues"))

dev.off()


# entropy
# axial z = 58 IPS, FEF, Motor
# sagittal x = 5, IPS, neg vmPFC