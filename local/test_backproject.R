
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


#' Plot fMRI data on an underlay image
#' 
#' @param underlay a 3D nifti image used for the image underlay (default b/w)
#' @param overlay a 4D nifti image used for plotting stats on underlay (color)
#' @param color_col a position in the 4th dim of overlay use to color plots
#' @param alpha_col a position in the 4th dim of overlay use to set alpha transparency of plots
#' @param underlay_colorscale A ggplot scale_fill_* function call used for coloration of underlay
#' @param overlay_colorscale A ggplot scale_fill_* function call used for coloration of overlay
ggbrain <- function(underlay=NULL, overlay=NULL, 
                    color_col=NULL, alpha_col=NULL,
                    underlay_colorscale=scale_fill_gradient(low="grey10", high="grey90"),
                    overlay_colorscale=scale_fill_gradient2(midpoint = 0, low = "blue", mid="grey90", high="red"),
                    axial_slices=c(.25, .50, .75),
                    sagittal_slices=c(.25, .50, .75),
                    coronal_slices=c(.25, .50, .75),
                    remove_null_space=TRUE
                    )

{
  checkmate::assert_file_exists(underlay)
  checkmate::assert_file_exists(overlay)
  
  underlay <- oro.nifti::readNIfTI(underlay, reorient = FALSE)
  overlay <- oro.nifti::readNIfTI(overlay, reorient = FALSE)
  
  #if overlay only has positive values, use
  #scale_fill_gradient(low = "grey90", high="red")
  
  #if overlay only has negative values, use
  #scale_fill_gradient(low = "grey90", high="blue")
  
  #verify that i,j,k (1,2,3) dimensions of underlay match dimensions of overlay
  stopifnot(identical(dim(underlay)[1:3], dim(overlay)[1:3]))
  
  if (isTRUE(remove_null_space)) {
    nzpos <- which(underlay != 0, arr.ind=TRUE)
    minx <- min(nzpos[,1]) - 1 #subtract one to give a small null space
    maxx <- max(nzpos[,1]) + 1
    miny <- min(nzpos[,2]) - 1
    maxy <- max(nzpos[,2]) + 1
    minz <- min(nzpos[,3]) - 1
    maxz <- max(nzpos[,3]) + 1
  
    xrange <- min(nzpos[,1]):max(nzpos[,1])
    yrange <- min(nzpos[,2]):max(nzpos[,2])
    zrange <- min(nzpos[,3]):max(nzpos[,3])
      
    #trim images
    underlay <- underlay[minx:maxx, miny:maxy, minz:maxz]
    overlay <- overlay[minx:maxx, miny:maxy, minz:maxz]
  } else {
    allpos <- which(underlay != Inf, arr.ind=TRUE)
    xrange <- min(allpos[,1]):max(allpos[,1])
    yrange <- min(allpos[,2]):max(allpos[,2])
    zrange <- min(allpos[,3]):max(allpos[,3])
  }
  
  #validate 0-1 bounds on arguments: axial_slices, sagittal_slices, coronal_slices
  anat_slices <- list(
    sagittal=aperm(underlay[unique(floor(quantile(xrange, sagittal_slices))),,], c(1,2,3)),
    coronal=aperm(underlay[,unique(floor(quantile(yrange, coronal_slices))),], c(2,1,3)),
    axial=aperm(underlay[,,unique(floor(quantile(zrange, axial_slices)))], c(3,1,2))
  )
  
  anat_slices <- dplyr::bind_rows(lapply(seq_along(anat_slices), function(ll) {
    df <- reshape2::melt(anat_slices[[ll]], varnames=c("slice", "dim1", "dim2"))
    df$type <- names(anat_slices)[ll]
    
    #uniquely identify slices by number and type
    df$slice <- paste0(substr(names(anat_slices)[ll],1,1), df$slice)
    return(df)
  }))
  
  slice_plots <- ggplot(anat_slices, aes(x=dim1, y=dim2, fill=value)) +
      geom_tile() + coord_fixed() + theme_void() +
      scale_fill_gradient(low="grey10", high="grey90") +
      facet_wrap(~slice) + guides(fill=FALSE)
  
  plot(slice_plots)
  
  #facet_wrap the slices, allow background color to be set as argument
  #also squish the panels together with no spacing to form contiguous image
  #remove panel headers
  #add annotation of x=<> (coordinate) on sagittal slices
  #add annotation of y=<> (coordinate) on coronal slices
  #add annotation of z=<> (coordinate) on axial slices
  #use oro.nifti::translateCoordinate() to lookup spatial locations within matrix to get annotations
  
  #add overlays based on structure of overlay image -- generate one ggplot per element in 4th dimension 
  #provide options for rendering to file (ggsave to pdf and so on)
  
}
  
