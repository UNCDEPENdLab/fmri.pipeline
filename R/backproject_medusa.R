#' @importFrom checkmate assert_data_frame assert_string assert_file_exists assert_character
#' @importFrom oro.nifti readNIfTI
#' @importFrom dplyr bind_rows filter select
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#' @importFrom iterators iter
backproject_medusa <- function(coef_df, brain_mask, plot_cols=NULL, parcel_col="roi", time_col="time", 
                               effect_col="term", output_dir=getwd(), ncpus=4L) {
  
  checkmate::assert_data_frame(coef_df)
  checkmate::assert_string(brain_mask)
  checkmate::assert_file_exists(brain_mask)
  checkmate::assert_character(plot_cols, null.ok = FALSE)
  checkmate::assert_string(parcel_col)
  checkmate::assert_string(time_col)
  checkmate::assert_string(effect_col)
  checkmate::assert_subset(parcel_col, names(coef_df))
  checkmate::assert_subset(time_col, names(coef_df))
  checkmate::assert_subset(effect_col, names(coef_df))
  
  #read mask
  mask_img <- oro.nifti::readNIfTI(brain_mask, reorient=FALSE)
  if (mask_img@dim_[1] != 3) { stop("mask_img does not appear to be a 3-D mask image.") }
  
  mvals <- unique(coef_df[[parcel_col]])
  imgvals <- sort(unique(as.vector(mask_img)))
  imgvals <- imgvals[imgvals != 0]
  
  #verify that roi mask values in the coefficients df exist in the mask
  checkmate::assert_subset(mvals, imgvals)
  
  imgpos <- dplyr::bind_rows(lapply(mvals, function(xx) {
    matches <- data.frame(which(mask_img==xx, arr.ind=TRUE))
    matches$roi <- xx
    return(matches)
  }))
  
  #unique effect coefficients
  effects <- sort(unique(coef_df[[effect_col]]))
  
  #parallelize over the combination of plot column and effect
  to_loop <- expand.grid(pc=plot_cols, ee=effects)
  
  if (ncpus > 1L) {
    cl <- makeCluster(ncpus)
    on.exit(try(stopCluster(cl)))
    registerDoParallel(cl)
  } else {
    registerDoSEQ()
  }
  
  ff <- foreach(ii=iter(to_loop, by="row"), .packages=c("dplyr", "oro.nifti")) %dopar% {
    
    if (!is.null(effect_col)) {
      this_df <- coef_df %>% dplyr::filter(!!sym(effect_col) == !!ii$ee)
    } else {
      this_df <- coef_df #just use all rows
    }
    
    timevals <- sort(unique(this_df[[time_col]]))
    n_t <- length(timevals)
    tr <- median(diff(timevals))
    
    results_img <- mask_img
    results_img@.Data <- array(0, dim=c(dim(mask_img)[1:3], n_t))
    results_img@dim_[1:5] <- c(4, results_img@dim_[2:4], n_t)
    results_img@datatype <- 16L #default to float
    results_img@bitpix <- 32L #floats require 32 bit storage
    results_img@pixdim[5] <- tr #set time step
    
    assign_mat <- list()
    
    for (mm in mvals) {
      for (tt in 1:n_t) {
        this_val <- this_df %>% dplyr::filter(!!sym(time_col) == !!timevals[tt] & !!sym(parcel_col) == !!mm) %>% pull(!!ii$pc)
        if (length(this_val) != 1L) { browser() }
        
        to_fill <- imgpos %>% dplyr::filter(roi == !!mm) %>% dplyr::select(-roi) %>% as.matrix()
        to_fill <- cbind(to_fill, tt, this_val) #4-D indices for this mask value and this timepoint plus the value to fill in col 5
        assign_mat[[paste(mm, tt, sep=".")]] <- to_fill
        
      }  
    }
    
    assign_mat <- do.call(rbind, assign_mat)
    results_img[assign_mat[,1:4]] <- assign_mat[,5] #fill all voxels at once with coef values
    results@cal_min <- min(as.vector(results_img)) #add min/max into nifti header
    results@cal_max <- max(as.vector(results_img))
    
    
    writeNIfTI(results_img, filename=file.path(output_dir, make.names(paste(ii$pc, ii$ee, sep="_"))))
    
    return(NULL)
  }
}

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
  
