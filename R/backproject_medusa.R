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
        this_val <- this_df %>%
          dplyr::filter(!!sym(time_col) == !!timevals[tt] & !!sym(parcel_col) == !!mm) %>%
          pull(!!ii$pc)
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
