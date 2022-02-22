#' Function to perform voxelwise deconvolution on an fMRI dataset using the fMRI model arguments object
#'
#' @param niftis A vector of processed fMRI timeseries images (4D files) to be deconvolved
#' @param add_metadata A data.frame with one row per value of \code{niftis}. Columns of this data.frame
#'   are added to the output file for identification.
#' @param out_dir Base output directory for deconvolved time series files. Default is \code{getwd()}.
#' @param out_file_expression Expression evaluated to resolve the filename for the deconvolved csv files. Default is
#'   to convert the \code{niftis} value for a given subject into a filename by replacing slashes with periods and adding
#'   the atlas image name. Note that the suffix \code{_deconvolved.csv.gz} or \code{_original.csv.gz} will be added
#'   to the expression, so don't pass this piece.
#' @param log_file Name (and path) of log file for any deconvolution errors or messages
#' @param TR the repetition time of the sequence in seconds. Required
#' @param time_offset The number of seconds that will be subtracted or added to the time field. Default: 0. Useful if some number of volumes
#'   have been dropped from the NIfTI data prior to deconvolution.
#' @param atlas an optional character vector specifying voxels used in deconvolution. If omitted, perform whole-brain deconvolution
#' @param mask an optional character string specifying a mask that should be used to constrain bounds of deconvolution.
#' @param nprocs The number of processors to use simultaneously for deconvolution
#' @param save_original_ts Whether to save the voxelwise BOLD data prior to deconvolution (for comparison/diagnosis). Default: TRUE
#' @param algorithm Which deconvolution algorithm to use for deconvolving voxelwise time series. Default: "bush2011". Alternative is
#'   "bush2015", which implements a resampling approach as well. If you use "bush2011", the function will try to call on a fast compiled
#'   version of the algorithm to support whole-brain processing.
#' @param decon_settings A list of settings passed to the deconvolution algorithm. If you have a compiled deconvolvefilter binary,
#'   pass it as \code{bush2011_binary}, which will be used in deconvolution.
#' @param afni_dir Full path to directory containing AFNI binaries (this function uses 3dMaskdump).
#'
#' @details The Bush 2011 algorithm is implemented in a compiled binary called deconvolvefilter
#'   (https://github.com/UNCDEPENdLab/deconvolution-filtering) that is much faster than the pure R (or original MATLAB) version.
#'   We recommend using this for whole-brain deconvolution. The package includes a binary for
#'   the Linux x86_64 architecture.
#'
#'   If you want to use subject metadata to name the output file, use \code{this_subj} in your \code{out_file_expression}, which will give you access to a one-row
#'   data.frame containing the metadata for the current subject in the loop.
#'
#' @examples
#'   \dontrun{
#'
#'     #name outputs according to subject metadata
#'     xx <- voxelwise_deconvolution(
#'       niftis="/proj/mnhallqlab/user/michael/test_nifti.nii.gz",
#'       out_dir="/proj/mnhallqlab/user/michael/decon_outputs",
#'       out_file_expression=expression(paste0(this_subj$subid, "_run", this_subj$run_num, "_", atlas_img_name))
#'     )
#'   }
#'
#' @return Nothing (invisible NULL).
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster
#' @importFrom data.table fread
#' @importFrom oro.nifti readNIfTI translateCoordinate
#' @importFrom checkmate assert_file_exists
#' @importFrom foreach foreach
#' @importFrom dplyr mutate mutate_at select left_join
#' @importFrom readr write_delim write_csv
#' @export
voxelwise_deconvolution <- function(
  niftis, add_metadata=NULL, out_dir=getwd(), out_file_expression=NULL, 
  log_file=file.path(out_dir, "deconvolve_errors"),
  TR = NULL, time_offset=0, atlas_files=NULL, mask=NULL, nprocs=20, save_original_ts=TRUE, algorithm="bush2011",
  decon_settings=list(nev_lr = .01, #neural events learning rate (default in algorithm)
    epsilon = .005, #convergence criterion (default)
    beta = 60, #best from Bush 2015 update
    kernel = spm_hrf(TR)$hrf, #canonical SPM difference of gammas
    n_resample = 25)) { #for Bush 2015 only

  sapply(niftis, checkmate::assert_file_exists)
  checkmate::assert_data_frame(add_metadata, nrows=length(niftis), null.ok=TRUE)
  sapply(atlas_files, checkmate::assert_file_exists)
  checkmate::assert_numeric(TR, lower=0.01)
  
  if (!is.null(decon_settings$kernel)) {
    hrf_pad <- length(decon_settings$kernel)
  } else {
    stop("At present no alternative way to set hrf_pad if $kernel not passed in decon_settings")
  }
  
  if (is.null(atlas_files)) {
    message("Performing whole-brain voxelwise deconvolution")
    if (is.null(mask)) {
      stop("Cannot conduct voxelwise deconvolution without a relevant mask")
    } else {
      checkmate::assert_file_exists(mask)
    }
    
    atlas_files <- mask
  }

  zero_thresh <- 1e-4 #for binarizing/indexing
  atlas_imgs <- lapply(atlas_files, readNIfTI, reorient=FALSE)
  if (!is.null(mask)) {
    checkmate::assert_file_exists(mask)
    mask <- readNIfTI(mask, reorient=FALSE) #go from mask as character to the mask image itself

    #ensure that mask is binary
    mask[mask < zero_thresh] <- 0.0
    mask[mask >= zero_thresh] <- 1.0
  }

  #setup cluster
  if (nprocs > 1) {
    cl <- makeCluster(nprocs)
    registerDoParallel(cl)
    on.exit(try(stopCluster(cl)))
  } else {
    registerDoSEQ()
  }

  #loop over atlas files
  for (ai in seq_along(atlas_files)) {
    cat("Working on atlas: ", atlas_files[ai], "\n")
    aimg <- atlas_imgs[[ai]]

    if (!is.null(mask)) { aimg <- aimg*mask } #multiply against binary mask to apply it

    #get indices of mask within matrix (ijk)
    a_indices <- which(aimg > zero_thresh, arr.ind=TRUE)

    #look up spatial coordinates of voxels in atlas (xyz)
    a_coordinates <- cbind(a_indices, t(apply(a_indices, 1, function(r) { translateCoordinate(i=r, nim=aimg, verbose=FALSE) })))
    a_coordinates <- as.data.frame(a_coordinates) %>%
      setNames(c("i", "j", "k", "x", "y", "z")) %>%
      dplyr::mutate(vnum = 1:n(), atlas_value = aimg[a_indices], atlas_name = basename(atlas_files[ai])) %>%
      mutate_at(vars(x, y, z), round, 2) %>%
      dplyr::select(vnum, atlas_value, everything())

    #setup output subdirectories for deconvolved files, named according to atlas
    atlas_img_name <- basename(sub(".nii(.gz)*", "", atlas_files[ai], perl=TRUE))
    dir.create(file.path(out_dir, atlas_img_name, "deconvolved"), showWarnings=FALSE, recursive=TRUE)
    if (isTRUE(save_original_ts)) { dir.create(file.path(out_dir, atlas_img_name, "original"), showWarnings=FALSE, recursive=TRUE) }

    #set a default filename based on the nifti -- hopefully this will not collide with other subjects/runs
    if (is.null(out_file_expression)) {
      out_file_expression <- expression(paste0(gsub("[/\\]", ".", niftis[si]), "_", atlas_img_name))
    }

    if (!is.null(add_metadata)) { add_metadata$.nifti <- NA_character_ } #initialize empty nifti string for population

    #loop over niftis in parallel
    ff <- foreach(si = seq_along(niftis), 
      .packages=c("dplyr", "readr", "data.table", "reshape2", "fmri.pipeline", "foreach", "iterators")) %dopar% {

      #get the si-th row of the metadata to match nifti, allow one to use this_subj in out_file_expression
      if (!is.null(add_metadata)) {
        this_subj <- add_metadata %>% dplyr::slice(si)
        add_metadata$.nifti[si] <- niftis[si]
      }

      out_name <- file.path(out_dir, atlas_img_name, "deconvolved", paste0(eval(out_file_expression), "_deconvolved.csv.gz"))

      if (file.exists(out_name)) {
        message("Deconvolved file already exists: ", out_name)
        return(NULL)
      }

      cat("  Deconvolving subject: ", niftis[si], "\n")
      dump_out <- tempfile()
      afnistat <- run_afni_command(paste0("3dmaskdump -mask ", atlas_files[ai], " -o ", dump_out, " ", niftis[si]))
      ts_out <- data.table::fread(dump_out) #read time series

      #to_deconvolve is a voxels x time matrix
      to_deconvolve <- as.matrix(ts_out[, -1:-3]) #remove ijk

      to_deconvolve <- t(apply(to_deconvolve, 1, scale)) #need to unit normalize for algorithm not to choke on near-constant 100-normed data

      # just demean, which will rescale to percent signal change around 0 (this matches Bush 2015)
      #to_deconvolve <- t(apply(to_deconvolve, 1, function(x) { scale(x, scale=FALSE) }))

      # pct signal change around 0
      #to_deconvolve <- t(apply(to_deconvolve, 1, function(x) { x/mean(x)*100 - 100 }))

      temp_i <- tempfile()
      temp_o <- tempfile()

      # to_deconvolve %>%  as_tibble() %>% write_delim(path=temp_i, col_names=FALSE)

      # zero pad tail end (based on various readings, but not original paper)
      # this was decided on because we see the deconvolved signal dropping to 0.5 for all voxels

      to_deconvolve %>%
        cbind(matrix(0, nrow = nrow(to_deconvolve), ncol = hrf_pad)) %>%
        as_tibble() %>%
        write_delim(path = temp_i, col_names = FALSE)

      #test1 <- deconvolve_nlreg(to_deconvolve[117,], kernel=decon_settings$kernel, nev_lr=decon_settings$nev_lr, epsilon=decon_settings$epsilon)
      #test2 <- deconvolve_nlreg(to_deconvolve[118,], kernel=decon_settings$kernel, nev_lr=decon_settings$nev_lr, epsilon=decon_settings$epsilon)

      if (algorithm == "bush2015") {
        #use R implementation of Bush 2015 algorithm
        alg_input <- as.matrix(data.table::fread(temp_i))
        deconv_mat <- foreach(vox_ts=iter(alg_input, by="row"), .combine="rbind", .packages=c("dependlab")) %do% {
          reg <- tryCatch(deconvolve_nlreg_resample(as.vector(vox_ts), kernel=decon_settings$kernel, nev_lr=decon_settings$nev_lr, epsilon=decon_settings$epsilon, n_resample=decon_settings$n_resample),
            error=function(e) { cat("Problem deconvolving: ", niftis[si], as.character(e), "\n", file=log_file, append=TRUE); return(rep(NA, length(vox_ts))) })

          if (is.list(reg)) { reg <- reg$NEVmean } #just keep the mean resampled events vector
          return(reg)
        }
      } else if (algorithm == "bush2011") {
        #this should use the new internal RcppArmadillo function
        alg_input <- as.matrix(t(data.table::fread(temp_i)))
        deconv_mat <- tryCatch(deconvolve_nlreg(BOLDobs = alg_input, kernel=decon_settings$kernel, nev_lr=decon_settings$nev_lr, epsilon=decon_settings$epsilon, beta=decon_settings$beta),
          error=function(e) {
            cat("Problem deconvolving: ", niftis[si], as.character(e), "\n", file=log_file, append=TRUE)
            return(matrix(NA, nrow=nrow(alg_input), ncol=ncol(alg_input)))
          })

        deconv_mat <- t(deconv_mat) #Rcpp function is time x voxels...
      } else if (algorithm == "bush2011_external") {
        #use C++ implementation of Bush 2011 algorithm, if possible
        decon_bin <- NULL #default to pure R algorithm
        if (is.null(decon_settings$bush2011_binary)) {
          ss <- Sys.info()
          if (ss["sysname"] == "Linux") {
            if (ss["machine"] == "x86_64") {
              #decon_bin <- system.file("bin", "linux_x64", "deconvolvefilter", package = "fmri.pipeline")
              decon_bin <- "/proj/mnhallqlab/users/michael/fmri.pipeline/inst/bin/linux_x64/deconvolvefilter"
            }
          }
        } else {
          decon_bin <- decon_settings$bush2011_binary
        }

        if (is.null(decon_bin)) {
          alg_input <- as.matrix(data.table::fread(temp_i))
          deconv_mat <- foreach(vox_ts=iter(alg_input, by="row"), .combine="rbind", .packages=c("dependlab")) %do% {
            reg <- tryCatch(deconvolve_nlreg(as.vector(vox_ts), kernel=decon_settings$kernel, nev_lr=decon_settings$nev_lr, epsilon=decon_settings$epsilon),
              error=function(e) { cat("Problem deconvolving: ", niftis[si], as.character(e), "\n", file=log_file, append=TRUE); return(rep(NA, length(vox_ts))) })
            return(reg)
          }
        } else {
          checkmate::assert_file_exists(decon_bin)

          #fo argument is 1/TR: https://github.com/UNCDEPENdLab/deconvolution-filtering/blob/f26df0ea1eb30f2019795f17f93e713517a220e4/ref/backup/deconvolve_filter.m
          #Looking at spm_hrf, it generates a vector of 33 values for the HRF for 1s TR. deconvolvefilter pads the time series at the beginning by this length
          #if you don't return a convolved result, it doesn't do the trimming for you...
          res <- system(paste0(decon_bin, " -i=", temp_i, " -o=", temp_o, " -convolved=0 -fo=", 1/TR, " -thread=2 >/dev/null 2>/dev/null"), intern=FALSE)
          if (res != 0) {
            cat("Problem deconvolving: ", niftis[si], "\n", file=log_file, append=TRUE)
            deconv_mat <- matrix(NA, nrow=nrow(to_deconvolve), ncol=ncol(to_deconvolve))
          } else {
            deconv_mat <- as.matrix(read.table(temp_o, header=FALSE)) %>% unname()  #remove names to avoid confusion in melt

            #NB. 17Apr2019. I modified the compiled C++ program to chop the leading zeros itself for all outputs (rather than leaving the leading hrf_pad)
            #deconv_mat <- deconv_mat[,c(-1*1:hrf_pad, seq(-ncol(deconv_mat), -ncol(deconv_mat)+hrf_pad-1))] #trim leading and trailing padding
          }
        }
      }

      #trim hrf end-padding for both C++ and R variants
      deconv_mat <- deconv_mat[,c(seq(-ncol(deconv_mat), -ncol(deconv_mat)+hrf_pad-1))] #trim trailing padding added above

      #melt this for combination
      deconv_melt <- reshape2::melt(deconv_mat, value.name="decon", varnames=c("vnum", "time"))
      deconv_melt$time <- (deconv_melt$time - 1)*TR + time_offset #convert back to seconds; first volume is time 0

      deconv_df <- deconv_melt %>% 
        dplyr::mutate(vnum=as.numeric(vnum)) %>%
        left_join(a_coordinates, by="vnum") %>%
        #mutate(nifti=niftis[si]) %>%
        dplyr::select(-i, -j, -k) #omitting i, j, k for now

      #add subject metadata, if relevant
      #if (!is.null(add_metadata)) { deconv_df <- deconv_df %>% cbind(this_subj) }

      readr::write_csv(deconv_df, file = out_name)

      #handle output of original time series (before deconvolution) -- mostly useful for debugging
      if (isTRUE(save_original_ts)) {
        to_deconvolve_melt <- reshape2::melt(to_deconvolve, value.name="BOLD_z", varnames=c("vnum", "time"))
        to_deconvolve_melt$time <- (to_deconvolve_melt$time - 1)*TR + time_offset #convert back to seconds; first volume is time 0

        orig_df <- to_deconvolve_melt %>% 
          dplyr::mutate(vnum=as.numeric(vnum)) %>%
          left_join(a_coordinates, by="vnum") %>%
          #mutate(nifti=niftis[si]) %>%
          dplyr::select(-i, -j, -k) #omitting i, j, k for now

        #if (!is.null(add_metadata)) { orig_df <- orig_df %>% cbind(this_subj) }

        #update output file for original
        out_name <- file.path(out_dir, atlas_img_name, "original", paste0(eval(out_file_expression), "_original.csv.gz"))
        readr::write_csv(orig_df, file = out_name)
      }
    }

    #write metadata as single data.frame that can be merged selectively, cutting down on storage demands in individual files
    if (!is.null(add_metadata)) {
      out_name <- file.path(out_dir, atlas_img_name, paste0(atlas_img_name, "_metadata.csv"))
      readr::write_csv(add_metadata, file = out_name)
    }

  }

  return(invisible(NULL))
}
