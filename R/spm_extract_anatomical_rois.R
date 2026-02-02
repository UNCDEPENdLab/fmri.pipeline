#' This is a wrapper around the spm_extract_anatomical_rois.m script in the inst directory
#'
#' @param l1spmdirs character vector of level 1 SPM directories containing SPM.mat files
#' @param masks character vector of NIfTI mask images for each anatomical ROI of interest
#' @param threshold p-value threshold applied to contrast within mask before extraction
#' @param threshdesc multiple comparisons correction on p-value. 'none' or 'FWE'
#' @param session which session (run) to use for extracting time series
#' @param extent exclude clusters having fewer than voxels than extent
#' @param adjust_F_index index of F-test in SPM.mat to adjust for all effects of interest
#' @param contrast_index index of t-test contrast in SPM.mat that is of interest
#' @param ncores number of cores to use in a parallel approach; function parallelizes over \code{l1spmdirs}
#' @param spm_path path to spm12 installation; added to MATLAB path at runtime
#' @param matlab_path location of MATLAB binary; used with matlabr for run_matlab_code()
#' 
#' @importFrom matlabr run_matlab_code
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach registerDoSEQ foreach %dopar%
#' @importFrom iterators iter
#' @importFrom parallel makeCluster stopCluster
#' @importFrom matlabr have_matlab run_matlab_code
#'
#' @author Michael Hallquist
#' 
#' @export
#' 
spm_extract_anatomical_rois <- function(l1spmdirs, masks, threshold=0.2, threshdesc='none', session=1, extent=0,
                                        adjust_F_index=1, contrast_index=NULL, ncores=1,
                                        spm_path="/gpfs/group/mnh5174/default/lab_resources/spm12",
                                        matlab_path="/opt/aci/sw/matlab/R2017b/bin") {
  
  if (missing(l1spmdirs)) { stop("Need to pass in a character vector of all level 1 directories containing SPM.mat files.") }
  stopifnot(all(dir.exists(l1spmdirs)))
  stopifnot(all(file.exists(masks)))
  if (is.null(names(masks))) {
    names(masks) <- sub("(\\.hdr|\\.nii|\\.img)(\\.gz)*$", "", basename(masks), perl=TRUE) #develop basic naming scheme
  }

  stopifnot(dir.exists(spm_path))
  stopifnot(dir.exists(matlab_path))

  #set the matlab path for matlabr
  options(matlab.path=matlab_path)

  if (!have_matlab()) { stop("Unable to find MATLAB installation") }
  
  if (ncores > 1) {
    cl <- makeCluster(ncores)
    on.exit(try(stopCluster(cl)))
    registerDoParallel(cl)
  } else {
    registerDoSEQ()
  }

  nifti_tmpdir <- tempdir()
  gzipped <- grepl(".nii.gz$", masks)
  dir.create(nifti_tmpdir, showWarnings=FALSE)
  masks[gzipped] <- sapply(masks[gzipped], function(x) {
    tmpout <- tempfile(fileext=".nii", tmpdir=nifti_tmpdir)
    system(paste0("gunzip -c ", x, " > ", tmpout))
    return(tmpout)
  })

  #add single quotes around mask strings if not present
  mask_names <- sub("^'?([^']+)'?", "'\\1'", names(masks), perl=TRUE)
  masks <- sub("^'?([^']+)'?", "'\\1'", masks, perl=TRUE)

  #TODO: support multi-session data

  spm_preamble <- c(
    ifelse(is.null(spm_path), "", paste0("addpath('", spm_path, "');")),
    "spm('defaults', 'fmri');",
    "spm_jobman('initcfg');",
    ""    
  )

  #contains matlab scripts for extraction
  matlab_scripts <- system.file("matlab", package = "dependlab")
  stopifnot(dir.exists(matlab_scripts))

  res <- foreach(dd=iter(l1spmdirs), .packages="matlabr") %dopar% {
    #set the matlab path for matlabr in each worker
    options(matlab.path=matlab_path)

    m_string <- c(spm_preamble,
      paste0("addpath('", matlab_scripts, "');"),
      "cfg = struct();",
      "masks = {",
      masks,
      "};",
      "names = {",
      mask_names,
      "};",
      "",
      paste0("threshold = ", threshold, ";"),
      paste0("threshdesc = ", sub("^'?([^']+)'?", "'\\1'", threshdesc, perl=TRUE), ";"),
      paste0("session = ", session, ";"),
      paste0("extent = ", extent, ";"),
      paste0("adjust_F_index = ", ifelse(is.null(adjust_F_index), "[]", adjust_F_index), "; %adjust time series for all effects of interest"),
      paste0("contrast_index = ", ifelse(is.null(contrast_index), "[]", contrast_index), "; %the contrast of interest"),
      "",
      "for jj = 1 : numel(masks)",
      paste0("  cfg(jj).target_dir = fullfile('", dd, "');"),
      "  cfg(jj).mask = masks{jj};",
      "  cfg(jj).adjust = adjust_F_index;",
      "  cfg(jj).session = session;",
      "  cfg(jj).name = names{jj};"
    )

    if (!is.null(contrast_index)) {
      m_string <- c(m_string,
      "  cfg(jj).contrast = contrast_index;",
      "  cfg(jj).threshold = threshold;",
      "  cfg(jj).threshdesc = threshdesc;",
      "  cfg(jj).extent = extent;"
      )
    }

    m_string <- c(m_string,
      "end",
      "spm_extract_anatomical_rois(cfg);"
    )      

    run_matlab_code(m_string, endlines = FALSE, verbose = TRUE, add_clear_all = FALSE)

    #this argument to run_matlab_code adds the paths at the end, not the beginning, which doesn't help
    #, paths_to_add = matlab_scripts)    
  }
}
