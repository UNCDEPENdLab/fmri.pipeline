#' This is a small helper function that returns the inputs provided in the feat_files field for
#' a set of .fsf files. 
#' 
#' @details 
#' One can also pass in .feat or .gfeat directories and the function will
#' use the design.fsf files within each of these to find the inputs
#' 
#' @param feat_list a vector of design.fsf filenames, and/or .feat/.gfeat directory names
#' @param recursive a boolean indicating whether to drill down and find lower-level inputs for .gfeat/.feat inputs
#' @author Michael Hallquist
#' @importFrom checkmate assert_file_exists assert_directory_exists
#' @export
read_feat_inputs <- function(feat_list, recursive=FALSE) {
  is_feat_dir <- grepl(".g?feat$", feat_list, perl=TRUE)
  
  if (any(!is_feat_dir)) {
    check_fsfs <- feat_list[!is_feat_dir]
    sapply(check_fsfs, checkmate::assert_file_exists)
    if (!all(goodfiles <- grepl(".fsf$", check_fsfs, perl=TRUE))) {
      error("The following inputs do not end in .fsf, .gfeat, or .feat: ", paste(check_fsfs[!goodfiles], collapse=", "))
    }
  }
  
  #add the /design.fsf suffix to directory inputs
  if (any(is_feat_dir)) { 
    sapply(feat_list[is_feat_dir], checkmate::assert_directory_exists)
    feat_list[is_feat_dir] <- file.path(feat_list[is_feat_dir], "design.fsf") 
  }
  
  inputs <- lapply(feat_list, function(ff) {
    fsf <- readLines(ff)
    infile <- grep("^set feat_files\\(\\d+\\)", fsf, perl=TRUE, value=TRUE)
    #input <- paste0(sub("set feat_files\\(d+\\) \"([^\"]+)\"", "\\1", nifti, perl=TRUE), ".nii.gz")
    fnumbers <- as.numeric(sub("set feat_files\\((\\d+)\\) \"[^\"]+\"", "\\1", infile, perl=TRUE)) #always return in sorted order
    infile <- sub("set feat_files\\(\\d+\\) \"([^\"]+)\"", "\\1", infile, perl=TRUE)
    return(infile[fnumbers])
  })
  
  inputs_are_featdirs <- grepl(".g?feat$", inputs, perl=TRUE)
  if (isTRUE(recursive) && any(inputs_are_featdirs)) {
    inputs <- lapply(1:length(inputs), function(ii) {
      if (isTRUE(inputs_are_featdirs[ii])) {
        return(read_feat_inputs(inputs[ii], recursive=recursive))
      } else {
        return(list(inputs[ii])) #single-element list
      }
    })
  }
  
  #need some sort of wrap-up function here in the recursive case
  return(inputs)
}
