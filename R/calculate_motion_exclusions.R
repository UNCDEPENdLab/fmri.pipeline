#' Function to calculate exclusions of runs and subjects based on excessive head movement
#' 
#' @param gpa a \code{glm_pipeline_arguments} object, as specified by glm_setup_pipeline
#' @return a modified version of the \code{gpa} object that contains a $motion_exclusions data.frame
#' @keywords internal
#' @author Michael Hallquist
calculate_motion_exclusions <- function(gpa) {

  
  cat("Excluding runs exceeding 10% frames with FD >= 0.9mm OR any movement > 5mm\n")
  # read each FD file, index it based on modeled fMRI data, then return FD statistics
  motexclude <- plyr::ldply(1:nrow(feat_runs), function(i) {
    fd <- read.table(feat_runs$fd_file[i], header = FALSE)$V1
    fd <- fd[(feat_runs$drop_volumes[i] + 1):feat_runs$trunc_lengths[i]] # only include volumes within run
    propSpikes_0p9 <- sum(as.integer(fd > 0.9)) / length(fd)
    spikeExclude <- if (propSpikes_0p9 > .10) 1 else 0
    maxFD <- max(fd)
    meanFD <- mean(fd)
    maxMotExclude <- if (maxFD > 5) 1 else 0
    data.frame(fd_file = feat_runs$fd_file[i], propSpikes_0p9, spikeExclude, meanFD, maxFD, maxMotExclude, stringsAsFactors = FALSE)
  })


}