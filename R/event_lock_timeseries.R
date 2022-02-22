#' function to get interpolated event locked data
#' @param fmri_obj An fmri_ts object containing time series data that are optionally keyed on one or more
#'   grouping variables (e.g., ROI). The fmri_obj must also contain trial-level data (in the $event_data field)
#'   for event-locking to proceed.
#' @param event Name of column in \code{fmri_obj$event_data} that identifies onset for event of interest
#' @param time_before How many seconds before the \code{event} do we want data. Default: -3
#' @param time_after How many seconds after the \code{event} do we want data. Default: 3
#' @param collide_before An optional vector of column names in \code{trial_df} that set a boundary on the earliest
#'   time point used in interpolation. This effectively truncates the time series for interpolation to a smaller window
#'   than specified by \code{time_before}.
#' @param collide_after An optional vector of column names in \code{trial_df} that set a boundary on the latest
#'   time point used in interpolation. This effectively truncates the time series for interpolation to a smaller window
#'   than specified by \code{time_after}.
#' @param pad_before Number of seconds to include in the epoch time window before the event of interest. Interpolation
#'   spans the window from \code{time_before} to \code{time_after}, but padding includes data points at the
#'   boundary that can help to have sufficient data to interpolate early and late times within the epoch.
#' @param pad_after Number of seconds to include in the epoch time window after the event of interest.
#' @param logfile Name of log file containing event-locking problems
#'
#' @importFrom dplyr filter pull
#' @importFrom data.table rbindlist
#' @importFrom checkmate assert_numeric assert_data_frame assert_string assert_subset
#'
#' @export
event_lock_ts <- function(fmri_obj, event=NULL, time_before=-3, time_after=3,
                          collide_before=NULL, collide_after=NULL, pad_before=-1, pad_after=1) {

  ts_data <- fmri_obj$get_ts(orig_names=TRUE) #so that vm pass-through to returned object works
  run_df <- fmri_obj$event_data
  vm <- fmri_obj$vm

  #input validation
  checkmate::assert_numeric(time_before, upper=0, max.len=1)
  checkmate::assert_numeric(time_after, lower=0, max.len=1)
  stopifnot(time_after >= time_before)
  checkmate::assert_data_frame(ts_data)
  checkmate::assert_data_frame(run_df)
  checkmate::assert_string(event)
  stopifnot(event %in% names(run_df))
  sapply(collide_before, checkmate::assert_string, null.ok=TRUE)
  sapply(collide_after, checkmate::assert_string, null.ok=TRUE)
  if (!is.null(collide_before)) { checkmate::assert_subset(collide_before, names(run_df)) }
  if (!is.null(collide_after)) { checkmate::assert_subset(collide_after, names(run_df)) }
  checkmate::assert_numeric(pad_before, upper=0, null.ok=FALSE)
  checkmate::assert_numeric(pad_after, lower=0, null.ok=FALSE)

  if (any(tt <- table(run_df[[ vm[["trial"]] ]]) > 1L)) {
    print(tt[tt > 1L])
    stop("run_df contains duplicated trials. Cannot figure out how to event-lock based on this.")
  }

  trials <- run_df[[ vm[["trial"]] ]] %>% sort() #sorted vector of trials

  #list containing event-locked data frames, one per trial
  tlist <- list()

  #loop over trials, retaining data within window around event of interest
  for (t in 1:length(trials)) {

    evt_onset <- run_df %>% filter(!!sym(vm[["trial"]]) == trials[t]) %>% pull(!!event)
    stopifnot(length(evt_onset) == 1L) #would be a bizarre and problematic result

    evt_before <- -Inf
    if (!is.null(collide_before)) {
      evts <- run_df %>% select(all_of(collide_before)) %>% unlist() %>% unname()
      if (any(evts < evt_onset)) { evt_before <- max(evts[evts < evt_onset]) } #most recent past event
    }

    evt_after <- Inf
    if (!is.null(collide_after)) {
      evts <- run_df %>% select(all_of(collide_after)) %>% unlist() %>% unname()
      if (any(evts > evt_onset)) { evt_after <- min(evts[evts > evt_onset]) } #nearest future event
    }

    if (is.na(evt_onset)) { next } #times are missing for some events on early trials (like rt_vmax_cum)
    ts_data$evt_time <- ts_data[[ vm[["time"]] ]] - evt_onset

    #add time on either side of interpolation grid to have preceding time points that inform linear interp
    #N.B. Occasionally, if the time series (fMRI run) ends before the late events, trial_ts will have 0 rows
    trial_ts <- ts_data %>% filter(evt_time >= time_before + pad_before & evt_time <= time_after + pad_after)

    #enforce colliding preceding event
    if (evt_before > -Inf) { trial_ts <- trial_ts %>% filter(!!sym(vm[["time"]]) > evt_before) }

    #enforce colliding subsequent event
    if (evt_after < Inf) { trial_ts <- trial_ts %>% filter(!!sym(vm[["time"]]) < evt_after) }

    #populate trial field
    trial_ts[[ vm[["trial"]] ]] <- if (nrow(trial_ts) > 0L) { trials[t] } else { numeric(0) } #allow for a zero-row df

    #add these data to trial list
    tlist[[t]] <- trial_ts

  }

  names(tlist) <- trials

  #combine rows to stack trial-locked time series for this atlas value
  tdf <- data.table::rbindlist(tlist)

  #convert to a keyed fmri_ts object for consistency
  res <- fmri_ts$new(ts_data=tdf, event_data=fmri_obj$event_data, vm=fmri_obj$vm, tr=fmri_obj$tr)
  res$replace_vm(time="evt_time")

  #add trial as keying variable
  res$add_keys(vm[["trial"]])

  return(res)

}

#' front-end function for taking a list of windowed time series by mask value, interpolating them onto a time grid,
#' and (optionally) averaging across voxels/units within a value to derive the mean interpolated time series
#'
#' @param fmri_obj an fmri_ts object containing a single run of data with corresponding events
#' @param event the event to which the fmri time series should be aligned (column in \code{fmri_obj$event_data})
#' @param time_before How many seconds before the \code{event} do we want data. Default: -3
#' @param time_after How many seconds after the \code{event} do we want data. Default: 3
#' @param collide_before An optional vector of column names in \code{trial_df} that set a boundary on the earliest
#'   time point used in interpolation. This effectively truncates the time series for interpolation to a smaller window
#'   than specified by \code{time_before}.
#' @param collide_after An optional vector of column names in \code{trial_df} that set a boundary on the latest
#'   time point used in interpolation. This effectively truncates the time series for interpolation to a smaller window
#'   than specified by \code{time_after}.
#' @param pad_before Number of seconds to include in the epoch time window before the event of interest. Interpolation
#'   spans the window from \code{time_before} to \code{time_after}, but padding includes data points at the
#'   boundary that can help to have sufficient data to interpolate early and late times within the epoch.
#' @param pad_after Number of seconds to include in the epoch time window after the event of interest.
#' @param output_resolution the sampling frequency (in seconds) of the interpolated data. Defaults to be
#'   the same as \code{tr}.
#' @param group_by return interpolated time series for each combination of group_by variables. Default is
#'   to provide one interpolated time series per trial.
#'
#' @author Michael Hallquist
#' @export
get_medusa_interpolated_ts <- function(fmri_obj, event=NULL, time_before=-3.0, time_after=3.0,
                                       collide_before=NULL, collide_after=NULL,
                                       pad_before=-1.5, pad_after=1.5, output_resolution=NULL,
                                       group_by="trial") {

  stopifnot(inherits(fmri_obj, "fmri_ts"))

  #make default output_resolution equal to fmri TR
  if (is.null(output_resolution)) { output_resolution <- fmri_obj$tr }

  # talign is an fmri_ts object keyed by trial (and other keying variables)
  talign <- event_lock_ts(fmri_obj,
    event = event, time_before = time_before, time_after = time_after,
    pad_before = pad_before, pad_after = pad_after, collide_before = collide_before, collide_after = collide_after
  )
  
  #need to interpolate by key variables
  interpolated_epochs <- interpolate_fmri_epochs(talign, time_before=time_before, time_after=time_after,
    output_resolution=output_resolution, group_by=group_by)  

  return(interpolated_epochs)
}

#' front-end function for taking a list of windowed time series by mask value, interpolating them onto a time grid,
#' and (optionally) averaging across voxels/units within a value to derive the mean interpolated time series
#'
#' @param fmri_obj an fmri_ts object containing a single run of data with corresponding events
#' @param event the event to which the fmri time series should be aligned (column in \code{fmri_obj$event_data})
#' @param time_before How many seconds before the \code{event} do we want data. Default: -3
#' @param time_after How many seconds after the \code{event} do we want data. Default: 3
#' @param collide_before An optional vector of column names in \code{trial_df} that set a boundary on the earliest
#'   time point used in interpolation. This effectively truncates the time series for interpolation to a smaller window
#'   than specified by \code{time_before}.
#' @param collide_after An optional vector of column names in \code{trial_df} that set a boundary on the latest
#'   time point used in interpolation. This effectively truncates the time series for interpolation to a smaller window
#'   than specified by \code{time_after}.
#' @param pad_before Number of seconds to include in the epoch time window before the event of interest. Interpolation
#'   spans the window from \code{time_before} to \code{time_after}, but padding includes data points at the
#'   boundary that can help to have sufficient data to interpolate early and late times within the epoch.
#' @param pad_after Number of seconds to include in the epoch time window after the event of interest.
#' @param output_resolution the sampling frequency (in seconds) of the interpolated data. Defaults to be
#'   the same as \code{tr}.
#' @param group_by return interpolated time series for each combination of group_by variables. Default is
#'   to provide one interpolated time series per trial.
#'
#' @importFrom checkmate assert_string assert_numeric assert_character
#' @importFrom data.table dcast melt
#' @author Michael Hallquist
#' @export
get_medusa_compression_score <- function(fmri_obj, event=NULL, time_before=-3, time_after=3,
                                         collide_before=NULL, collide_after=NULL, group_by=NULL) {

  #input validation
  checkmate::assert_string(event)
  checkmate::assert_numeric(time_before, upper=0, max.len=1)
  checkmate::assert_numeric(time_after, lower=0, max.len=1)
  checkmate::assert_character(collide_before, null.ok=TRUE, unique=TRUE)
  checkmate::assert_character(collide_after, null.ok=TRUE, unique=TRUE)

  #get event-locked windowed data
  stopifnot(inherits(fmri_obj, "fmri_ts"))

  #talign is an aligned fmri_ts object keyed by trial (and other keying variables)
  talign <- event_lock_ts(fmri_obj, event=event, time_before=time_before, time_after=time_after,
    pad_before=0, pad_after=0, collide_before=collide_before, collide_after=collide_after)

  interp_dt <- talign$get_ts(orig_names=TRUE) %>%
    melt(measure.vars=talign$vm$value, variable.name=".vkey", values.name="value")

  ## yy <- dcast(xx, evt_time ~ vnum, value.var="value")
  ## yy[,evt_time:=NULL]
  ## zz <- as.matrix(yy)
  ## compress_mts_pca(zz, pexp_target=0.9, scale_columns=TRUE)

  #small wrapper function to cast a given matrix into 
  calculate_compression <- function(dt, keys=NULL, ...) {
    #convert to wide format with multiple signals (e.g., from many voxels) as columns
    ff <- as.formula(paste0("evt_time ~ ", paste(keys, collapse=" + ")))
    thismat <- data.table::dcast(dt, ff, value.var="value")
    thismat[, evt_time := NULL] #omit time column itself from matrix
    thismat <- as.matrix(thismat)
    compress_mts_pca(thismat, ...)
  }

  #find the key variables that are not in the group_by
  other_keys <- setdiff(talign$vm$key, group_by)

  #so, on the columns of the dcast, it's all keys that are not in the group_by... making them simultaneous columns
  #xyz <- interp_dt[, compress_mts_pca(.SD, pexp_target=0.9, scale_columns=TRUE), by=c(".vkey", group_by)]
  comp_calc <- interp_dt[, calculate_compression(.SD, keys=other_keys, pexp_target=0.9, scale_columns=TRUE), by=c(".vkey", group_by)]

  #reshape multiple signals to wide format again
  output_dt <- data.table::dcast(comp_calc, formula = ... ~ .vkey, value.var=grep("a*compress_.*", names(comp_calc), value=TRUE), sep="|")

  #make names more readable: like decon_compress_0.9, not compress_0.9_decon
  setnames(output_dt, names(output_dt), sub("(a*compress_[^|]+)\\|(.+)", "\\2_\\1", names(output_dt), perl=TRUE))

  return(output_dt)

}

#' Interpolate an aligned fmri_ts object onto a consistent time grid with respect to a target event
#'
#' @param a_obj an fmri_ts object that has data aligned to an event of interest
#' @param evt_time column name in \code{a_obj} event data to which time series should be aligned and interpolated
#' @param time_before The earliest time point (in seconds) relative to the event to be output in interpolation
#' @param time_after The latest time point (in seconds) relative to the event to be output in interpolation
#' @param output_resolution The timestep (in seconds) used for interpolation
#' @param group_by a character vector of keying variables used for aggregation of data prior to interpolation
#'
#' @importFrom data.table dcast melt setnames
#' @importFrom dplyr summarize group_by across
#' @author Michael Hallquist
#' @keywords internal
interpolate_fmri_epochs <- function(a_obj, evt_time="evt_time", time_before=-3, time_after=3,
                                    output_resolution=1.0, group_by=NULL) {

  #basically need to split align_obj on keying units, then interpolate within keys

  stopifnot(inherits(a_obj, "fmri_ts"))

  #worker to apply a set of summary functions to each timepoint within a given cluster/unit
  interpolate_worker <- function(to_interpolate, funs=list(mean=mean, median=median, sd=sd)) {
    #For now, we only support linear interpolation. Because of this, interpolation of a mean time
    #  series is the same is the mean of interpolated time series. But the former is faster to compute.
    #  https://math.stackexchange.com/questions/15596/mean-of-interpolated-data-or-interpolation-of-means-in-geostatistics

    interp_agg <- to_interpolate %>% 
      group_by(across(evt_time)) %>%
      dplyr::summarize(across("value", funs, na.rm=TRUE), .groups="drop")

    #rare, but if we have no data at tail end of run, we may not be able to interpolate
    if (nrow(interp_agg) < 2) {
      #cat("For subject:", trial_df[[idcol]][1], ", insufficient interpolation data for run:",
      #  trial_df[[runcol]][1], ", trial:", trials[t], "\n", file=logfile, append=TRUE)
      return(NULL)
    }

    vcols <- grep("value_.*", names(interp_agg), value=TRUE)
    tout <- seq(time_before, time_after, by=output_resolution)
    df <- data.frame(evt_time=tout, sapply(vcols, function(cname) {
      if (all(is.na(interp_agg[[cname]]))) {
        rep(NA_real_, length(tout)) #if all inputs are NA, do not try to interpolate
      } else {
        approx(x=interp_agg[[evt_time]], y=interp_agg[[cname]], xout=tout)$y
      }
    }, simplify=FALSE))

    return(df)
  }

  #we need to ensure that multivariate outcomes are stored as a key-value pair of columns, rather than
  # in wide format. This allows .SD to refer to the entire data.table per grouped chunk

  interp_dt <- a_obj$get_ts(orig_names=TRUE) %>%
    melt(measure.vars=a_obj$vm$value, variable.name=".vkey", values.name="value")

  interp_dt <- interp_dt[, interpolate_worker(.SD), by=c(".vkey", group_by)]

  #reshape multiple signals to wide format again
  output_dt <- data.table::dcast(interp_dt, formula = ... ~ .vkey, value.var=grep("value_.*", names(interp_dt), value=TRUE), sep="|")

  #make names more readable: like decon_mean, not value_mean_decon
  setnames(output_dt, names(output_dt), sub("value_([^|]+)\\|(.+)", "\\2_\\1", names(output_dt), perl=TRUE))

  return(output_dt)
}


#' Core function for computing multivariate time series compression scores by principal components
#'
#' @param mts multivariate time series structured time x signals
#' @param pexp_target Proportion of variance explained by principal components. This can be a vector,
#'   in which case compression is calculated at different thresholds.
#' @param scale_columns whether to z-score the time series prior to eigendecomposition (recommended)
#'
#' @details This function accepts a time x signals (e.g., voxels) time series matrix. It computes the
#'   singular value decomposition (SVD) and then examines how many eigenvectors are needed to explain at least
#'   \code{pexp_target} proportion of variance.
#'
#' Compression scores are normalized 0 -- 1.0 by the equation: 1 - (n_components / n_timeseries). Thus,
#'   if 6 components explain 91% of the variance in a 12-time series matrix, then the compression score is 0.5.
#'
#' Given that the number of eigenvectors is an integer, a linear approximation to the exact proportion of variance
#'   explained is also calculated. For example, if 3 components explain 84% of variance, while 4 components explain
#'   92%, we 'overshoot' our target of 90% slightly with 4 components. We can then calculated how much compression would
#'   get us to precisely 90% by linear interpolation.
#'
#' @return a list containing compression estimates of the matrix. For each pexp_target value, two values are included,
#'  one representing the compression calculated using integer
#' @importFrom checkmate assert_matrix assert_numeric assert_logical
#' @author Michael Hallquist
#' @export
compress_mts_pca <- function(mts, pexp_target=0.9, scale_columns=TRUE) {
  checkmate::assert_matrix(mts)
  checkmate::assert_numeric(pexp_target, lower=1e-2, upper=1.0, unique=TRUE)
  checkmate::assert_logical(scale_columns, max.len=1L)
  
  #orig <- mts
  mts <- mts[, !apply(mts, 2, function(x) all(is.na(x))), drop=FALSE] #drop columns/time series that are all NA (edge of brain)
  mts <- na.omit(mts) #now drop rows that have NAs (event censoring)

  #screen for lousy time series (no variation or extreme means)
  mms <- apply(mts, 2, mean)
  sds <- apply(mts, 2, sd)

  bad <- sds < 1e-3 | mms < 1e-6 #very small SD or super-low mean
  mts <- mts[, !bad, drop=FALSE]

  if (nrow(mts) <= 3 || ncol(mts) <= 1) {
    warning("mts matrix cannot be compressed. Dims: ", paste(dim(mts), collapse=", "))
    compress_ret <- list() #setup dummy list
    for (pp in pexp_target) {
      compress_ret[[paste0("compress_", pp)]] <- NA_real_
      compress_ret[[paste0("acompress_", pp)]] <- NA_real_
    }

    compress_ret[["compress_nt"]] <- NA_integer_
  } else {
    #if (isTRUE(scale_columns)) { mts <- scale(mts) } #z-score columns of matrix before compression (standard practice)
    #for some reason, scale is blowing up svd (infinite/missing values -- just use manual z-scoring)
    if (isTRUE(scale_columns)) { 
      mts <- apply(mts, 2, function(x) {  (x - mean(x))/mean(x) })
    }

    dd <- svd(mts)
    pexp <- cumsum(dd$d^2)/sum(dd$d^2) #proportion of variance explained

    compress_ret <- list()
    for (pp in pexp_target) {
      nexceed <- min(which(pexp > pp))

      #linear approximation of rank, allowing non-integer values
      napprox <- approx(y=1:length(dd$d), x=pexp, xout=pp, ties="ordered")$y
      compress_ret[[paste0("compress_", pp)]] <- 1-(nexceed/length(dd$d))
      compress_ret[[paste0("acompress_", pp)]] <- 1-(napprox/length(dd$d))
    }

    compress_ret[["compress_nt"]] <- length(dd$d) #number of timepoints

  }

  return(compress_ret)
}


#' Align a set of deconvolved time series files from \code{voxelwise_deconvolution} to a set of events in a task-related design
#' @param atlas_files A character vector containing filenames of atlases that used for aggregating voxel-level deconvolved time series
#' @param decon_dir The directory containing the deconvolved data from \code{voxelwise_deconvolution}. Should be the same as
#'   \code{out_dir} from \code{voxelwise_deconvolution}.
#' @param alignments A list of alignment specifications, each of which will be computed by this function.
#'   See \code{get_medusa_interpolated_ts} for help with arguments, and see Details.
#' @param nbins For atlases with continuous values, how many bins should be used to discretize the mask values, leading to aggregation
#'   by bin.
#' @param overwrite if TRUE, overwrite existing output files
#' @param tr The repetition time of the scan in seconds.
#' @param ncpus The number of cores to use for each parallel job
#' @param walltime The time requested for each event alignment job. Default: 1:00:00 (1 hour).
#'
#' @details
#'   This function will create a separate R batch job on the scheduler for each combination of atlas files and alignments. Each of these
#'   jobs can be allocated a number of cores, which will be used to loop over the deconvolved time series files for event alignment. The
#'   number of cores requested for each alignment job is specified using \code{ncpus}, which defaults to 8.
#'
#'   The alignments argument is a list containing information about how to event-align the deconvolved time series. Elements include:
#'
#' \itemize{
#'   \item \code{$evt_col}: (Required) The name of the column in \code{trial_df} containing the event of interest.
#'   \item \code{$time_before}: (Required) The number of seconds before the event of interest to include in the interpolation window.
#'   \item \code{$time_after}: (Required) The number of seconds after the event of interest to include in the interpolation window.
#'   \item \code{$pad_before}: (Optional) The number of seconds before the earliest event to include in the interpolation window
#'      to avoid missing values at the edge of the interpolation window. Default: -1.5
#'   \item \code{$pad_after}: (Optional) The number of seconds after the latest event to include in the interpolation window
#'      to avoid missing values at the edge of the interpolation window. Default: 1.5
#' }
#' 
#' @examples 
#'
#' \dontrun{
#'   atlas_files <- c(
#'     "/proj/mnhallqlab/projects/clock_analysis/fmri/ph_da_striatum/masks/bilateral_striatum_tight_7Networks_2.3mm.nii.gz",
#'     "/proj/mnhallqlab/projects/clock_analysis/fmri/ph_da_striatum/masks/pauli_combined_integermask_2.3mm.nii.gz"
#'   )
#'
#'   decon_dir <- "/proj/mnhallqlab/users/michael/sceptic_decon" # has the outputs of voxelwise_deconvolution for these atlases
#'   trial_df <- get_trial_data(repo_directory = "/proj/mnhallqlab/projects/clock_analysis", dataset = "mmclock_fmri", groupfixed = TRUE) 
#'
#'   alignments <- list(
#'     clock_long = list(
#'       evt_col = "clock_onset",
#'       time_before = -5,
#'       time_after = 10,
#'       collide_before = "iti_onset", # censor data if we bump into the end of the prior trial
#'       collide_after = "clock_onset" # censor data it we hit the next trial
#'     ),
#'     rt_long = list(
#'       evt_col = "rt_time",
#'       time_before = -4,
#'       time_after = 7,
#'       collide_before = "iti_onset",
#'       collide_after = "clock_onset"
#'     )
#'   )
#'
#'   # run all atlases and alignments as separate slurm jobs (2 hours each, 8 cpus) 
#'   run_decon_alignment(atlas_files, decon_dir, trial_df, alignments, 
#'     overwrite = TRUE, tr = 1.0, ncpus = 8, walltime = "2:00:00", scheduler = "slurm")
#'
#' }
#' @export
run_decon_alignment <- function(atlas_files, decon_dir, trial_df, alignments = list(), nbins = 12, overwrite = FALSE, tr = NULL,
                               ncpus = 8, walltime = "1:00:00", scheduler = "slurm") {
  checkmate::assert_file_exists(atlas_files)
  checkmate::assert_directory_exists(decon_dir)
  checkmate::assert_data_frame(trial_df)
  checkmate::assert_list(alignments, names = "named") # require named list input
  checkmate::assert_integerish(nbins, lower = 1, upper = 1e4, len = 1L)
  checkmate::assert_logical(overwrite, len = 1L)
  checkmate::assert_number(tr, lower = 0.01, upper = 1000)
  checkmate::assert_integerish(ncpus, lower = 1, upper = 1e3, len = 1L)
  walltime <- validate_dhms(walltime) # validate and correct ambiguities
  checkmate::assert_string(scheduler)
  checkmate::assert_subset(scheduler, c("local", "slurm", "torque"))

  # validate padding and alignment settings, input defaults
  for (aa in seq_along(alignments)) {
    if (!"pad_before" %in% names(alignments[[aa]])) {
      message(glue("For alignment {names(alignments)[aa]}, defaulting pad_before to -1.5"))
      alignments[[aa]]$pad_before <- -1.5
    } else {
      checkmate::assert_number(alignments[[aa]]$pad_before, upper = 0)
    }

    if (!"pad_after" %in% names(alignments[[aa]])) {
      message(glue("For alignment {names(alignments)[aa]}, defaulting pad_after to 1.5"))
      alignments[[aa]]$pad_after <- 1.5
    } else {
      checkmate::assert_number(alignments[[aa]]$pad_after, lower = 0)
    }

    if (!"output_resolution" %in% names(alignments[[aa]])) {
      message(glue("For alignment {names(alignments)[aa]}, defaulting output_resolution (sampling rate) of 1.0"))
      alignments[[aa]]$output_resolution <- 1.0
    } else {
      checkmate::assert_number(alignments[[aa]]$output_resolution, lower = 0.01, upper = 1e2)
    }

    has_ingredients <- c("pad_before", "pad_after", "evt_col", "time_before", "time_after", "output_resolution") %in% names(alignments[[aa]])
    if (!all(has_ingredients)) {
      stop(glue("For alignment {names(alignments)[aa]}, missing required fields: {paste(names(alignments[[aa]]), collapse=', ')}"))
    }

    stopifnot(alignments[[aa]]$evt_col %in% names(trial_df))
  }

  # loop over atlases and alignments
  for (af in atlas_files) {
    if (!checkmate::test_file_exists(af)) {
      message("Cannot find input atlas file: {af}. Skipping alignment for this atlas.")
      next
    }

    aname <- file_sans_ext(basename(af)) # atlas name

    mask <- RNifti::readNifti(af)
    mi <- which(mask > 0, arr.ind = TRUE)
    maskvals <- unique(mask[mi])

    if (isTRUE(checkmate::test_integerish(maskvals))) {
      continuous <- FALSE
      message("Atlas: ", aname, " has integer values. Using unique mask values for aggregation in event-aligning.")
      atlas_cuts <- NULL
    } else {
      continuous <- TRUE
      message("Atlas: ", aname, "has continuous values. Dividing into ", nbins, "bins")
      atlas_cuts <- seq(min(mask[mi]) - 1e-5, max(mask[mi]) + 1e-5, length.out = nbins + 1)
    }

    adir <- file.path(decon_dir, aname, "deconvolved")
    if (!dir.exists(adir)) {
      message(glue("Cannot find expected deconvolved output directory: {adir}. Have you run voxelwise_deconvolution yet?"))
      next
    }

    d_files <- list.files(path = adir, pattern = "_deconvolved\\.csv\\.gz", full.names = TRUE)
    if (length(d_files) == 0L) {
      message(glue("No deconvolved files found in directory: {adir}"))
      next
    }

    # loop over event alignments within this atlas
    for (ee in names(alignments)) {
      out_file <- file.path(decon_dir, aname, glue("{aname}_{ee}_decon_aligned.csv.gz"))
      this_alignment <- alignments[[ee]]
      this_alignment$aname <- aname # populate atlas name so that it can be added in downstream alignment function

      if (file.exists(out_file) && isFALSE(overwrite)) {
        message("Output file already exists: ", out_file)
        next
      }

      # mask <- list(continuous = continuous, atlas_cuts = atlas_cuts, nifti = mask)

      d_batch <- R_batch_job$new(
        job_name = glue("evtalign_{aname}_{ee}"), n_cpus = ncpus, mem_per_cpu = "4g",
        wall_time = walltime, scheduler = scheduler,
        input_objects = named_list(d_files, trial_df, this_alignment, tr, atlas_cuts, out_file, ncpus), # pass relevant vars to the batch
        r_packages = "fmri.pipeline",
        r_code = c(
          "evt_align_decon_files(d_files, trial_df, this_alignment, tr, atlas_cuts, out_file, ncpus)"
        )
      )

      d_batch$submit()
    }
  }
}

#' Align a set of deconvolved time series files to an event of interest function.
#' @details
#'   This is intended to be used internally by \code{run_decon_alignment}, which accepts a set of mask/atlas files
#'   and alignments, then processes these in parallel.
#' @param d_files A vector of deconvolved .csv.gz files created by \code{voxelwise_deconvolution}.
#' @param trial_df The trial-level data.frame containing id and run for each subject represented in \code{d_files}.
#' @param alignment A list containing alignment details passed to get_medusa_interpolated_ts
#' @param tr The repetition time of the sequence in seconds.
#' @param atlas_cuts For a continuous-valued atlas (e.g., containing a gradient of interest), a vector of cut points for binning values
#' @param mask A list of mask-related information, including ...
#' @param output_dir The output directory for the event-aligned csv file. If NULL, nothing it output
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach registerDoSEQ
#' @export
evt_align_decon_files <- function(d_files, trial_df, alignment = list(), tr = NULL, atlas_cuts, out_file = NULL, ncpus = 8) {
  checkmate::assert_character(d_files)
  checkmate::assert_file_exists(d_files)
  checkmate::assert_data_frame(trial_df)
  checkmate::assert_list(alignment)
  checkmate::assert_number(tr, lower = 0.01, upper = 1000)
  checkmate::assert_numeric(atlas_cuts, null.ok = TRUE)
  checkmate::assert_string(out_file, null.ok = TRUE)

  if (ncpus > 1) {
    clusterobj <- parallel::makeCluster(ncpus)
    on.exit(try(stopCluster(clusterobj)))
    registerDoParallel(clusterobj)
  } else {
    registerDoSEQ()
  }

  elist <- foreach(fname = iter(d_files), .packages = c("dplyr", "readr", "data.table", "fmri.pipeline")) %dopar% {

    # add sub and run for now since I screwed this up in the outputs...
    id <- as.numeric(sub("^.*/sub(\\d+)_.*", "\\1", fname))
    run <- as.numeric(sub("^.*/sub\\d+_run(\\d+).*", "\\1", fname))
    d <- read_csv(fname) %>% dplyr::select(-atlas_name, -x, -y, -z)
    if (all(is.na(d$decon))) {
      cat(glue("For file {fname}, all decon values are NA, suggesting a failue in deconvolution. Returning NULL."))
      return(NULL)
    }

    # discretize atlas value into bins (continuous) or unique values (integer mask)
    if (!is.null(atlas_cuts)) {
      d <- d %>% dplyr::mutate(atlas_value = cut(atlas_value, atlas_cuts))
    } else {
      # d <- d %>% mutate(atlas_value=factor(atlas_value))
      d <- d %>% dplyr::mutate(atlas_value = as.numeric(atlas_value))
    }

    # run data
    subj_df <- trial_df %>% filter(id == !!id & run == !!run)

    if (nrow(subj_df) == 0L) {
      cat("Cannot find trial-level data for id {id}, run {run}.\n")
      return(NULL)
    }

    tsobj <- fmri_ts$new(
      ts_data = d, event_data = subj_df, tr = tr,
      vm = list(value = c("decon"), key = c("vnum", "atlas_value"))
    )

    subj_align <- tryCatch(
      get_medusa_interpolated_ts(tsobj,
        event = alignment$evt_col,
        time_before = alignment$time_before, time_after = alignment$time_after,
        collide_before = alignment$collide_before, collide_after = alignment$collide_after,
        pad_before = alignment$pad_before, pad_after = alignment$pad_after, output_resolution = alignment$output_resolution,
        group_by = c("atlas_value", "trial")
      ), # one time series per region and trial
      error = function(err) {
        cat("Problems with event aligning ", fname, " for event: ", e, "\n  ",
          as.character(err), "\n\n",
          file = "evt_align_errors.txt", append = TRUE
        )
        return(NULL)
      }
    )

    if (!is.null(subj_align)) {
      # tack on run and id to interpolated time series object
      subj_align[, id := id]
      subj_align[, run := run]
    }

    return(list(ts = subj_align))
  }

  all_e <- dplyr::bind_rows(lapply(elist, "[[", "ts"))
  all_e$atlas <- alignment$aname
  if (!is.null(out_file)) {
    message("Writing output: ", out_file)
    readr::write_csv(all_e, file = out_file)
  }

  return(all_e)
}
