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
#' @importFrom dplyr filter pull bind_rows
#' @importFrom checkmate assert_numeric assert_data_frame assert_string assert_subset
#' 
#' @export
event_lock_ts <- function(fmri_obj, event=NULL, time_before=-3, time_after=3,
                          collide_before=NULL, collide_after=NULL, pad_before=-1, pad_after=1) {

  ts_data <- fmri_obj$get_ts(orig_names=TRUE) #so that vm pass-through to returned object works
  run_df <- fmri_obj$event_data
  vm <- fmri_obj$vm
  
  #input validation
  checkmate::assert_numeric(time_before, upper=-1e-2, max.len=1)
  checkmate::assert_numeric(time_after, lower=1e-2, max.len=1)
  stopifnot(time_after > time_before)
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
    trial_ts <- ts_data %>% filter(evt_time > time_before + pad_before & evt_time < time_after + pad_after)

    #enforce colliding preceding event
    if (evt_before > -Inf) { trial_ts <- trial_ts %>% filter(!!sym(vm[["time"]]) > evt_before) }
    
    #enforce colliding subsequent event
    if (evt_after < Inf) { trial_ts <- trial_ts %>% filter(!!sym(vm[["time"]]) < evt_after) }

    #populate trial field
    trial_ts[[ vm[["trial"]] ]] <- trials[t]
    
    #add these data to trial list
    tlist[[t]] <- trial_ts

  }

  names(tlist) <- trials
  
  #combine rows to stack trial-locked time series for this atlas value
  tdf <- bind_rows(tlist)

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
  
  #talign is an fmri_ts object keyed by trial (and other keying variables)
  talign <- event_lock_ts(fmri_obj, event=event, time_before=time_before, time_after=time_after,
    pad_before=pad_before, pad_after=pad_after, collide_before=collide_before, collide_after=collide_after)
  
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
  checkmate::assert_numeric(time_before, upper=-1e-2, max.len=1)
  checkmate::assert_numeric(time_after, lower=1e-2, max.len=1)
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
    thismat <- dcast(dt, ff, value.var="value")
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
  output_dt <- dcast(comp_calc, formula = ... ~ .vkey, value.var=grep("a*compress_.*", names(comp_calc), value=TRUE), sep="|")

  #make names more readable: like decon_compress_0.9, not compress_0.9_decon
  setnames(output_dt, names(output_dt), sub("(a*compress_[^|]+)\\|(.+)", "\\2_\\1", names(output_dt), perl=TRUE))

  return(output_dt)
  
}

#' Interpolate an aligned fmri_ts object onto a consistent time grid with respect to a target event
#'
#' @param a_obj an fmri_ts object that has data aligned to an event of interest
#' @param evt_time
#' @param time_before The earliest time point (in seconds) relative to the event to be output in interpolation
#' @param time_after The latest time point (in seconds) relative to the event to be output in interpolation
#' @param output_resolution The timestep (in seconds) used for interpolation
#' @param group_by a character vector of keying variables used for aggregation of data prior to interpolation
#'
#' @importFrom data.table dcast melt
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

    interp_agg <- to_interpolate %>% group_by(across(evt_time)) %>%
      summarize(across("value", funs, na.rm=TRUE), .groups="drop")

    #rare, but if we have no data at tail end of run, we may not be able to interpolate
    if (nrow(interp_agg) < 2) {
      #cat("For subject:", trial_df[[idcol]][1], ", insufficient interpolation data for run:",
      #  trial_df[[runcol]][1], ", trial:", trials[t], "\n", file=logfile, append=TRUE)
      next
    }
    
    vcols <- grep("value_.*", names(interp_agg), value=TRUE)
    tout <- seq(time_before, time_after, by=output_resolution)
    df <- data.frame(evt_time=tout, sapply(vcols, function(cname) {
      approx(x=interp_agg[[evt_time]], y=interp_agg[[cname]], xout=tout)$y
    }, simplify=FALSE))

    return(df)
  }

  #we need to ensure that multivariate outcomes are stored as a key-value pair of columns, rather than
  # in wide format. This allows .SD to refer to the entire data.table per grouped chunk

  interp_dt <- a_obj$get_ts(orig_names=TRUE) %>%
    melt(measure.vars=a_obj$vm$value, variable.name=".vkey", values.name="value")

  interp_dt <- interp_dt[, interpolate_worker(.SD), by=c(".vkey", group_by)]

  #reshape multiple signals to wide format again
  output_dt <- dcast(interp_dt, formula = ... ~ .vkey, value.var=grep("value_.*", names(interp_dt), value=TRUE), sep="|")

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
  mts <- mts[,!apply(mts, 2, function(x) all(is.na(x)))] #drop columns/time series that are all NA (edge of brain)
  mts <- na.omit(mts) #now drop rows that have NAs (event censoring)

  #screen for lousy time series (no variation or extreme means)
  mms <- apply(mts, 2, mean)
  sds <- apply(mts, 2, sd)

  bad <- sds < 1e-3 | mms < 1e-6 #very small SD or super-low mean
  mts <- mts[,!bad]
  
  #if (isTRUE(scale_columns)) { mts <- scale(mts) } #z-score columns of matrix before compression (standard practice)
  if (isTRUE(scale_columns)) { #for some reason, scale is blowing up svd (infinite/missing values -- just use manual z-scoring)
    mts <- apply(mts, 2, function(x) {  (x - mean(x))/mean(x) })
  }
    
  dd <- svd(mts)
  pexp <- cumsum(dd$d^2)/sum(dd$d^2) #proportion of variance explained

  compress_ret <- list()
  for (pp in pexp_target) {
    nexceed <- min(which(pexp > pp))
    napprox <- approx(y=1:length(dd$d), x=pexp, xout=pp, ties="ordered")$y #linear approximation of rank, allowing non-integer values
    compress_ret[[paste0("compress_", pp)]] <- 1-(nexceed/length(dd$d))
    compress_ret[[paste0("acompress_", pp)]] <- 1-(napprox/length(dd$d))
  }

  compress_ret[["compress_nt"]] <- length(dd$d) #number of timepoints
  return(compress_ret)
}
