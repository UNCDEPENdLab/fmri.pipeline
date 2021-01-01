g#' function to get interpolated event locked data
#' @param t_by_value List of time series data where each element is a bin (cutting continuous variable) or
#'   unique mask value (integer mask)
#' @param trial_df Trial-level data frame containing timing of events, used for event-locking time series.
#'   This should be a data frame only for the time series passed in for this subject/run (i.e., \code{t_by_value}).
#' @param idcol Column name of the unique participant identifier in \code{trial_df}.
#' @param runcol Column name of the fMRI run identifier (for multi-run datasets) in \code{trial_df}. If
#'   NULL is passed, a run column with a value of 1 is added for consistency in outputs.
#' @param event Name of column in \code{trial_df} that identifies onset for event of interest
#' @param time_before How many seconds before the \code{event} do we want data. Default: -3
#' @param time_after How many seconds after the \code{event} do we want data. Default: 3
#' @param collide_before An optional vector of column names in \code{trial_df} that set a boundary on the earliest
#'   time point used in interpolation. This effectively truncates the time series for interpolation to a smaller window
#'   than specified by \code{time_before}.
#' @param collide_after An optional vector of column names in \code{trial_df} that set a boundary on the latest
#'   time point used in interpolation. This effectively truncates the time series for interpolation to a smaller window
#'   than specified by \code{time_after}.
#' @param tr the repetition time of the scanning sequence, in seconds.
#' @param output_resolution the sampling frequency (in seconds) of the interpolated data. Defaults to be
#'   the same as \code{tr}.
#' @param logfile Name of log file containing event-locking problems
#' 
#' @importFrom dplyr filter pull bind_rows
#' 
#' @export
event_lock_decon <- function(fmri_obj, event=NULL, time_before=-3, time_after=3,
                             collide_before=NULL, collide_after=NULL, pad_before=-1, pad_after=1) {

  ts_data <- fmri_obj$get_ts()
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
    ts_data$time_map <- ts_data[[ vm[["time"]] ]] - evt_onset

    #add time on either side of interpolation grid to have preceding time points that inform linear interp
    trial_ts <- ts_data %>% filter(time_map > time_before - pad_before & time_map < time_after + pad_after)

    #enforce colliding preceding event
    if (evt_before > -Inf) { trial_ts <- trial_ts %>% filter(!!sym(vm[["time"]]) > evt_before) }
    
    #enforce colliding subsequent event
    if (evt_after < Inf) { trial_ts <- trial_ts %>% filter(!!sym(vm[["time"]]) < evt_after) }
    
    #add these data to trial list
    tlist[[t]] <- trial_ts
  }

  names(tlist) <- trials
  
  #combine rows to stack trial-locked time series for this atlas value
  #tdf <- bind_rows(tlist)
  return(tlist)

}

#' front-end function for taking a list of windowed time series by mask value, interpolating them onto a time grid,
#' and (optionally) averaging across voxels/units within a value to derive the mean interpolated time series
#event_lock_decon <- function(ts_data, run_df,
#                             vm=c(id="id", run="run", trial="trial", run_trial="trial", time="time", atlas_value="atlas_value"),
#                             event=NULL, time_before=-3, time_after=3, collide_before=NULL, collide_after=NULL,
#                             pad_before=-1.5*tr, pad_after=1.5*tr, tr=1.0) {

get_medusa_interpolated_ts <- function(fmri_obj, event=NULL, time_before=-3.0, time_after=3.0, collide_before=NULL, collide_after=NULL,
                                       pad_before=-4.5, pad_after=4.5, output_resolution=1.0,
                                       group_by="trial", logfile="evtlockerrors.txt") {

  stopifnot(inherits(fmri_obj, "fmri_ts"))

  #divide ts_data into usable chunks that can be passed to group_modify

  #align to event of interest
  talign <- event_lock_decon(fmri_obj, event=event, time_before=time_before, time_after=time_after,
    pad_before=pad_before, pad_after=pad_after, collide_before=collide_before, collide_after=collide_after)

  #talign is list of aligned data by trial

  #need to interpolate by cluster variables
  decon_res <- interpolate_decon(talign, time_before=time_before, time_after=time_after, output_resolution=output_resolution)
  
  res <- ts_data %>% group_by(grouping) %>%
    group_modify(interp(.x))
      
  windowed_ts <- event_lock_decon(ts_data
    subj_df, event = evt_col, time_before=time_before, time_after=time_after),

  #interpolate and aggregate
  
}

get_medusa_compression_score <- function() {
  #get event-locked windowed data
  windowed_ts <- event_lock_decon(t_by_value)

  compress_output <- lapply(windowed_ts, function(byval) {
    byval <- byval %>% group_by() %>%
      
  })
}


interpolate_decon <- function(trial_list, time_before=-3, time_after=3, output_resolution=1.0, aggregate=TRUE) {

  ##yy <- yy %>% group_by(across(matches("cluster\\d+")))
  interpolated_result <- lapply(trial_list, function(bytrial) {
    
    #add dummy cluster column if not clustered data
    if (!any(grepl("cluster\\d+", names(bytrial), perl=TRUE))) { bytrial[, cluster1 := 1] }
    
    #convert multivariate value columns to long form
    bytrial <- bytrial %>% pivot_longer(cols=matches("value\\d+"), names_to="key", values_to="value")

    group_string <- c(grep("cluster\\d+", names(bytrial), perl=TRUE, value=TRUE), "key", "time_map")

    to_interpolate <- bytrial %>% 
      select(key, value, time_map, matches("cluster\\d+")) %>% group_by(across(group_string)) %>%
      summarize(vox_sd=sd(value, na.rm=TRUE), decon=mean(value, na.rm=TRUE), .groups="drop") #, decon_med=median(decon, na.rm=TRUE))
    
    #for checking heterogeneity
    #ggplot(to_interpolate, aes(x=time_map, y=decon, color=factor(vnum))) + geom_line()

    #rare, but if we have no data at tail end of run, we may not be able to interpolate
    if (nrow(to_interpolate) < 2) {
      cat("For subject:", trial_df[[idcol]][1], ", insufficient interpolation data for run:",
        trial_df[[runcol]][1], ", trial:", trials[t], "\n", file=logfile, append=TRUE)
      next
    }
    
    #put everything onto same time grid (colliding events generate NAs, which is good)
    interp <- approx(x=to_interpolate$time_map, y=to_interpolate$decon, xout=seq(time_before, time_after, by=output_resolution)) %>%
      as.data.frame() %>% setNames(c("evt_time", "decon_interp"))
    interp_sd <- approx(x=to_interpolate$time_map, y=to_interpolate$vox_sd, xout=seq(time_before, time_after, by=output_resolution)) %>%
      as.data.frame() %>% setNames(c("evt_time", "sd_interp"))


    
  })
  
}

mean


      
      
      interp_all <- interp %>% left_join(interp_sd, by="evt_time") %>% 
        mutate("{trialcol}" := trials[t], atlas_value=this_data$atlas_value[1])
      
      tlist[[t]] <- interp_all
    }
    
    tdf <- bind_rows(tlist)
    return(tdf)
  })

  results <- bind_rows(results) %>%
    mutate("{idcol}" := trial_df[[idcol]][1], #tack on identifying columns
      "{runcol}" := trial_df[[runcol]][1]) %>%
    arrange(!!sym(idcol), !!sym(runcol), atlas_value, !!sym(trialcol), evt_time) %>%
    mutate(event=!!event) %>% select(!!idcol, !!runcol, atlas_value, !!trialcol, event, evt_time, everything())
  
  return(results)
}

#' helper function for computing time series compression
#'
compress_mts_pca <- function(mts, pexp_target=0.9, scale_columns=TRUE) {
  if (isTRUE(scale_columns)) { mts <- scale(mts) } #z-score columns of matrix before compression (standard practice)
  dd <- svd(mts)
  pexp <- cumsum(dd$d^2)/sum(dd$d^2) #proportion of variance explained

  compress_ret <- c()
  for (pp in pexp_target) {
    nexceed <- min(which(pexp > pp))
    napprox <- approx(y=1:length(dd$d), x=pexp, xout=pp)$y #linear approximation of rank, allowing non-integer values

    compress_ret[paste0("compress_", pp)] <- 1-(nexceed/length(dd$d))
    compress_ret[paste0("acompress_", pp)] <- 1-(napprox/length(dd$d))
  }

  return(compress_ret)
}
