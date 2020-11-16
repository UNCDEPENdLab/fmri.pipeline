#' function to get interpolated event locked data
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
#' @param tr the repetition time of the scanning sequence, in seconds.
#' @param output_resolution the sampling frequency (in seconds) of the interpolateed data. Defaults to be
#'   the same as \code{tr}.
#' @param logfile Name of log file containing event-locking problems
#' 
#' @importFrom dplyr filter pull bind_rows
#' 
#' @export
event_lock_decon <- function(t_by_value, trial_df, idcol="id", runcol="run", trialcol="run_trial",
                             event="feedback_onset", time_before=-3, time_after=3,
                             tr=1.0, output_resolution=tr, logfile="evtlockerrors.txt") {
  
  stopifnot(time_after > time_before)
  stopifnot(time_before <= 0)
  stopifnot(time_after >= 0)
  checkmate::assert_data_frame(trial_df)
  checkmate::assert_string(event)
  
  if (is.null(runcol)) {
    runcol <- "run"
    trial_df[[runcol]] <- 1
  }

  stopifnot(idcol %in% names(trial_df))
  stopifnot(trialcol %in% names(trial_df))

  if (any(tt <- table(trial_df[[trialcol]]) > 1L)) {
    print(tt[tt > 1L])
    stop("trial_df contains duplicated trials. Cannot figure out how to event-lock based on this.")
  }

  results <- lapply(t_by_value, function(this_data) {
    #subj_df <- trial_df %>% filter(!!sym(idcol)==this_data[[idcol]][1] & !!sym(runcol)==this_data[[runcol]][1])
    trials <- trial_df %>% pull(!!trialcol) %>% sort()

    tlist <- list()
    for (t in trials) {
      evt_onset <- trial_df %>% filter(!!sym(trialcol) == t) %>% pull(!!event)
      stopifnot(length(evt_onset) == 1L) #would be a bizarre and problematic result
      
      if (is.na(evt_onset)) { next } #times are missing for some events on early trials (like rt_vmax_cum)
      this_data$time_map <- this_data$time - evt_onset
      
      #add time on either side of interpolation grid to have preceding time points that inform linear interp
      to_interpolate <- this_data %>% filter(time_map > time_before - 1.5*tr & time_map < time_after + 1.5*tr) %>% 
        select(time_map, decon, vnum) %>% group_by(time_map) %>%
        summarize(vox_sd=sd(decon, na.rm=TRUE), decon=mean(decon, na.rm=TRUE), .groups="drop") #, decon_med=median(decon, na.rm=TRUE))
      
      #for checking heterogeneity
      #ggplot(to_interpolate, aes(x=time_map, y=decon, color=factor(vnum))) + geom_line()

      #rare, but if we have no data at tail end of run, we may not be able to interpolate
      if (nrow(to_interpolate) < 2) {
        cat("For subject:", trial_df[[idcol]][1], ", insufficient interpolation data for run:",
          trial_df[[runcol]][1], ", trial:", t, "\n", file=logfile, append=TRUE)
        next
      }
      
      #put everything onto same time grid
      interp <- approx(x=to_interpolate$time_map, y=to_interpolate$decon, xout=seq(time_before, time_after, by=output_resolution)) %>%
        as.data.frame() %>% setNames(c("evt_time", "decon_interp"))
      interp_sd <- approx(x=to_interpolate$time_map, y=to_interpolate$vox_sd, xout=seq(time_before, time_after, by=output_resolution)) %>%
        as.data.frame() %>% setNames(c("evt_time", "sd_interp"))
      
      interp_all <- interp %>% left_join(interp_sd, by="evt_time") %>% 
        mutate("{trialcol}" := t, atlas_value=this_data$atlas_value[1])
      
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
