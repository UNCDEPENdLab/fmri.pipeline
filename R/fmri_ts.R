library(R6)

#essential structure
# - time data: time, signal, unit (voxel/region) ([voxels + time] x signals)
# - metadata: id, trial, events (trials x signals)
# - time -> unit mapping: vnum, atlas_value, region etc. (units x grouping)





#ts_data
# $time: time index
# $signal: only used if there is more than one value column (multivariate ts)
# $value: DV value of time series
# $cluster: 

#' R6 class representing a multivariate time series object for fMRI analysis
#'
#' @importFrom R6 R6Class
#' @export
fmri_ts <- R6::R6Class("fmri_ts",
  public=list(
    #' @field ts_data time x signals data.table
    ts_data=NULL,

    #' @field ts_keys RLE-encoded clustering/keying variables
    ts_keys=NULL,

    #' @field event_data trial x event data, used for aligning time series with events
    event_data=NULL,

    #' @vm a list of variable mappings between internal constructs and input variable names
    vm=NULL,

    #' @description Create a new fmri_ts object
    #' @param ts_data a data.frame or data.table containing time series
    #' @param event_data a data.frame containing trial-level events that occurred in the time period represented by \code{ts_data}
    #' @param vm a list of variable names used in \code{ts_data} and \code{event_data} that map onto internal constructs
    initialize = function(ts_data=NULL, event_data=NULL, vm=NULL) {
      default_vm <- list(id="id", run="run", trial="trial", run_trial="trial", time="time", value="value", cluster="cluster")
      for (nn in names(default_vm)) { #populate default variable mappings if not provided in input vector
        if (!nn %in% names(vm)) { vm[nn] <- default_vm[nn] }
      }

      checkmate::assert_data_frame(ts_data)
      checkmate::assert_data_frame(event_data)
      #stopifnot(all(vm[c("time", "signal")] %in% names(ts_data)))

      sapply(vm[c("time", "value", "cluster")], function(x) { stopifnot(all(x %in% names(ts_data))) } )
      
      #always make ts_data long/tidy? this is nice, but a) double or quadruples RAM demand, and b) adds compute time
      #if (length(vm[["signal"]] > 1L)) { ts_data <- ts_data %>% tidyr::pivot_longer(cols=vm[["signal"]], names_to="signal", values_to="value") }
      #if (length(vm[["cluster"]] > 1L)) { ts_data <- ts_data2 %>% tidyr::pivot_longer(cols=vm[["cluster"]], names_to="cluster_var", values_to="cluster") }

      #we could convert to data.table, then split. If we drop the key columns, we can get the RAM pressure down 20-30%
      #but, it generates a complicated data structure and slows down per-group calculations
      #dd <- split(xx, by=vm[["cluster"]], keep.by=T)

      #conclusion: for internal object storage, sort by clustering variables, then use RLE encoding of keys to compress object
      #  add method $get_ts that returns the rehydrated data (with inverse.rle)

      ts_data <- data.table(ts_data)
      setorderv(ts_data, vm[["cluster"]])
      ts_keys <- sapply(vm[["cluster"]], function(x) { rle(ts_data[[x]]) }, simplify=FALSE)
      ts_data[,vm[["cluster"]]:=NULL] #drop key columns

      self$ts_data <- ts_data
      self$ts_keys <- ts_keys
      self$event_data <- event_data
      self$vm <- vm
    },

    #helper function to get rehydrated object with key values
    get_ts = function() {
      tsd <- self$ts_data
      for (kk in 1:length(self$ts_keys)) { tsd[, names(self$ts_keys)[kk] := inverse.rle(self$ts_keys[[kk]])] }
      setcolorder(tsd, self$vm[["cluster"]]) #put clustering variables first in object
    },
    export = function(filename) {
      
    }
  )
)
