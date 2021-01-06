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
  private=list(
    cvars=NULL, #internal names of cluster variables
    vvars=NULL, #internal names of value variables
    vmvec=NULL  #mapping between input names and internal names
  ),
  public=list(
    #' @field ts_data time x signals data.table
    ts_data=NULL,

    #' @field ts_keys RLE-encoded clustering/keying variables
    ts_keys=NULL,

    #' @field event_data trial x event data, used for aligning time series with events
    event_data=NULL,

    #' @field vm a list of variable mappings between internal constructs and input variable names
    vm=NULL,

    #' @field tr the repetition time (TR) of the fMRI sequence in seconds
    tr=NULL,

    #' @description Create a new fmri_ts object
    #' @param ts_data a data.frame or data.table containing time series
    #' @param event_data a data.frame containing trial-level events that occurred in the time period represented by \code{ts_data}
    #' @param vm a list of variable names used in \code{ts_data} and \code{event_data} that map onto internal constructs
    initialize = function(ts_data=NULL, event_data=NULL, vm=NULL, tr=NULL) {
      if (is.null(tr)) { stop("tr must be provided at object initialization") }
      checkmate::assert_numeric(tr, lower=1e-2, null.ok=FALSE)
      
      default_vm <- list(id="id", run="run", trial="trial", run_trial="run_trial", time="time", value="value")
      if ("cluster" %in% names(ts_data)) { default_vm[["cluster"]] <- "cluster" } #add default cluster mapping if clustered
      for (nn in names(default_vm)) { #populate default variable mappings if not provided in input vector
        if (!nn %in% names(vm)) { vm[nn] <- default_vm[nn] }
      }

      checkmate::assert_data_frame(ts_data)
      checkmate::assert_data_frame(event_data, null.ok=TRUE)

      #setup standardized naming
      private$vmvec <- unlist(vm) #yields cluster1, cluster2, etc.
      if ("cluster" %in% names(vm)) { private$cvars <- paste0("cluster", 1:length(vm$cluster)) }
      if ("value" %in% names(vm)) { private$vvars <- paste0("value", 1:length(vm$value)) }
      
      if (!is.null(event_data)) {
        if (!is.null(vm$id)) {
          stopifnot(vm$id %in% names(event_data))
          if (length(unique(event_data[[vm$id]])) > 1L) {
            stop("fmri_ts objects only support single runs of data for single IDs. You can combine fmri_ts objects using combine_ts()")
          }
        }
        if (!is.null(vm$run)) {
          stopifnot(vm$run %in% names(event_data))
          if (length(unique(event_data[[vm$run]])) > 1L) {
            stop("fmri_ts objects only support single runs of data for single IDs. You can combine fmri_ts objects using combine_ts()")
          }
        }
        
        checkmate::assert_string(vm$trial) #singular string
        stopifnot(vm$trial %in% names(event_data))

        event_data <- data.table(event_data)
        
        #handle internal renaming to make programming with these objects easy
        setnames(event_data, private$vmvec, names(private$vmvec), skip_absent=TRUE)

        setorderv(event_data, vm$trial) #order by trial
      }
      
      #verify presence of required columns      
      sapply(vm[c("time", "value", "cluster")], function(x) { stopifnot(all(x %in% names(ts_data))) } )
      
      #always make ts_data long/tidy? this is nice, but a) double or quadruples RAM demand, and b) adds compute time
      #if (length(vm[["signal"]] > 1L)) { ts_data <- ts_data %>% tidyr::pivot_longer(cols=vm[["signal"]], names_to="signal", values_to="value") }
      #if (length(vm[["cluster"]] > 1L)) { ts_data <- ts_data2 %>% tidyr::pivot_longer(cols=vm[["cluster"]], names_to="cluster_var", values_to="cluster") }

      #we could convert to data.table, then split. If we drop the key columns, we can get the RAM pressure down 20-30%
      #but, it generates a complicated data structure and slows down per-group calculations
      #dd <- split(xx, by=vm[["cluster"]], keep.by=T)

      #conclusion: for internal object storage, sort by clustering variables, then use RLE encoding of keys to compress object
      #  add method $get_ts that returns the rehydrated data (with inverse.rle)

      #convert to data table for speed, memory management
      ts_data <- data.table(ts_data)

      #handle internal renaming to make programming with these objects easy
      setnames(ts_data, private$vmvec, names(private$vmvec), skip_absent=TRUE)

      setorderv(ts_data, private$cvars)
      ts_keys <- sapply(private$cvars, function(x) { rle(ts_data[[x]]) }, simplify=FALSE)
      ts_data[, private$cvars := NULL] #drop key columns

      self$ts_data <- ts_data
      self$ts_keys <- ts_keys
      self$event_data <- event_data
      self$vm <- vm
      self$tr <- tr
    },

    #' @description method to get rehydrated time series object with key values
    #' @param orig_names boolean indicating whether to return data.table with original naming scheme. Default: FALSE
    get_ts = function(orig_names=FALSE) {
      tsd <- data.table::copy(self$ts_data) #ensure that we copy the object to avoid altering $ts_data
      for (kk in 1:length(self$ts_keys)) { tsd[, names(self$ts_keys)[kk] := inverse.rle(self$ts_keys[[kk]])] }
      setcolorder(tsd, private$cvars) #put clustering variables first in object
      if (isTRUE(orig_names)) { setnames(tsd, names(private$vmvec), private$vmvec, skip_absent=TRUE) }
      return(tsd)
    },

    #' @description method to add a variable in ts_data to the set of clustering variables for further use
    add_clusters = function(cv) {
      stopifnot(all(cv %in% names(self$ts_data)))
      for (vname in cv) {
        nclus <- length(private$cvars)
        newvar <- paste0("cluster", nclus+1)
        private$cvars <- c(private$cvars, newvar)
        self$ts_keys[[newvar]] <- rle(self$ts_data[[vname]])      
        
      }

      self$vm[["cluster"]] <- c(self$vm[["cluster"]], cv)
      private$vmvec <- unlist(self$vm) #yields cluster1, cluster2, etc.

      self$ts_data[, (cv) := NULL] #drop new keys
    },
    
    get_cvars = function() { #simple get method to allow access to cluster variables
      return(private$cvars)
    },
    export = function(filename) {
      
    }
  )
)


