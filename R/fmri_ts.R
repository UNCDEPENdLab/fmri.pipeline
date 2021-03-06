library(R6)

#essential structure
# - time data: time, signal, unit (voxel/region) ([voxels + time] x signals)
# - metadata: id, trial, events (trials x signals)
# - time -> unit mapping: vnum, atlas_value, region etc. (units x grouping)

#ts_data
# $time: time index
# $signal: only used if there is more than one value column (multivariate ts)
# $value: DV value of time series
# $key: 

#' R6 class representing a multivariate time series object for fMRI analysis
#'
#' @importFrom R6 R6Class
#' @export
fmri_ts <- R6::R6Class("fmri_ts",
  private=list(
    kvars=NULL, #internal names of key variables
    vvars=NULL, #internal names of value variables
    vmvec=NULL,  #mapping between input names and internal names
    names_to_original=function(dt) {
      setnames(dt, names(private$vmvec), private$vmvec, skip_absent=TRUE)
    },    
    names_to_internal=function(dt) {
      #look for naming collisions
      poss_conf <- private$vmvec[ intersect(names(private$vmvec), names(dt)) ]
      poss_conf <- poss_conf[poss_conf != names(poss_conf)]
      
      if (length(poss_conf) > 0L) {
        new_names <- if(length(poss_conf) > 1L) { paste0(".aux", 1:length(poss_conf)) } else { ".aux" }
        cat("To avoid naming conflicts, renaming these columns:", paste(names(poss_conf), collapse=", "),
          "to:", paste0(new_names, collapse=", "), "\n")
        self$vm[[".aux"]] <- names(poss_conf)
        setnames(dt, self$vm[[".aux"]], new_names, skip_absent=FALSE)
        private$vmvec <- unlist(self$vm) #yields key1, key2, etc.          
      }

      #handle internal renaming to make programming with these objects easy
      setnames(dt, private$vmvec, names(private$vmvec), skip_absent=TRUE)
    }
    
  ),
  public=list(
    #' @field ts_data time x signals data.table
    ts_data=NULL,

    #' @field ts_keys RLE-encoded keying variables
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
      if ("key" %in% names(ts_data)) { default_vm[["key"]] <- "key" } #add default key mapping if keyed
      for (nn in names(default_vm)) { #populate default variable mappings if not provided in input vector
        if (!nn %in% names(vm)) { vm[nn] <- default_vm[nn] }
      }

      checkmate::assert_data_frame(ts_data)
      checkmate::assert_data_frame(event_data, null.ok=TRUE)

      #setup standardized naming
      private$vmvec <- unlist(vm) #yields key1, key2, etc.
      if ("key" %in% names(vm)) { private$kvars <- paste0("key", 1:length(vm$key)) }
      if ("value" %in% names(vm)) { private$vvars <- paste0("value", 1:length(vm$value)) }
      
      if (!is.null(event_data)) {
        if (!is.null(vm$id)) {
          checkmate::assert_string(vm$id)
          stopifnot(vm$id %in% names(event_data))
          if (length(unique(event_data[[vm$id]])) > 1L) {
            stop("fmri_ts objects only support single runs of data for single IDs. You can combine fmri_ts objects using combine_ts()")
          }
        }
        if (!is.null(vm$run)) {
          checkmate::assert_string(vm$run)
          stopifnot(vm$run %in% names(event_data))
          if (length(unique(event_data[[vm$run]])) > 1L) {
            stop("fmri_ts objects only support single runs of data for single IDs. You can combine fmri_ts objects using combine_ts()")
          }
        }
        
        checkmate::assert_string(vm$trial) #singular string
        stopifnot(vm$trial %in% names(event_data))

        event_data <- data.table(event_data)

        #handle internal renaming to make programming with these objects easy
        private$names_to_internal(event_data)
        
        setorderv(event_data, vm$trial) #order by trial
      }
      
      #verify presence of required columns      
      sapply(vm[c("time", "value", "key")], function(x) { stopifnot(all(x %in% names(ts_data))) } )

      #convert to data table for speed, memory management
      ts_data <- data.table(ts_data)

      #handle internal renaming to make programming with these objects easy
      private$names_to_internal(ts_data)

      setorderv(ts_data, private$kvars)
      ts_keys <- sapply(private$kvars, function(x) { rle(ts_data[[x]]) }, simplify=FALSE)
      ts_data[, private$kvars := NULL] #drop key columns

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
      setcolorder(tsd, private$kvars) #put keying variables first in object
      if (isTRUE(orig_names)) { private$names_to_original(tsd) }
      return(tsd)
    },

    #' @description method to add a variable in ts_data to the set of keying variables for further use
    add_keys = function(kv) {
      stopifnot(all(kv %in% names(self$ts_data)))
      for (vname in kv) {
        nkeys <- length(private$kvars)
        newvar <- paste0("key", nkeys+1)
        private$kvars <- c(private$kvars, newvar)
        self$ts_keys[[newvar]] <- rle(self$ts_data[[vname]]) # RLE-compress the new key        
      }

      self$vm[["key"]] <- c(self$vm[["key"]], kv)
      private$vmvec <- unlist(self$vm) #yields key1, key2, etc.

      self$ts_data[, (kv) := NULL] #drop new keys from rectangular data
    },

    #' @description method to replace one or more variable mappings in the object
    #  TODO: create a private validate_vm method and run both initialize and replace through it
    replace_vm = function(...) {
      replist <- list(...)
      repfields <- names(replist)
      stopifnot(all(repfields %in% self$vm)) #must be replacement

      #revert to original names before we modify
      private$names_to_original(self$ts_data)
      
      for (rr in 1:length(replist)) {
        this_field <- repfields[rr]
        self$vm[[this_field]] <- replist[[rr]]
      }

      #update internal renaming scheme
      private$vmvec <- unlist(self$vm) #yields key1, key2, etc.

      #convert to internal naming scheme
      private$names_to_internal(self$ts_data)
    },
    get_kvars = function() { #simple get method to allow access to key variables
      return(private$kvars)
    },
    get_vmvec = function() {
      return(private$vmvec)
    },
    export = function(filename) {
      
    }
  )
)


