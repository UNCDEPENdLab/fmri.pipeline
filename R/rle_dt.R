#' R6 class for a keyed data.table object that uses run-length encoding to reduce storage size
#' 
#' @importFrom R6 R6Class
#' @importFrom data.table data.table
#' @importFrom gtools permutations
#' @importFrom checkmate assert_data_frame assert_logical
#' @export
rle_dt <- R6::R6Class("rle_dt",
  private=list(
    #' @field keyed data.table object
    data=NULL,

    #' @field keys RLE-encoded clustering/keying variables
    keys=NULL,

    key_metadata=NULL, #at present, levels of factor and whether ordered
  ),
  public=list(
    #' @description Create an RLE-encoded copy of the data
    #' @param data A data.frame or data.table object containing original data
    #' @keys A character vector of keys within \code{data} that should be RLE-encoded
    #' @optimize_order A boolean indicating whether to search for the smallest encoding scheme.
    #'    The default it to search in data with 4 or fewer keys, but not to search for 5+.
    initialize=function(data, keys=NULL, optimize_order=(length(keys) <= 4)) {
      checkmate::assert_data_frame(data)
      checkmate::assert_character(keys, null.ok=TRUE)
      checkmate::assert_logical(optimize_order, max.len=1)
      
      #convert to data table for speed, memory management
      dt <- data.table(data)

      if (!is.null(keys)) {        
        #determine nesting/ordering
        stopifnot(all(keys %in% names(data)))

        costs <- c()
        klist <- list()
        
        if (isTRUE(optimize_order)) {
          perms <- gtools::permutations(n=length(keys), r=length(keys), v=keys, repeats.allowed=FALSE)
        } else {
          perms <- matrix(data=keys, nrow=1)
        }

        key_factors <- sapply(keys, function(x) { is.factor(dt[[x]]) })
        if (any(key_factors)) {
          to_convert <- keys[key_factors]
          for (nn in to_convert) { private$key_metadata[[nn]] <- list(levels=levels(dt[[nn]]), ordered=is.ordered(dt[[nn]])) }
          dt[, (to_convert) := lapply(.SD, as.character), .SDcols=to_convert]
        }

        #search for best RLE encoding of keys
        for(ii in 1:nrow(perms)) {
          setorderv(dt, perms[ii,])
          klist[[ii]] <- sapply(keys, function(x) { rle(dt[[x]]) }, simplify=FALSE)
          costs[ii] <- object.size(klist[[ii]])
        }

        best_order <- which.min(costs)
        setorderv(dt, perms[best_order,])
        
        dt[, (keys) := NULL] #drop key columns
        private$keys <- klist[[best_order]] #set keys field
      }

      #set data field
      private$data <- dt
    },

    #' @description Simple method to return the data.table with all columns in their original form.
    #' @detail Note that the data are modified slightly in that the keys columns are placed first,
    #'   and the data are ordered in the order of the keys (as originally provided, left-to-right)
    get = function() {
      dd <- data.table::copy(private$data) #ensure that we copy the object to avoid altering $data

      #rehydrate key columns
      if (!is.null(private$keys)) {
        for (kk in names(private$keys)) {
          dd[, (kk) := inverse.rle(private$keys[[kk]])]

          #handle factors
          if (!is.null(meta <- private$key_metadata[[kk]])) {
            if (isTRUE(meta$ordered)) {
              set(dd, j=kk, value=ordered(dd[[kk]], levels=meta$levels))
            } else {
              set(dd, j=kk, value=factor(dd[[kk]], levels=meta$levels))
            }
          }
        }
        setcolorder(dd, names(private$keys)) #put clustering variables first in object
        setorderv(dd, names(private$keys)) #order by original key inputs
      }
      return(dd)
    }
  )
)
