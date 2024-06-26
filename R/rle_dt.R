#' R6 class for a keyed data.table object that uses run-length encoding to reduce storage size
#' 
#' @importFrom R6 R6Class
#' @importFrom data.table data.table
#' @importFrom checkmate assert_data_frame assert_logical
#' @export
rle_dt <- R6::R6Class("rle_dt",
  private=list(
    # data: keyed data.table object
    data=NULL,

    # keys: RLE-encoded clustering/keying variables
    keys=NULL,

    # key_metadata: at present, levels of factor and whether ordered
    key_metadata=NULL,

    # compression_ratio: ratio of uncompressed to compressed object representing savings in RAM
    compression_ratio=NULL
  ),
  public=list(
    #' @description Create an RLE-encoded copy of the data
    #' @param data A data.frame or data.table object containing original data
    #' @param keys A character vector of keys within \code{data} that should be RLE-encoded
    #' @param optimize_order A boolean indicating whether to search for the smallest encoding scheme.
    #'    The default it to search in data with 4 or fewer keys, but not to search for 5+.
    #'    Optionally, if a positive integer, randomly test ordering over this number of permutations.
    initialize=function(data, keys = NULL, optimize_order = (length(keys) <= 4)) {
      checkmate::assert_data_frame(data)
      checkmate::assert_character(keys, null.ok=TRUE)
      checkmate::assert_integerish(as.numeric(optimize_order), max.len=1, lower=0)

      dsize <- object.size(data)

      #convert to data table for speed, memory management
      dt <- data.table::data.table(data)
      rm(data)

      if (!is.null(keys)) {
        #determine nesting/ordering
        stopifnot(all(keys %in% names(dt)))

        data.table::setkeyv(dt, keys) #start with key and sort order provided by user

        costs <- c()
        klist <- list()
        optimize_order <- as.integer(optimize_order)

        if (optimize_order == 1L && length(keys) > 10) {
          message("More than 10 keys specified. Will disable order optimization because number of permutations is huge!")
          perms <- matrix(data = keys, nrow = 1) # take keys as stated
        } else if (optimize_order != 0L) { #TRUE or numeric
          perms <- compute_permutations(keys = keys) # compute all possible permutations of keys
          if (optimize_order==1L) {
            perms <- perms #test all permutations
          } else {
            perms <- perms[sample(seq_len(nrow(perms)), min(optimize_order, nrow(perms))), ] #random sample
          }
        } else {
          perms <- matrix(data=keys, nrow=1) #take keys as stated
        }

        key_factors <- sapply(keys, function(x) { is.factor(dt[[x]]) })
        if (any(key_factors)) {
          to_convert <- keys[key_factors]
          for (nn in to_convert) { private$key_metadata[[nn]] <- list(levels=levels(dt[[nn]]), ordered=is.ordered(dt[[nn]])) }
          dt[, (to_convert) := lapply(.SD, as.character), .SDcols=to_convert]
        }

        #search for best RLE encoding of keys
        for (ii in seq_len(nrow(perms))) {
          data.table::setorderv(dt, perms[ii, ])
          klist[[ii]] <- sapply(keys, function(x) { rle(dt[[x]]) }, simplify=FALSE)
          costs[ii] <- object.size(klist[[ii]])
        }

        best_order <- which.min(costs)
        data.table::setorderv(dt, perms[best_order, ])

        dt[, (keys) := NULL] #drop key columns
        private$keys <- klist[[best_order]] #set keys field
      }

      #set data field
      private$data <- dt

      savings <- 100 * (1 - (object.size(private$data) + object.size(private$keys))/dsize)
      message("RLE-encoding of keys savings: ", round(savings, 2), "%")
      private$compression_ratio <- savings
    },

    #' @description Simple method to return the data.table with all columns in their original form.
    #' @details Note that the data are modified slightly in that the keys columns are placed first,
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
              data.table::set(dd, j=kk, value=ordered(dd[[kk]], levels=meta$levels))
            } else {
              data.table::set(dd, j=kk, value=factor(dd[[kk]], levels=meta$levels))
            }
          }
        }

        data.table::setkeyv(dd, names(private$keys)) #add keys back to object
        data.table::setcolorder(dd, names(private$keys)) #put clustering variables first in object
        data.table::setorderv(dd, names(private$keys)) #order by original key input order (rather than optimized order)
      }
      return(dd)
    }
  )
)
