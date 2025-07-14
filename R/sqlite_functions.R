#' helper function to insert a keyed data.frame into the sqlite storage database
#'
#' @param id the id of the subject to whom these data belong
#' @param session the session of these data
#' @param run_number the run_number of these data
#' @param data A \code{data.frame} containing the data to be inserted into the sqlite db
#' @param table A character string of the table name to be modified
#' @param delete_extant Whether to delete any existing records for this id + session + run_number combination
#' @param append Whether to append records to the table (passed through to dbWriteTable)
#' @param overwrite Whether to overwrite the existing table (passed through to dbWriteTable)
#' @param immediate Whether to open unique connection, commit transaction, then close the connection.
#'   This should be useful for SQLite concurrency issues in a parallel compute environment, but at present
#'   we are still getting errors even with the immediate approach.
#'
#' @return a TRUE/FALSE indicating whether the record was successfully inserted
#' @importFrom checkmate assert_integerish test_null assert_data_frame assert_string
#' @importFrom DBI dbDataType dbConnect dbDisconnect dbIsValid dbCommit dbRollback dbBegin
#' @importFrom glue glue_sql
insert_df_sqlite <- function(gpa = NULL, id = NULL, session = NULL, run_number = NULL, data = NULL,
                             table = NULL, delete_extant = TRUE, append = TRUE, overwrite = FALSE, immediate=FALSE) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  if (checkmate::test_null(id)) {
    stop("insert_df_sqlite requires a specific id for keying data")
  }
  checkmate::assert_integerish(session, lower = 1, null.ok = TRUE)
  if (is.null(session)) session <- 1
  checkmate::assert_integerish(run_number, lower = 1, null.ok = TRUE)
  checkmate::assert_data_frame(data, null.ok = FALSE)
  checkmate::assert_string(table, null.ok = FALSE)
  checkmate::assert_logical(delete_extant, null.ok = FALSE, len=1L)
  checkmate::assert_logical(append, null.ok = FALSE, len = 1L)
  checkmate::assert_logical(overwrite, null.ok = FALSE, len = 1L)
  
  # open connection if needed
  if (isTRUE(immediate)) {
    # cf. https://blog.r-hub.io/2021/03/13/rsqlite-parallel/
    con <- DBI::dbConnect(RSQLite::SQLite(), gpa$output_locations$sqlite_db)
    RSQLite::sqliteSetBusyHandler(con, get_default_sqlite_busy_timeout()) # retry write operations several times
    on.exit(try(DBI::dbDisconnect(con)))
  } else if (is.null(gpa$sqlite_con) || !DBI::dbIsValid(gpa$sqlite_con)) {
    con <- DBI::dbConnect(RSQLite::SQLite(), gpa$output_locations$sqlite_db)
    on.exit(try(DBI::dbDisconnect(con)))
  } else {
    con <- gpa$sqlite_con # recycle connection
  }
  
  # handle columns in data that are not in table
  has_table <- DBI::dbExistsTable(con, table)
  if (isTRUE(has_table)) {
    table_names <- DBI::dbListFields(con, table)
    uniq_df <- setdiff(names(data), table_names)
    if (length(uniq_df) > 0L) {
      DBI::dbBegin(con) # begin transaction
      alter_failed <- FALSE
      for (nn in uniq_df) {
        dtype <- DBI::dbDataType(con, data[[nn]])
        query <- glue::glue_sql("ALTER TABLE {table} ADD COLUMN {nn} {dtype};", .con = con)
        q_result <- tryCatch(DBI::dbExecute(con, query), error = function(e) {
          message("Error with query: ", query)
          message(as.character(e))
          DBI::dbRollback(con)
          return(FALSE)
        })
        if (isFALSE(q_result)) {
          alter_failed <- TRUE
          break # end loop
        }
      }
      if (!alter_failed) DBI::dbCommit(con) # commit transaction
    }
  }
  
  
  # treat the delete and append as a single transaction so that if either fails, the table is unchanged
  DBI::dbBegin(con)
  transaction_failed <- FALSE
  
  # delete any existing record
  if (isTRUE(delete_extant) && isTRUE(has_table)) {
    query <- glue::glue_sql(
      "DELETE FROM {`table`}",
      "WHERE id = {id} AND session = {session}",
      ifelse(is.null(run_number), "", "AND run_number = {run_number}"),
      .con = con, .sep = " "
    )
    q_result <- tryCatch(DBI::dbExecute(con, query), error = function(e) {
      message("Problem with query: ", query)
      message(as.character(e))
      DBI::dbRollback(con)
      return(FALSE)
    })
    if (isFALSE(q_result)) { transaction_failed <- TRUE }
  }
  
  # add record -- include keying fields for lookup
  data$id <- id
  data$session <- session
  if (!is.null(run_number)) data$run_number <- run_number
  
  q_result <- tryCatch(
    DBI::dbWriteTable(conn = con, name = table, value = data, append = append, overwrite = overwrite),
    error = function(e) {
      print(as.character(e))
      return(FALSE)
    }
  )
  
  if (isFALSE(q_result)) {
    transaction_failed <- TRUE
  }
  
  # commit delete and insert if no errors in subcomponents
  if (isFALSE(transaction_failed)) { DBI::dbCommit(con) }
  
  return(invisible(NULL))
}


#' helper function to lookup a keyed data.frame from the sqlite storage database
#'
#' @param gpa A \code{glm_pipeline_arguments} object used to lookup location of SQLite database for this analysis
#' @param db_file An optional string specifying the SQLite database from which to read
#' @param id the id of the subject to whom these data belong
#' @param session the session of these data
#' @param run_number the run_number of these data
#' @param table A character string of the table name from which to read
#' @param drop_keys whether to drop identifying metatdata columns from data before returning the object
#' @param quiet a logical indicating whether to issue a warning if the table is not found
#'
#' @return a data.frame containing the requested data. Will return NULL if not found
#' @importFrom checkmate assert_integerish test_null assert_data_frame assert_string
#' @importFrom glue glue_sql
#' @importFrom DBI dbConnect dbDisconnect dbGetQuery dbExistsTable
#' @keywords internal
read_df_sqlite <- function(gpa = NULL, db_file=NULL, id = NULL, session = NULL, run_number = NULL, table = NULL, drop_keys=TRUE, quiet=TRUE) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments", null.ok = TRUE)
  if (is.null(gpa)) {
    checkmate::assert_string(db_file)
    checkmate::assert_file_exists(db_file)
    extant_con <- NULL
  } else {
    db_file <- gpa$output_locations$sqlite_db
    extant_con <- gpa$sqlite_con
  }
  if (checkmate::test_null(id)) stop("read_df_sqlite requires a specific id for lookup")
  checkmate::assert_integerish(session, lower = 1, null.ok = TRUE)
  if (is.null(session)) session <- 1
  checkmate::assert_integerish(run_number, lower = 1, null.ok = TRUE)
  checkmate::assert_string(table, null.ok = FALSE)
  checkmate::assert_logical(drop_keys, len = 1L)
  
  # open connection if needed
  if (is.null(extant_con) || !DBI::dbIsValid(extant_con)) {
    con <- DBI::dbConnect(RSQLite::SQLite(), db_file)
    on.exit(try(DBI::dbDisconnect(con)))
  } else {
    con <- extant_con # recycle connection
  }
  
  # if table does not exist, then query is invalid (just return NULL)
  if (!DBI::dbExistsTable(con, table)) {
    if (isFALSE(quiet)) warning(sprintf("Cannot find SQLite table %s in file %s.", table, db_file))
    return(NULL)
  }
  
  # lookup any existing record
  query <- glue::glue_sql(
    "SELECT * FROM {`table`}",
    "WHERE id = {id} AND session = {session}",
    ifelse(is.null(run_number), "", "AND run_number = {run_number}"),
    .con = con, .sep = " "
  )
  
  data <- tryCatch(DBI::dbGetQuery(con, query), error = function(e) {
    if (isFALSE(quiet)) warning("Failed to obtain records for query: ", query)
    return(data.frame())
  })
  
  if (nrow(data) > 0L && isTRUE(drop_keys)) {
    data <- data %>% dplyr::select(-id, -session)
    if (!is.null(run_number)) {
      data <- data %>% dplyr::select(-run_number)
    }
  }
  
  # return NULL in case of zero matches
  if (nrow(data) == 0L) data <- NULL
  
  return(data)
}


#' helper function to establish sqlite connection and submit query
#' 
#' @param sqlite_db Character path to SQLite database
#' @param str Character query statement to execute
#' @param params Optional list of parameters to pass to statement
#' @param busy_timeout Time (in s) after which to retry write operations; default is 10 s
#' @param return_result Logical. If TRUE submits DBI::dbGetQuery instead of DBI::dbExecute;
#'                        Only use if expecting something in return for your query
#' 
#' @importFrom DBI dbConnect dbExecute dbDisconnect
#' @importFrom RSQLite sqliteSetBusyHandler
#' 
#' @keywords internal
submit_sqlite_query <- function(str = NULL, sqlite_db = NULL, param = NULL, 
                                busy_timeout = NULL, return_result = FALSE) {
  checkmate::assert_logical(return_result)
  if(is.null(str) | is.null(sqlite_db)) return(invisible(NULL))
  if(is.null(busy_timeout)) busy_timeout <- get_default_sqlite_busy_timeout()
  
  con <- dbConnect(RSQLite::SQLite(), sqlite_db) # establish connection
  sqliteSetBusyHandler(con, busy_timeout) # busy_timeout arg in seconds * 1000 ms
  
  if (isTRUE(return_result)) {
    res <- dbGetQuery(con, str, param = param) # execute query and return result
  } else {
    res <- dbExecute(con, str, param = param) # execute query
  }
  
  dbDisconnect(con) # disconnect
  
  return(invisible(res))
}


#' helper function to check if a table exists in an SQLite database
#' 
#' @param sqlite_db Character path to SQLite database
#' @param table Character string name of table
#' 
#' @importFrom DBI dbConnect dbExistsTable dbDisconnect
#' @importFrom RSQLite sqliteSetBusyHandler
#' 
#' @keywords internal
sqlite_table_exists <- function(sqlite_db = NULL, table = NULL, busy_timeout = NULL) {
  
  if(is.null(sqlite_db) | is.null(table)) return(invisible(FALSE))
  if(is.null(busy_timeout)) busy_timeout <- get_default_sqlite_busy_timeout()
  
  con <- dbConnect(RSQLite::SQLite(), sqlite_db) # establish connection
  sqliteSetBusyHandler(con, busy_timeout) # busy_timeout arg in seconds * 1000 ms
  table_exists <- dbExistsTable(con, table)
  dbDisconnect(con)
  
  return(table_exists)
  
}

#' helper function to get generic busy timeout
#' 
#' @param ms Logical. Return in milliseconds? If FALSE, returns in seconds
#' 
#' @keywords internal
get_default_sqlite_busy_timeout <- function(ms = T) {
  
  default_seconds <- 10
  
  if (isTRUE(ms)) {
    return(default_seconds * 1000)
  } else {
    return(default_seconds)
  }
  
}

