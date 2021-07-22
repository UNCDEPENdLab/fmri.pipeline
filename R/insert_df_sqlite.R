#' helper function to insert a keyed data.frame into the sqlite storage database
#'
#' @param id the id of the subject to whom these data belong
#' @param sesssion the session of these data
#' @param run_number the run_number of these data
#' @param data A \code{data.frame} containing the data to be inserted into the sqlite db
#' @param table A character string of the table name to be modified
#' @param delete_extant Whether to delete any existing records for this id + session + run_number combination
#' @param append Whether to append records to the table (passed through to dbWriteTable)
#' @param overwrite Whether to overwrite the existing table (passed through to dbWriteTable)
#' 
#' @return a TRUE/FALSE indicating whether the record was successfully inserted
#' @importFrom checkmate assert_integerish test_null assert_data_frame assert_string
insert_df_sqlite <- function(gpa = NULL, id = NULL, session = NULL, run_number = NULL, data = NULL,
                             table = NULL, delete_extant = TRUE, append = TRUE, overwrite = FALSE) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  if (checkmate::test_null(id)) {
    stop("insert_df_sqlite requires a specific id for keying data")
  }
  checkmate::assert_integerish(session, lower = 1, null.ok = TRUE)
  if (is.null(session)) session <- 1
  checkmate::assert_integerish(run_number, lower = 1, null.ok = TRUE)
  checkmate::assert_data_frame(data, null.ok = FALSE)
  checkmate::assert_string(table, null.ok = FALSE)

  # open connection if needed
  if (is.null(gpa$sqlite_con) || !DBI::dbIsValid(gpa$sqlite_con)) {
    con <- DBI::dbConnect(RSQLite::SQLite(), gpa$sqlite_db)
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

  if (isFALSE(transaction_failed)) { DBI::dbCommit(con) }

  return(invisible(NULL))
}