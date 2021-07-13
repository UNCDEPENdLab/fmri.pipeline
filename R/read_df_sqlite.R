#' helper function to lookup a keyed data.frame from the sqlite storage database
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
#' @importFrom glue glue_sql
#' @importFrom DBI dbConnect dbDisconnect dbGetQuery dbExistsTable
read_df_sqlite <- function(gpa = NULL, id = NULL, session = NULL, run_number = NULL, table = NULL, drop_keys=TRUE) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  if (checkmate::test_null(id)) stop("read_df_sqlite requires a specific id for lookup")
  checkmate::assert_integerish(session, lower = 1, null.ok = TRUE)
  if (is.null(session)) session <- 1
  checkmate::assert_integerish(run_number, lower = 1, null.ok = TRUE)
  checkmate::assert_string(table, null.ok = FALSE)
  checkmate::assert_logical(drop_keys, len = 1L)

  # open connection if needed
  if (is.null(gpa$sqlite_con) || !DBI::dbIsValid(gpa$sqlite_con)) {
    con <- DBI::dbConnect(RSQLite::SQLite(), gpa$sqlite_db)
    on.exit(try(DBI::dbDisconnect(con)))
  } else {
    con <- gpa$sqlite_con # recycle connection
  }

  # if table does not exist, then query is invalid (just return NULL)
  if (!DBI::dbExistsTable(con, table)) {
    return(NULL)
  }

  # delete any existing record
  query <- glue::glue_sql(
    "SELECT * FROM {`table`}",
    "WHERE id = {id} AND session = {session}",
    ifelse(is.null(run_number), "", "AND run_number = {run_number}"),
    .con = con, .sep = " "
  )

  data <- tryCatch(DBI::dbGetQuery(con, query), error = function(e) {
    message("Failed to obtain records for query: ", query)
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