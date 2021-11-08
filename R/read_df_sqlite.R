#' helper function to lookup a keyed data.frame from the sqlite storage database
#'
#' @param gpa A \code{glm_pipeline_arguments} object used to lookup location of SQLite database for this analysis
#' @param db_file An optional string specifying the SQLite database from which to read
#' @param id the id of the subject to whom these data belong
#' @param session the session of these data
#' @param run_number the run_number of these data
#' @param table A character string of the table name from which to read
#' @param drop_keys whether to drop identifying metatdata columns from data before returning the object
#'
#' @return a data.frame containing the requested data. Will return NULL if not found
#' @importFrom checkmate assert_integerish test_null assert_data_frame assert_string
#' @importFrom glue glue_sql
#' @importFrom DBI dbConnect dbDisconnect dbGetQuery dbExistsTable
read_df_sqlite <- function(gpa = NULL, db_file=NULL, id = NULL, session = NULL, run_number = NULL, table = NULL, drop_keys=TRUE) {
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
    warning(sprintf("Cannot find SQLite table %s in file %s.", table, db_file))
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