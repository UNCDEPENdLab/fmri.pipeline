#' Creates a SQLite database connection
#' @param db_path Path to SQLite database
#' @param sqlite_con May pass existing DB connection to check
#' 
#' @return DB connection
#' @importFrom DBI dbDataType dbConnect dbDisconnect dbIsValid dbCommit dbRollback dbBegin
#' @export
get_sqlite_conn <- function(db_path, sqlite_con=NULL) {
    if (is.null(sqlite_con) || !DBI::dbIsValid(sqlite_con)) {
        sqlite_con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
    }
    return(sqlite_con)
}