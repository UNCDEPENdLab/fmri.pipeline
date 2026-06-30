# helper function to establish sqlite connection and submit query

helper function to establish sqlite connection and submit query

## Usage

``` r
submit_sqlite_query(
  str = NULL,
  sqlite_db = NULL,
  param = NULL,
  busy_timeout = NULL,
  return_result = FALSE
)
```

## Arguments

- str:

  Character query statement to execute

- sqlite_db:

  Character path to SQLite database

- param:

  Optional list of parameters to pass to statement

- busy_timeout:

  Time (in s) after which to retry write operations; default is 10 s

- return_result:

  Logical. If TRUE submits DBI::dbGetQuery instead of DBI::dbExecute;
  Only use if expecting something in return for your query
