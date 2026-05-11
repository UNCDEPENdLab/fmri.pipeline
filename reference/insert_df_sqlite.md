# helper function to insert a keyed data.frame into the sqlite storage database

helper function to insert a keyed data.frame into the sqlite storage
database

## Usage

``` r
insert_df_sqlite(
  gpa = NULL,
  id = NULL,
  session = NULL,
  run_number = NULL,
  data = NULL,
  table = NULL,
  delete_extant = TRUE,
  append = TRUE,
  overwrite = FALSE,
  immediate = FALSE
)
```

## Arguments

- id:

  the id of the subject to whom these data belong

- session:

  the session of these data

- run_number:

  the run_number of these data

- data:

  A `data.frame` containing the data to be inserted into the sqlite db

- table:

  A character string of the table name to be modified

- delete_extant:

  Whether to delete any existing records for this id + session +
  run_number combination

- append:

  Whether to append records to the table (passed through to
  dbWriteTable)

- overwrite:

  Whether to overwrite the existing table (passed through to
  dbWriteTable)

- immediate:

  Whether to open unique connection, commit transaction, then close the
  connection. This should be useful for SQLite concurrency issues in a
  parallel compute environment, but at present we are still getting
  errors even with the immediate approach.

## Value

a TRUE/FALSE indicating whether the record was successfully inserted
