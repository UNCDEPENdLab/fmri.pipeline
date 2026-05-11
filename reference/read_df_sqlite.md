# helper function to lookup a keyed data.frame from the sqlite storage database

helper function to lookup a keyed data.frame from the sqlite storage
database

## Usage

``` r
read_df_sqlite(
  gpa = NULL,
  db_file = NULL,
  id = NULL,
  session = NULL,
  run_number = NULL,
  table = NULL,
  drop_keys = TRUE,
  quiet = TRUE
)
```

## Arguments

- gpa:

  A `glm_pipeline_arguments` object used to lookup location of SQLite
  database for this analysis

- db_file:

  An optional string specifying the SQLite database from which to read

- id:

  the id of the subject to whom these data belong

- session:

  the session of these data

- run_number:

  the run_number of these data

- table:

  A character string of the table name from which to read

- drop_keys:

  whether to drop identifying metatdata columns from data before
  returning the object

- quiet:

  a logical indicating whether to issue a warning if the table is not
  found

## Value

a data.frame containing the requested data. Will return NULL if not
found
