# Function for quick viewing of log files in the command line

Function for quick viewing of log files in the command line

## Usage

``` r
view_log(input, from_end = TRUE, lines = 15L, level = NULL, id = NULL)
```

## Arguments

- input:

  A character path to a project folder or a gpa object.

- from_end:

  Logical. If TRUE, view lines from end of file instead of beginning.

- lines:

  An integer. How many lines from beginnning/end to view.

- level:

  An integer equal to 1, 2 or 3 (model level).

- id:

  An integer or NULL. The ID number if wishing to view individual-level
  log.

## Value

invisible NULL

## Author

Zach Vig
