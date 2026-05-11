# Export collinearity diagnostics to CSV files

Writes the collinearity summary data to CSV files for external review or
documentation purposes.

## Usage

``` r
export_collinearity_to_csv(
  collin_summary,
  output_dir,
  prefix = "collinearity_"
)
```

## Arguments

- collin_summary:

  An l1_collinearity_summary object from extract_l1_collinearity

- output_dir:

  Directory where CSV files should be written. If NULL, uses the gpa
  output directory.

- prefix:

  Prefix for output file names. Default is "collinearity\_".

## Value

Invisibly returns a character vector of file paths that were written.
