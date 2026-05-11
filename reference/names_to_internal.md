# helper function to rename columns of input data.frame to internal nomenclature based on the variable mapping (vm) vector

helper function to rename columns of input data.frame to internal
nomenclature based on the variable mapping (vm) vector

## Usage

``` r
names_to_internal(df, vm)
```

## Arguments

- df:

  a data.frame containing columns to be renamed to internal standards

- vm:

  a named vector of columns in `df` that identify internal constructs
  such as id, session, and run.

## Value

a modified version of `df` with column names modified to use internal
names
