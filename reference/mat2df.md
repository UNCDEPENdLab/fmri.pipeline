# Fast conversion of a 2D numeric matrix to a 3-column data.frame

Converts a 2D numeric matrix into a 3-column data.frame

## Usage

``` r
mat2df(
  mat,
  na_zeros = FALSE,
  varnames = as.character(c("dim1", "dim2")),
  value_name = "value"
)
```

## Arguments

- mat:

  A `matrix` to convert to data.frame

- na_zeros:

  Whether near-zero values should be converted to `NA`

- varnames:

  Names for the two dimension key columns

- value_name:

  Name for the value column

## Value

A 3-column data.frame with two dimension key columns and one value
column

## Details

This function is a fast matrix-to-long-data-frame converter for the
simple 2D case.

## Author

Michael Hallquist
