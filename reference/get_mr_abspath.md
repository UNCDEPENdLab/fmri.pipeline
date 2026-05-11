# small helper function to pull absolute paths to a given column in run_data or trial_data

small helper function to pull absolute paths to a given column in
run_data or trial_data

## Usage

``` r
get_mr_abspath(mr_df, col = "run_nifti")
```

## Arguments

- mr_df:

  a data.frame from a gpa object, following standard variable mapping
  nomenclature

- col:

  a character string denoting the column in `mr_df` to be used for
  looking up absolute paths

## Value

Absolute paths to all files in the specified `col`

## Details

Note that if a given value of the requested column is an absolute path,
it will not be combined with \$mr_dir to generate the combined path.
Thus, the `col` in `mr_df` can contain a mixture of relative an absolute
paths. The relative paths will be combined with
