# Build concatenated L1 confounds for multi-run analyses

Internal helper that uses stored confound column metadata to align
shared columns and union spike regressors across runs with zero padding.

## Usage

``` r
read_l1_confound_columns(gpa, id, session, run_number, lg = NULL)
```
