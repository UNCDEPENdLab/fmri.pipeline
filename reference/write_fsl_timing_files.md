# Write FSL 3-column timing files

Writes timing files in FSL's 3-column format (onset, duration, value).
Optionally centers values and removes zero-value events.

## Usage

``` r
write_fsl_timing_files(
  dmat,
  output_directory,
  runs_to_output,
  center_values = TRUE
)
```

## Arguments

- dmat:

  2D array of design matrices (runs x regressors)

- output_directory:

  path to output directory

- runs_to_output:

  numeric vector of run numbers for naming

- center_values:

  logical, whether to mean-center values

## Value

matrix of file paths (runs x regressors)
