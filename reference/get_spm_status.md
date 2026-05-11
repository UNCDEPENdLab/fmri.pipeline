# Helper to check whether expected SPM outputs exist

Helper to check whether expected SPM outputs exist

## Usage

``` r
get_spm_status(spm_dir, lg = NULL, prefix = NULL)
```

## Arguments

- spm_dir:

  Directory containing SPM outputs

- lg:

  Optional logger

- prefix:

  Optional prefix for column names

## Value

A data.frame with SPM status columns including completion status, timing
information, and output file existence checks.
