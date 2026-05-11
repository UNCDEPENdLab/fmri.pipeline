# helper function to build level 1 signals

helper function to build level 1 signals

## Usage

``` r
bl1_build_signals(
  l1_model_set,
  trial_data,
  block_data = NULL,
  subtrial_data = NULL,
  ppi_data = NULL,
  lg = NULL
)
```

## Arguments

- l1_model_set:

  An `l1_model_set` object whose signals should be created or modified

- trial_data:

  A data.frame containing trial-level signal information

- block_data:

  An optional data.frame containing block-level signal information

- subtrial_data:

  An optional data.frame containing subtrial-level signal information

- ppi_data:

  An optional data.frame containing physiological signals for PPI
  analysis

- lg:

  an lgr logger used for status and warning messages

## Value

a modified version of `l1_model_set` with updated `$signals`
