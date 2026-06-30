# helper function to manipulate confounds

helper function to manipulate confounds

## Usage

``` r
confound_manipulations(
  gpa,
  all_confounds,
  expected_l1_confound_file,
  run_df,
  demean,
  lg,
  confound_path = NA_character_,
  motion_path = NA_character_
)
```

## Arguments

- gpa:

  A `glm_pipeline_arguments` object.

- all_confounds:

  Combined confound data.

- expected_l1_confound_file:

  Expected output confound filename.

- run_df:

  One-row run metadata.

- demean:

  Whether to demean confounds.

- lg:

  Logger object.

- confound_path:

  Path to the raw confound file.

- motion_path:

  Path to the raw motion file.
