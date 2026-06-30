# helper function to read confounds and motion parameters and combine them

helper function to read confounds and motion parameters and combine them

## Usage

``` r
read_confounds_motion_parameters(
  gpa,
  id,
  session,
  run_number,
  run_df,
  expected_l1_confound_file,
  lg,
  confound_path = NA_character_,
  motion_path = NA_character_
)
```

## Arguments

- gpa:

  A `glm_pipeline_arguments` object.

- id:

  subject id.

- session:

  session number.

- run_number:

  run number.

- run_df:

  One-row run metadata.

- expected_l1_confound_file:

  Expected output confound filename.

- lg:

  Logger object.

- confound_path:

  Path to the raw confound file.

- motion_path:

  Path to the raw motion file.
