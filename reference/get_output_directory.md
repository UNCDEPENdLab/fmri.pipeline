# small helper function to return the location of an l1 directory based on id, session, and run number

small helper function to return the location of an l1 directory based on
id, session, and run number

## Usage

``` r
get_output_directory(
  id = NULL,
  session = NULL,
  run_number = NULL,
  l1_model = NULL,
  l2_model = NULL,
  l3_model = NULL,
  l1_contrast = NULL,
  l2_contrast = NULL,
  l3_contrast = NULL,
  what = "l1",
  gpa,
  glm_software = "fsl",
  create_if_missing = FALSE
)
```

## Arguments

- id:

  The id of a participant

- session:

  The session number to lookup

- run_number:

  The run number to lookup

- gpa:

  A `glm_pipeline_arguments` object

- glm_software:

  which software is being used for the analysis (since directories may
  vary)

- create_if_missing:

  whether to create the directory if it does not exist
