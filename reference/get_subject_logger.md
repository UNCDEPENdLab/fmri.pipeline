# Create or retrieve a subject level logger

Create or retrieve a subject level logger

## Usage

``` r
get_subject_logger(base_logger, id, gpa, log_prefix = "setup_l1_models")
```

## Arguments

- base_logger:

  Name of parent logger (e.g., "glm_pipeline/l1_setup")

- id:

  Subject identifier

- gpa:

  glm_pipeline_arguments object containing logging settings

- log_prefix:

  Optional filename prefix for subject log files
