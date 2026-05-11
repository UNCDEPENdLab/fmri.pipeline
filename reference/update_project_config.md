# helper function to create/update project configuration json file in gpa output directory

helper function to create/update project configuration json file in gpa
output directory

## Usage

``` r
update_project_config(gpa, job_sequence, sequence_id, batch_directory = NULL)
```

## Arguments

- gpa:

  a glm_pipeline_arguments object

- job_sequence:

  a list with named logicals indicating whether that part of the
  pipeline was submitted or not (currently finalize, l1, l2, l3,
  cleanup)

- sequence_id:

  the unique identifer used for this batch sequence submission

- batch_directory:

  a path to the batch directory
