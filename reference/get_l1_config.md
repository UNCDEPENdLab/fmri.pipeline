# Function to output l1_models into '.yaml' object

Function to output l1_models into '.yaml' object

## Usage

``` r
get_l1_config(gpa)
```

## Arguments

- gpa:

  the `gpa` object

## Value

a yaml document of l1_models for an object

## Details

This function exports the configuration of level 1 (run-level) GLM
settings in a glm_pipeline_arguments (gpa) object to a YAML file. This
allows the user to then setup their models from this configuration file,
rather than having to go through the \`build_l1_models\` step manually.

## Examples

``` r
if (FALSE) { # \dontrun{
  gpa <- setup_glm_pipeline(analysis_name="flanker_test",
    output_directory="/proj/mnhallqlab/no_backup/flanker_gnomes",
    n_expected_runs = 2, tr = 2.0,
    run_data = run_df, subject_data = subject_df, trial_data = trial_df,
    l1_models=NULL, l2_models=NULL, l3_models=NULL
  )

 # setup level 1 models through menu system
 gpa <- build_l1_models(gpa)

 # get level 1 model settings as a list
 l1_settings <- get_l1_config(gpa)

 # we should typically use export_glm_config to write to file, but we could convert manually
 l1_yaml <- as.yaml(l1_settings)
 writeLines(l1_yaml, "my_l1_config.yaml")

 # in future, we can populate the models using
 gpa <- build_l1_models(gpa, from_spec_file="my_l1_config.yaml")
} # }
```
