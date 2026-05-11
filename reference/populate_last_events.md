# Internal function to add last_onset and last_offset columns to gpa\$run_data based on events in l1_models

Internal function to add last_onset and last_offset columns to
gpa\$run_data based on events in l1_models

## Usage

``` r
populate_last_events(gpa, lg)
```

## Arguments

- gpa:

  a `glm_pipeline_arguments` object containing valid \$run_data and
  \$l1_models objects

- lg:

  a Logger object for logging results

## Value

a modified copy of gpa with \$run_data populated with last_offset and
last_onset columns (times in seconds)

## Details

The last_onset and last_offset columns are calculated for each run based
on the timing of all events in the \$l1_models\$events list. These are
then used to facilitate run truncation if the user requests truncation
after a final onset or offset.
