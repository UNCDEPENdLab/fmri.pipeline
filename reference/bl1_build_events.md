# helper function to build events consisting of onsets and durations

helper function to build events consisting of onsets and durations

## Usage

``` r
bl1_build_events(l1_model_set, spec_list = NULL, trial_data, lg = NULL)
```

## Arguments

- l1_model_set:

  an l1_model_set object that may have extant events in it

- spec_list:

  an optional parsed YAML/JSON model specification used to populate
  events non-interactively

- trial_data:

  the trial_data object from the `gpa` object

- lg:

  an lgr logger used for status and warning messages

## Value

a modified copy of l1_model_set with events added/updated
