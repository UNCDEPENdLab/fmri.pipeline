# Refresh GLM backend status for a given analysis level

Refresh GLM backend status for a given analysis level

## Usage

``` r
refresh_glm_status(gpa, level = 1L, lg = NULL, glm_software = NULL)
```

## Arguments

- gpa:

  a `glm_pipeline_arguments` object

- level:

  the level of analysis to be refreshed (1, 2, or 3)

- lg:

  an optional lgr logger object used for logging

- glm_software:

  optional backend(s) to refresh. If NULL, refreshes all backends in the
  gpa object.

## Value

a modified copy of `gpa` with backend status columns refreshed
