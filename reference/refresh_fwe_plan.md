# Refresh the state of a compiled FWE plan

Recheck method inputs and expected artifacts without rerunning a
correction. This is especially useful after scheduler-submitted jobs
have finished.

## Usage

``` r
refresh_fwe_plan(plan)
```

## Arguments

- plan:

  an `fwe_plan` or path to an RDS plan written by
  [`write_fwe_plan()`](https://hallquistlab.github.io/fmri.pipeline/reference/write_fwe_plan.md).

## Value

a validated `fwe_plan` with refreshed task and artifact state.
