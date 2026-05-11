# Derive a session-level data.frame from run- and subject-level inputs

Derive a session-level data.frame from run- and subject-level inputs

## Usage

``` r
derive_session_data(run_data, subject_data, lg)
```

## Arguments

- run_data:

  A run-level data.frame with id/session/run_number columns.

- subject_data:

  A subject-level data.frame with unique id/session rows.

- lg:

  Logger used for diagnostics.

## Value

A data.frame with one row per id/session.
