# Helper to trace back to highest level parent from a given reference id

Helper to trace back to highest level parent from a given reference id

## Usage

``` r
get_top_parent(tracking_df, ref_id)
```

## Arguments

- tracking_df:

  a data.frame object extracted from the job tracking SQLite database

- ref_id:

  the id of the job from which to trace back
