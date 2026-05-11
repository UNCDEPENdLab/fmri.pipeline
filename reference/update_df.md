# this helper function replaces matching rows of a current data.frame with rows in a new data.frame

this helper function replaces matching rows of a current data.frame with
rows in a new data.frame

## Usage

``` r
update_df(current = NULL, new = NULL, id_cols = NULL, sort = TRUE)
```

## Arguments

- current:

  The current data.frame that may contains rows matching the new
  data.frame

- new:

  A data.frame containing new data that should be added to the current
  data.frame

- id_cols:

  A character vector of columns names that must exist in both the
  current and new data.frames. These are used to determine which rows
  match in the two datasets.

## Details

The goal here is to keep records in the current data.frame that don't
overlap with the new data.frame and replace overlapping records with
those in the new data.frame. This is helpful when you want to update a
master data.frame with new records that may overlap current records or
may be truly new. If there is no overlap in the datasets, this function
basically just binds them together using dplyr::bind_rows.
