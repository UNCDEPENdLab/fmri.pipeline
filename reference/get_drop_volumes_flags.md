# Get internal processing flags for drop_volumes handling

Returns a list of internal flags that control how drop_volumes is
applied to different data types. These are not currently exposed to
users.

## Usage

``` r
get_drop_volumes_flags()
```

## Value

list with the following logical flags:

- `shift_nifti`: whether to truncate NIfTI files (FALSE, not supported)

- `shift_timing`: whether to shift event timing by drop_volumes\*tr
  (TRUE)

- `shorten_additional`: whether to drop volumes from additional
  regressors (TRUE)

- `shorten_ts`: whether to drop volumes from ts_multipliers (TRUE)
