# Lookup FSL FEAT output images by model and contrast

Build a tidy lookup table that maps the human-readable pipeline model
and contrast names onto the FEAT folder structure for level 1, 2, and/or
3 FSL outputs.

## Usage

``` r
lookup_feat_outputs(
  gpa,
  level = c(1L, 2L, 3L),
  what = c("cope", "varcope", "zstat", "tstat"),
  include_missing = TRUE,
  include_internal = FALSE,
  source = c("auto", "setup", "cache", "filesystem"),
  cache_dir = NULL,
  refresh_status = FALSE,
  lg = NULL
)
```

## Arguments

- gpa:

  a `glm_pipeline_arguments` object. The function is most informative
  when populated FSL setup tables are present, but `source = "auto"` can
  also search scheduler caches and crawl FEAT folders.

- level:

  integer vector selecting FEAT levels to lookup. Valid values are `1`,
  `2`, and `3`. Defaults to all levels.

- what:

  statistic image types to include. Valid values are `"cope"`,
  `"varcope"`, `"zstat"`, and `"tstat"`.

- include_missing:

  if `TRUE`, include expected output paths even when the image does not
  exist. If `FALSE`, return only existing images.

- include_internal:

  if `TRUE`, retain setup/debug columns such as FEAT execution
  timestamps and L2 input mode. The default `FALSE` returns a compact,
  user-facing lookup table.

- source:

  where to look for output metadata. `"setup"` uses only the
  `gpa$l*_model_setup$fsl` tables, `"cache"` searches scheduler
  `run_pipeline_cache*.RData` files, `"filesystem"` crawls FEAT folders,
  and `"auto"` tries these in that order.

- cache_dir:

  optional directory containing scheduler batch caches. If `NULL`,
  caches are searched below `gpa$output_directory` and
  `gpa$output_locations$scheduler_scripts`.

- refresh_status:

  if `TRUE`, refresh the FEAT status columns in the setup table before
  building the lookup.

- lg:

  optional
  [`lgr::Logger`](https://s-fleck.github.io/lgr/reference/Logger.html)
  object.

## Value

a data.frame with one row per model/contrast/statistic image.
