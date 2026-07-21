# helper function to look up core stats outputs from a .feat folder

helper function to look up core stats outputs from a .feat folder

## Usage

``` r
read_feat_dir(
  feat_dir,
  what = "all",
  statistics = c("cope", "varcope", "zstat", "tstat"),
  include_missing = FALSE
)
```

## Arguments

- feat_dir:

  a .feat folder containing the outputs of an FSL analysis

- what:

  which parts of the .feat folder should be parsed. If "all", extract
  everything.

- statistics:

  statistic image families to discover. The returned object retains all
  existing statistic fields; unrequested fields contain aligned `NA`
  values.

- include_missing:

  if `TRUE`, include expected statistic paths when the corresponding
  image does not yet exist.

## Value

a list containing sorted vectors of each stat output
