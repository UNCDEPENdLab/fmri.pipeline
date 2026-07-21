# helper function to look up core stats outputs from a .gfeat folder

helper function to look up core stats outputs from a .gfeat folder

## Usage

``` r
read_gfeat_dir(
  gfeat_dir,
  what = "all",
  statistics = c("cope", "varcope", "zstat", "tstat"),
  include_missing = FALSE
)
```

## Arguments

- gfeat_dir:

  a .gfeat folder containing the outputs of an FSL analysis

- what:

  What to parse in each folder. Currently just passed through to
  read_feat_dir

- statistics:

  statistic image families to discover in each cope directory.

- include_missing:

  if `TRUE`, include expected statistic paths when the corresponding
  image does not yet exist.

## Value

a list containing sorted vectors of each stat output
