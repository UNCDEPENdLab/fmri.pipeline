# helper function to look up core stats outputs from a .feat folder

helper function to look up core stats outputs from a .feat folder

## Usage

``` r
read_feat_dir(feat_dir, what = "all")
```

## Arguments

- feat_dir:

  a .feat folder containing the outputs of an FSL analysis

- what:

  which parts of the .feat folder should be parsed. If "all", extract
  everything.

## Value

a list containing sorted vectors of each stat output
