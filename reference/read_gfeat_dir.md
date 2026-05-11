# helper function to look up core stats outputs from a .gfeat folder

helper function to look up core stats outputs from a .gfeat folder

## Usage

``` r
read_gfeat_dir(gfeat_dir, what = "all")
```

## Arguments

- gfeat_dir:

  a .gfeat folder containing the outputs of an FSL analysis

- what:

  What to parse in each folder. Currently just passed through to
  read_feat_dir

## Value

a list containing sorted vectors of each stat output
