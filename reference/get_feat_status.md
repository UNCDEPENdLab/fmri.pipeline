# helper function to look at whether feat ingredients exist for a .feat/.gfeat directory

helper function to look at whether feat ingredients exist for a
.feat/.gfeat directory

## Usage

``` r
get_feat_status(feat_dir, feat_fsf = NULL, lg = NULL, prefix = NULL)
```

## Arguments

- feat_dir:

  an expected .feat/.gfeat directory for an analysis

- feat_fsf:

  an expected .fsf file corresponding to the feat analysis

- lg:

  an optional lgr logger object to be used for logging. If not passed,
  the root logger will be used

## Value

a data.frame containing information about whether the .feat analysis is
complete and whether various ingredients are present
