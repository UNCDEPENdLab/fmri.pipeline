# helper function to generate a contrast matrix from an lm() object and a set of user-specified contrasts using emmeans

helper function to generate a contrast matrix from an lm() object and a
set of user-specified contrasts using emmeans

## Usage

``` r
get_contrasts_from_spec(mobj, lmfit = NULL)
```

## Arguments

- mobj:

  a model object created by build_l\<X\>\_models

- lmfit:

  an optional lm() object used for emmeans calculations. If provided,
  this object will be used instead of mobj\$lmfit (the parent lm on the
  overall dataset).

## Value

a modified copy of the model object `mobj` with \$contrast_list and
\$contrasts fully populated
