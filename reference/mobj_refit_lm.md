# Helper function to re-calculate the l3 design matrix for a model based on available data

Helper function to re-calculate the l3 design matrix for a model based
on available data

## Usage

``` r
mobj_refit_lm(mobj, new_data)
```

## Arguments

- mobj:

  a `hi_model_spec` object containing the L3 GLM model to run

- new_data:

  The run-level data frame containing data for all ids and sessions.
  This will be split into individual chunks

## Value

a modified copy of `mobj` where the \$by_subject field has been added

## Details

The function adds the \$by_subject field, which contains the design
matrices and contrasts for each subject and session in `data` based on
the available data for that session. For example, if a subject is
missing a few runs (or these are dropped from analysis), then some
contrasts may change or drop out of the model.

The \$by_subject field is a keyed data.table object containing list
elements for the cope_list (mapping cope numbers to contrast names), the
contrasts, and the design matrix for each session.
