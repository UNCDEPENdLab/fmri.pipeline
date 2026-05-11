# Recalculate an L2 model by grouping units derived from available run data

Recalculate an L2 model by grouping units derived from available run
data

## Usage

``` r
respecify_l2_models_by_subject(
  mobj,
  data,
  split_on = c("id", "session"),
  aggregated_session = 0L
)
```

## Arguments

- mobj:

  an `l1_model_spec` or `hi_model_spec` object containing the GLM model
  to run

- data:

  The run-level data frame containing data for all ids and sessions.
  This is split into chunks by `split_on`.

- split_on:

  Character vector defining grouping columns used to respecify the
  model. Must include `"id"` and may include `"session"` (default:
  `c("id", "session")`).

- aggregated_session:

  Integer session value assigned when grouping across sessions (i.e.,
  when `split_on = "id"`). This keeps output schemas consistent for
  downstream joins and cope metadata. Default is `0L`.

## Value

a modified copy of `mobj` where the `$by_subject` field has been added

## Details

The function adds the `$by_subject` field, which contains the design
matrices and contrasts for each grouping unit in `data`, based on the
runs that remain available for that unit. For example, if a subject is
missing runs (or they are dropped from analysis), some contrasts may
change or drop out of the model.

The `$by_subject` field is a keyed data.table containing list elements
for `cope_list` (mapping cope numbers to contrast names), `contrasts`,
and the design matrix for each unit. It also includes `grouping_scope`
(`"id_session"` or `"id"`).
