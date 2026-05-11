# Internal function to walk through l2/l3 model setup, including populated from a specification (YAML) file

Internal function to walk through l2/l3 model setup, including populated
from a specification (YAML) file

## Usage

``` r
create_new_hi_model(
  data,
  to_modify = NULL,
  level = NULL,
  cur_model_names = NULL,
  spec_list = NULL,
  lg = NULL
)
```

## Arguments

- data:

  a data.frame containing data that can be included in the model

- to_modify:

  an existing \`hi_model_spec\` object to be modified

- level:

  the model level: 2 or 3

- cur_model_names:

  a character vector of the names of all current models in the model
  list

- spec_list:

  a list containing model specifications read in from a YAML/JSON file

- lg:

  the current logger

## Value

a \`hi_model_spec\` object containing the L2/L3 model to be added to
gpa\$l2_models or gpa\$l3_models

## Details

if \`spec_list\` is passed in, the resulting model object will be
populated from the fields of the specification list, rather than
prompting the user for input
