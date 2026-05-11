# Internal helper function to update tracker_args object

Internal helper function to update tracker_args object

## Usage

``` r
populate_list_arg(list_to_populate, arg_name, new_value = NULL, append = FALSE)
```

## Arguments

- list_to_populate:

  The list whose argument will be populated/updated

- arg_name:

  The named list element to update

- new_value:

  The new value to update the element with

- append:

  If TRUE, appends the new value to the current value using the paste
  function. Default: FALSE
