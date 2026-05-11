# Obtain user input from the console

Obtain user input from the console

## Usage

``` r
prompt_input(
  prompt = "",
  prompt_eol = ">",
  instruct = NULL,
  type = "character",
  lower = -Inf,
  upper = Inf,
  len = NULL,
  min.len = NULL,
  max.len = NULL,
  split = NULL,
  among = NULL,
  required = TRUE,
  uniq = FALSE,
  default = NULL
)
```

## Arguments

- prompt:

  The character string to place on the line preceding the user input
  prompt. For example, "Enter location"

- prompt_eol:

  The character string to place at the end of the prompt line. For
  example, "\>"

- instruct:

  The instructions to display above the prompt.

- lower:

  For numeric inputs, the lowest valid value

- upper:

  For numeric inputs, the highest valid value

- len:

  The number of expected values to be returned. If NULL, the user can
  enter any number of values.

- min.len:

  The minimum number of values to be returned. If NULL, the user can
  enter any number of values.

- max.len:

  The maximum number of values to be returned. If NULL, the user can
  enter any number of values.

- split:

  The character(s) to split the input string into multiple values. Only
  relevant if len \> 1.

- among:

  A vector of valid values for the input. If NULL, any value is
  accepted.

- required:

  If TRUE, the user must provide a value. If FALSE, the user can skip
  the input by pressing Enter.

- uniq:

  If TRUE, all entries must be unique.

## Value

The user input, converted to the appropriate type (numeric, integer, or
character).

## Details

The function will keep prompting the user until valid input is provided.
It will also display instructions and feedback about the expected input
format.

## Note

This function is intended for interactive use and may not work as
expected in non-interactive environments (e.g., R scripts run in batch
mode).
