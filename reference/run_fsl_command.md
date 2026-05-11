# Wrapper for running an FSL command safely within R

Wrapper for running an FSL command safely within R

## Usage

``` r
run_fsl_command(
  args,
  fsldir = NULL,
  stdout = NULL,
  stderr = NULL,
  echo = TRUE,
  ...
)
```

## Arguments

- args:

  FSL command string to be run

- fsldir:

  Location of FSL installation. If NULL, the function will search the
  environment for FSLDIR or FSL commands in the PATH

- stdout:

  File target for redirecting stdout. If NULL, stdout will not be
  captured

- stderr:

  File target for redirecting stderr. If NULL, stderr will not be
  captured

- echo:

  Whether to print FSL command to the screen. Default: TRUE

- ...:

  Arguments passed through to the \`system\` command.

## Value

The exit status of the executed FSL command. 0 for success, non-zero for
failure

## Details

This command ensures that FSL command are run in an environment with FSL
setup correctly.

## Author

Michael Hallquist

## Examples

``` r
if (FALSE) { # \dontrun{
run_fsl_command("fslmaths test_data copy_data")
} # }
```
