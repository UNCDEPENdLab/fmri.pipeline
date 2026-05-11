# This function convolves a stimulus vector with the double-gamma hrf

This function convolves a stimulus vector with the double-gamma hrf

## Arguments

- x:

  A vector of volume numbers used to evaluate the function at each value

- a1:

  The a1 parameter of the double gamma

- a2:

  The a2 parameter of the double gamma

- b1:

  The b1 parameter of the double gamma

- b2:

  The b2 parameter of the double gamma

- cc:

  The cc parameter of the double gamma

## Value

A vector of the double-gamma HRF at each value of `x`

## Details

This is an internal function that is used by convolve_hrf

## Author

Michael Hallquist
