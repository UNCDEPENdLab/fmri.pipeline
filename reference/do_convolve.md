# Internal port of fsl do_convolve

Internal port of fsl do_convolve

## Arguments

- input:

  the input

- kernel:

  the kernel to convolve

- phase:

  an integer for phase-shifting the HRF

- renorm:

  boolean indicating whether to renormalize the output by the sum

## Value

A vector containing the convolution of input and kernel

## Details

This is an internal function used for testing alignment with FSL HRF
convolution
