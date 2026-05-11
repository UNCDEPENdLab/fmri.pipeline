# C++ port of Bush and Cisler 2013, Magnetic Resonance Imaging Adapted from the original provided by Keith Bush as well as C++ code from Jiang Bian

C++ port of Bush and Cisler 2013, Magnetic Resonance Imaging Adapted
from the original provided by Keith Bush as well as C++ code from Jiang
Bian

## Arguments

- BOLDobs:

  matrix of observed BOLD timeseries (n_timepoints x n_signals)

- kernel:

  assumed kernel of the BOLD signal (e.g., from spm_hrf)

- nev_lr:

  learning rate for the assignment of neural events. Default: .01

- epsilon:

  relative error change (termination condition). Default: .005

- beta:

  slope of the sigmoid transfer function (higher = more nonlinear)

- normalize:

  whether to unit-normalize (z-score) `BOLDobs` before deconvolution.
  Default: TRUE

- trim_kernel:

  whether to remove the first K time points from the deconvolved vector,
  corresponding to kernel leftovers from convolution. Default: TRUE

## Value

A time series of the same length containing reconstructed neural events

## Details

This function deconvolves the BOLD signal using Bush 2011 method

Author: Keith Bush, PhD Institution: University of Arkansas at Little
Rock Date: Aug. 9, 2013

The original code did not unit normalize the BOLD signal in advance, but
in my testing, this proves useful in many cases (unless you want to mess
with the learning rate a lot), especially when the time series has a
non-zero mean (e.g., mean 100).

## Author

Michael Hallquist
