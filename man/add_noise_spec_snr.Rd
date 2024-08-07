% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{add_noise_spec_snr}
\alias{add_noise_spec_snr}
\title{Add noise to an mrs_data object to match a given SNR.}
\usage{
add_noise_spec_snr(
  mrs_data,
  target_snr,
  sig_region = c(4, 0.5),
  ref_data = NULL
)
}
\arguments{
\item{mrs_data}{data to add noise to.}

\item{target_snr}{desired spectral SNR, note this assumes the input data is
noise-free, eg simulated data. Note the SNR is estimated from the first
scan in the dataset and the same noise level is added to all spectra.}

\item{sig_region}{spectral limits to search for the strongest spectral data
point.}

\item{ref_data}{measure the signal from the first scan in this reference data
and apply the same target noise level to mrs_data.}
}
\value{
mrs_data object with additive normally distributed noise.
}
\description{
Add noise to an mrs_data object to match a given SNR.
}
