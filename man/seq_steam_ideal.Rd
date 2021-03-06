% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pulse_sequences.R
\name{seq_steam_ideal}
\alias{seq_steam_ideal}
\title{STEAM sequence with ideal pulses.}
\usage{
seq_steam_ideal(spin_params, ft, ref, TE = 0.03, TM = 0.02)
}
\arguments{
\item{spin_params}{spin system definition.}

\item{ft}{transmitter frequency in Hz.}

\item{ref}{reference value for ppm scale.}

\item{TE}{sequence parameter in seconds.}

\item{TM}{sequence parameter in seconds.}
}
\value{
list of resonance amplitudes and frequencies.
}
\description{
STEAM sequence with ideal pulses.
}
