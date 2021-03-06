% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_2d_psf}
\alias{get_2d_psf}
\title{Get the point spread function (PSF) for a 2D phase encoded MRSI scan.}
\usage{
get_2d_psf(FOV = 160, mat_size = 16, sampling = "circ", hamming = FALSE)
}
\arguments{
\item{FOV}{field of view in mm.}

\item{mat_size}{acquisition matrix size (not interpolated).}

\item{sampling}{can be either "circ" for circular or "rect" for rectangular.}

\item{hamming}{should Hamming k-space weighting be applied (default FALSE).}
}
\value{
A matrix of the PSF with 1mm resolution.
}
\description{
Get the point spread function (PSF) for a 2D phase encoded MRSI scan.
}
