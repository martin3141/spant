% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{crop_td_pts_pot}
\alias{crop_td_pts_pot}
\title{Crop \code{mrs_data} object data points in the time-domain rounding down to
the next smallest power of two (pot). Data that already has a pot length will
not be changed.}
\usage{
crop_td_pts_pot(mrs_data)
}
\arguments{
\item{mrs_data}{MRS data.}
}
\value{
cropped \code{mrs_data} object.
}
\description{
Crop \code{mrs_data} object data points in the time-domain rounding down to
the next smallest power of two (pot). Data that already has a pot length will
not be changed.
}
