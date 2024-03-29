% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_read_ima.R
\name{read_ima_dyn_dir}
\alias{read_ima_dyn_dir}
\title{Read a directory containing Siemens MRS IMA files and combine along the
dynamic dimension. Note that the coil ID is inferred from the sorted file
name and should be checked when consistency is required.}
\usage{
read_ima_dyn_dir(dir, extra = NULL, verbose = FALSE)
}
\arguments{
\item{dir}{data directory path.}

\item{extra}{an optional data frame to provide additional variables for use
in subsequent analysis steps, eg id or grouping variables.}

\item{verbose}{output extra information to the console.}
}
\value{
mrs_data object.
}
\description{
Read a directory containing Siemens MRS IMA files and combine along the
dynamic dimension. Note that the coil ID is inferred from the sorted file
name and should be checked when consistency is required.
}
