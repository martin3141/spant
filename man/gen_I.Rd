% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qm_simulation.R
\name{gen_I}
\alias{gen_I}
\title{Generate the I product operator for a single spin.}
\usage{
gen_I(n, spin_num, op)
}
\arguments{
\item{n}{spin index number for the required operator.}

\item{spin_num}{vector of spin numbers in the system.}

\item{op}{operator, one of "x", "y", "z", "p", "m".}
}
\value{
I product operator matrix.
}
\description{
Generate the I product operator for a single spin.
}
