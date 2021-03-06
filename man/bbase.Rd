% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abfit.R
\name{bbase}
\alias{bbase}
\title{Generate a spline basis, slightly adapted from : "Splines, knots, and
penalties", Eilers 2010.}
\usage{
bbase(N, number, deg = 3)
}
\arguments{
\item{N}{number of data points.}

\item{number}{number of spline functions.}

\item{deg}{spline degree : deg = 1 linear, deg = 2 quadratic, deg = 3 cubic.}
}
\value{
spline basis as a matrix.
}
\description{
Generate a spline basis, slightly adapted from : "Splines, knots, and
penalties", Eilers 2010.
}
