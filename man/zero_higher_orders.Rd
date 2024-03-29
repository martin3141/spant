% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qm_simulation.R
\name{zero_higher_orders}
\alias{zero_higher_orders}
\title{Zero all coherences including and above a given order.}
\usage{
zero_higher_orders(sys, rho, order)
}
\arguments{
\item{sys}{spin system object.}

\item{rho}{density matrix.}

\item{order}{states higher than or equal to this argument will be set to
zero.}
}
\value{
density matrix.
}
\description{
Zero all coherences including and above a given order.
}
