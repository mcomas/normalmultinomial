% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coda.R
\name{ilr_coordinates}
\alias{ilr_coordinates}
\title{ilr coordinates}
\usage{
ilr_coordinates(X)
}
\arguments{
\item{X}{compositional dataset. Either a matrix, a data.frame or a vector}
}
\value{
coordinates with respect basis ilr_basis
}
\description{
Calculate the ilr coordinates with respect basis given by function
ilr_basis
}
\examples{
X = rnormalmultinomial(100, 100, rep(0,4), diag(4))
nm_fit(X, verbose = T)
}

