% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitting.R
\name{nm_fit}
\alias{nm_fit}
\title{Log-ratio normal-multinomial parameters estimation. The parameters mu and sigma are expressed with respect
basis given by function @ilrBase.}
\usage{
nm_fit(X, eps = 0.001, nsim = 1000, parallel.cluster = NULL,
  max.em.iter = 100, verbose = FALSE)
}
\arguments{
\item{X}{count data set}

\item{eps}{maximum error accepted on the last EM-step}

\item{nsim}{number of simulations used in the E-step}

\item{parallel.cluster}{parallel Socket Cluster created with function @makeCluster}

\item{max.em.iter}{maximum number of steps allowed in the EM-algorithm}

\item{verbose}{show information during estimation}
}
\value{
A list with parameters mu and sigma and the number of iterations before convergence
}
\description{
Log-ratio normal-multinomial parameters estimation. The parameters mu and sigma are expressed with respect
basis given by function @ilrBase.
}
\examples{
X = rnormalmultinomial(100, 100, rep(0,4), diag(4))
nm_fit(X, verbose = T)
}
