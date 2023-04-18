\name{causal}
\alias{causal.tlseModel}
\alias{causal.tlseFit}
\alias{causal}
\title{
Causal Effect based on thresholding Least Squares models
}
\description{
These are the main functions to estimate the causal effects using the
semiparametric thresholding least squares method.
}
\usage{
\method{causal}{tlseModel}(object, seType=c("analytical", "lm"),
       causal = c("ALL","ACT","ACE","ACN"), vcov.=NULL, ...)
\method{causal}{tlseFit}(object, seType=c("analytical", "lm"),
       causal = c("ALL","ACT","ACE","ACN"), vcov.=NULL, ...)
}
\arguments{
  \item{object}{A model of class \code{tlseModel} created by the
    \code{setModel} function or \code{tlseFit} object.}
  \item{causal}{What causal effect should we compute.}
  \item{seType}{Should the causal effect standard errors be computed
    based on an analytical expression or based on least squares
    asymptotics.}
  \item{vcov.}{An alternative function to compute the covariance matrix
    of the least squares estimators. The default is the \code{vcov}
    method for \code{lm} objects.}
  \item{...}{Additional arguments to pass to \code{vcov.}}
}

\examples{
data(simDat3)
mod1 <- setModel(Y~Z|~X1*X2, data=simDat3)
CE1 <- causal(mod1, vcov.=sandwich::vcovHC, type="HC3")
fit1 <- estModel(mod1)
CE2 <- causal(fit1)
}

\keyword{causal effects}
