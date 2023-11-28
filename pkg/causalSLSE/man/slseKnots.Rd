\name{slseKnots}
\alias{slseKnots}
\title{
Knots creator for basis functions
}
\description{
The functions create an object of class \code{slseKnots} for
semiparametric least squares (SLSE) models. It returns an object of 
class \code{slseKnots}.
}
\usage{
slseKnots(form, data, X, nbasis = function(n) n^0.3, 
          knots)
}
\arguments{
  \item{form}{A formula that determines the covariates to include in
    the regression. It is a regular formula and only the right hand side
    is considered.} 
  \item{data}{A \code{data.frame} with all variables included in
    \code{form}. It is required when \code{X} and/or \code{Z} are missing.}
  \item{nbasis}{A function to determined the number of bases. It has to
    be a function of one argument, the sample size. To fix it, simply
    input a function that returns an integer.}
  \item{knots}{An optional list of knots. It is a  list of A list of
    knots, one for each covariates. The choice for each covariate is
    NULL for no knots or numeric to provide specific values. If missing
    (the default), the knots are automatically generated.}
  \item{X}{An optional matrix of covariates. When provided, \code{form}
    is not needed.} 
}

\examples{
data(simDat3)
k <- slseKnots(Y~X1*X2, data=simDat3)
k
}

\keyword{causal effects}
