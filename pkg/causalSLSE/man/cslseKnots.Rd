\name{cslseKnots}
\alias{cslseKnots}
\alias{slseKnots}
\title{
Knots creator for basis functions
}
\description{
The functions create an object of class \code{slseKnots} for
semiparametric least squares (SLSE) models or \code{cslseKnots} for
causal SLSE models. For the SLSE, it contains a list of knots, one for
each covariates. For the causal SLSE, it is a list of \code{slseKnots}
objects, one for each treatment group.
}
\usage{
slseKnots(form, data, X, nbasis = function(n) n^0.3, 
          knots)
cslseKnots(form, data, X, Z, nbasis = function(n) n^0.3, 
           knots)
}
\arguments{
  \item{form}{A formula that determines the covariates to include in
    the regression. For \code{slseKnots}, it is a regular formula and
    only the right hand side is considered. For \code{cslseKnots}, it is
    a formula for the outcome versus treatment indicator and a formula
    for the covariates to include in the regression, separated by \code{|}.}
  \item{data}{A \code{data.frame} with all variables included in
    \code{form}. It is required when \code{X} and/or \code{Z} are missing.}
  \item{nbasis}{A function to determined the number of bases. It has to
    be a function of one argument, the sample size. To fix it, simply
    input a function that returns an integer.}
  \item{knots}{An optional list of knots. For \code{slseKnots}, it is a
    list of A list of knots, one for each covariates. For
    \code{cslseKnots}, it is one list of knots per treatment group. Each
    element of the list is also a list  of length  equal to the number
    of covariates. The choice for each covariate is NULL for no knots or 
    numeric to provide specific values. If missing (the default), the
    knots are automatically generated.}
  \item{X, Z}{An optional matrix of covariates and vector of treatment
    indicator. When provided, \code{form} and \code{data} are not needed.}
}

\examples{
data(simDat3)
## knots for SLSE models
k1 <- slseKnots(Y~X1*X2, data=simDat3)
## knots for causal SLSE models
k2 <- cslseKnots(Y~Z|~X1*X2, data=simDat3)
}

\keyword{causal effects}
