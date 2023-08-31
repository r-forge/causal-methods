\name{causalSLSE}
\alias{causalSLSE}
\alias{causalSLSE.slseModel}
\alias{causalSLSE.slseFit}
\alias{causalSLSE.formula}
\title{
Causal Effect based on semiparametric Least Squares models
}
\description{
This is the main method to estimate the causal effects using the semiparametric
least squares method. It returns an object of class
\code{causalslse} and is registered for objects of class
\code{slseModel} and \code{slseFit}.
}
\usage{
\method{causalSLSE}{slseModel}(object,
                               selType=c("SLSE","BSLSE","FSLSE"),
                               selCrit = c("AIC", "BIC", "PVT"),
                               causal = c("ALL","ACT","ACE","ACN"),
                               seType=c("analytic", "lm"),
                               pvalT = function(p) 1/log(p),
                               vcov.=vcovHC, ...)
\method{causalSLSE}{slseFit}(object, seType=c("analytic", "lm"),
       causal = c("ALL","ACT","ACE","ACN"), vcov.=vcovHC, ...)
\method{causalSLSE}{formula}(object, data, nbasis=function(n) n^0.3,
                             knots0, knots1, 
                             selType=c("SLSE","BSLSE","FSLSE"),
                             selCrit = c("AIC", "BIC", "PVT"),
                             causal = c("ALL","ACT","ACE","ACN"),
                             seType=c("analytic", "lm"),
                             pvalT = function(p) 1/log(p),
                             vcov.=vcovHC, ...)
}
\arguments{
  \item{object}{An object of class \code{slseModel} created by the
    \code{setModel} function, \code{slseFit} created by \code{estModel}
    or \code{formula}.}
  \item{data}{A \code{data.frame} with all variables included in \code{form}.}
  \item{nbasis}{A function to determined the number of bases. It has to
    be a function of one argument, the sample size. To fix it, simply
    input a function that returns an integer.}
  \item{knots0}{A list of knots for the control group. The length must
    be equal to the number of covariates. The number of knots is 0 if an
    element of the list is set to \code{NULL}. If missing, the knots are
  automatically generated.}
  \item{knots1}{A list of knots for the control group. See \code{knots0}.}  
  \item{selType}{The method for selecting the piecewise polynomial
    knots. By default (SLSE), all knots from the model are kept.}
  \item{selCrit}{The criterion to select the piecewise polynomial knots.}
  \item{causal}{What causal effect should we compute.}
  \item{seType}{Should the causal effect standard errors be computed
    based on an analytic expression or based on least squares
    asymptotics (\code{vcov.lm}).}
  \item{pvalT}{A function to determine the p-value threshold for the
    significance of the coefficients. It has to be a function of one
    parameter, which is the average number of knots in the model.}
  \item{vcov.}{An alternative function to compute the covariance matrix
    of the least squares estimators. The default is the \code{vcovHC}
    method for \code{lm} objects.}
  \item{...}{Additional arguments to pass to \code{vcov.}}
}

\examples{
data(simDat3)
mod1 <- setModel(Y~Z|~X1*X2, data=simDat3)
fit1 <- causalSLSE(mod1, selType="BSLSE", type="HC1")
causalSLSE(Y~Z|~X1*X2, data=simDat3, selType="BSLSE", type="HC1")
}

\keyword{causal effects}
