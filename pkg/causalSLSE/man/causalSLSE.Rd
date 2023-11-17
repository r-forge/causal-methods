\name{causalSLSE}
\alias{causalSLSE}
\alias{causalSLSE.cslseModel}
\alias{causalSLSE.slseFit}
\alias{causalSLSE.formula}
\title{
Causal Effect based on semiparametric Least Squares models
}
\description{
This is the main method to estimate the causal effects using the semiparametric
least squares method. It returns an object of class
\code{cslse} and is registered for objects of class
\code{cslseModel} and \code{slseFit}.
}
\usage{
\method{causalSLSE}{cslseModel}(object,
                               selType=c("SLSE","BLSE","FLSE"),
                               selCrit = c("AIC", "BIC", "PVT"),
                               causal = c("ALL","ACT","ACE","ACN"),
                               pvalT = function(p) 1/log(p),
                               vcov.=vcovHC, reSelect=FALSE, ...)
\method{causalSLSE}{slseFit}(object, causal = c("ALL","ACT","ACE","ACN"),
                             vcov.=vcovHC, ...)
\method{causalSLSE}{formula}(object, data, nbasis=function(n) n^0.3,
                             knots, 
                             selType=c("SLSE","BLSE","FLSE"),
                             selCrit = c("AIC", "BIC", "PVT"),
                             causal = c("ALL","ACT","ACE","ACN"),
                             pvalT = function(p) 1/log(p),
                             vcov.=vcovHC, reSelect=FALSE, ...)
}
\arguments{
  \item{object}{An object of class \code{cslseModel} created by the
    \code{cslseModel} function, \code{slseFit} created by \code{estSLSE}
    or \code{formula}.}
  \item{data}{A \code{data.frame} with all variables included in \code{form}.}
  \item{nbasis}{A function to determined the number of bases. It has to
    be a function of one argument, the sample size. To fix it, simply
    input a function that returns an integer.}
  \item{knots}{A list of knots for the treated and nontreated
    groups. The list must be named using the group names. Each element
    of the list is also a list  of length  equal to the number of
    covariates. The choice for each covariate is NULL for no knots or
    numeric to provide specific values. If missing, the knots are 
    automatically generated.}
  \item{selType}{The method for selecting the piecewise polynomial
    knots. By default (SLSE), all knots from the model are kept.}
  \item{selCrit}{The criterion to select the piecewise polynomial knots.}
  \item{causal}{What causal effect should we compute.}
  \item{pvalT}{A function to determine the p-value threshold for the
    significance of the coefficients. It has to be a function of one
    parameter, which is the average number of knots in the model.}
  \item{vcov.}{An alternative function to compute the covariance matrix
    of the least squares estimators. The default is the \code{vcovHC}
    method for \code{lm} objects.}
  \item{reSelect}{The previous selections are stored in the model
    object. By default, it is the previous selections that is used. If
    \code{reSelect} is set to \code{TRUE}, the selection is re-computed.}  
  \item{...}{Additional arguments to pass to \code{vcov.}}
}

\examples{
data(simDat3)
mod1 <- cslseModel(Y~Z|~X1*X2, data=simDat3)
fit1 <- causalSLSE(mod1, selType="BLSE", type="HC1")
causalSLSE(Y~Z|~X1*X2, data=simDat3, selType="BLSE", type="HC1")
}

\keyword{causal effects}
