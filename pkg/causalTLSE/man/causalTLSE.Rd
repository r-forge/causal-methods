\name{causalTLSE}
\alias{causalTLSE}
\alias{causalTLSE.tlseModel}
\alias{causalTLSE.tlseFit}
\alias{causalTLSE.formula}
\title{
Causal Effect based on thresholding Least Squares models
}
\description{
This is the main method to estimate the causal effects using the semiparametric
thresholding least squares method. It returns an object of class
\code{causaltlse} and is registered for objects of class
\code{tlseModel} and \code{tlseFit}.
}
\usage{
\method{causalTLSE}{tlseModel}(object,
                               selType=c("SLSE","BTLSE","FTLSE"),
                               selCrit = c("AIC", "BIC", "ASY"),
                               causal = c("ALL","ACT","ACE","ACN"),
                               seType=c("analytic", "lm"),
                               pvalT = function(p) 1/log(p),
                               vcov.=vcovHC, ...)
\method{causalTLSE}{tlseFit}(object, seType=c("analytic", "lm"),
       causal = c("ALL","ACT","ACE","ACN"), vcov.=vcovHC, ...)
\method{causalTLSE}{formula}(object, data, nbasis=function(n) n^0.3,
                             knots0, knots1, 
                             selType=c("SLSE","BTLSE","FTLSE"),
                             selCrit = c("AIC", "BIC", "ASY"),
                             causal = c("ALL","ACT","ACE","ACN"),
                             seType=c("analytic", "lm"),
                             pvalT = function(p) 1/log(p),
                             vcov.=vcovHC, ...)
}
\arguments{
  \item{object}{An object of class \code{tlseModel} created by the
    \code{setModel} function, \code{tlseFit} created by \code{estModel}
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
fit1 <- causalTLSE(mod1, selType="BTLSE", type="HC1")
causalTLSE(Y~Z|~X1*X2, data=simDat3, selType="BTLSE", type="HC1")
}

\keyword{causal effects}
