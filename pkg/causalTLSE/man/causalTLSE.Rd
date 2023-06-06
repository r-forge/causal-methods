\name{causalTLSE}
\alias{causalTLSE}
\alias{causalTLSE.tlseModel}
\alias{causalTLSE.tlseFit}
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
                               vcov.=NULL, ...)
\method{causalTLSE}{tlseFit}(object, seType=c("analytic", "lm"),
       causal = c("ALL","ACT","ACE","ACN"), vcov.=NULL, ...)
}
\arguments{
  \item{object}{A model of class \code{tlseModel} created by the
    \code{setModel} function.}
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
    of the least squares estimators. The default is the \code{vcov}
    method for \code{lm} objects.}
  \item{...}{Additional arguments to pass to \code{vcov.}}
}

\examples{
data(simDat3)
mod1 <- setModel(Y~Z|~X1*X2, data=simDat3)
fit1 <- causalTLSE(mod1, selType="BTLSE", vcov.=sandwich::vcovHC, type="HC3")
}

\keyword{causal effects}
