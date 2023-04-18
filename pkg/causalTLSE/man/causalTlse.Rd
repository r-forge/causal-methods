\name{causalTlse}
\alias{causalTlse}
\title{
Causal Effect based on thresholding Least Squares with knots selection
}
\description{
These are the main functions to select the knots using different methods
and estimate different causal effects using the semiparametric
thresholding least squares method.
}
\usage{
causalTlse(model, selType=c("SLSE","BTLSE","FTLSE"),
     selCrit = c("ASY", "AIC", "BIC"),
     causal = c("ALL","ACT","ACE","ACN"),
     seType=c("analytical", "lm"),
     minPV = function(p) 1/log(p), vcov.=NULL, ...)
}
\arguments{
  \item{model}{A model of class \code{tlseModel} created by the
    \code{getModel} function.}
  \item{selType}{The method for selecting the piecewise polynomial knots.}
  \item{selCrit}{The criterion to select the piecewise polynomial knots.}
  \item{causal}{What causal effect should we compute.}
  \item{seType}{Should the causal effect standard errors be computed
    based on an analytical expression or based on least squares
    asymptotics.}
  \item{minPV}{A function to determine the threshold for the
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
fit1 <- causalTlse(mod1, selType="BTLSE", vcov.=sandwich::vcovHC, type="HC3")
}

\keyword{causal effects}
