\name{selTLSE}
\alias{selTLSE}
\title{
Thresholding Least Squares estimator
}
\description{
This is the main function to estimate the causal effects using the
semiparametric thresholding least squares method.
}
\usage{
selTLSE(model, method=c("FTLSE", "BTLSE"),
        crit = c("ASY", "AIC", "BIC"), 
        minPV = function(p) 1/(p * log(p)), vcov.=NULL, ...)
}
\arguments{
  \item{model}{A model of class \code{tlseModel} created by
    the \code{getModel} function.}
  \item{method}{The selection method: backward for forward.}
  \item{crit}{The criterion to select the piecewise polynomial knots.}
  \item{minPV}{A function to determine the threshold for the
    significance of the coefficients. It has to be a function of one
    parameter, which is the average number of knots in the model.}
  \item{vcov.}{An alternative function to compute the covariance matrix
    of the least squares estimators. The default is the \code{vcov}
    method for \code{lm} objects.}
  \item{...}{Additional arguments to pass to \code{vcov.}}
}

\examples{
data(simData)
mod1 <- setModel(Y~Z|~X1*X2, data=simData)
selTLSE(mod1, "BTLSE")
}

\keyword{causal effects}
