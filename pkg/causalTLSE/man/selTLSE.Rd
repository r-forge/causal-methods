\name{selTLSE}
\alias{selTLSE}
\title{
Knot selection method
}
\description{
This is the main function to select the knots in \code{tlseModel} objects. It returns
a model with non significant knots removed.
}
\usage{
selTLSE(model, selType=c("BTLSE", "FTLSE"),
        selCrit = c("AIC", "BIC", "PVT"), 
        pvalT = function(p) 1/log(p), vcov.=vcovHC, ...)
}
\arguments{
  \item{model}{A model of class \code{tlseModel} created by
    the \code{getModel} function.}
  \item{selType}{The selection method: backward for forward.}
  \item{selCrit}{The criterion to select the piecewise polynomial knots.}
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
selTLSE(mod1, "FTLSE")
}

\keyword{causal effects}
