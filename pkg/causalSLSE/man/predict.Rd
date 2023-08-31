\name{predict}
\alias{predict.slseFit}
\title{Plot method}
\description{
The method computes the predicted outcome for each group with standard
errors and confidence intervals.
}
\usage{

\method{predict}{slseFit}(object, interval=c("none","confidence"),
                       se.fit=FALSE, newdata=NULL, level=0.95,
                       vcov.=vcovHC, ...)
}
\arguments{
  \item{object}{Object of class \code{slseFit}.}
  \item{interval}{If set to \code{interval}, it returns the predicted 
    values along with the lower and upper bounds of
    the confidence interval.}
  \item{se.fit}{Shoud the function returns the standard error of the
    predicted values?}
  \item{level}{The confidence interval level if interval is set to
    \code{confidence}.}
  \item{newdata}{A \code{data.frame} of new data. It must include values
    for all covariates.}
  \item{vcov.}{An alternative function to compute the covariance
    matrix of the least squares estimates. The default is the
    \code{vcovHC} method for \code{lm} objects.}
  \item{...}{Additional argument to pass to the \code{vcov.} function.}  
}

\examples{
data(simDat3)
mod <- setModel(Y~Z|~X1+X2, data=simDat3)
fit <- causalSLSE(mod, seType="lm")
pr <- predict(fit, interval="confidence")
}
