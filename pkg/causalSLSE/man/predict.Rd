\name{predict}
\alias{predict.cslseFit}
\alias{predict.slseFit}
\title{Outcome Prediction}
\description{
The method computes the predicted outcome for each group with standard
errors and confidence intervals.
}
\usage{

\method{predict}{cslseFit}(object, interval=c("none","confidence"),
                       se.fit=FALSE, newdata=NULL, level=0.95,
                       vcov.=vcovHC, ...)
\method{predict}{slseFit}(object, interval=c("none","confidence"),
                       se.fit=FALSE, newdata=NULL, level=0.95,
                       vcov.=vcovHC, ...)
}
\arguments{
  \item{object}{Object of class \code{cslseFit} or  \code{slseFit}
    created by \code{\link{estSLSE}}.}
  \item{interval}{If set to \code{"confidence"}, it returns the
    predicted values along with the lower and upper bounds of the
    confidence interval.}   
  \item{se.fit}{Should the function return the standard errors of the 
    predicted values?}
  \item{level}{The confidence interval level if interval is set to
    \code{"confidence"}.}
  \item{newdata}{A \code{data.frame} of new data. It must include values
    for all covariates, and for the treatment indicator in the case of
    \code{cslseFit} objects.}
  \item{vcov.}{An alternative function to compute the covariance
    matrix of the least squares estimates. The default is the
    \code{\link[sandwich]{vcovHC}}.}
  \item{...}{Additional argument to pass to the \code{vcov.} function.}  
}

\value{ For \code{slseFit} objects, it returns the predicted outcome if
  \code{se.fit} is \code{FALSE} or a list of the following two elements
  otherwise:

  \item{fit}{The predicted outcome.}

  \item{se.fit}{The standard errors of the predicted outcomes.}

  If the argument \code{confidence} is set to \code{"interval"}, the
  predicted outcome is a matrix with the predicted outcome, and the
  lower and upper bounds of the confidence intervals.

  For objects of class 'cslseFit', the same is returned for each
  treatment group in a list. The elements of the list are \code{treated}
  and \code{nontreated} (until the package allows for more than one
  treatment).  
  
}

\examples{
data(simDat3)
mod <- cslseModel(Y ~ Z | ~ X1 + X2, data = simDat3)
fit <- causalSLSE(mod)

## Predicting outcome for all observations

pr <- predict(fit, interval = "confidence")

## Predicting outcome with new data

ndat <- data.frame(X1 = c(-2, 1, 2, 3), X2 = c(-4, -2, 0, 1),
                   Z = c(1, 1, 0, 0))
predict(fit, newdata = ndat)

}
