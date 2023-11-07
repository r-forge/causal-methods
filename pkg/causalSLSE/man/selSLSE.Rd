\name{selSLSE}
\alias{selSLSE}
\alias{selSLSE.cslseModel}
\title{
Knot selection method
}
\description{
This is the main function to select the knots in \code{cslseModel}
objects. It returns a model with non significant knots removed. 
}
\usage{
\method{selSLSE}{cslseModel}(model, selType=c("BLSE", "FLSE"),
        selCrit = c("AIC", "BIC", "PVT"), 
        pvalT = function(p) 1/log(p), vcov.=vcovHC, fortran=FALSE, ...)
}
\arguments{
  \item{model}{A model of class \code{cslseModel} created by
    the \code{cslseModel} function.}
  \item{selType}{The selection method: backward for forward.}
  \item{selCrit}{The criterion to select the piecewise polynomial knots.}
  \item{pvalT}{A function to determine the p-value threshold for the
    significance of the coefficients. It has to be a function of one
    parameter, which is the average number of knots in the model.}
  \item{vcov.}{An alternative function to compute the covariance matrix
    of the least squares estimators. The default is the \code{vcovHC}
    method for \code{lm} objects.}
  \item{fortran}{Should the selection be done using the experimental
    Fortran subroutines? If yes, the function returns a list of 3
    models, one for each \code{selCrit}. This alternative is faster for
    larger models and largewr sample sizes.}
  \item{...}{Additional arguments to pass to \code{vcov.}}
}

\examples{
data(simDat3)
mod1 <- cslseModel(Y~Z|~X1*X2, data=simDat3)
selSLSE(mod1, "FLSE")
}

\keyword{causal effects}
