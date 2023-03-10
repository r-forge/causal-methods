\name{polSelect}
\alias{selASY}
\alias{selIC}
\title{
Polynomial selection methods.
}
\description{
This is a collection of methods to select the piecewise polynomials.
}
\usage{
selASY(form, data, pFact = 0.3, splineMet = c("manual", "bs"),
       HCtype="HC", mZeroProp=0.1, knots0=NA, knots1=NA,
       minPV=function(p) 1/(p*log(p)))
selIC(form, data, pFact = 0.3, type=c("AIC", "BIC", "CV"),
      splineMet = c("manual", "bs"), HCtype="HC",
      mZeroProp=0.1, knots0=NA, knots1=NA) 
}
\arguments{
  \item{form}{A formula to identify the outcome, the treatment
    indicator and the covariate vectors to be converted into
    spline matrices. See the example below}
  \item{data}{A \code{data.frame} that contains all variables from the
    formula \code{form}}
  \item{pFact}{The maximum number of knots when the argument \code{knots} is set to
    \code{NA} if \code{n^pFact}, where n is the length of \code{X}.}    
  \item{type}{The type of informatrion criterion. "CV" stands for cross-validation.}
  \item{splineMet}{Should the method be homemade (manual) of based on the
    \code{bs} function from the splines package?}
  \item{HCtype}{This is the type of HCCM matrix to use to for
    inference. The default, "HC" is HC0 (See \code{\link{vcovHC}})}
  \item{mZeroProp}{If the proportion of zeros is greater than this
    value, the knots and spline matrices are based on the non-zero
    observations. This is particularly useful when binary variables are
    interacted with continuous variables.}
  \item{minPV}{A function to determine the cuttoff point for the
    significance. the argument of the function is the total number of
    basis.}
  \item{knots0}{The piecewise polynomial knots for the control group.
    If set to \code{NA}, the knots are set to the \code{p} equally
    spaced quantiles of \code{X}. If \code{NULL}, the function returns
    \code{X}.}
  \item{knots1}{The piecewise polynomial knots for the treated group.
    If set to \code{NA}, the knots are set to the \code{p} equally
    spaced quantiles of \code{X}. If \code{NULL}, the function returns
    \code{X}.}    
}

\keyword{selection, polynomial}

\examples{
data(simData)
## Y is the outcome, Z the treatment indicator
## The covariate vectors are X1, X2 and X1:X2
sel <- selASY(Y~Z|~X1*X2, simData, .3, "manual", "HC0")
sel$knots0
sel$knots1
}

