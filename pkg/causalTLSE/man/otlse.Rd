\name{otlse}
\alias{otlse}
\title{
Optimal Thresholding Least Squares
}
\description{
This is the main function to estimate the causal effects using the
optimal thresholding least squares method.
}
\usage{
otlse(form, data, crit = c("ASY", "AIC", "BIC", "CV", "NONE"),
      pFact=0.3, splineMet=c("manual","bs"), HCtype="HC",
      mZeroProp=0.1, knots0=NA, knots1=NA, ...)
}
\arguments{
  \item{form}{A formula to identify the outcome, the treatment
    indicator and the covariate vectors to be converted into
    spline matrices. See the example below}
  \item{data}{A \code{data.frame} that contains all variables from the
    formula \code{form}}
  \item{crit}{The method to select the piecewise polynomial knots. If
    set to "NONE", the all bases are kept. It is included for
    experimentation purpose.}  
  \item{pFact}{The maximum number of knots when the argument \code{knots} is set to
    \code{NA} if \code{n^pFact}, where n is the length of \code{X}.}    
  \item{splineMet}{Should the method be homemade (manual) of based on the
    \code{bs} function from the splines package?}
  \item{HCtype}{This is the type of HCCM matrix to use to for
    inference. The default, "HC" is HC0 (See \code{\link{vcovHC}})}
  \item{mZeroProp}{If the proportion of zeros is greater than this
    value, the knots and spline matrices are based on the non-zero
    observations. This is particularly useful when binary variables are
    interacted with continuous variables.}
  \item{knots0}{The piecewise polynomial knots for the control group.
    If set to \code{NA}, the knots are set to the \code{p} equally
    spaced quantiles of \code{X}. If \code{NULL}, the function returns
    \code{X}.}
  \item{knots1}{The piecewise polynomial knots for the treated group.
    If set to \code{NA}, the knots are set to the \code{p} equally
    spaced quantiles of \code{X}. If \code{NULL}, the function returns
    \code{X}.}      
  \item{...}{Other arguments to pass to the basis selection function.}
}

\examples{
data(simData)
## Y is the outcome, Z the treatment indicator
## The covariate vectors are X1, X2 and X1:X2
fit <- otlse(Y~Z|~X1*X2*male, simData, crit="AIC", HCtype="HC3")
fit
}

\keyword{causal effects}
