\name{selSLSE}
\alias{selSLSE}
\alias{selSLSE.cslseModel}
\alias{selSLSE.slseModel}
\title{
Knot selection method
}
\description{
This is the main function to select the knots in \code{cslseModel} or
\code{slseModel} objects. It returns a model with non significant knots
removed.  
}
\usage{
\method{selSLSE}{cslseModel}(model, selType=c("BLSE", "FLSE"),
        selCrit = c("AIC", "BIC", "PVT"), 
        pvalT = function(p) 1/log(p),
        vcovType = c("HC0", "Classical", "HC1", "HC2", "HC3"),
        reSelect=FALSE, ...)
\method{selSLSE}{slseModel}(model, selType=c("BLSE", "FLSE"),
        selCrit = c("AIC", "BIC", "PVT"), 
        pvalT = function(p) 1/log(p),
        vcovType = c("HC0", "Classical", "HC1", "HC2", "HC3"),
        reSelect=FALSE, ...)
}
\arguments{
  \item{model}{A model of class \code{cslseModel} created by
    the \code{cslseModel} function or \code{slseModel} created by
    the \code{slseModel} function.}
  \item{selType}{The selection method: backward for forward.}
  \item{selCrit}{The criterion to select the piecewise polynomial knots.}
  \item{pvalT}{A function to determine the p-value threshold for the
    significance of the coefficients. It has to be a function of one
    parameter, which is the average number of knots in the model.}
  \item{vcovType}{The type of least squares covariance matrix used to
    compute the p-values needed for the selection.}
  \item{reSelect}{The previous selections are stored in the model
    object. By default, it is the previous selections that is used. If
    \code{reSelect} is set to \code{TRUE}, the selection is re-computed.}
  \item{...}{Additional arguments to pass to other methods. Currently
    not used.}
}

\examples{
data(simDat3)
mod1 <- cslseModel(Y~Z|~X1*X2, data=simDat3)
selSLSE(mod1, "FLSE")
}

\keyword{causal effects}
