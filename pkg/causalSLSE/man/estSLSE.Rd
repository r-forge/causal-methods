\name{estSLSE}
\alias{estSLSE}
\alias{estSLSE.cslseModel}
\alias{estSLSE.slseModel}
\title{
Least squares estimate of \code{cslseModel} or \code{slseModel} objects
}
\description{
This is the main function to estimate \code{cslseModel} or
\code{slseModel} objects. It generates the basis functions based on the
knots specified in the model and estimates it by least squares.
}
\usage{
\method{estSLSE}{cslseModel}(model, selKnots, ...)
}
\arguments{
  \item{model}{A model of class \code{cslseModel} created by the
    \code{cslseModel} function.}
  \item{selKnots}{An optional list of integers to select the knots from
    the list of knots specified by the model. If the model is a
    \code{cslseModel} object, it  must be a named list with the names
    begin either "treated",  "nontreated" or both. By default, the knots
    from the model are used.} 
  \item{...}{Additional arguments to pass to other methods. Currently
    not used.}
}

\examples{
data(simDat3)
mod1 <- cslseModel(Y~Z|~X1*X2, data=simDat3)
fit1 <- estSLSE(mod1)
mod2 <- slseModel(Y~X1*X2, data=simDat3)
fit2 <- estSLSE(mod2)
}

\keyword{spline regression}
