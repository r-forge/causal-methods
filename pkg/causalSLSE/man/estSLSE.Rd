\name{estSLSE}
\alias{estSLSE}
\alias{estSLSE.cslseModel}
\title{
Least squares estimate of \code{cslseModel} objects
}
\description{
This is the main function to estimate the \code{cslseModel} objects. It
generates the basis functions based on the knots specified in the model
and estimates it by least squares.
}
\usage{
\method{estSLSE}{cslseModel}(model, selKnots, ...)
}
\arguments{
  \item{model}{A model of class \code{cslseModel} created by the
    \code{cslseModel} function.}
  \item{selKnots}{An optional list of integers to select the knots from
    the list of knots specified by the model. It 
    must be a named list with the names begin either "treated",
    "nontreated" or both. By default, the knots from the model are used.}
  \item{...}{Additional arguments to pass to other methods. Currently
    not used.}
}

\examples{
data(simDat3)
mod1 <- cslseModel(Y~Z|~X1*X2, data=simDat3)
fit <- estSLSE(mod1)
}

\keyword{spline regression}
