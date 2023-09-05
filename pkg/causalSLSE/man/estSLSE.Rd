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
\method{estSLSE}{cslseModel}(model, w0=NULL, w1=NULL, ...)
}
\arguments{
  \item{model}{A model of class \code{cslseModel} created by the
    \code{cslseModel} function.}
  \item{w0}{An optional list of integers to select the knots for the
    control group from the list of knots specified by the model. By
    default, the knots from the model are used.}
  \item{w1}{An optional list of integers to select the knots for the
    treated group from the list of knots specified by the model. By
    default, the knots from the model are used.}
  \item{...}{Additional arguments to pass to other methods. Currently
    not used.}
}

\examples{
data(simDat3)
mod1 <- cslseModel(Y~Z|~X1*X2, data=simDat3)
fit <- estSLSE(mod1)
}

\keyword{spline regression}
