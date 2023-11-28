\name{llSplines}
\alias{llSplines}
\alias{llSplines.slseModel}
\alias{llSplines.cslseModel}
\title{
Local linear splines generator for model objects.
}
\description{
It generates a matrix of basis functions using local linear splines. The
number of basis functions and the breaking points are determined by the
knot specifications of the \code{slseModel} or \code{cslseModel} model.
}
\usage{
\method{llSplines}{slseModel}(object, ...)
\method{llSplines}{cslseModel}(object, ...)
}
\arguments{
  \item{object}{A model of class \code{slseModel}, created by the
    \code{slseModel} function, or \code{cslseModel}, created by the
    \code{cslseModel} function.}
  \item{...}{Additional arguments to pass to other methods. Currently
    not used.}
}

\examples{
data(simDat3)
mod1 <- slseModel(Y~X1*X2, data=simDat3)
U <- llSplines(mod1)
fit <- lm(Y~U, data=simDat3)
}

\keyword{local linear spline function}
