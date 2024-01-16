\name{llSplines}
\alias{llSplines}
\alias{llSplines.slseModel}
\alias{llSplines.cslseModel}
\title{
Local Linear Splines Generator for Model Objects
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
  \item{object}{A model of class \code{\link{slseModel}} or
    \code{\link{cslseModel}}.}  
  \item{...}{Additional arguments to pass to other methods. Currently
    not used.}
}

\value{ The function returns a matrix of basis functions used to
estimate the semiparametric model. It is used directly as regressor in
\code{\link{lm}} as shown in the example below.}

\examples{

data(simDat3)

## We manually estimate the semiparametric model

mod1 <- slseModel(Y ~ X1 * X2, data = simDat3)
U <- llSplines(mod1)
fit1 <- lm(Y ~ U, data = simDat3)

## We use estSLSE instead (results are identical)

fit2 <- estSLSE(mod1)

}

\keyword{local linear splines}
