\name{estModel}
\alias{estModel}
\title{
Least squares estimate of \code{tlseModel} objects
}
\description{
This is the main function to estimate the \code{tlseModel} objects. It
generates the basis functions based on the knots specified in the model
and estimates it by least squares.
}
\usage{
estModel(model, w0=NULL, w1=NULL)
}
\arguments{
  \item{model}{A model of class \code{tlseModel} created by the
    \code{setModel} function.}
  \item{w0}{An optional list of integers to select the knots for the
    control group from the list of knots specified by the model. By
    default, the knots from the model are used.}
  \item{w1}{An optional list of integers to select the knots for the
    treated group from the list of knots specified by the model. By
    default, the knots from the model are used.}  
}

\examples{
data(simDat3)
mod1 <- setModel(Y~Z|~X1*X2, data=simDat3)
fit <- estModel(mod1)
}

\keyword{spline regression}
