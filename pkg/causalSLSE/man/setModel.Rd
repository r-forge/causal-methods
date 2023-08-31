\name{setModel}
\alias{setModel}
\title{
Thresholding Least Squares estimator model
}
\description{
This function creates an object of class \code{slseModel}. The model
contains all specifications about the model being estimated including
the starting knots for the control and treated group for the piecewise
polynomial.
}
\usage{
setModel(form, data, nbasis = function(n) n^0.3, 
         knots0, knots1)
}
\arguments{
  \item{form}{A formula for the regression (outcome versus treatment
    indicator) and a formula for the covariates to include in the
    regression, separated by \code{|}. See details.}
  \item{data}{A \code{data.frame} with all variables included in \code{form}.}
  \item{nbasis}{A function to determined the number of bases. It has to
    be a function of one argument, the sample size. To fix it, simply
    input a function that returns an integer.}
  \item{knots0}{A list of knots for the control group. The length must
    be equal to the number of covariates. The number of knots is 0 if an
    element of the list is set to \code{NULL}. If missing, the knots are
  automatically generated.}
  \item{knots1}{A list of knots for the control group. See \code{knots0}.}  
}

\examples{
data(simDat3)
mod1 <- setModel(Y~Z|~X1*X2, data=simDat3)
mod1
}

\keyword{causal effects}
