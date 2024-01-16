\name{slseKnots}
\alias{slseKnots}
\title{
Knots Creator for Basis Functions
}
\description{
The function creates an object of class \code{slseKnots} for
semiparametric least squares models. It returns an object of 
class \code{slseKnots}.
}
\usage{
slseKnots(form, data, X, nbasis = function(n) n^0.3, 
          knots)
}
\arguments{
  \item{form}{A formula that determines the covariates to include in
    the regression. It is a regular formula and only the right hand side
    is considered.} 
  \item{data}{A \code{data.frame} with all variables included in
    \code{form}. It is required when \code{X} is missing.}
  \item{nbasis}{A function to determine the number of basis
    functions. It has to be a function of one argument, the sample
    size.} 
  \item{knots}{An optional list of knots. Each element of the list (one
    for each covariate) is either a vector of numeric knots or it is
    \code{NULL}. The latter implies that the number of knots is set to
    0. If missing (the default), the knots are automatically generated
    (see the details below).}
  \item{X}{An optional matrix of covariates. When provided, \code{form}
    is not needed.} 
}

\details{The automatic selection is as follows. It is applied to each
  variable created by the \code{\link{model.matrix}} of \code{form} except for
  the intercept:

  The number of knots is the ceiling of what the \code{nbasis} function
  returns minus 1. Let \code{p-1} be the number of knots. Then, we
  compute the \code{p+1} empirical quantiles of the variable for equally
  spaced probabilities going from 0 to 1 and drop the first and last
  ones. This is done using the function \code{\link{quantile}} with
  \code{type=1}. We then remove the duplicated values and the ones equal
  to either the \code{min} or the \code{max} of the variable. if the number
  of remaining knots is equal to 0, the set of knots for this variable
  is set to \code{NULL}.

  For manual selection, see the vignette.
}

\value{

  It returns an object of class \code{slseKnots}. It is a list for which
  the length and names are respectively the number of columns and the
  column names of the \code{\link{model.matrix}} of \code{form} after
  the intercept has been removed. Each element is a numeric vector of
  knots, unless the number of knots is set to 0, in  which case it is
  \code{NULL}. 

}

\examples{
data(simDat3)
k <- slseKnots(Y ~ X1 * X2, data = simDat3)
k

## We can extract the set for one variable

k$X1

}

\keyword{causal effects}
