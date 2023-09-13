\name{cslseModel}
\alias{cslseModel}
\title{
Semiparametric Least Squares Estimator model
}
\description{
This function creates an object of class \code{slseModel}. The model
contains all specifications about the model being estimated including
the starting knots for the control and treated group for the piecewise
polynomial.
}
\usage{
cslseModel(form, data, nbasis = function(n) n^0.3, 
         knots)
}
\arguments{
  \item{form}{A formula for the regression (outcome versus treatment
    indicator) and a formula for the covariates to include in the
    regression, separated by \code{|}. See details.}
  \item{data}{A \code{data.frame} with all variables included in \code{form}.}
  \item{nbasis}{A function to determined the number of bases. It has to
    be a function of one argument, the sample size. To fix it, simply
    input a function that returns an integer.}
  \item{knots}{A list of knots for the treated and nontreated
    groups. The list must be named using the group names. Each element
    of the list is also a list  of length  equal to the number of
    covariates. The choice for each covariate is NULL for no knots or
    numeric to provide specific values. If missing, the knots are 
    automatically generated.}
}

\examples{
data(simDat3)
mod1 <- cslseModel(Y~Z|~X1*X2, data=simDat3)
mod1
}

\keyword{causal effects}
