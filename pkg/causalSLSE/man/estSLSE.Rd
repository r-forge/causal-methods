\name{estSLSE}
\alias{estSLSE}
\alias{estSLSE.cslseModel}
\alias{estSLSE.slseModel}
\title{
Least Squares Estimate of \code{cslseModel} or \code{slseModel} Objects
}
\description{
This is the main function to estimate \code{cslseModel} or
\code{slseModel} objects. It generates the basis functions based on the
knots specified in the model and estimates it by least squares.
}
\usage{
\method{estSLSE}{cslseModel}(model, selKnots, ...)

\method{estSLSE}{slseModel}(model, selKnots, ...)
}
\arguments{
  \item{model}{A model of class \code{\link{cslseModel}}  or
    \code{\link{slseModel}}.} 
  \item{selKnots}{An optional list of integers to select the knots from
    the list of knots specified by the model. If the model is a
    \code{cslseModel} object, it  must be a named list with the names
    being either "treated",  "nontreated" or both. By default, the knots
    from the model stored in the element \code{knots} are used.} 
  \item{...}{Additional arguments to pass to other methods. Currently
    not used.}
}

\details{

  The method for \code{slseModel} objects generates the matrix
of basis functions implied by the set of knots included in the model and 
estimate the model by the least squares. Let \eqn{Y}{Y} be the
outcome and \eqn{U}{U} be the matrix of basis functions. Then, the
function estimates the model using the code \code{lm(Y~U)}.

For \code{cslseModel}, we could estimate the model using
\code{lm(Y~Z+I(Z-1)+I(U0*(1-Z))+I(U1*Z))}, where \code{Z} is a binary
variable equal to 1 for the treated and 0 for the nontreated, and
\code{U0} and \code{U1} are the matrices of basis functions for the
nontreated and treated, but the model is estimated separately for each
group. Therefore, the function \code{estSLSE.cslseModel} calls the
function \code{estSLSE.slseModel} for each \code{slseModel} objects
included in the  \code{cslseModel} object.
}

\value{

It returns an object of class \code{slseFit} or \code{cslseFit}
  depending on which method is called. An object of class \code{slseFit}
  is a list with the following elements: 

\item{LSE}{This is the least squares estimate of the semiparametric
  model. It is an object of class \code{\link[stats]{lm}}.}

\item{model}{An object of class \code{slseModel}. It is the model that
  is being estimated. The model may have been created by
  \code{\link{slseModel}} or modified by \code{\link{selSLSE}} or
  \code{\link[=update.slseModel]{update}}.}

An object of class \code{cslseFit} is a list of \code{slseFit} objects,
one for each treatment group. It also contains the following additional
attributes:

\item{treatedVar}{The name of the variable in the dataset that represents
  the treatment indicator.}

\item{groupInd}{A named vector with the value of the treatment indicator
corresponding do each treatment group.}
 
}

\examples{
data(simDat3)

## Estimating a causal semiparametric model

mod1 <- cslseModel(Y ~ Z | ~ X1 * X2, data = simDat3)
fit1 <- estSLSE(mod1)

## Estimating a semiparametric model

mod2 <- slseModel(Y ~ X1 * X2, data = simDat3)
fit2 <- estSLSE(mod2)
}

\keyword{spline regression}
