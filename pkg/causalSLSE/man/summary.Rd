\name{summary}
\alias{summary.cslse}
\alias{summary.cslseFit}
\alias{summary.slseFit}
\title{Summary Method for Fitted Models}
\description{
The method computes summary statistics for estimated semiparametric and
causal models. 
}
\usage{
\method{summary}{cslse}(object, ...)
\method{summary}{cslseFit}(object, vcov.=vcovHC, ...)
\method{summary}{slseFit}(object, vcov.=vcovHC, ...)
}
\arguments{
  \item{object}{A model estimated by \code{\link{estSLSE}} or
    \code{\link{causalSLSE}}.} 
  \item{vcov.}{A function to compute the strandard error of the
    least squares coefficients. The default is
    \code{\link[sandwich]{vcovHC}}.}  
\item{...}{Argument for other types of objects}
}

\value{

The function \code{summary.slseFit} returns the summary statistics of
the least squares estimation of \code{\link{slseModel}} objects
estimated by \code{\link{estSLSE}}. It is an object of class
\code{summary.slseFit}, which is a list with the following elements:

\item{model}{The \code{\link{slseModel}} model being estimated.}

\item{lseSum}{An object of class \code{\link[stats]{summary.lm}}.}

The function \code{summary.cslseFit} returns an object of class
\code{summary.cslseFit}, which is a list of \code{summary.slseFit} objects,
one for each treatment group. 

The function \code{summary.cslse} returns the summary statistics of the
causal effect estimates. It is an object of class \code{summary.cslse},
which is a list with the following elements:

\item{causal}{A matrix with the causal effect estimates in the first
  column, their standard errors in the second, their t-ratios in the
  third and p-values in the fourth.}

\item{lse}{A list of \code{summary.slseFit} objects, one for each
  treatment group.}

}


\examples{

## A causal semiparametric model with causal effect estimates 

data(simDat3)
mod1 <- cslseModel(Y ~ Z | ~ X1 * X2, data = simDat3)
CE <- causalSLSE(mod1, selType = "BLSE")
summary(CE)

## A causal semiparametric model

fit <- estSLSE(mod1)
summary(fit)

## A semiparametric model

mod2 <- slseModel(Y ~ X1 * X2, data = simDat3)
fit2 <- estSLSE(mod2)
summary(fit2)
}

