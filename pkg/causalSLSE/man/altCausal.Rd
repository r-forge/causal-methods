\name{altCausal}
\alias{altCausal}
\alias{matching}
\alias{LLmatching}
\alias{ipw}
	
\title{Alternative Causal Effect estimation methods.}

\description{This documentation file presents a collection of popular
methods used to estimate the average causal effect, the causal effect on
the treated and the causal effect on the nontreated.
} 

\usage{

matching(form,  balm, data, type=c("ACE","ACT","ACN"), M=4,
         psForm=NULL, bcForm=NULL, vJ=4)

LLmatching(form, psForm, data, type=c("ACE","ACT","ACN"),
           kern=c("Gaussian","Epanechnikov"),tol=1e-4,
           h=NULL, from=.00001, to=5, ngrid=10, maxit=100,
           hMethod=c("Brent","Grid"))

ipw(form, psForm, data, type=c("ACE","ACT","ACN"),
    normalized=FALSE, tolPS=0, ...)
}
\arguments{

  \item{form}{A formula that links the outcome variable to the treatment
    indicator. For the moment, only one treatment group is allowed.}
  
  \item{balm}{A formula or a matrix with balancing covariates to be
    matched.}

  \item{data}{A data.frame or a matrix with column names.}

  \item{type}{The type of causal effect to compute. \code{ACE} stands
    for average causal effect, \code{ACT} for causal effect on the treated
    and \code{ACN} for causal effect on the control or nontreated.}

  \item{M}{The minimum number of matches}
  
  \item{psForm}{It is the \code{formula} argument passed to
    \code{\link{glm}} to estimate the propensity scores by a logistic
    regression.}
  
  \item{bcForm}{A formula that represents the right hand side in the
    regression used for the bias correction.}

  \item{kern}{The type of kernel to use in the local linear regression
    method.}
    
  \item{tol}{The tolerance level for the stopping rule used to compute
    the optimal bandwidth.}

  \item{tolPS}{Observations for which the propensity scores are not
  between \code{tolPS}  and \code{(1-tolPS)} are removed. It
  is meant to satisfy the assumption that the propensity score is
  bounded away from, 0 and 1. The default is 0, in which case
  observations with negative or greater than 1 propensity scores are
  dropped.}

  \item{h}{A fixed bandwidth. By default, the optimal bandwidth is found
    by minimizing a cross-validation.}
  
  \item{from}{The lower bound for the search of the optimal bandwidth.}

  \item{to}{The upper bound for the search of the optimal bandwidth.}

  \item{ngrid}{The number of grid points if the optimal bandwidth is
    obtained by grid search.}

  \item{maxit}{The maximum number of iterations for the minimization of
    the cross-validation.}

  \item{hMethod}{The method used to find the optimal bandwidth.} 

  \item{normalized}{Should the weights be normalized. If set to
    \code{GPE}, the GPE method is used.}

  \item{vJ}{The minimum number of matches for the standard error estimation}
  
  \item{...}{Additional arguments to pass to \code{\link{optim}} when
    the \code{GPE} method is used.}
	    
}

\value{
  All methods return an object of classes \code{"causalfit"}. It is a
  list with the following elements:

  \item{estim}{The causal effect estimate.}
  \item{se}{The standard error estimate of the causal effect.}
  \item{type}{The type of causal effect. It is the value of the
    argument \code{type}.}
  \item{method}{A description of the method.}
  \item{form}{A list of formulas used for the estimation.}
  \item{details}{A more detailed description of the method.}
  \item{info}{Convergence information if relevant.}  
 }

\examples{
data(simDat3)

balm <-  ~ X1 * X2
g <- Y ~ Z
ps <- Z~X1*X2

## A print and summary methods exist for this class
           
fit1 <- matching(form=g, balm=balm, data=simDat3, type="ACE", bcForm=balm)
fit1

fit2 <- LLmatching(form=g, psForm=ps, data=simDat3, type="ACE")
fit2

fit3 <- ipw(form=g, psForm=ps, data=simDat3, type="ACE", normalized=TRUE)
fit3

## A print and summary methods exist for this class

summary(fit1)

}

