\name{causalSLSE}
\alias{causalSLSE}
\alias{causalSLSE.cslseModel}
\alias{causalSLSE.cslseFit}
\alias{causalSLSE.formula}
\title{
Causal Effect based on semiparametric Least Squares models
}
\description{
This is the main method to estimate the causal effects using the semiparametric
least squares method. It returns an object of class
\code{cslse} and is registered for objects of class
\code{\link{cslseModel}} and \code{\link[=estSLSE]{cslseFit}}.
}
\usage{
\method{causalSLSE}{cslseModel}(object,
                               selType=c("SLSE","BLSE","FLSE"),
                               selCrit = c("AIC", "BIC", "PVT"),
                               selVcov = c("HC0", "Classical", "HC1", "HC2", "HC3"),    
                               causal = c("ALL","ACT","ACE","ACN"),
                               pvalT = function(p) 1/log(p),
                               vcov.=vcovHC, reSelect=FALSE, ...)
\method{causalSLSE}{cslseFit}(object, causal = c("ALL","ACT","ACE","ACN"),
                             vcov.=vcovHC, ...)
\method{causalSLSE}{formula}(object, data, nbasis=function(n) n^0.3,
                             knots, 
                             selType=c("SLSE","BLSE","FLSE"),
                             selCrit = c("AIC", "BIC", "PVT"),
                             selVcov = c("HC0", "Classical", "HC1", "HC2", "HC3"),        
                             causal = c("ALL","ACT","ACE","ACN"),
                             pvalT = function(p) 1/log(p),
                             vcov.=vcovHC, reSelect=FALSE, ...)
}
\arguments{
  \item{object}{An object of class \code{cslseModel} created by the
    \code{\link{cslseModel}} function, \code{cslseFit} created by \code{\link{estSLSE}}
    or \code{formula}.}
  \item{data}{A \code{data.frame} with all variables included in \code{form}.}
  \item{nbasis}{A function to determined the number of basis
    functions. It has to  be a function of one argument, the sample
    size.} 
  \item{knots}{A list of knots for the treated and nontreated
    groups. The list must be named using the group names. Each element
    of the list is also a list  of length  equal to the number of
    confounders. The choice for each confounders is NULL for no knots or
    numeric for specific values. If missing, the knots are 
    automatically generated.}
  \item{selType}{The method for selecting the knots. By default (SLSE),
    all knots from the model are used.} 
  \item{selCrit}{The criterion to select the knots.}
  \item{causal}{What causal effect should we compute.}
  \item{pvalT}{A function to determine the p-value threshold for the
    significance of the coefficients. It has to be a function of one
    parameter, which is the average number of knots in the model. This
    value may differ across treatment group.}
  \item{selVcov}{The type of least squares covariance matrix used to
    compute the p-values needed for the selection.}  
  \item{vcov.}{An alternative function to compute the covariance matrix
    of the least squares estimators. The default is
    \code{\link[=vcovHC.default]{vcovHC}}. This function is used to
    compute the standard errors of the causal effect estimators and the
    SLSE coefficients.}   
  \item{reSelect}{By default, the stored selections are used. If
    \code{reSelect} is set to \code{TRUE}, the selection is re-computed.}  
  \item{...}{Additional arguments to pass to \code{vcov.}}
}

\value{It returns an object of class \code{cslse}, which inherits from
  the class \code{cslseFit}. It is a list with the following elements:

  \item{treated, nontreated}{They are objects of class \code{slseFit}
  obtained by \code{\link{estSLSE}}.}

  \item{ACE, ACT, ACN}{Estimates of the average causal effect, the
  causal effect on the treated and the causal effect on the
  nontreated. Each of them is a vector of two elements: the estimate and
  its estimated standard error. All three are included only if the
  argument \code{causal} is set to \code{"ALL"}.}

Also, the object contains the following additional attributes:

\item{treatedVar}{The name of the variable in the dataset that represents
  the treatment indicator.}

\item{groupInd}{A named vector with the value of the treatment indicator
corresponding do each treatment group.}
  
}

\seealso{\code{\link{estSLSE}} for the  estimation of the model,
  \code{\link{slseKnots}} for the format of knots, and \code{\link{selSLSE}}
  and \code{\link[=update.cslseModel]{update}} for the knots selection
  and to understand how stored selections are used.}

\examples{
data(simDat3)

## A causal SLSE model with the outcome Y
## the treatment indicator Z and the confounders X1, X2 and X1:X2

mod1 <- cslseModel(Y ~ Z | ~ X1 * X2, data = simDat3)

## The causal effects are estimated using the backward method and AIC criterion
## The HC1 type is used for the least squares covariance matrix

fit1 <- causalSLSE(mod1, selType = "BLSE", type = "HC1")

## This is the same for formula objects

fit2 <- causalSLSE(Y ~ Z | ~ X1 * X2, data = simDat3, selType = "BLSE", type = "HC1")
}

\keyword{causal effects}
