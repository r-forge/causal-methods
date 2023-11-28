\name{cslseModel}
\alias{cslseModel}
\alias{slseModel}
\title{
Semiparametric Least Squares Estimator model
}
\description{
This function creates an object of class  \code{cslseModel} or
\code{slseModel}. The model contains all specifications about the model
being estimated including the starting knots for each group, if any, for
the local linear polynomial. 
}
\usage{
cslseModel(form, data, nbasis = function(n) n^0.3, 
           knots, groupn=c(treated=1, nontreated=0))
slseModel(form, data, nbasis = function(n) n^0.3, 
           knots)
}
\arguments{
  \item{form}{A formula for the regression in  the case of a simple spline
    regression, or two formulas separated by \code{|}, one for the
    outcome versus the treatment indicator and one of the
    confounders. See vignette for details.} 
  \item{data}{A \code{data.frame} with all variables included in \code{form}.}
  \item{nbasis}{A function to determined the number of bases. It has to
    be a function of one argument, the sample size. To fix it, simply
    input a function that returns an integer.}
  \item{knots}{An optional list of knots. For \code{cslseModel}, it must 
    be a named list using the group names (treated and nontreated). Each element  
    of the list is also a list  of length  equal to the number of
    covariates. The choice for each covariate is NULL for no knots or
    numeric to provide specific values. If missing, the knots are 
    automatically generated.}
  \item{groupn}{A named vector with the group names and values. By
    default, the treatment indicator defined in \code{form} is equal to
    1 for the treated group and 0 for the nontreated. For now, these are
    the only allowed names, but the values can differ. In particular,
    the treatment indicator may be a factor with character values.}
}  
\examples{
data(simDat3)
mod1 <- cslseModel(Y~Z|~X1*X2, data=simDat3)
mod1
mod2 <- slseModel(Y~X1*X2, data=simDat3)
mod2
}

\keyword{causal effects}
