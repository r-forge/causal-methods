\name{extract,cslse-method}
\alias{extract,cslse-method}
\alias{extract.cslse}
\title{\code{\link{extract}} Method for \code{cslse} Objects}
\usage{
\S4method{extract}{cslse}(
    model,
    include.nobs = TRUE,
    include.nknots = TRUE,
    include.numcov = TRUE,
    include.rsquared = TRUE,
    include.adjrs=TRUE,
    separated.rsquared = FALSE,    
    which = c("ALL","ACE","ACT","ACN","ACE-ACT","ACE-ACN","ACT-ACN"),
    ...)
}
\arguments{
\item{model}{A causal effect object created by
  \code{\link{causalSLSE}}.} 

\item{include.nobs}{Report the number of obervations?}

\item{include.nknots}{Report the total number of knots for each group?}

\item{include.numcov}{Report the total number of covariates (including
  interactions is any) per group?}

\item{include.rsquared}{Report the R-squared of the final regression
  estimtion}

\item{include.adjrs}{Report the adjusted R-squared of the final
  regression estimtion}

\item{separated.rsquared}{Should we print the R-squared separately for
  each group? This applies as well to the adjusted R-squared. }

\item{which}{Which causal effect measures should be printed?}

\item{...}{Custom parameters, which are handed over to subroutines. Currently
not in use.}
}
\description{
\code{\link{extract}} method for \code{cslse} objects created by the
\code{\link[causalSLSE]{causalSLSE}} function.
}

\value{
It returns an object of class \code{\link{texreg}}.  
}
