\name{extract}
\alias{extract,slseFit-method}
\alias{extract.slseFit}
\alias{extract,cslse-method}
\alias{extract.cslse}
\alias{extract,altCausal-method}
\title{\code{\link{extract}} methods for some objects.}
\usage{
\S4method{extract}{slseFit}(
  model,
  include.rsquared = TRUE,
  include.adjrs = TRUE, 
  include.nobs = TRUE,
  include.fstatistic = FALSE,
  include.rmse = FALSE,
  ...)

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

\S4method{extract}{altCausal}(
    model,
    include.nobs = TRUE,
    ...)

}
\arguments{
\item{model}{A model object.} 

\item{include.nobs}{Report the number of obervations?}

\item{include.nknots}{Report the total number of knots for each group?}

\item{include.numcov}{Report the total number of covariates (including
  interactions is any) per group?}

\item{include.rsquared}{Report the R-squared of the final regression
  estimation}

\item{include.adjrs}{Report the adjusted R-squared of the final
  regression estimtion}

\item{separated.rsquared}{Should we print the R-squared separately for
  each group? This applies as well to the adjusted R-squared. }

\item{which}{Which causal effect measures should be printed?}

\item{include.fstatistic}{Report the F-statistics?}

\item{include.rmse}{Report the RMSE?}

\item{...}{Custom parameters, which are handed over to subroutines. Currently
not in use.}

}
\description{
\code{\link{extract}} method for \code{slseFit} objects created by the
\code{\link[causalSLSE]{estSLSE}} function.

\code{\link{extract}} method for \code{cslse} objects created by the
\code{\link[causalSLSE]{causalSLSE}} function.

\code{\link{extract}} method for \code{altCausal} objects created by the
functions in the list of alternative methods: \code{\link{altCausal}}.
}

\value{
It returns an object of class \code{\link{texreg}}.  
}
