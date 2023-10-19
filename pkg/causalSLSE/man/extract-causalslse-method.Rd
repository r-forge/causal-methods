\name{extract,cslse-method}
\alias{extract,cslse-method}
\alias{extract.cslse}
\title{\code{\link{extract}} method for \code{cslse} objects}
\usage{
\S4method{extract}{cslse}(
  model,
  include.nobs = TRUE,
  include.nknots = TRUE,
  include.numcov = TRUE,
  include.rsquared = TRUE,
  include.adjrs=TRUE,
  which = c("ALL","ACE","ACT","ACN","ACE-ACT","ACE-ACN","ACT-ACN"),
  ...)
}
\arguments{
\item{model}{A statistical model object.}

\item{include.nobs}{Report the number of obervations?}

\item{include.nknots}{Report the total number of knots for each group?}

\item{include.numcov}{Report the total number of covariates (including
  interactions is any) per group?}

\item{include.rsquared}{Report the R-squared of the final regression
  estimtion}

\item{include.adjrs}{Report the adjusted R-squared of the final
  regression estimtion}

\item{which}{Which causal effect measures should be printed?}

\item{...}{Custom parameters, which are handed over to subroutines. Currently
not in use.}
}
\description{
\code{\link{extract}} method for \code{cslse} objects created by the
\code{\link[causalSLSE]{causalSLSE}} function.
}
