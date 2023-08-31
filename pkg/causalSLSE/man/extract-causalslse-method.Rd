\name{extract,causalslse-method}
\alias{extract,causalslse-method}
\alias{extract.causalslse}
\title{\code{\link{extract}} method for \code{causalslse} objects}
\usage{
\S4method{extract}{causalslse}(
  model,
  include.nobs = TRUE,
  include.nknots = TRUE,
  include.numcov = TRUE,
  include.rsquared = TRUE,
  include.adjrsquared=TRUE,
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

\item{include.adjrsquared}{Report the adjusted R-squared of the final
  regression estimtion}

\item{which}{Which causal effect measures should be printed?}

\item{...}{Custom parameters, which are handed over to subroutines. Currently
not in use.}
}
\description{
\code{\link{extract}} method for \code{causalslse} objects created by the
\code{\link[causalSLSE]{causalSLSE}} function.
}
