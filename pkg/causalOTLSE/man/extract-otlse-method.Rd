\name{extract,otlse-method}
\alias{extract,otlse-method}
\alias{extract.otlse}
\title{\code{\link{extract}} method for \code{otlse} objects}
\usage{
\S4method{extract}{otlse}(
  model,
  include.nobs = TRUE,
  include.nknots = TRUE,
  include.numcov = TRUE,
  which = c("ALL","ACE","ACT","ACN","ACE-ACT","ACE-ACN","ACT-ACN"),
  ...)
}
\arguments{
\item{model}{A statistical model object.}

\item{include.nobs}{Report the number of obervations?}

\item{include.nknots}{Report the total number of knots for each group?}

\item{include.numcov}{Report the total number of covariates (including
  interactions is any) per group?}

\item{which}{Which causal effect measures should be printed?}

\item{...}{Custom parameters, which are handed over to subroutines. Currently
not in use.}
}
\description{
\code{\link{extract}} method for \code{otlse} objects created by the
\code{\link[causalOTLSE]{otlse}} function.
}
