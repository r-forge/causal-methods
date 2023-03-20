\name{extract,tlse-method}
\alias{extract,tlse-method}
\alias{extract.tlse}
\title{\code{\link{extract}} method for \code{tlse} objects}
\usage{
\S4method{extract}{tlse}(
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
\code{\link{extract}} method for \code{tlse} objects created by the
\code{\link[causalTLSE]{tlse}} function.
}
