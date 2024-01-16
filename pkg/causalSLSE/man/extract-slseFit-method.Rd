\name{extract,slseFit-method}
\alias{extract,slseFit-method}
\alias{extract.slseFit}
\title{\code{\link{extract}} Method for \code{slseFit} Objects}
\usage{
\S4method{extract}{slseFit}(
  model,
  include.rsquared = TRUE,
  include.adjrs = TRUE, 
  include.nobs = TRUE,
  include.fstatistic = FALSE,
  include.rmse = FALSE,
  ...)

}
\arguments{
\item{model}{A model estimated by \code{\link{estSLSE}}.} 

\item{include.rsquared}{Report the R-squared of the final regression
  estimation}

\item{include.adjrs}{Report the adjusted R-squared of the final
  regression estimtion}

\item{include.nobs}{Report the number of obervations?}

\item{include.fstatistic}{Report the F-statistics?}

\item{include.rmse}{Report the RMSE?}

\item{...}{Custom parameters, which are handed over to subroutines. Currently
not in use.}

}
\description{
\code{\link{extract}} method for \code{slseFit} objects created by the
\code{\link[causalSLSE]{estSLSE}} function.
}

\value{
It returns an object of class \code{\link{texreg}}.  
}
