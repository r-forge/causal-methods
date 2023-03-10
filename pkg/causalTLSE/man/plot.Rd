\name{plot}
\alias{plot.otlse}
\title{Plot methods}
\description{
Some plot utilities for object of class \code{otlse}.
}
\usage{
\method{plot}{otlse}(x, y, which=y, addInterval=FALSE, level=0.95,
                    leg.pos="topright", ...)
}
\arguments{
  \item{x}{Object of class \code{otlse}.}
  \item{y}{alias for \code{which} for compatibility with \code{plot}}
  \item{which}{Which covariate to plot against the outcome variable. It
    could be an integer of a character}
  \item{addInterval}{Should we add a confidence interval?}
  \item{level}{The confidence interval level if included.}
  \item{leg.pos}{The legend position.}
  \item{...}{Currently not used.}  
}

\examples{
data(simData)
fit <- otlse(Y~Z|~X1*X2*male, simData, crit="AIC", HCtype="HC3")
plot(fit,"X1")
}












