\name{plot}
\alias{plot.tlse}
\title{Plot methods}
\description{
Some plot utilities for object of class \code{otlse}.
}
\usage{
\method{plot}{tlse}(x, y, which=y, interval=c("none","confidence"),
                    level=0.95, legendPos="topright",
                    vcov.=NULL, ...)
}
\arguments{
  \item{x}{Object of class \code{tlse}.}
  \item{y}{alias for \code{which} for compatibility with \code{plot}}
  \item{which}{Which covariate to plot against the outcome variable. It
    could be an integer or a character}
  \item{interval}{The type of confidence interval. The default is \code{none}.}
  \item{level}{The confidence interval level if included.}
  \item{legendPos}{The legend position.}
  \item{vcov.}{An alternative function to compute the covariance
    matrix of the least squares estimates. The default \code{vcov}
    method for \code{lm} objects.}
  \item{...}{Additional argument to pass to the \code{vcov.} function.}  
}

\examples{
data(simData)
mod <- setModel(Y~Z|~X1+X2, data=simData)
fit <- tlse(mod, seType="lm", vcov.=sandwich::vcovHC, type="HC")
plot(fit, "X1")
}
