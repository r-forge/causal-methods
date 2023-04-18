\name{plot}
\alias{plot.tlseFit}
\title{Plot methods}
\description{
Some plot utilities for object of class \code{tlseFit}.
}
\usage{
\method{plot}{tlseFit}(x, y, which = y, interval = c("none", "confidence"), 
                    level = 0.95, newdata = NULL, legendPos = "topright",
                    vcov. = NULL, col0=2, col1=5, lty0=1, lty1=2,
                    add.=FALSE, addToLegend=NULL,
                    cex=1, ylim.=NULL, addPoints=FALSE, ...) 
}
\arguments{
  \item{x}{Object of class \code{tlseFit}.}
  \item{y}{alias for \code{which} for compatibility with \code{plot}}
  \item{which}{Which covariate to plot against the outcome variable. It
    could be an integer or a character}
  \item{interval}{The type of confidence interval. The default is \code{none}.}
  \item{level}{The confidence interval level if included.}
  \item{newdata}{A named vector to fix some covariates. By default,
    covariates not selected by \code{which} are set to their sample means.}
  \item{legendPos}{The legend position.}
  \item{vcov.}{An alternative function to compute the covariance
    matrix of the least squares estimates. The default \code{vcov}
    method for \code{lm} objects.}
  \item{col0, col1}{The line color for the control group and teated
    groups.}
  \item{lty0, lty1}{The type of line for the control group and teated
    groups.}  
  \item{add.}{Should the plot be added to an existing one?}
  \item{addToLegend}{A character string to add to the legend next to
    treated and control.}
  \item{cex}{The font size for the legend.}
  \item{ylim.}{An optional range for the y-axis.}
  \item{addPoints}{Should we add the scatterplot of the outcome and
    covariate on the graph?}
  \item{...}{Additional argument to pass to the \code{vcov.} function.}  
}

\examples{
data(simDat3)
mod <- setModel(Y~Z|~X1+X2, data=simDat3)
fit <- causalTlse(mod, seType="lm", vcov.=sandwich::vcovHC, type="HC")
plot(fit, "X1")
}
