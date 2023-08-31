\name{plot}
\alias{plot.slseFit}
\title{Plot methods}
\description{
Some plot utilities for object of class \code{slseFit}.
}
\usage{
\method{plot}{slseFit}(x, y, which = y, interval = c("none", "confidence"), 
                   counterfactual = FALSE, level = 0.95, fixedCov0 = NULL,
                   fixedCov1 = fixedCov0, legendPos = "topright", 
                   vcov. = vcovHC, col0 = 1, col1 = 2, lty0 = 1, lty1 = 2, add. = FALSE, 
                   addToLegend = NULL, cex = 1, ylim. = NULL, xlim. = NULL, 
                   addPoints = FALSE, FUN = mean, main = NULL, plot=TRUE, ...)
}
\arguments{
  \item{x}{Object of class \code{slseFit}.}
  \item{y}{alias for \code{which} for compatibility with \code{plot}}
  \item{which}{Which covariate to plot against the outcome variable. It
    could be an integer or a character}
  \item{interval}{The type of confidence interval. The default is
    \code{none}.}
  \item{counterfactual}{Should we plot the counterfactual prediction?}
  \item{level}{The confidence interval level if included.}
  \item{fixedCov0, fixedCov1}{Named lists of covariates to fix to
    specific values for the nontreated and treated groups. By default,
    covariates not selected by \code{which} are set to their group
    specific sample  means. If only \code{fixedCov0} is set,
    \code{fixedCov1} is set to the same list.}
  \item{legendPos}{The legend position.}
  \item{vcov.}{An alternative function to compute the covariance
    matrix of the least squares estimates. The default is the
    \code{vcovHC} method for \code{lm} objects.}
  \item{col0, col1}{The line color for the control group and teated
    groups.}
  \item{lty0, lty1}{The type of line for the control group and teated
    groups.}  
  \item{add.}{Should the plot be added to an existing one?}
  \item{addToLegend}{A character string to add to the legend next to
    treated and control.}
  \item{cex}{The font size for the legend.}
  \item{ylim., xlim.}{An optional range for the y-axis and x-axis.}
  \item{addPoints}{Should we add the scatterplot of the outcome and
    covariate on the graph?}
  \item{FUN}{The function to determine the fixed value for the
    covariates not fixed by \code{fixedCov0} and \code{fixedCov1}. The
    default is \code{mean}.} 
  \item{main}{An optional title to replace the default one.}
  \item{plot}{If set to \code{FALSE}, a \code{data.frame} is returned
    for each group with the covariate selected by \code{which} and the
    prediction.} 
  \item{...}{Additional argument to pass to the \code{vcov.} function.}  
}

\examples{
data(simDat3)
mod <- cslseModel(Y~Z|~X1+X2, data=simDat3)
fit <- causalSLSE(mod, seType="lm") 
plot(fit, "X1")
}
