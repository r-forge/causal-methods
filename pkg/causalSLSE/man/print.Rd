\name{print}
\alias{print.cslseModel}
\alias{print.slseModel}
\alias{print.cslse}
\alias{print.cslseFit}
\alias{print.slseFit}
\alias{print.slseKnots}
\alias{print.summary.cslse}
\alias{print.summary.cslseFit}
\alias{print.summary.slseFit}
\title{Print Methods}
\description{
Print methods for different objects from the package.
}
\usage{
\method{print}{cslseModel}(x, which=c("Model", "selKnots", "Pvalues"),
                           digits = max(3L, getOption("digits") - 3L), ...)
\method{print}{slseModel}(x, which=c("Model", "selKnots", "Pvalues"),
                           digits = max(3L, getOption("digits") - 3L), ...)
\method{print}{cslse}(x, digits = max(3L, getOption("digits") - 3L), ...)
\method{print}{cslseFit}(x, digits = max(3L, getOption("digits") - 3L), ...)
\method{print}{slseFit}(x, digits = max(3L, getOption("digits") - 3L), ...)
\method{print}{summary.cslse}(x, digits = max(3L, getOption("digits") - 3L),
                              signif.stars = getOption("show.signif.stars"),
                              beta=FALSE, knots = FALSE,...)
\method{print}{summary.cslseFit}(x, groups, digits = max(3L, getOption("digits") - 3L),
                                signif.stars = getOption("show.signif.stars"),
                                ...)
\method{print}{summary.slseFit}(x, digits = max(3L, getOption("digits") - 3L),
                                signif.stars = getOption("show.signif.stars"),
                                ...)
\method{print}{slseKnots}(x, header=c("None", "All", "Select"),
                          digits = max(3L, getOption("digits") - 3L), ...)
}
\arguments{
  \item{x}{A model or a fit object.}
  \item{digits}{The number of digits to print.}
  \item{signif.stars}{Should we print the significant stars?}
  \item{beta}{Should we print the coefficient matrix of the basis
    functions?} 
  \item{knots}{Should we print the set of knots?}
  \item{which}{What info we should print? It prints a summary of the
    model when set to \code{"Model"}, the current knots when set to
    \code{"selKnots"} and the p-values and initial knots when set to
    \code{"Pvalues"}. For the latter, it is only available when a knots
    selection method has been applied to the object.} 
  \item{header}{What description of the object should be printed.}
  \item{groups}{The names of the group to be printed. By default,
    they are all printed.}
\item{...}{Argument for other types of objects.}
}

\value{
All methods return \code{NULL} invisibly.
}

\examples{

## For cslseModel objects

data(simDat3)
mod1 <- cslseModel(Y ~ Z | ~ X1 * X2, data = simDat3)
print(mod1)

## For cslse objects

causal <- causalSLSE(mod1)
print(causal)

## For summary.cslse objects

s <- summary(causal)
print(s)

## For cslseFit objects

fit <- estSLSE(mod1)
print(fit)

## For summary.cslseFit objects

summary(fit)
}












