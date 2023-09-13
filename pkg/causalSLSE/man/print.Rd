\name{print}
\alias{print.cslseModel}
\alias{print.cslse}
\alias{print.slseFit}
\alias{print.slseKnots}
\alias{print.cslseKnots}
\alias{print.summary.cslse}
\alias{print.summary.slseFit}
\title{Print methods}
\description{
Print methods for different objects from the package.
}
\usage{
\method{print}{cslseModel}(x, which=c("Model", "selKnots", "Pvalues"),
                           digits = max(3L, getOption("digits") - 3L), ...)
\method{print}{cslse}(x, digits = max(3L, getOption("digits") - 3L), ...)
\method{print}{slseFit}(x, digits = max(3L, getOption("digits") - 3L), ...)
\method{print}{summary.cslse}(x, digits = max(3L, getOption("digits") - 3L),
                              signif.stars = getOption("show.signif.stars"),
                              beta=FALSE, knots = FALSE,...)
\method{print}{summary.slseFit}(x, digits = max(3L, getOption("digits") - 3L),
                                signif.stars = getOption("show.signif.stars"),
    ...)
\method{print}{slseKnots}(x, which=c("selKnots", "Pvalues"),
                          header=c("None", "All", "Select"),
                          digits = max(3L, getOption("digits") - 3L), ...)
\method{print}{cslseKnots}(x, which=c("selKnots", "Pvalues"),
                           header = c("None", "All", "Select"),
                           digits = max(3L, getOption("digits") - 3L), ...)
}
\arguments{
  \item{x}{Model of fit objects.}
  \item{digits}{The number of digits to print.}
  \item{signif.stars}{Should we print the significant stars?}
  \item{beta}{Should we print the coefficient matrix of the piecewise
    polynomials?}
  \item{knots}{Should we print the knots of the piecewise polynomial?}
  \item{which}{What info from the model should be printed? It prints a
    summary of the model when set to \code{Model}, the current knots
    when set to \code{selKnots} and the p-values and initial knots when
    set to \code{Pvalues}. For the latter, it is only available when a
    knots selection method has been applied to the object.}
  \item{header}{What description of the object should be printed.}
\item{...}{Argument for other type of objects}
}

\examples{
data(simDat3)
mod1 <- cslseModel(Y~Z|~X1*X2, data=simDat3)
print(mod1)
causal <- causalSLSE(mod1)
print(causal)
s <- summary(causal)
print(s)
fit <- estSLSE(mod1)
print(fit)
summary(fit)
## the following is a cslseKnots object
print(mod1$knots, "selKnots")
## the following is as slseKnots object
print(mod1$knots$treated, "selKnots")
}












