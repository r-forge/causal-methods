\name{print}
\alias{print.cslseModel}
\alias{print.cslse}
\alias{print.slseFit}
\alias{print.summary.cslse}
\alias{print.summary.slseFit}
\title{Print methods}
\description{
Print methods for different objects from the package.
}
\usage{
\method{print}{cslseModel}(x, which=c("Model", "selKnots", "Pvalues"), 
                            ...)
\method{print}{cslse}(x, ...)
\method{print}{slseFit}(x, ...)
\method{print}{summary.cslse}(x, digits = 4,
                              signif.stars = getOption("show.signif.stars"),
                              beta=FALSE, knots = FALSE,...)
\method{print}{summary.slseFit}(x, digits = 4,
                                signif.stars = getOption("show.signif.stars"),
                                ...)
}
\arguments{
  \item{x}{Model of fit objects.}
  \item{digits}{The number of digits to print.}
  \item{signif.stars}{Should we print the significant stars?}
  \item{beta}{Should we print the coefficient matrix of the piecewise
    polynomials?}
  \item{knots}{Should we print the knots of the piecewise polynomial?}
  \item{which}{What info form the model should be printed?}
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
}












