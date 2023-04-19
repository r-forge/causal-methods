\name{print}
\alias{print.tlseModel}
\alias{print.causaltlse}
\alias{print.tlseFit}
\alias{print.summary.causaltlse}
\alias{print.summary.tlseFit}
\title{Print methods}
\description{
Print methods for different objects from the package.
}
\usage{
\method{print}{tlseModel}(x, knots=FALSE, ...)
\method{print}{causaltlse}(x, ...)
\method{print}{tlseFit}(x, ...)
\method{print}{summary.causaltlse}(x, digits = 4,
                              signif.stars = getOption("show.signif.stars"),
                              beta=FALSE, knots = FALSE,...)
\method{print}{summary.tlseFit}(x, digits = 4,
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
\item{...}{Argument for other type of objects}
}

\examples{
data(simDat3)
mod1 <- setModel(Y~Z|~X1*X2, data=simDat3)
print(mod1)
causal <- causalTLSE(mod1)
print(causal)
s <- summary(causal)
print(s)
fit <- estModel(mod1)
print(fit)
summary(fit)
}












