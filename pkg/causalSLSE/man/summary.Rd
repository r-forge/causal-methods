\name{summary}
\alias{summary.cslse}
\alias{summary.cslseFit}
\alias{summary.slseFit}
\title{Summary for objects of class oslse}
\description{
The method builds the coefficient matrices. 
}
\usage{
\method{summary}{cslse}(object, ...)
\method{summary}{cslseFit}(object, vcov.=vcovHC, ...)
\method{summary}{slseFit}(object, vcov.=vcovHC, ...)
}
\arguments{
  \item{object}{Any fitted models from the package.}
  \item{vcov.}{A function to compute the strandard error of the
    coefficients. The default is \code{vcovHC} for \code{lm} objects.}
\item{...}{Argument for other type of objects}
}

\examples{
data(simDat3)
mod1 <- cslseModel(Y~Z|~X1*X2, data=simDat3)
CE <- causalSLSE(mod1, selType="BLSE")
summary(CE)
fit <- estSLSE(mod1)
summary(fit)
mod2 <- slseModel(Y~X1*X2, data=simDat3)
fit2 <- estSLSE(mod2)
summary(fit2)
}

