\name{summary}
\alias{summary.causalslse}
\alias{summary.slseFit}
\title{Summary for objects of class oslse}
\description{
The method builds the coefficient matrices. 
}
\usage{
\method{summary}{causalslse}(object, ...)
\method{summary}{slseFit}(object, vcov.=vcovHC, ...)
}
\arguments{
  \item{object}{An object of class \code{causalSlse}}
  \item{vcov.}{A function to compute the strandard error of the
    coefficients. The default is \code{vcovHC} for \code{lm} objects.}
\item{...}{Argument for other type of objects}
}

\examples{
data(simDat3)
mod1 <- setModel(Y~Z|~X1*X2, data=simDat3)
CE <- causalSLSE(mod1, selType="BSLSE")
summary(CE)
fit <- estModel(mod1)
summary(fit)
}

