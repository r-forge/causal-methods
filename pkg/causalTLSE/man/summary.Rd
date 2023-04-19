\name{summary}
\alias{summary.causaltlse}
\alias{summary.tlseFit}
\title{Summary for objects of class otlse}
\description{
The method builds the coefficient matrices. 
}
\usage{
\method{summary}{causaltlse}(object, ...)
\method{summary}{tlseFit}(object, vcov.=NULL, ...)
}
\arguments{
  \item{object}{An object of class \code{causalTlse}}
  \item{vcov.}{A function to compute the strandard error of the
    coefficients. The default is \code{vcov.lm}}
\item{...}{Argument for other type of objects}
}

\examples{
data(simDat3)
mod1 <- setModel(Y~Z|~X1*X2, data=simDat3)
CE <- causalTLSE(mod1, selType="BTLSE", vcov.=sandwich::vcovHC, type="HC3")
summary(CE)
fit <- estModel(mod1)
summary(fit)
}

