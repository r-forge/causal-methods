\name{summary}
\alias{summary.tlse}
\title{Summary for objects of class otlse}
\description{
The method builds the coefficient matrices. 
}
\usage{
\method{summary}{tlse}(object, ...)
}
\arguments{
 \item{object}{An object of class \code{otlse}}
\item{...}{Argument for other type of objects}
}

\examples{
data(simData)
mod1 <- setModel(Y~Z|~X1*X2, data=simData)
fit <- tlse(mod1, selType="BTLSE", vcov.=sandwich::vcovHC, type="HC3")
summary(fit)
}

