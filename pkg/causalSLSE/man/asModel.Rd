\name{as.model}
\alias{as.model}
\alias{as.model.slseFit}
\alias{as.model.cslseFit}
\alias{as.model.cslse}
\title{
Extract the model from objects.
}
\description{
When the information about a model is available, it reconstructs it and returns a valid model object. 
}

\usage{
\method{as.model}{slseFit}(x, ...)

\method{as.model}{cslseFit}(x, ...)

\method{as.model}{cslse}(x, ...)

}

\arguments{
\item{x}{An object containing the model to extract.}
\item{...}{Other arguments to pass to other methods. Currently not used.}
}


\examples{
data(simDat3)
mod <- cslseModel(Y~Z|~X1*X2, data=simDat3)
fit <- estSLSE(mod)
as.model(fit)
cs <- causalSLSE(mod)
as.model(cs)
}
