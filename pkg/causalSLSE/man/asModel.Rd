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

\value{
The method returns an object of class \code{\link{slseModel}} or
\code{\link{cslseModel}}. 
}

\examples{
data(simDat3)
mod <- cslseModel(Y ~ Z | ~ X1 * X2, data = simDat3)
fit <- estSLSE(mod)

## Extract the model from a cslseModel object

as.model(fit)

## Extract the model from a cslse object

cs <- causalSLSE(mod)
as.model(cs)
}
