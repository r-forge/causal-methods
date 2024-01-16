\name{update}
\alias{update.cslseModel}
\alias{update.slseModel}
\alias{update.slseKnots}
\title{Update Methods}
\description{
The method updates an object by modifying its specification. Currently,
it is used to change the set of knots by either selecting them manually
or by specifying the selection methods.
}
\usage{

\method{update}{cslseModel}(object, selType, selCrit="AIC",
                            selKnots, ...)

\method{update}{slseModel}(object, selType, selCrit="AIC",
                           selKnots, ...)

\method{update}{slseKnots}(object, selKnots, ...)

}

\arguments{
  \item{object}{An object to be modified.}
  \item{selKnots}{An optional list of integers to select the knots from
    the original list of knots. If missing, the current knots are kept.}  
  \item{selType}{The selection method: \code{"BLSE"} for the backward method, 
  \code{"FLSE"} for the forward method or \code{"none"} to recover the
  originnal set of knots.}  
  \item{selCrit}{The criterion to select the optimal set of knots.}
  \item{...}{Argument for other types of objects. Currently not used.}
}

\details{
  The method for \code{slseKnots} is explained in the vignette, but it
  is mostly used internally. For the model objects, the method can be used
  to choose a set of knots already stored in the object. It avoids
  having to re-compute them which can be computationally intensive for
  large samples. It returns an error message if the type of selection
  requested does not exist. In that case, you need to run the
  \code{selSLSE} method. If the \code{selType} argument is set to
  \code{None}, the method returns the  original model.  
}

\value{
  It returns an object of the same class, but with a different set of knots.
}

\seealso{\code{\link{selSLSE}} for more details on how stored knots can
  be selected from a model, \code{\link{slseModel}} for model
  description and \code{\link{slseKnots}} for the format of knots}

\examples{
data(simDat3)
mod1 <- cslseModel(Y ~ Z | ~ X1 * X2, data = simDat3)
mod2 <- selSLSE(mod1, "BLSE", "AIC")

## We changed the knots to the BLSE-BIC selection
## already stored in the model object

update(mod2, "BLSE", "BIC")

## We recover the original set of knots

update(mod2, "None")
}












