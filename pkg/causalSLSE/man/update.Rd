\name{update}
\alias{update.cslseModel}
\alias{update.slseModel}
\alias{update.slseKnots}
\title{Update methods}
\description{
Method to modify objects from the package. 
}
\usage{

\method{update}{cslseModel}(object, selType, selCrit="AIC",
                            selKnots, ...)

\method{update}{slseModel}(object, selType, selCrit="AIC",
                           selKnots, ...)

\method{update}{slseKnots}(object, selKnots, ...)

}

\arguments{
  \item{object}{An object for be modified.}
  \item{selKnots}{An optional list of integers to select the knots from
    the list of knots specified by the model. It missing, the current
    knots are kept.} 
  \item{selType}{The selection method: backward for forward.}
  \item{selCrit}{The criterion to select the piecewise polynomial knots.}
  \item{...}{Argument for other type of objects. Currently not used.}
}

\details{
  The method for \code{slseKnots} is explained in the vignette, but it
  is mostly used internally. For the model objects, the method can be used
  to choose a set of knots already stored in the object. It avoids
  having to re-select them which can be computationally intensive for
  large samples. It returns an error message if the type of selection
  requested does not exist. In that case, you need to run the
  \code{selSLSE} method. If the \code{selType} argument is set to
  \code{None}, the method returns the  original model.  
}

\examples{
data(simDat3)
mod1 <- cslseModel(Y~Z|~X1*X2, data=simDat3)
mod2 <- selSLSE(mod1, "BLSE", "AIC")
update(mod2, "BLSE", "BIC")
update(mod2, "None")
}












