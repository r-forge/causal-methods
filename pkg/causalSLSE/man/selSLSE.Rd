\name{selSLSE}
\alias{selSLSE}
\alias{selSLSE.cslseModel}
\alias{selSLSE.slseModel}
\title{
Knots Selection Method
}
\description{
This is the main function to select the knots in \code{cslseModel} or
\code{slseModel} objects. It returns a model with an optimal set of
knots. 
}
\usage{
\method{selSLSE}{cslseModel}(model, selType=c("BLSE", "FLSE"),
        selCrit = c("AIC", "BIC", "PVT"), 
        pvalT = function(p) 1/log(p),
        vcovType = c("HC0", "Classical", "HC1", "HC2", "HC3"),
        reSelect=FALSE, ...)
\method{selSLSE}{slseModel}(model, selType=c("BLSE", "FLSE"),
        selCrit = c("AIC", "BIC", "PVT"), 
        pvalT = function(p) 1/log(p),
        vcovType = c("HC0", "Classical", "HC1", "HC2", "HC3"),
        reSelect=FALSE, ...)
}
\arguments{
  \item{model}{A model of class \code{\link{cslseModel}} or
    \code{\link{slseModel}}.} 
  \item{selType}{The selection method: backward (\code{"BLSE"}) or
    forward (\code{"FLSE"}).} 
  \item{selCrit}{The criterion to select the piecewise polynomial knots.}
  \item{pvalT}{A function to determine the p-value threshold for the
    significance of the coefficients. It has to be a function of one
    parameter, which is the average number of knots in the model.}
  \item{vcovType}{The type of least squares covariance matrix used to
    compute the p-values needed for the selection.}
  \item{reSelect}{By default, the stored selections are used. If
    \code{reSelect} is set to \code{TRUE}, the selection is re-computed.}
  \item{...}{Additional arguments to pass to other methods. Currently
    not used.}
}

\details{

It selects the knots using one of the two methods, FLSE or BLSE, with
either the AIC, BIC or a p-value threshold (see the vignette for more
details). Any of these selection methods requires several least squares
estimations and it is performed only if the method has not been applied
yet and \code{reSelect} is set to \code{TRUE}. This is possible because
any new knots selection is saved into the returned model. A model may
have more than one selection saved into it. The active knots (the ones
used when we estimate the model) is stored into the element \code{knots}
of the model and the saved selections are stored into the element
\code{selections}. See below for what is included in this element.  

Note that the selections for the three criteria AIC, BIC and PVT are
computed and saved automatically in the returned model when
\code{selCrit} is set to either \code{"AIC"} or \code{"BIC"}, because it
does not require many more operations to select them all once we do it
for AIC or BIC. However, it is only computed for PVT when \code{selCrit}
is set to \code{"PVT"}. 

The knots are selected jointly for all treatment groups in
\code{cslseModel} objects. However, the active knots and all saved
selections are stored separately for each treatment group. For example,
the active knots for the treated in the \code{cslseModel} object
\code{mod} are stored in \code{mod$treated$knots}. See the Vignette for
more details.

}

\value{

  The method returns an object of class \code{\link{slseModel}} or
\code{\link{cslseModel}} depending on which object it is applied
to. When it does not already exist, the element \code{selections} is added
to the \code{slseModel} object (or to each \code{slseModel} object in
\code{cslseModel} objects). The element \code{selections} is a list with
one or more of the following elements:

\item{originalKnots}{The original knots as selected initially by
  \code{\link{slseModel}} or \code{\link{cslseModel}}. This element is
  always included.}

\item{FLSE, BLSE}{This is where selections based on the forward (FLSE) and
  the backward (BLSE) methods are stored.}

Finally, \code{BLSE} and \code{FLSE} are lists that may contain the following elements:   

\item{AIC, BIC}{A list of integer vectors, one for each covariate in the
  \code{slseModel} model. The vectors of integers indicate which sets the
  original set of knots minimizes the AIC or the BIC.}

\item{PVT}{Same as the AIC and BIC, but the selection is based on a p-value
  threshold.}

\item{JAIC, JBIC}{This criteria is for \code{cslseModel} models. It is
  like the AIC and BIC, but it indicates that the AIC or BIC was
  computes by estimating the models of all treatment groups
  jointly. This selection is not the same as a selection group by group.}

\item{Threshold}{The p-value threshold used for the PVT criterion.}

\item{pval}{A list of p-values, one for each original knots. See
  vignette for a definition of the p-values.}

}

\seealso{\code{\link{slseModel}} and \code{\link{cslseModel}} for model
  objects description and \code{\link[=update.slseModel]{update}} for
  ways of selecting stored selections}

\examples{
data(simDat3)
mod1 <- cslseModel(Y~Z|~X1*X2, data=simDat3)
mod1 <- selSLSE(mod1, selType="FLSE", selCrit="AIC")

## The following does not require additional computation
## because the selection is stored in mod1

mod1 <- selSLSE(mod1, selType="FLSE", selCrit="BIC")


## But the following does

mod1 <- selSLSE(mod1, selType="BLSE", selCrit="BIC")

## See one selection:

mod1$treated$selections$BLSE$JBIC

}

\keyword{causal effects}
