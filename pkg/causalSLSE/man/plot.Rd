\name{plot}
\alias{plot.cslseFit}
\alias{plot.slseFit}
\title{Plot of Predicted Outcome}
\description{
For objects of class \code{slseFit}, the method plots the predicted
outcome with respect to a given covariate. It is the same for objects of class
\code{cslseFit}, but the predicted outcome is plotted for each treatment
group separately.   
}
\usage{
\method{plot}{cslseFit}(x, y, which = y, interval = c("none", "confidence"), 
                        level = 0.95, fixedCov = list(),
                        vcov. = vcovHC,  add = FALSE, addToLegend = NULL,
                        addPoints = FALSE, FUN = mean,
                        plot=TRUE, graphPar=list(), ...)
\method{plot}{slseFit}(x, y, which = y, interval = c("none", "confidence"), 
                       level = 0.95, fixedCov = NULL,
                       vcov. = vcovHC,  add = FALSE, 
                       addPoints = FALSE, FUN = mean,
                       plot=TRUE, graphPar=list(), ...)
}
\arguments{
  \item{x}{Object of class \code{cslseFit} or \code{slseFit} created by
    \code{\link{estSLSE}}.}
  \item{y}{alias for \code{which} for compatibility with \code{plot}.}
  \item{which}{Which covariate to plot against the outcome variable. It
    could be an integer or a character.}
  \item{interval}{The type of confidence interval. The default is
    \code{none}.}
  \item{level}{The confidence interval level if included.}
  \item{fixedCov}{List of covariates to fix to
    specific values for the nontreated and treated groups. By default,
    covariates not selected by \code{which} are set to their group
    specific sample means.}
  \item{vcov.}{An alternative function to compute the covariance
    matrix of the least squares estimates. The default is the
    \code{vcovHC} method for \code{lm} objects.}
  \item{add}{Should the plot be added to an existing one?}
  \item{addToLegend}{A character string to add to the legend next to
    treated and control.}
  \item{addPoints}{Should we add the scatterplot of the outcome and
    covariate on the graph?}
  \item{FUN}{The function to determine the fixed value for the
    covariates not fixed by \code{fixedCov}. The
    default is \code{mean}.} 
  \item{plot}{If set to \code{FALSE}, a \code{data.frame} is returned
    for each group with the covariate selected by \code{which} and the
    prediction.}
  \item{graphPar}{A list of graphical parameters. See Details.}
  \item{...}{Additional argument to pass to the \code{vcov.} function.}  
}

\details{

The default set of parameters can be obtained by running the command
\code{causalSLSE:::.initPar()}. It returns a list of four elements:
\code{treated} for the parameter of the lines or points of the treated,
\code{nontreated} for the parameters of the nontreated, \code{common}
for the common parameters not specific to a group like the main title or
the axis labels, and \code{legend} for the legend parameters. The
elements \code{treated} and \code{nontreated} are lists of two elements:
\code{points} and \code{lines}, which are lists of graphical parameters
for the scatterplot (when \code{addPoints} is \code{TRUE}) and the
lines. The argument \code{graphPar} can be use to modify existing
parameters or to add new ones. It must respect the structure of the
initial list. See the example below.

To fix covariates to the same values for both groups, \code{fixedCov} is
just a named list with values for the covariates associated with the
names. To fix the covariates to different values for the treated and
nontreated, \code{fixedCov} is a named list of at least 1 element with
names being either treated, nontreated or both. If only one group is
specified, the covariates for the other groups are determined by
\code{FUN}.

}

\value{

It returns \code{NULL} invisibly if the argument \code{plot} is
\code{TRUE}. Otherwise, it returns a \code{data.frame} with the
necessary variables to manually create the plot. For \code{slseModel}
objects, it is a list with the following elements:

\item{Outcome}{The outcome variable from the model. The name of this
  element is the actual name of the outcome variable in the model.}

\item{which}{The covariate selected by the argument \code{which}. The
  name of this element is the name of the selected covariate.}

\item{fit}{The predicted outcome from the model fit.}

\item{lower, upper}{The lower and upper bounds of the confidence
  interval. It is only available when the argument \code{interval} is
  set to \code{"confidence"}.

  }

Note that all returned variables are ordered with respect to the selected
covariates. See the example below.

For \code{cslseModel}, the above list of variables is returned
separately for each treatment group. 
}

\examples{
data(simDat3)

##  For cslse objects

mod <- cslseModel(Y ~ Z | ~ X1 + X2, data = simDat3)
fit <- causalSLSE(mod) 
plot(fit, "X1")

## Let's change the type of points for the treated and lines for the
## nontreated

gpar <- list(treated = list(points = list(pch = 24, col = 5)),
             nontreated = list(lines = list(lty = 5, col =
             "darkgreen")), common = list(xlab = "New X",  main =
             "Plot with changed parameters"))  

plot(fit, "X1", addPoints = TRUE, graphPar = gpar)

## For slseModel objects:

mod2 <- slseModel(Y ~ X1 + X2, data = simDat3)
fit2 <- estSLSE(mod2)
plot(fit2, "X1", interval = "confidence", addPoints = TRUE)

##  The same graph produced manually

p2 <- plot(fit2, "X1", interval = "confidence", plot = FALSE)
plot(p2$X1, p2$Y, pch = 21, main = "Y against X1", xlab = "X1", ylab =
     "Y") 
lines(p2$X1, p2$fit, lwd = 2)
lines(p2$X1, p2$lower, lty = 3, lwd = 2)
lines(p2$X1, p2$upper, lty = 3, lwd = 2)
}
