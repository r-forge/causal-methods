\name{simDat2}
\alias{simDat2}
\docType{data}
\title{
Simulated data.
}
\description{
This dataset is used in several documentation files to illustrate the
different functionalities of the package.
}
\usage{data("simDat2")}
\format{
  A data frame with 300 observations on the following 11 variables. 
  \describe{
    \item{\code{X}}{Continuous Covariate}
    \item{\code{Z}}{Treatment indicator.}
    \item{\code{Y}}{Observed outcome.}
    \item{\code{Y1}}{Outcome for the treated group (not observed in
      practice).}      
    \item{\code{Y0}}{Outcome for the control group (not observed in
      practice).}
    \item{\code{U01},\code{U02},\code{U03}}{The true bases of the covariates
      for the control group.}
    \item{\code{U11},\code{U12},\code{U13}}{The true bases of the covariates
      for the treated group.}    
  }
}

\keyword{datasets}
