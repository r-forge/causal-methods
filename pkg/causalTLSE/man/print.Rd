\name{print}
\alias{print.tlseModel}
\alias{print.tlse}
\alias{print.summary.tlse}
\title{Print methods}
\description{
Print results of different methods applied to \code{tlse} ands
\code{tlseModel} Objects.
}
\usage{
\method{print}{tlseModel}(x, knots=FALSE, ...)
\method{print}{tlse}(x, ...)
\method{print}{summary.tlse}(x, digits = 4,
                              signif.stars = getOption("show.signif.stars"),
                              beta=FALSE, knots = FALSE,...)
}
\arguments{
  \item{x}{Onject of class \code{tlse} or \code{summary.tlse}.}
  \item{digits}{The number of digits to print.}
  \item{signif.stars}{Should we print the significant stars?}
  \item{beta}{Should we print the coefficient matrix of the piecewise
    polynomials?}
  \item{knots}{Should we print the knots of the piecewise polynomial?}
\item{...}{Argument for other type of objects}
}

\examples{
data(simData)
mod1 <- setModel(Y~Z|~X1*X2, data=simData)
print(mod1)
fit <- tlse(mod1)
print(fit)
s <- summary(fit)
print(s)
}












