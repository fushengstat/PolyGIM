\name{polygim_v}
\alias{polygim_v}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  PolyGIM
}
\description{Integrating Summary Statistics and Individual-level Data with Fixed V}
\usage{
polygim_v(formula, data, models, ncase, nctrl, V)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A formula.}
  \item{data}{A internal individual-level data frame containing all variables that are specified in \code{formula} and \code{models}. }
  \item{models}{A list of external model and summary statistics from multiple external studies.}
  \item{ncase}{A matrix specifying the number of cases shared in datasets that are used to fit the working models given in \code{models}.}
  \item{nctrl}{A matrix specifying the number of controls shared in datasets that are used to fit the working models given in \code{models}.}
  \item{V}{The variance-covariance matrix of external estimate summary data in \code{models}.}
}

\examples{
  library("PolyGIM")
  data(data, package = "PolyGIM")
  formula = "y~score"
  V = diag(length(models))
  fit = polygim_v(formula, int, models, ncase, nctrl, V)
  # estimate
  fit$theta
  # the corresponding standard errors
  fit$se
}
