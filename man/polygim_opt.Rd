\name{polygim_opt}
\alias{polygim_opt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  PolyGIM
}
\description{Integrating Summary Statistics and Individual-level Data Using the Optimal V}
\usage{
 polygim_opt(formula, data, models, ncase, nctrl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A formula.}
  \item{data}{A internal individual-level data frame containing all variables that are specified in \code{formula} and \code{models}. }
  \item{models}{A list of external model and summary statistics from multiple external studies.}
  \item{ncase}{A matrix specifying the number of cases shared in datasets that are used to fit the working models given in \code{models}.}
  \item{nctrl}{A matrix specifying the number of controls shared in datasets that are used to fit the working models given in \code{models}.}
}

\examples{
  data(data, package="PolyGIM")
  formula = "y~score"
  fit = polygim_opt(formula, int, models, ncase, nctrl)
  # estimate
  fit$theta
  # the corresponding standard errors
  fit$se
}
