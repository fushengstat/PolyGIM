\name{data}
\alias{data}
\docType{data}
\title{
Data for example in \code{\link{PolyGIM}}
}
\description{
\code{data} is a data frame used in the example of \code{\link{PolyGIM}}.
}
\usage{data("data")}
\format{
  A data frame with 4000 observations on the following 4 variables.
  \describe{
    \item{\code{int}}{a data frame with 21 covairates (X1,..,X21), a PRS score covariate (score) and polytoumous outcomes (y)}
    \item{\code{moedls}}{a list of 147 summary data and models from marginal binary logistic regression models}
    \item{\code{nctrl}}{a matrix specifying the number of cases shared in datasets that are used to fit the working models given in \code{models}}
    \item{\code{nctrl}}{a matrix specifying the number of controls shared in datasets that are used to fit the working models given in \code{models}}
  }
}
\details{
This is a dataset from which internal and external data are extracted for the example.
}

\keyword{datasets}
