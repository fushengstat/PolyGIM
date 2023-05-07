\name{PolyGIM-package}
\alias{PolyGIM-package}
\docType{package}
\title{
Improve the model of disease subtype heterogeneity by leveraging external summary data
}
\description{
  This package proposes a generalized integration procedure that utilizes individual-level data and summary data from multiple binary logistic
  regression models obtained from external sources, in order to fit the polytomous logistic regression (PLR) model for more efficient statistical
  inference.
}
\details{
  Researchers are often interested in understanding the disease subtype heterogeneity by
  testing whether a risk exposure has the same level of effect on different disease subtypes. The
  polytomous logistic regression (PLR) model provides a flexible tool for such an evaluation.
  Disease subtype heterogeneity can also be investigated with a case-only study that uses a
  case-case comparison procedure to directly assess the difference between risk effects on two disease subtypes.
  We develop PolyGIM, a procedure to fit the PLR model by integrating individual-level data with summary
  data extracted from multiple studies under different designs.
  The summary data consist of coefficient estimates from working logistic regression models established by external studies.
}
\author{
Sheng Fu, Kai Yu
}
\references{
 Fu, S., Purdue M. P., Zhang H., Qin J., Song L., Berndt S. I., and Yu K. (2023)
 Improve the model of disease subtype heterogeneity by leveraging external summary data.
 Submitted

  Zhang, H., Deng, L., Schiffman, M., Qin, J., Yu, K. (2020)
  Generalized integration model for improved statistical inference by leveraging external summary data.\cr
  Biometrika. asaa014, https://doi.org/10.1093/biomet/asaa014
}

\keyword{Case-control studies; Disease subtype heterogeneity; Empirical likelihood; Genome-wide association studies;
Polygenic risk score; Polytomous logistic regression model}
