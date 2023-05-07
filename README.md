# PolyGIM: A procedure to integrate individual-level data and summary data for polytomous outcome

# Overview
**PolyGIM** is a novel method that integrates individual-level data and summary statistics with polytomous outcome

# Installation
**PolyGIM** PACKAGE can be installed via Github. To install the latest version of PolyGIM package via Github, run following commands in R:
```{r, include = FALSE}
library(devtools)
devtools::install_github("fushengstat/PolyGIM")
```

# Usage
This is a good start point for using **PolyGIM** package.
```{r,include = FALSE}
data(data, package = "PolyGIM")
formula = "y~score"
V = diag(length(models))
fit = polygim_v(formula, int, models, ncase, nctrl, V)

# estimate
fit$theta
# the corresponidng standard errors
fit$se
```

Please see the [PolyGIM user manual](https://github.com/fushengstat/PolyGIM/blob/main/doc/PolyGIM.pdf) for more details. 



# References
Fu S.,Purdue M. P., Zhang H., Qin J., Song L., Berndt S. I., & Yu K. Improve the model of disease subtype heterogeneity by leveraging external summary data. Manuscript, 2023. \
Zhang, H., Deng, L., Schiffman, M., Qin, J., & Yu, K. (2020). Generalized integration model for improved statistical inference by leveraging external summary data. Biometrika, 107(3), 689-703. 
