# PolyGIM: A procedure to integrate individual-level data and summary data for polytomous outcomes

# Overview
**PolyGIM** is a novel method that integrates individual-level data and summary statistics with polytomous outcomes

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
# the corresponding standard errors
fit$se
```
**PolyGIM** employs a user interface similar to that of GIM. For more information on using the package, please consult the [user manual](https://github.com/fushengstat/PolyGIM/blob/main/doc/PolyGIM.pdf) and the [GIM vignette](https://cran.rstudio.com/web/packages/gim/vignettes/gim.html).



# References
Fu S., Purdue M. P., Zhang H., Wheeler, W., Qin J., Song L., Berndt S. I., & Yu K. (2023). Improve the model of disease subtype heterogeneity by leveraging external summary data. PLOS Computational Biology, Accepted. \
Fu, S., Deng, L., Zhang, H., Qin, J., & Yu, K. (2023). Integrative analysis of individual-level data and high-dimensional summary statistics. Bioinformatics, 39(4). \
Zhang, H., Deng, L., Wheeler, W., Qin, J., & Yu, K. (2022). Integrative analysis of multiple case‐control studies. Biometrics, 78(3), 1080-1091. \
Zhang, H., Deng, L., Schiffman, M., Qin, J., & Yu, K. (2020). Generalized integration model for improved statistical inference by leveraging external summary data. Biometrika, 107(3), 689-703. 
