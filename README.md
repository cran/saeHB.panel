
<!-- README.md is generated from README.Rmd. Please edit that file -->

# saeHB.panel

We designed this package to provide several functions for area level of
small area estimation using hierarchical Bayesian (HB) method. This
package provides model using Rao-Yu Model for variable interest.This
package also provides a dataset produced by a data generation. The
“rjags” package is employed to obtain parameter estimates. Model-based
estimators involves the HB estimators which include the mean and the
variation of mean. For the reference, see Rao and Molina (2015) and
Torabi and Shokoohi (2012).

## Author

Velia Tri Marliana, Azka Ubaidillah

## Maintaner

Velia Tri Marliana <221810642@stis.ac.id>

## Function

-   `RaoYuAr1()` This function gives estimation of y using Hierarchical
    Bayesian under Rao Yu Model
-   `Panel()` This function gives estimation of y using Hierarchical
    Bayesian under Rao Yu Model when rho = 0

## Installation

You can install the development version of saeHB.panel from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Veliatrimarliana/saeHB.panel")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(saeHB.panel)
data(dataAr1)
formula = ydi ~ xdi1 + xdi2
area = max(dataAr1[, "area"])
period = max(dataAr1[,"period"])
vardir = dataAr1[,4]
result <- Raoyu.Ar1(formula, area, period,  vardir = vardir, data = dataAr1)
```

Extract area mean estimation

``` r
result$Est
```

Extract coefficient estimation

``` r
result$coefficient
```

Extract area random effect variance

``` r
result$refVar
```

##References \* Rao, J.N.K & Molina. (2015). Small Area Estimation 2nd
Edition. New York: John Wiley and Sons, Inc. \* Torabi, M., & Shokoohi,
F. (2012). Likelihood inference in small area estimation by combining
time-series and cross-sectional data. Journal of Multivariate Analysis,
111, 213–221. <https://doi.org/10.1016/j.jmva.2012.05.016>
