
<!-- README.md is generated from README.Rmd. Please edit that file -->

# coglasso - Collaborative Graphical Lasso

<!-- badges: start -->

[![R-CMD-check](https://github.com/DrQuestion/coglasso/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/DrQuestion/coglasso/actions/workflows/R-CMD-check.yaml)

[![codecov](https://codecov.io/gh/DrQuestion/coglasso/graph/badge.svg?token=Q370RQ1CAD)](https://app.codecov.io/gh/DrQuestion/coglasso)

[![downloads](https://cranlogs.r-pkg.org/badges/coglasso)](https://cran.r-project.org/package=coglasso)
<!-- badges: end -->

Coglasso implements *collaborative graphical lasso*, an algorithm for
network reconstruction from multi-omics data sets ([Albanese, Kohlen and
Behrouzi, 2024](#references)). Our algorithm joins the principles of the
*graphical lasso* by Friedman, Hastie and Tibshirani
([2008](#references)) and *collaborative regression* by Gross and
Tibshirani ([2015](#references)).

## Installing coglasso

You can install the CRAN release of coglasso with:

``` r
install.packages("coglasso")
```

## Installing the development version

To install the development version of coglasso from
[GitHub](https://github.com/) you need to make sure to install devtools
with:

``` r
if (!require("devtools")) {
  install.packages("devtools")
}
```

You can then install the development version with:

``` r
devtools::install_github("DrQuestion/coglasso")
```

## Usage

Here follows an example on how to reconstruct a multi-omics network with
*collaborative graphical lasso*. For a more exhaustive example we refer
the user to the vignette `vignette("coglasso")`. The package provides
example multi-omics data sets of different dimensions, here we will use
`multi_omics_sd_small`. Please notice that the current version of the
coglasso package expects multi-omics data sets with *two* “omic” layers,
where the single layers are grouped by column. For example, in
`multi_omics_sd_small` the first 14 columns represent transcript
abundances, and the other 5 columns represent metabolite abundances. To
default usage of `coglasso()` only needs the input dataset and the
dimension of the first “omic” layer.

``` r
library(coglasso)

cg <- coglasso(multi_omics_sd_small, pX = 14)
```

`coglasso()` explores several combinations of the hyperparameters
characterizing *collaborative graphical lasso*. To select the
combination yielding the best network according to the chosen model
selection method, the package provides the function `select_coglasso()`.
Among others, this function implements *eXtended Efficient StARS*
(*XEStARS*), a significantly faster and memory-efficient version of
*eXtended StARS* (*XStARS*, [Albanese, Kohlen and Behrouzi,
2024](#ref)). These are coglasso-adapted versions of the *StARS*
selection algorithm ([Liu, Roeder and Wasserman, 2010](#references))
selecting the hyperparameter combination that yields the most stable,
yet sparse network. *XEStARS* is the default option for the parameter
`method`, so it is enough to call:

``` r
sel_cg <- select_coglasso(cg)
```

## References

Albanese, A., Kohlen, W., & Behrouzi, P. (2024). Collaborative graphical
lasso (arXiv:2403.18602). *arXiv*
<https://doi.org/10.48550/arXiv.2403.18602>

Friedman, J., Hastie, T., & Tibshirani, R. (2008). Sparse inverse
covariance estimation with the graphical lasso. *Biostatistics*, 9(3),
432–441. <https://doi.org/10.1093/biostatistics/kxm045>

Gross, S. M., & Tibshirani, R. (2015). Collaborative regression.
*Biostatistics*, 16(2), 326–338.
<https://doi.org/10.1093/biostatistics/kxu047>

Liu, H., Roeder, K., & Wasserman, L. (2010). Stability Approach to
Regularization Selection (StARS) for High Dimensional Graphical Models
(arXiv:1006.3316). *arXiv* <https://doi.org/10.48550/arXiv.1006.3316>
