# PROregTMB

Reimplementation of [PROreg](https://cran.r-project.org/package=PROreg) beta-binomial regression models with [TMB](https://github.com/kaskr/adcomp).

## Status

| Feature | Status |
|---------|--------|
| `dBB` / `rBB` | done |
| `BBest` (MM + MLE) | done |
| `BBreg` (marginal BB logistic) | done (TMB) |
| `BBmm` (mixed effects + Laplace) | done (TMB) |
| Multidimensional `BBmm` (`nDim`) | done (basic) |

## Install

```r
# install.packages(c("TMB", "remotes"))
remotes::install_github("idaejin/PROregTMB")
```

For local development without installing:

```r
devtools::load_all("/path/to/PROregTMB")
```

First call to `BBreg()` / `BBmm()` compiles the TMB template (needs a C++ toolchain).

## Quick example

```r
library(PROregTMB)
set.seed(18)
k <- 500; m <- 10
x <- rnorm(k, 5, 3)
p <- 1/(1+exp(-(-10 + 2*x)))
y <- rBB(k, m, p, phi = 1.2)
fit <- BBreg(y ~ x, m)
summary(fit)
```

## Parity check vs PROreg

```bash
Rscript scripts/parity_BBreg.R
Rscript scripts/parity_BBmm.R
```

**Note on `BBmm`:** when the random-effect signal is weak (few groups / large \(\phi\)), ML Laplace can shrink \(\sigma\to 0\) while PROreg's adjusted profile h-likelihood keeps a small positive value. With identifiable designs (many groups, smaller \(\phi\)), estimates match closely and TMB is typically faster after the one-time compile.