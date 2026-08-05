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
| `BBjm` one-stage joint (BB + Weibull) | done (TMB) |

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

## One-stage joint model

```r
fit <- BBjm(long, surv, m = 24)  # long: id,time,y; surv: id,time,status
summary(fit)
```

Demo vs two-stage:

```bash
Rscript scripts/demo_BBjm.R
```

Pilot comparing TSBB stage-1 engines (`PROreg::BBmm` vs `PROregTMB::BBmm`):

```bash
Rscript scripts/sim_tsbb_galan.R --quick --nsim=8 --N=80
# full paper grid (slow): omit --quick; default N=250
```

Results: `scripts/bench_out/tsbb_galan_summary.csv`

**Note on `BBmm`:** when the random-effect signal is weak (few groups / large \(\phi\)), ML Laplace can shrink \(\sigma\to 0\) while PROreg's adjusted profile h-likelihood keeps a small positive value. With identifiable designs (many groups, smaller \(\phi\)), estimates match closely and TMB is typically faster after the one-time compile.