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

First call to `BBreg()` / `BBmm()` / `BBjm()` compiles the TMB template (needs a C++ toolchain).

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

### Mixed model (`BBmm`)

Random-intercept beta-binomial mixed model:

```r
set.seed(42)
n_g <- 30; n_per <- 6; m <- 10
x <- rnorm(n_g * n_per)
z <- factor(rep(seq_len(n_g), each = n_per))
u <- rnorm(n_g, 0, 0.8)
eta <- 0.5 - 0.4 * x + u[z]
y <- rBB(length(x), m, 1/(1+exp(-eta)), phi = 0.15)
dat <- data.frame(y, x, z)

fit_mm <- BBmm(y ~ x, random.formula = ~ z, m = m, data = dat)
summary(fit_mm)
```

For multiple crossed factors use e.g. `random.formula = ~ site + doc`. For a custom `Z`, pass `Z` and `nRandComp` instead of `random.formula`.

### One-stage joint model (`BBjm`)

Shared RI/RS beta-binomial longitudinal model + Weibull PH survival with
current-value association \(\alpha\, m_i(t)\).

```r
# long: id, time, y (BB score); surv: id, time, status (1 = event)
fit_jm <- BBjm(long, surv, m = 30)
summary(fit_jm)
```

**COPD application (St George ACTIVITY + incident fall).** Time in days since
baseline; ACTIVITY binned to 0–24 (`m = 30`); exclude baseline falls.
On N = 506 (110 events): \(\hat\alpha \approx 0.089\) (SE 0.019, p < 0.001) —
higher predicted activity limitation is associated with higher fall hazard.
Two-stage TSBB on the same data gave \(\hat\alpha \approx 0.085\).

Sketch of the data prep used in that analysis:

```r
# COPD has columns id, date, time (visit), ACTIVITY, fall (1 = fall, 2 = no)
# y <- activity_bin(ACTIVITY)          # 0..24
# long <- data.frame(id, time = days_since_baseline, y)
# surv <- one row per id: time to first post-baseline fall (or censor), status
fit_jm <- BBjm(long, surv, m = 30)
# compare: stage-1 BBmm + survival::coxph counting process (TSBB)
```

## Benchmarks and vignette

Precomputed pilot results (multi-RE, TSBB, BBjm vs TSBB, multidimensional BBMM):

```r
browseVignettes("PROregTMB")
# or from source:
# rmarkdown::render("vignettes/PROregTMB-results.Rmd")
```

| Study | Script | Pilot scale | Output |
|-------|--------|-------------|--------|
| Multi-RE `BBmm` | `scripts/bench_multi_RE.R` | crossed `~site+doc` | `bench_out/bench_multi_RE.csv` |
| TSBB stage-1 engines | `scripts/sim_tsbb_galan.R` | N=80, nsim=8 | `tsbb_galan_summary.csv` |
| BBjm vs TSBB | `scripts/sim_BBjm_vs_TSBB.R` | N=100, nsim=40 | `bbjm_vs_tsbb_summary.csv` |
| Multidimensional BBMM | `scripts/sim_multi_BB.R` | nsim=30, σ∈{0.5,1} | `multi_BB_summary.csv` |

```bash
Rscript scripts/sim_tsbb_galan.R --quick --nsim=8 --N=80
Rscript scripts/sim_BBjm_vs_TSBB.R --quick --nsim=40 --N=100
Rscript scripts/sim_multi_BB.R --quick --nsim=30
Rscript scripts/demo_BBjm.R
# full paper grids (slow): omit --quick; see script headers
```

**Pilot takeaways**

- Multi-RE / multi-outcome `BBmm`: estimates close to PROreg; TMB typically **~3–30×** faster after one-time compile.
- TSBB with TMB stage-1: \(\hat\alpha\) matches PROreg stage-1 (~30–47× faster on the pilot).
- BBjm vs TSBB (TMB stage-1): similar bias/CP for \(\alpha\); TSBB faster; both 100% conv on the pilot.

**Note on `BBmm`:** when the random-effect signal is weak (few groups / large \(\phi\)), ML Laplace can shrink \(\sigma\to 0\) while PROreg's adjusted profile h-likelihood keeps a small positive value. With identifiable designs (many groups, smaller \(\phi\)), estimates match closely.
