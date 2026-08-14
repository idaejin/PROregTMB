# Demo: BBreg / BBmm with method = "mle" | "bayes" (tmbstan)
# Requires tmbstan + rstan. On R 4.4:
#   install.packages("https://cran.r-project.org/src/contrib/Archive/tmbstan/tmbstan_1.0.91.tar.gz",
#                    repos = NULL, type = "source")
#
#   Rscript scripts/demo_tmbstan_BBreg.R

suppressPackageStartupMessages({
  library(rstan)
  library(tmbstan)
})

args <- commandArgs(trailingOnly = FALSE)
f <- grep("^--file=", args, value = TRUE)
pkg_root <- if (length(f)) {
  normalizePath(file.path(dirname(sub("^--file=", "", f)), ".."))
} else normalizePath(".")
devtools::load_all(pkg_root, quiet = TRUE)

set.seed(7)
k <- 200L; m <- 10L
x <- rnorm(k)
y <- rBB(k, m, 1/(1+exp(-(-0.5 + 0.8*x))), phi = 0.25)
dat <- data.frame(y, x)

cat("=== BBreg method='bayes' ===\n")
fit <- BBreg(y ~ x, m = m, data = dat, method = "bayes",
             chains = 2, iter = 1000, warmup = 400, seed = 1)
print(fit)

# Small mixed model
set.seed(42)
n_g <- 20; n_per <- 5
z <- factor(rep(seq_len(n_g), each = n_per))
u <- rnorm(n_g, 0, 0.6)
xx <- rnorm(n_g * n_per)
yy <- rBB(length(xx), m, 1/(1+exp(-(0.2 - 0.5*xx + u[z]))), phi = 0.2)
dat_mm <- data.frame(y = yy, x = xx, z = z)

cat("\n=== BBmm method='bayes' (laplace_bayes=TRUE) ===\n")
fit_mm <- BBmm(y ~ x, random.formula = ~ z, m = m, data = dat_mm,
               method = "bayes", chains = 2, iter = 600, warmup = 200,
               seed = 1, laplace_bayes = TRUE)
print(fit_mm)

out_dir <- file.path(pkg_root, "scripts", "bench_out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(fit$posterior, file.path(out_dir, "tmbstan_BBreg_demo.csv"),
          row.names = FALSE)
write.csv(fit_mm$posterior, file.path(out_dir, "tmbstan_BBmm_demo.csv"),
          row.names = FALSE)
cat("Saved scripts/bench_out/tmbstan_BBreg_demo.csv and tmbstan_BBmm_demo.csv\n")
