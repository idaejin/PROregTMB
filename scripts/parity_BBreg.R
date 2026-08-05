# Parity check: PROregTMB::BBreg vs PROreg::BBreg
# Run from package root: Rscript scripts/parity_BBreg.R

suppressPackageStartupMessages({
  library(PROreg)
})

pkg_root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Install devtools to run this script with load_all")
}
devtools::load_all(pkg_root, quiet = TRUE)

set.seed(18)
k <- 1000
m <- 10
x <- rnorm(k, 5, 3)
beta_true <- c(-10, 2)
p <- 1 / (1 + exp(-(beta_true[1] + beta_true[2] * x)))
phi_true <- 1.2
y <- PROregTMB::rBB(k, m, p, phi_true)

cat("=== Compiling / fitting TMB ===\n")
fit_tmb <- PROregTMB::BBreg(y ~ x, m)
cat("=== Fitting PROreg ===\n")
fit_old <- PROreg::BBreg(y ~ x, m)

cmp <- data.frame(
  param = c(paste0("beta", 0:1), "phi"),
  true = c(beta_true, phi_true),
  TMB = c(as.numeric(fit_tmb$coefficients), fit_tmb$phi),
  PROreg = c(as.numeric(fit_old$coefficients), fit_old$phi)
)
cmp$diff <- cmp$TMB - cmp$PROreg
print(cmp, digits = 4)

cat("\nTMB nll:", fit_tmb$nll, "  iters:", fit_tmb$iter,
    "  conv:", fit_tmb$conv, "\n")
cat("PROreg iters:", fit_old$iter, "  conv:", fit_old$conv, "\n")
