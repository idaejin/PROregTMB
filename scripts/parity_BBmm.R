# Parity: PROregTMB::BBmm vs PROreg::BBmm
# Rscript scripts/parity_BBmm.R

suppressPackageStartupMessages(library(PROreg))

args <- commandArgs(trailingOnly = FALSE)
f <- grep("^--file=", args, value = TRUE)
pkg_root <- if (length(f)) {
  normalizePath(file.path(dirname(sub("^--file=", "", f)), ".."))
} else {
  normalizePath(".")
}
devtools::load_all(pkg_root, quiet = TRUE)

compare_fit <- function(dat, m, label) {
  cat("\n######## ", label, " ########\n", sep = "")
  t0 <- proc.time()
  fit_old <- PROreg::BBmm(
    fixed.formula = y ~ x, random.formula = ~ z, m = m, data = dat
  )
  t_old <- (proc.time() - t0)[["elapsed"]]

  t1 <- proc.time()
  fit_tmb <- PROregTMB::BBmm(
    fixed.formula = y ~ x, random.formula = ~ z, m = m, data = dat
  )
  t_tmb <- (proc.time() - t1)[["elapsed"]]

  cmp <- data.frame(
    param = c(paste0("beta", 0:1), "phi", "sigma"),
    TMB = c(as.numeric(fit_tmb$fixed.coef), fit_tmb$phi.coef,
            as.numeric(fit_tmb$sigma.coef)),
    PROreg = c(as.numeric(fit_old$fixed.coef), fit_old$phi.coef,
               as.numeric(fit_old$sigma.coef))
  )
  cmp$diff <- cmp$TMB - cmp$PROreg
  print(cmp, digits = 4)
  cat("Cor(u_TMB, u_PROreg):",
      round(cor(fit_tmb$random.coef, fit_old$random.coef), 6), "\n")
  cat("Time PROreg:", round(t_old, 3), "s | TMB:", round(t_tmb, 3),
      "s | ratio TMB/PROreg:", round(t_tmb / max(t_old, 1e-6), 2), "\n")
  cat("TMB nll:", round(fit_tmb$nll, 4), " conv:", fit_tmb$conv,
      "| PROreg conv:", fit_old$conv, "\n")
  invisible(list(tmb = fit_tmb, old = fit_old, cmp = cmp))
}

# --- Warm-up compile (excluded from timings below) ---
cat("Warm-up compile...\n")
set.seed(1)
dat0 <- data.frame(y = rBB(40, 5, 0.4, 0.2), x = rnorm(40),
                   z = factor(rep(1:8, each = 5)))
invisible(PROregTMB::BBmm(y ~ x, random.formula = ~ z, m = 5, data = dat0))

# A) Man-page style (weak RE / few levels) — sigma may differ
set.seed(7)
k <- 100
m <- 10
phi <- 0.5
beta <- c(1.5, -1.1)
sigma <- 0.5
x <- runif(k, 0, 10)
z <- as.factor(PROreg::rBI(k, 4, 0.5, 2))
Z <- model.matrix(~ z - 1)
u <- rnorm(nlevels(z), 0, sigma)
eta <- beta[1] + beta[2] * x + as.numeric(Z %*% u)
p <- 1 / (1 + exp(-eta))
y <- PROregTMB::rBB(k, m, p, phi)
dat_A <- data.frame(y, x, z)
cat("\nTrue (A man-page): beta=", beta, " phi=", phi, " sigma=", sigma,
    " nlev=", nlevels(z), "\n", sep = " ")
res_A <- compare_fit(dat_A, m, "A) PROreg man-page simulation")

# B) Identifiable RE: many groups, smaller phi
set.seed(42)
n_g <- 40
n_per <- 8
k <- n_g * n_per
m <- 10
phi <- 0.15
beta <- c(0.5, -0.4)
sigma <- 0.8
x <- rnorm(k)
z <- factor(rep(seq_len(n_g), each = n_per))
Z <- model.matrix(~ z - 1)
u <- rnorm(n_g, 0, sigma)
eta <- beta[1] + beta[2] * x + as.numeric(Z %*% u)
p <- 1 / (1 + exp(-eta))
y <- PROregTMB::rBB(k, m, p, phi)
dat_B <- data.frame(y, x, z)
cat("\nTrue (B strong RE): beta=", beta, " phi=", phi, " sigma=", sigma,
    " nlev=", n_g, "\n", sep = " ")
res_B <- compare_fit(dat_B, m, "B) Strong random-effect signal")
