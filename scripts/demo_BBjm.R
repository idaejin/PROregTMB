# Demo: 1-stage joint model (BBjm) vs two-stage TSBB (BBmm + coxph)
# Rscript scripts/demo_BBjm.R

suppressPackageStartupMessages(library(survival))

args <- commandArgs(trailingOnly = FALSE)
f <- grep("^--file=", args, value = TRUE)
pkg_root <- if (length(f)) {
  normalizePath(file.path(dirname(sub("^--file=", "", f)), ".."))
} else normalizePath(".")
devtools::load_all(pkg_root, quiet = TRUE)

invlogit <- function(x) 1 / (1 + exp(-x))
mi_fun <- function(t, eta0, eta1, m) m * invlogit(eta0 + eta1 * t)

cumhaz_tv <- function(t, eta0, eta1, alpha, m, w1 = 0.1, w2 = 1.6) {
  if (t <= 0) return(0)
  f <- function(s) {
    s <- pmax(s, 1e-12)
    w1 * w2 * s^(w2 - 1) * exp(alpha * mi_fun(s, eta0, eta1, m))
  }
  integrate(f, 0, t, rel.tol = 1e-4)$value
}

rweibull_tv <- function(eta0, eta1, alpha, m, w1 = 0.1, w2 = 1.6, t_max = 30) {
  u <- runif(1)
  target <- -log(u)
  if (cumhaz_tv(t_max, eta0, eta1, alpha, m, w1, w2) < target) return(t_max)
  uniroot(function(tt) cumhaz_tv(tt, eta0, eta1, alpha, m, w1, w2) - target,
          c(1e-6, t_max))$root
}

simulate_jm <- function(N, beta, sigma, m, phi, alpha,
                        visits = c(0, 1, 2, 5), seed = 1) {
  set.seed(seed)
  u0 <- rnorm(N, 0, sigma[1])
  u1 <- rnorm(N, 0, sigma[2])
  Ai <- max(visits)
  long <- list()
  surv <- data.frame(id = seq_len(N), time = NA_real_, status = NA_integer_)
  for (i in seq_len(N)) {
    eta0 <- beta[1] + u0[i]
    eta1 <- beta[2] + u1[i]
    Tstar <- rweibull_tv(eta0, eta1, alpha, m)
    Ci <- if (runif(1) < 0.1) runif(1, 0, Ai) else Inf
    Ti <- min(Tstar, Ai, Ci)
    delta <- as.integer(abs(Ti - Tstar) < 1e-10)
    surv$time[i] <- Ti
    surv$status[i] <- delta
    tj <- visits[visits <= Ti + 1e-10]
    if (!length(tj)) tj <- 0
    y <- rBB(length(tj), m, invlogit(eta0 + eta1 * tj), phi)
    long[[i]] <- data.frame(id = i, time = tj, y = as.numeric(y))
  }
  list(long = do.call(rbind, long), surv = surv, m = m,
       truth = list(beta = beta, sigma = sigma, phi = phi, alpha = alpha,
                    w1 = 0.1, w2 = 1.6, u0 = u0, u1 = u1))
}

# Two-stage TSBB with PROregTMB::BBmm
fit_tsbb <- function(dat) {
  long <- dat$long
  N <- nrow(dat$surv)
  idf <- factor(long$id)
  Z0 <- model.matrix(~ idf - 1)
  Z <- cbind(Z0, Z0 * long$time)
  X <- model.matrix(~ time, long)
  t0 <- proc.time()[["elapsed"]]
  bb <- BBmm(X = X, y = long$y, Z = Z, nRandComp = c(N, N), m = dat$m)
  t_bb <- proc.time()[["elapsed"]] - t0
  if (!identical(bb$conv, "yes")) return(list(ok = FALSE, time = t_bb))
  beta <- as.numeric(bb$fixed.coef)
  u <- as.numeric(bb$random.coef)
  u0 <- u[seq_len(N)]
  u1 <- u[N + seq_len(N)]
  # coxph counting process
  rows <- list()
  for (i in seq_len(N)) {
    Ti <- dat$surv$time[i]
    delta <- dat$surv$status[i]
    eta0 <- beta[1] + u0[i]
    eta1 <- beta[2] + u1[i]
    grid <- seq(0, Ti, by = 0.25)
    if (tail(grid, 1) < Ti - 1e-8) grid <- c(grid, Ti)
    if (length(grid) < 2) grid <- c(0, Ti)
    for (k in seq_len(length(grid) - 1)) {
      rows[[length(rows) + 1]] <- data.frame(
        tstart = grid[k], tstop = grid[k + 1],
        event = as.integer(k == length(grid) - 1 & delta == 1),
        yhat = mi_fun(grid[k], eta0, eta1, dat$m)
      )
    }
  }
  dd <- do.call(rbind, rows)
  t1 <- proc.time()[["elapsed"]]
  cx <- coxph(Surv(tstart, tstop, event) ~ yhat, data = dd)
  t_cox <- proc.time()[["elapsed"]] - t1
  list(
    ok = TRUE,
    alpha = unname(coef(cx)["yhat"]),
    se = summary(cx)$coefficients["yhat", "se(coef)"],
    beta = beta,
    phi = as.numeric(bb$phi.coef)[1],
    sigma = as.numeric(bb$sigma.coef)[1:2],
    time = t_bb + t_cox,
    time_bbmm = t_bb
  )
}

# ---- run ----
cat("Simulating Scenario-1-like data...\n")
dat <- simulate_jm(
  N = 120,
  beta = c(-0.19, 0.03),
  sigma = c(1.2, 0.05),
  m = 24,
  phi = 0.5,
  alpha = 0.05,
  seed = 42
)
cat("n obs:", nrow(dat$long), " events:", sum(dat$surv$status), "\n")

cat("\n=== Two-stage TSBB (BBmm + coxph) ===\n")
ts <- fit_tsbb(dat)
print(ts[c("ok", "alpha", "se", "beta", "phi", "sigma", "time")])

cat("\n=== One-stage BBjm (TMB) ===\n")
# force recompile if template changed
ensure_tmb_dll("bb_jm", force = TRUE)
t0 <- proc.time()[["elapsed"]]
jm <- BBjm(dat$long, dat$surv, m = dat$m, n_quad = 40)
t_jm <- proc.time()[["elapsed"]] - t0
print(jm)
cat("Wall time BBjm:", round(t_jm, 3), "s\n")

cmp <- data.frame(
  param = c("beta0", "beta1", "phi", "sigma0", "sigma1", "alpha", "w1", "w2"),
  true = c(dat$truth$beta, dat$truth$phi, dat$truth$sigma,
           dat$truth$alpha, dat$truth$w1, dat$truth$w2),
  TSBB = c(ts$beta, ts$phi, ts$sigma, ts$alpha, NA, NA),
  BBjm = c(jm$fixed.coef, jm$phi, as.numeric(jm$sigma),
           jm$alpha, jm$w1, jm$w2)
)
cmp$diff_jm_ts <- cmp$BBjm - cmp$TSBB
cat("\n=== Comparison ===\n")
print(cmp, digits = 4)
cat("\nTime TSBB:", round(ts$time, 3), "s | BBjm:", round(t_jm, 3), "s\n")

out <- file.path(pkg_root, "scripts", "bench_out", "demo_BBjm.csv")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
utils::write.csv(cmp, out, row.names = FALSE)
cat("Wrote", out, "\n")
