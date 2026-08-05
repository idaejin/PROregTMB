# Monte Carlo: one-stage BBjm vs two-stage TSBB (PROregTMB::BBmm + coxph)
# Metrics for alpha (and longitudinal params): %Bias, ESD, ASD, CP95, time, conv
#
# Design: Galán-Arcicollar et al. SORT 2024 Sec. 6
#
# Usage:
#   Rscript scripts/sim_BBjm_vs_TSBB.R --quick --nsim=40 --N=100
#   Rscript scripts/sim_BBjm_vs_TSBB.R --nsim=100 --N=250
#   Rscript scripts/sim_BBjm_vs_TSBB.R --scenario=1 --phi=0.5 --alpha=0.05 --nsim=50

suppressPackageStartupMessages(library(survival))

args <- commandArgs(trailingOnly = FALSE)
f <- grep("^--file=", args, value = TRUE)
pkg_root <- if (length(f)) {
  normalizePath(file.path(dirname(sub("^--file=", "", f)), ".."))
} else normalizePath(".")

user_args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), user_args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}
has_flag <- function(name) any(user_args == paste0("--", name))

NSIM <- as.integer(get_arg("nsim", "40"))
N_SUBJ <- as.integer(get_arg("N", "250"))
QUICK <- has_flag("quick")
ONLY_SC <- get_arg("scenario", NULL)
ONLY_PHI <- get_arg("phi", NULL)
ONLY_ALPHA <- get_arg("alpha", NULL)
if (QUICK && is.null(get_arg("N", NULL))) N_SUBJ <- 100L

devtools::load_all(pkg_root, quiet = TRUE)

# ---- DGP (paper) -------------------------------------------------------------

invlogit <- function(x) 1 / (1 + exp(-x))
mi_fun <- function(t, eta0, eta1, m) m * invlogit(eta0 + eta1 * t)

cumhaz_tv <- function(t, eta0, eta1, alpha, m, w1 = 0.1, w2 = 1.6) {
  if (t <= 0) return(0)
  f <- function(s) {
    s <- pmax(s, .Machine$double.eps)
    w1 * w2 * s^(w2 - 1) * exp(alpha * mi_fun(s, eta0, eta1, m))
  }
  tryCatch(
    stats::integrate(f, 0, t, rel.tol = 1e-4)$value,
    error = function(e) {
      gg <- seq(1e-6, t, length.out = 200)
      sum(f(gg) * diff(c(0, gg)))
    }
  )
}

rweibull_tv <- function(eta0, eta1, alpha, m, w1 = 0.1, w2 = 1.6, t_max = 30) {
  u <- runif(1)
  target <- -log(u)
  if (!is.finite(cumhaz_tv(t_max, eta0, eta1, alpha, m, w1, w2)) ||
      cumhaz_tv(t_max, eta0, eta1, alpha, m, w1, w2) < target) {
    return(t_max)
  }
  uniroot(function(tt) cumhaz_tv(tt, eta0, eta1, alpha, m, w1, w2) - target,
          c(1e-6, t_max), tol = 1e-4)$root
}

simulate_one <- function(N, beta, sigma, m, phi, alpha,
                         visit_times = c(0, 1, 2, 5),
                         w1 = 0.1, w2 = 1.6, dropout_rate = 0.10,
                         seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  u0 <- rnorm(N, 0, sigma[1])
  u1 <- rnorm(N, 0, sigma[2])
  Ai <- max(visit_times)
  long <- vector("list", N)
  surv <- data.frame(id = seq_len(N), time = NA_real_, status = NA_integer_)
  for (i in seq_len(N)) {
    eta0 <- beta[1] + u0[i]
    eta1 <- beta[2] + u1[i]
    Tstar <- rweibull_tv(eta0, eta1, alpha, m, w1, w2)
    Ci <- if (runif(1) < dropout_rate) runif(1, min(visit_times), Ai) else Inf
    Ti <- min(Tstar, Ai, Ci)
    delta <- as.integer(abs(Ti - Tstar) < 1e-10)
    surv$time[i] <- Ti
    surv$status[i] <- delta
    tj <- visit_times[visit_times <= Ti + 1e-10]
    if (!length(tj)) tj <- 0
    y <- rBB(length(tj), m, invlogit(eta0 + eta1 * tj), phi)
    long[[i]] <- data.frame(id = i, time = tj, y = as.numeric(y))
  }
  list(long = do.call(rbind, long), surv = surv, m = m)
}

# ---- Fits --------------------------------------------------------------------

fit_tsbb <- function(dat) {
  long <- dat$long
  N <- nrow(dat$surv)
  idf <- factor(long$id, levels = seq_len(N))
  Z0 <- model.matrix(~ idf - 1)
  Z <- cbind(Z0, Z0 * long$time)
  X <- model.matrix(~ time, long)
  t0 <- proc.time()[["elapsed"]]
  bb <- tryCatch(
    BBmm(X = X, y = long$y, Z = Z, nRandComp = c(N, N), m = dat$m, silent = TRUE),
    error = function(e) e
  )
  t_bb <- proc.time()[["elapsed"]] - t0
  if (inherits(bb, "error") || !identical(bb$conv, "yes")) {
    return(list(ok = FALSE, time = t_bb, time_bb = t_bb))
  }
  beta <- as.numeric(bb$fixed.coef)
  u <- as.numeric(bb$random.coef)
  if (length(u) != 2 * N) return(list(ok = FALSE, time = t_bb, time_bb = t_bb))
  u0 <- u[seq_len(N)]
  u1 <- u[N + seq_len(N)]

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
  cx <- tryCatch(
    coxph(Surv(tstart, tstop, event) ~ yhat, data = dd),
    error = function(e) e
  )
  t_cox <- proc.time()[["elapsed"]] - t1
  if (inherits(cx, "error")) {
    return(list(ok = FALSE, time = t_bb + t_cox, time_bb = t_bb))
  }
  list(
    ok = TRUE,
    alpha = unname(coef(cx)["yhat"]),
    se = unname(summary(cx)$coefficients["yhat", "se(coef)"]),
    beta0 = beta[1], beta1 = beta[2],
    phi = as.numeric(bb$phi.coef)[1],
    sigma0 = as.numeric(bb$sigma.coef)[1],
    sigma1 = as.numeric(bb$sigma.coef)[2],
    time = t_bb + t_cox,
    time_bb = t_bb
  )
}

fit_bbjm <- function(dat) {
  t0 <- proc.time()[["elapsed"]]
  jm <- tryCatch(
    BBjm(dat$long, dat$surv, m = dat$m, n_quad = 40L,
         use_bbmm_start = FALSE, silent = TRUE, maxiter = 200),
    error = function(e) e
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  if (inherits(jm, "error") || !identical(jm$conv, "yes")) {
    return(list(ok = FALSE, time = elapsed))
  }
  list(
    ok = TRUE,
    alpha = jm$alpha,
    se = jm$alpha.se,
    beta0 = unname(jm$fixed.coef[1]),
    beta1 = unname(jm$fixed.coef[2]),
    phi = jm$phi,
    sigma0 = unname(jm$sigma[1]),
    sigma1 = unname(jm$sigma[2]),
    w1 = jm$w1,
    w2 = jm$w2,
    time = elapsed,
    nll = jm$nll
  )
}

summarize_alpha <- function(est, se, true, times, nsim) {
  ok <- which(is.finite(est) & is.finite(se))
  if (!length(ok)) {
    return(data.frame(
      n_ok = 0, pct_conv = 0, mean = NA_real_, bias = NA_real_,
      pct_bias = NA_real_, RMSE = NA_real_, ESD = NA_real_, ASD = NA_real_,
      CP95 = NA_real_, time_mean = mean(times, na.rm = TRUE)
    ))
  }
  a <- est[ok]
  s <- se[ok]
  cover <- mean(a - 1.96 * s <= true & a + 1.96 * s >= true)
  data.frame(
    n_ok = length(ok),
    pct_conv = 100 * length(ok) / nsim,
    mean = mean(a),
    bias = mean(a) - true,
    pct_bias = (mean(a) - true) / true,
    RMSE = sqrt(mean((a - true)^2)),
    ESD = stats::sd(a),
    ASD = mean(s),
    CP95 = 100 * cover,
    time_mean = mean(times, na.rm = TRUE)
  )
}

# ---- Grid --------------------------------------------------------------------

sc_main <- list(
  list(id = 1L, beta = c(-0.19, 0.03), sigma = c(1.2, 0.05), m = 24L,
       alphas = c(0.01, 0.05, 0.10)),
  list(id = 2L, beta = c(0.40, -0.15), sigma = c(1.5, 0.30), m = 8L,
       alphas = c(-0.05, -0.10, -0.15))
)
phis <- c(0.05, 0.5, 1)

grid <- list()
for (sc in sc_main) {
  for (ph in phis) {
    for (al in sc$alphas) {
      grid[[length(grid) + 1]] <- list(
        scenario = sc$id, beta = sc$beta, sigma = sc$sigma,
        m = sc$m, phi = ph, alpha = al
      )
    }
  }
}

if (QUICK) {
  grid <- list(
    list(scenario = 1L, beta = c(-0.19, 0.03), sigma = c(1.2, 0.05),
         m = 24L, phi = 0.05, alpha = 0.01),
    list(scenario = 1L, beta = c(-0.19, 0.03), sigma = c(1.2, 0.05),
         m = 24L, phi = 0.5, alpha = 0.05),
    list(scenario = 1L, beta = c(-0.19, 0.03), sigma = c(1.2, 0.05),
         m = 24L, phi = 0.5, alpha = 0.10),
    list(scenario = 2L, beta = c(0.40, -0.15), sigma = c(1.5, 0.30),
         m = 8L, phi = 0.5, alpha = -0.10),
    list(scenario = 2L, beta = c(0.40, -0.15), sigma = c(1.5, 0.30),
         m = 8L, phi = 1, alpha = -0.15)
  )
}

if (!is.null(ONLY_SC)) {
  grid <- Filter(function(g) g$scenario == as.integer(ONLY_SC), grid)
}
if (!is.null(ONLY_PHI)) {
  grid <- Filter(function(g) abs(g$phi - as.numeric(ONLY_PHI)) < 1e-12, grid)
}
if (!is.null(ONLY_ALPHA)) {
  grid <- Filter(function(g) abs(g$alpha - as.numeric(ONLY_ALPHA)) < 1e-12, grid)
}

out_dir <- file.path(pkg_root, "scripts", "bench_out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
rep_path <- file.path(out_dir, "bbjm_vs_tsbb_reps.csv")
sum_path <- file.path(out_dir, "bbjm_vs_tsbb_summary.csv")

cat("BBjm vs TSBB Monte Carlo\n")
cat("nsim =", NSIM, "| N =", N_SUBJ, "| cells =", length(grid),
    if (QUICK) "(quick)" else "", "\n")
flush.console()

cat("Warm-up TMB DLLs...\n")
ensure_tmb_dll("bb_mm")
ensure_tmb_dll("bb_jm")
wu <- simulate_one(40, c(-0.19, 0.03), c(1.2, 0.05), 24, 0.5, 0.05, seed = 0)
invisible(BBjm(wu$long, wu$surv, m = 24, use_bbmm_start = FALSE, silent = TRUE))

rep_rows <- list()
sum_rows <- list()

for (g in grid) {
  tag <- sprintf("Sc%d_phi%g_alpha%g", g$scenario, g$phi, g$alpha)
  cat("\n====", tag, "====\n")
  flush.console()

  a_ts <- se_ts <- t_ts <- rep(NA_real_, NSIM)
  a_jm <- se_jm <- t_jm <- rep(NA_real_, NSIM)
  ok_ts <- ok_jm <- 0L

  for (s in seq_len(NSIM)) {
    seed <- 7000L * g$scenario + as.integer(g$phi * 100) +
      as.integer(abs(g$alpha) * 1000) + s
    dat <- simulate_one(
      N_SUBJ, g$beta, g$sigma, g$m, g$phi, g$alpha, seed = seed
    )
    ts <- fit_tsbb(dat)
    jm <- fit_bbjm(dat)

    if (isTRUE(ts$ok)) {
      ok_ts <- ok_ts + 1L
      a_ts[s] <- ts$alpha
      se_ts[s] <- ts$se
    }
    t_ts[s] <- ts$time
    if (isTRUE(jm$ok)) {
      ok_jm <- ok_jm + 1L
      a_jm[s] <- jm$alpha
      se_jm[s] <- jm$se
    }
    t_jm[s] <- jm$time

    row <- data.frame(
      tag = tag, scenario = g$scenario, phi = g$phi, alpha_true = g$alpha,
      N = N_SUBJ, sim = s,
      ok_tsbb = isTRUE(ts$ok), ok_bbjm = isTRUE(jm$ok),
      alpha_tsbb = if (isTRUE(ts$ok)) ts$alpha else NA_real_,
      se_tsbb = if (isTRUE(ts$ok)) ts$se else NA_real_,
      alpha_bbjm = if (isTRUE(jm$ok)) jm$alpha else NA_real_,
      se_bbjm = if (isTRUE(jm$ok)) jm$se else NA_real_,
      beta0_tsbb = if (isTRUE(ts$ok)) ts$beta0 else NA_real_,
      beta0_bbjm = if (isTRUE(jm$ok)) jm$beta0 else NA_real_,
      beta1_tsbb = if (isTRUE(ts$ok)) ts$beta1 else NA_real_,
      beta1_bbjm = if (isTRUE(jm$ok)) jm$beta1 else NA_real_,
      phi_tsbb = if (isTRUE(ts$ok)) ts$phi else NA_real_,
      phi_bbjm = if (isTRUE(jm$ok)) jm$phi else NA_real_,
      time_tsbb = ts$time,
      time_bbjm = jm$time,
      stringsAsFactors = FALSE
    )
    rep_rows[[length(rep_rows) + 1]] <- row
    utils::write.csv(do.call(rbind, rep_rows), rep_path, row.names = FALSE)

    if (s %% 5 == 0 || s == NSIM) {
      cat(sprintf(
        "  sim %d/%d | conv TS:%d JM:%d | aTS=%s aJM=%s | tTS=%.1fs tJM=%.1fs\n",
        s, NSIM, ok_ts, ok_jm,
        if (isTRUE(ts$ok)) sprintf("%.4f", ts$alpha) else "NA",
        if (isTRUE(jm$ok)) sprintf("%.4f", jm$alpha) else "NA",
        ts$time, jm$time
      ))
      flush.console()
    }
  }

  st <- summarize_alpha(a_ts, se_ts, g$alpha, t_ts, NSIM)
  sj <- summarize_alpha(a_jm, se_jm, g$alpha, t_jm, NSIM)

  sum_rows[[length(sum_rows) + 1]] <- data.frame(
    tag = tag, scenario = g$scenario, phi = g$phi, alpha_true = g$alpha,
    method = c("TSBB", "BBjm"),
    n_ok = c(st$n_ok, sj$n_ok),
    pct_conv = c(st$pct_conv, sj$pct_conv),
    mean_alpha = c(st$mean, sj$mean),
    bias = c(st$bias, sj$bias),
    pct_bias = c(st$pct_bias, sj$pct_bias),
    RMSE = c(st$RMSE, sj$RMSE),
    ESD = c(st$ESD, sj$ESD),
    ASD = c(st$ASD, sj$ASD),
    CP95 = c(st$CP95, sj$CP95),
    time_mean = c(st$time_mean, sj$time_mean),
    speedup_vs_tsbb = c(NA_real_, st$time_mean / max(sj$time_mean, 1e-8)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(do.call(rbind, sum_rows), sum_path, row.names = FALSE)

  cat(sprintf(
    "  TSBB  %%Bias=%+.3f RMSE=%.4f CP=%.1f time=%.2fs (conv %.0f%%)\n",
    st$pct_bias, st$RMSE, st$CP95, st$time_mean, st$pct_conv
  ))
  cat(sprintf(
    "  BBjm  %%Bias=%+.3f RMSE=%.4f CP=%.1f time=%.2fs (conv %.0f%%) speedup=%.2fx\n",
    sj$pct_bias, sj$RMSE, sj$CP95, sj$time_mean, sj$pct_conv,
    st$time_mean / max(sj$time_mean, 1e-8)
  ))
}

sum_df <- do.call(rbind, sum_rows)
cat("\n=== Summary ===\n")
print(sum_df, digits = 4)
cat("\nWrote:\n ", rep_path, "\n ", sum_path, "\n")
