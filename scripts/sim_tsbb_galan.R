# Two-stage joint model simulation (Galán-Arcicollar et al., SORT 2024)
# Compare TSBB stage-1 engine: PROreg::BBmm vs PROregTMB::BBmm
# Stage-2 Cox identical (survival::coxph).
#
# Paper settings (Section 6):
#   N=250, max 4 visits at times 0,1,2,5 from entry, Weibull(w1=0.1,w2=1.6)
#   Sc1: beta=(-0.19,0.03), sigma=(1.2,0.05), m=24, alpha in {0.01,0.05,0.10}
#   Sc2: beta=(0.40,-0.15), sigma=(1.5,0.30), m=8,  alpha in {-0.05,-0.10,-0.15}
#   phi in {0.05, 0.5, 1}
#
# Usage:
#   Rscript scripts/sim_tsbb_galan.R
#   Rscript scripts/sim_tsbb_galan.R --nsim=20 --quick
#   Rscript scripts/sim_tsbb_galan.R --nsim=50 --scenario=1 --phi=0.5 --alpha=0.05

suppressPackageStartupMessages({
  library(PROreg)
  library(survival)
})

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

NSIM <- as.integer(get_arg("nsim", "25"))
N_SUBJ <- as.integer(get_arg("N", "250"))
QUICK <- has_flag("quick")
ONLY_SC <- get_arg("scenario", NULL)
ONLY_PHI <- get_arg("phi", NULL)
ONLY_ALPHA <- get_arg("alpha", NULL)
if (QUICK && is.null(get_arg("N", NULL))) N_SUBJ <- 100L

devtools::load_all(pkg_root, quiet = TRUE)

# ---- helpers -----------------------------------------------------------------

invlogit <- function(x) 1 / (1 + exp(-x))

# True subject mean score m_i(t)
mi_fun <- function(t, eta0, eta1, m) {
  m * invlogit(eta0 + eta1 * t)
}

# Cumulative hazard under Weibull baseline + association with m_i(t)
cumhaz_tv <- function(t, eta0, eta1, alpha, m, w1 = 0.1, w2 = 1.6) {
  if (t <= 0) return(0)
  f <- function(s) {
    s <- pmax(s, .Machine$double.eps)
    w1 * w2 * s^(w2 - 1) * exp(alpha * mi_fun(s, eta0, eta1, m))
  }
  tryCatch(
    stats::integrate(f, lower = 0, upper = t, rel.tol = 1e-4)$value,
    error = function(e) {
      # fallback quadrature
      gg <- seq(1e-6, t, length.out = 200)
      sum(f(gg) * diff(c(0, gg)))
    }
  )
}

# Invert H(T*) = -log(U)
rweibull_tv <- function(eta0, eta1, alpha, m, w1 = 0.1, w2 = 1.6,
                        t_max = 20) {
  u <- stats::runif(1)
  target <- -log(u)
  Hmax <- cumhaz_tv(t_max, eta0, eta1, alpha, m, w1, w2)
  if (!is.finite(Hmax) || Hmax < target) return(t_max)
  uniroot(
    function(tt) cumhaz_tv(tt, eta0, eta1, alpha, m, w1, w2) - target,
    interval = c(1e-6, t_max),
    tol = 1e-4
  )$root
}

build_Z_ri_rs <- function(id, time) {
  id <- factor(id)
  Z0 <- model.matrix(~ id - 1)
  Z1 <- Z0 * as.numeric(time)
  colnames(Z0) <- paste0("ri_", levels(id))
  colnames(Z1) <- paste0("rs_", levels(id))
  list(Z = cbind(Z0, Z1), nRandComp = c(ncol(Z0), ncol(Z1)))
}

# Simulate one dataset (paper design, time from entry)
simulate_one <- function(N = 250, beta, sigma, m, phi, alpha,
                         visit_times = c(0, 1, 2, 5),
                         w1 = 0.1, w2 = 1.6,
                         dropout_rate = 0.10,
                         seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  u0 <- rnorm(N, 0, sigma[1])
  u1 <- rnorm(N, 0, sigma[2])
  Ai <- max(visit_times)

  long_list <- vector("list", N)
  surv <- data.frame(
    id = seq_len(N), time = NA_real_, status = NA_integer_,
    u0 = u0, u1 = u1
  )

  for (i in seq_len(N)) {
    eta0 <- beta[1] + u0[i]
    eta1 <- beta[2] + u1[i]
    Tstar <- rweibull_tv(eta0, eta1, alpha, m, w1, w2, t_max = 30)
    Ci <- if (runif(1) < dropout_rate) {
      runif(1, min(visit_times), Ai)
    } else {
      Inf
    }
    Ti <- min(Tstar, Ai, Ci)
    delta <- as.integer(abs(Ti - Tstar) < 1e-10)
    surv$time[i] <- Ti
    surv$status[i] <- delta

    tj <- visit_times[visit_times <= Ti + 1e-10]
    if (!length(tj)) tj <- 0
    pij <- invlogit(eta0 + eta1 * tj)
    yij <- PROregTMB::rBB(length(tj), m, pij, phi)
    long_list[[i]] <- data.frame(
      id = i, time = tj, y = as.numeric(yij),
      stringsAsFactors = FALSE
    )
  }

  long <- do.call(rbind, long_list)
  list(long = long, surv = surv, m = m, visit_times = visit_times)
}

# Stage-1 BBmm (random intercept + slope) via Z API
fit_bbmm_engine <- function(long, m, engine = c("PROreg", "PROregTMB")) {
  engine <- match.arg(engine)
  zz <- build_Z_ri_rs(long$id, long$time)
  X <- model.matrix(~ time, data = long)
  y <- long$y

  if (engine == "PROreg") {
    # PROreg prints iteration spam to stdout; silence it
    zznull <- file(nullfile(), open = "wt")
    sink(zznull)
    sink(zznull, type = "message")
    fit <- tryCatch(
      PROreg::BBmm(
        X = X, y = y, Z = zz$Z, nRandComp = zz$nRandComp, m = m,
        maxiter = 50
      ),
      error = function(e) e
    )
    sink(type = "message")
    sink()
    close(zznull)
  } else {
    fit <- tryCatch(
      PROregTMB::BBmm(
        X = X, y = y, Z = zz$Z, nRandComp = zz$nRandComp, m = m,
        maxiter = 100, silent = TRUE
      ),
      error = function(e) e
    )
  }
  fit
}

extract_bbmm <- function(fit, n_subjects) {
  if (inherits(fit, "error") || is.null(fit$conv) || fit$conv != "yes") {
    return(NULL)
  }
  beta <- as.numeric(fit$fixed.coef)
  # Random effects order: first n intercepts, then n slopes
  u <- as.numeric(fit$random.coef)
  if (length(u) != 2 * n_subjects) return(NULL)
  list(
    beta = beta,
    u0 = u[seq_len(n_subjects)],
    u1 = u[n_subjects + seq_len(n_subjects)],
    phi = as.numeric(fit$phi.coef)[1],
    sigma = as.numeric(fit$sigma.coef)[1:2]
  )
}

# Stage-2 Cox with counting-process & predicted score ŷ_i(t)
fit_cox_stage2 <- function(pars, surv, m, grid_step = 0.25) {
  rows <- list()
  for (i in seq_len(nrow(surv))) {
    Ti <- surv$time[i]
    delta <- surv$status[i]
    eta0 <- pars$beta[1] + pars$u0[i]
    eta1 <- pars$beta[2] + pars$u1[i]
    grid <- seq(0, Ti, by = grid_step)
    if (tail(grid, 1) < Ti - 1e-10) grid <- c(grid, Ti)
    if (length(grid) < 2) grid <- c(0, Ti)
    for (k in seq_len(length(grid) - 1)) {
      t0 <- grid[k]
      t1 <- grid[k + 1]
      yhat <- mi_fun(t0, eta0, eta1, m)
      rows[[length(rows) + 1]] <- data.frame(
        id = surv$id[i], tstart = t0, tstop = t1,
        event = as.integer(k == length(grid) - 1 & delta == 1),
        yhat = yhat
      )
    }
  }
  dd <- do.call(rbind, rows)
  fit <- tryCatch(
    survival::coxph(Surv(tstart, tstop, event) ~ yhat, data = dd),
    error = function(e) e
  )
  if (inherits(fit, "error")) return(NULL)
  s <- summary(fit)$coefficients
  list(
    alpha = unname(coef(fit)["yhat"]),
    se = unname(s["yhat", "se(coef)"]),
    conv = fit$iter < 100 # coxph doesn't expose conv cleanly
  )
}

run_one_rep <- function(dat, engine) {
  t0 <- proc.time()[["elapsed"]]
  bb <- fit_bbmm_engine(dat$long, dat$m, engine)
  t_bb <- proc.time()[["elapsed"]] - t0
  N <- nrow(dat$surv)
  pars <- extract_bbmm(bb, N)
  if (is.null(pars)) {
    return(list(ok = FALSE, time_bbmm = t_bb, time_cox = NA_real_,
                time_total = t_bb, engine = engine))
  }
  t1 <- proc.time()[["elapsed"]]
  cox <- fit_cox_stage2(pars, dat$surv, dat$m)
  t_cox <- proc.time()[["elapsed"]] - t1
  if (is.null(cox)) {
    return(list(ok = FALSE, time_bbmm = t_bb, time_cox = t_cox,
                time_total = t_bb + t_cox, engine = engine,
                beta = pars$beta, phi = pars$phi, sigma = pars$sigma))
  }
  list(
    ok = TRUE,
    engine = engine,
    alpha = cox$alpha,
    se = cox$se,
    beta = pars$beta,
    phi = pars$phi,
    sigma = pars$sigma,
    time_bbmm = t_bb,
    time_cox = t_cox,
    time_total = t_bb + t_cox
  )
}

# ---- scenario grid -----------------------------------------------------------

sc_main <- list(
  list(id = 1, beta = c(-0.19, 0.03), sigma = c(1.2, 0.05), m = 24L,
       alphas = c(0.01, 0.05, 0.10)),
  list(id = 2, beta = c(0.40, -0.15), sigma = c(1.5, 0.30), m = 8L,
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
  # 4 representative cells from the paper grid
  grid <- list(
    list(scenario = 1, beta = c(-0.19, 0.03), sigma = c(1.2, 0.05),
         m = 24L, phi = 0.05, alpha = 0.01),
    list(scenario = 1, beta = c(-0.19, 0.03), sigma = c(1.2, 0.05),
         m = 24L, phi = 0.5, alpha = 0.05),
    list(scenario = 2, beta = c(0.40, -0.15), sigma = c(1.5, 0.30),
         m = 8L, phi = 0.5, alpha = -0.10),
    list(scenario = 2, beta = c(0.40, -0.15), sigma = c(1.5, 0.30),
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

cat("TSBB simulation: PROreg vs PROregTMB\n")
cat("nsim =", NSIM, "| N =", N_SUBJ, "| cells =", length(grid),
    if (QUICK) "(quick subset)" else "", "\n")
flush.console()

# Warm-up compile
cat("Warm-up TMB compile...\n")
wu <- simulate_one(N = 40, beta = c(-0.19, 0.03), sigma = c(1.2, 0.05),
                   m = 24, phi = 0.5, alpha = 0.05, seed = 0)
invisible(fit_bbmm_engine(wu$long, wu$m, "PROregTMB"))

out_dir <- file.path(pkg_root, "scripts", "bench_out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rep_rows <- list()
sum_rows <- list()

for (g in grid) {
  tag <- sprintf("Sc%d_phi%g_alpha%g", g$scenario, g$phi, g$alpha)
  cat("\n====", tag, "====\n")
  alphas_p <- alphas_t <- ses_p <- ses_t <- numeric(0)
  times_p <- times_t <- numeric(0)
  ok_p <- ok_t <- 0L

  for (s in seq_len(NSIM)) {
    dat <- simulate_one(
      N = N_SUBJ, beta = g$beta, sigma = g$sigma, m = g$m,
      phi = g$phi, alpha = g$alpha, seed = 1000 * g$scenario +
        as.integer(g$phi * 100) + as.integer(abs(g$alpha) * 1000) + s
    )
    rp <- run_one_rep(dat, "PROreg")
    rt <- run_one_rep(dat, "PROregTMB")

    if (isTRUE(rp$ok)) {
      ok_p <- ok_p + 1L
      alphas_p <- c(alphas_p, rp$alpha)
      ses_p <- c(ses_p, rp$se)
      times_p <- c(times_p, rp$time_total)
    }
    if (isTRUE(rt$ok)) {
      ok_t <- ok_t + 1L
      alphas_t <- c(alphas_t, rt$alpha)
      ses_t <- c(ses_t, rt$se)
      times_t <- c(times_t, rt$time_total)
    }

    row <- data.frame(
      tag = tag, scenario = g$scenario, phi = g$phi, alpha_true = g$alpha,
      N = N_SUBJ, sim = s,
      ok_proreg = isTRUE(rp$ok), ok_tmb = isTRUE(rt$ok),
      alpha_proreg = if (isTRUE(rp$ok)) rp$alpha else NA_real_,
      alpha_tmb = if (isTRUE(rt$ok)) rt$alpha else NA_real_,
      se_proreg = if (isTRUE(rp$ok)) rp$se else NA_real_,
      se_tmb = if (isTRUE(rt$ok)) rt$se else NA_real_,
      time_proreg = rp$time_total,
      time_tmb = rt$time_total,
      time_bbmm_proreg = rp$time_bbmm,
      time_bbmm_tmb = rt$time_bbmm,
      stringsAsFactors = FALSE
    )
    rep_rows[[length(rep_rows) + 1]] <- row
    # incremental checkpoint
    utils::write.csv(do.call(rbind, rep_rows),
                     file.path(out_dir, "tsbb_galan_reps.csv"),
                     row.names = FALSE)

    cat(sprintf("  sim %d/%d | conv P:%d T:%d | times P:%.1fs T:%.1fs | aP=%s aT=%s\n",
                s, NSIM, ok_p, ok_t,
                rp$time_total, rt$time_total,
                if (isTRUE(rp$ok)) sprintf("%.4f", rp$alpha) else "NA",
                if (isTRUE(rt$ok)) sprintf("%.4f", rt$alpha) else "NA"))
    flush.console()
  }

  summarize_engine <- function(a, se, times, true_a, ok, nsim) {
    if (!length(a)) {
      return(list(
        n_ok = 0, pct_conv = 0, mean = NA, bias = NA, pct_bias = NA,
        esd = NA, asd = NA, cp = NA, time_mean = mean(times, na.rm = TRUE)
      ))
    }
    cover <- mean(a - 1.96 * se < true_a & a + 1.96 * se > true_a, na.rm = TRUE)
    list(
      n_ok = length(a),
      pct_conv = 100 * length(a) / nsim,
      mean = mean(a),
      bias = mean(a) - true_a,
      pct_bias = (mean(a) - true_a) / true_a,
      esd = stats::sd(a),
      asd = mean(se),
      cp = 100 * cover,
      time_mean = mean(times, na.rm = TRUE)
    )
  }

  sp <- summarize_engine(alphas_p, ses_p, times_p, g$alpha, ok_p, NSIM)
  st <- summarize_engine(alphas_t, ses_t, times_t, g$alpha, ok_t, NSIM)

  sum_rows[[length(sum_rows) + 1]] <- data.frame(
    tag = tag, scenario = g$scenario, phi = g$phi, alpha_true = g$alpha,
    engine = c("PROreg", "PROregTMB"),
    n_ok = c(sp$n_ok, st$n_ok),
    pct_conv = c(sp$pct_conv, st$pct_conv),
    mean_alpha = c(sp$mean, st$mean),
    bias = c(sp$bias, st$bias),
    pct_bias = c(sp$pct_bias, st$pct_bias),
    ESD = c(sp$esd, st$esd),
    ASD = c(sp$asd, st$asd),
    CP95 = c(sp$cp, st$cp),
    time_mean = c(sp$time_mean, st$time_mean),
    speedup = c(NA_real_, sp$time_mean / st$time_mean),
    stringsAsFactors = FALSE
  )

  cat(sprintf(
    "  PROreg  %%Bias=%.3f ESD=%.4f CP=%.1f time=%.2fs (conv %.0f%%)\n",
    sp$pct_bias, sp$esd, sp$cp, sp$time_mean, sp$pct_conv
  ))
  cat(sprintf(
    "  TMB     %%Bias=%.3f ESD=%.4f CP=%.1f time=%.2fs (conv %.0f%%) speedup=%.2fx\n",
    st$pct_bias, st$esd, st$cp, st$time_mean, st$pct_conv,
    sp$time_mean / max(st$time_mean, 1e-8)
  ))
}

rep_df <- do.call(rbind, rep_rows)
sum_df <- do.call(rbind, sum_rows)

rep_path <- file.path(out_dir, "tsbb_galan_reps.csv")
sum_path <- file.path(out_dir, "tsbb_galan_summary.csv")
utils::write.csv(rep_df, rep_path, row.names = FALSE)
utils::write.csv(sum_df, sum_path, row.names = FALSE)

cat("\n=== Summary ===\n")
print(sum_df, digits = 4)
cat("\nWrote:\n ", rep_path, "\n ", sum_path, "\n")
