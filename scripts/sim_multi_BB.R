# Multidimensional BBMM simulation (Najera-Zuloaga et al., Statistical Modelling 2024)
# Shared random-effects joint analysis of L=3 BB outcomes.
# Compare PROreg::BBmm vs PROregTMB::BBmm (nDim=3).
#
# Paper design (Sec. 3.1 / Table 1):
#   n=100, nsim=100, x~U(0,1), u~N(0, sigma_u^2)
#   Dim1: gamma=(1,-1.5),  phi=0.5, m=20
#   Dim2: gamma=(-2,2.5),  phi=1.0, m=10
#   Dim3: gamma=(0.5,0.5), phi=1.5, m=4
#   sigma_u in {0, 0.5, 1, 1.5}
#
# Usage:
#   Rscript scripts/sim_multi_BB.R --quick --nsim=30
#   Rscript scripts/sim_multi_BB.R --nsim=100 --sigma=1

suppressPackageStartupMessages({
  library(PROreg)
  library(Matrix)
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

NSIM <- as.integer(get_arg("nsim", "50"))
N <- as.integer(get_arg("N", "100"))
QUICK <- has_flag("quick")
ONLY_SIGMA <- get_arg("sigma", NULL)
if (QUICK && is.null(get_arg("nsim", NULL))) NSIM <- 30L

devtools::load_all(pkg_root, quiet = TRUE)

# ---- Truth (Table 1) ---------------------------------------------------------

truth <- list(
  gamma0 = c(1, -2, 0.5),
  gamma1 = c(-1.5, 2.5, 0.5),
  phi = c(0.5, 1.0, 1.5),
  m = c(20L, 10L, 4L)
)
L <- 3L
param_names <- c(
  paste0("g0_", 1:L), paste0("g1_", 1:L),
  paste0("phi_", 1:L), "sigma"
)

# ---- Simulate one replicate --------------------------------------------------

simulate_multi <- function(n, sigma_u, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x <- runif(n)
  u <- if (sigma_u <= 0) rep(0, n) else rnorm(n, 0, sigma_u)
  y_list <- vector("list", L)
  for (l in seq_len(L)) {
    eta <- truth$gamma0[l] + truth$gamma1[l] * x + u
    p <- 1 / (1 + exp(-eta))
    y_list[[l]] <- PROregTMB::rBB(n, truth$m[l], p, truth$phi[l])
  }
  list(x = x, y = y_list, u = u, n = n, sigma_u = sigma_u)
}

# Design matrices for multidimensional BBmm
build_multi_design <- function(sim) {
  n <- sim$n
  x <- sim$x
  Xl <- lapply(seq_len(L), function(l) cbind(1, x))
  X <- as.matrix(bdiag(Xl))
  Z <- kronecker(matrix(1, L, 1), diag(n))
  y <- unlist(sim$y)
  m <- unlist(lapply(seq_len(L), function(l) rep(truth$m[l], n)))
  list(X = X, y = y, Z = Z, m = m, nRandComp = n, nDim = L)
}

# ---- Fit engines -------------------------------------------------------------

fit_proreg <- function(des) {
  t0 <- proc.time()[["elapsed"]]
  zznull <- file(nullfile(), open = "wt")
  sink(zznull)
  sink(zznull, type = "message")
  fit <- tryCatch(
    PROreg::BBmm(
      X = des$X, y = des$y, Z = des$Z,
      nRandComp = des$nRandComp, m = des$m, nDim = des$nDim,
      maxiter = 50
    ),
    error = function(e) e
  )
  sink(type = "message")
  sink()
  close(zznull)
  elapsed <- proc.time()[["elapsed"]] - t0
  list(fit = fit, time = elapsed)
}

fit_tmb <- function(des) {
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    PROregTMB::BBmm(
      X = des$X, y = des$y, Z = des$Z,
      nRandComp = des$nRandComp, m = des$m, nDim = des$nDim,
      maxiter = 100, silent = TRUE
    ),
    error = function(e) e
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  list(fit = fit, time = elapsed)
}

extract_est <- function(fit) {
  if (inherits(fit, "error") || is.null(fit$conv) || fit$conv != "yes") {
    return(NULL)
  }
  beta <- as.numeric(fit$fixed.coef)
  if (length(beta) != 2 * L) return(NULL)
  # Order from bdiag(X1,X2,X3): (g01,g11, g02,g12, g03,g13)
  g0 <- beta[c(1, 3, 5)]
  g1 <- beta[c(2, 4, 6)]
  phi <- as.numeric(fit$phi.coef)
  if (length(phi) == 1L) phi <- rep(phi, L)
  if (length(phi) < L) return(NULL)
  sigma <- as.numeric(fit$sigma.coef)[1]
  # SE for fixed effects
  vcov <- fit$fixed.vcov
  se_beta <- if (!is.null(vcov) && all(dim(vcov) == c(6, 6))) {
    sqrt(diag(vcov))
  } else {
    rep(NA_real_, 6)
  }
  se_g0 <- se_beta[c(1, 3, 5)]
  se_g1 <- se_beta[c(2, 4, 6)]
  list(
    g0 = g0, g1 = g1, phi = phi[1:L], sigma = sigma,
    se_g0 = se_g0, se_g1 = se_g1
  )
}

# ---- Metrics -----------------------------------------------------------------

summarize_param <- function(est, se, true, nsim) {
  ok <- which(is.finite(est))
  if (!length(ok)) {
    return(c(n_ok = 0, mean = NA, bias = NA, pct_bias = NA,
             RMSE = NA, ESD = NA, ASD = NA, CP95 = NA))
  }
  e <- est[ok]
  s <- se[ok]
  cover <- if (all(is.finite(s))) {
    mean(e - 1.96 * s <= true & e + 1.96 * s >= true)
  } else NA_real_
  c(
    n_ok = length(ok),
    mean = mean(e),
    bias = mean(e) - true,
    pct_bias = (mean(e) - true) / ifelse(abs(true) < 1e-8, NA_real_, true),
    RMSE = sqrt(mean((e - true)^2)),
    ESD = stats::sd(e),
    ASD = if (all(is.finite(s))) mean(s) else NA_real_,
    CP95 = if (is.finite(cover)) 100 * cover else NA_real_
  )
}

# ---- Grid --------------------------------------------------------------------

sigmas <- c(0.5, 1.0, 1.5)
if (QUICK) sigmas <- c(0.5, 1.0)
if (!is.null(ONLY_SIGMA)) sigmas <- as.numeric(ONLY_SIGMA)

out_dir <- file.path(pkg_root, "scripts", "bench_out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
rep_path <- file.path(out_dir, "multi_BB_reps.csv")
sum_path <- file.path(out_dir, "multi_BB_summary.csv")

cat("Multidimensional BBMM: PROreg vs PROregTMB\n")
cat("Najera-Zuloaga et al. Stat Modelling 2024, Table 1\n")
cat("nsim =", NSIM, "| n =", N, "| sigma_u =", paste(sigmas, collapse = ", "), "\n")
flush.console()

cat("Warm-up TMB...\n")
ensure_tmb_dll("bb_mm")
wu <- simulate_multi(40, 0.5, seed = 0)
invisible(fit_tmb(build_multi_design(wu)))

rep_rows <- list()
sum_rows <- list()

for (su in sigmas) {
  cat("\n==== sigma_u =", su, "====\n")
  flush.console()

  # store estimates: nsim x params
  store <- function() {
    list(
      g0 = matrix(NA_real_, NSIM, L),
      g1 = matrix(NA_real_, NSIM, L),
      phi = matrix(NA_real_, NSIM, L),
      sigma = rep(NA_real_, NSIM),
      se_g0 = matrix(NA_real_, NSIM, L),
      se_g1 = matrix(NA_real_, NSIM, L),
      time = rep(NA_real_, NSIM),
      ok = rep(FALSE, NSIM)
    )
  }
  P <- store()
  Tm <- store()

  for (s in seq_len(NSIM)) {
    sim <- simulate_multi(N, su, seed = 10000L + as.integer(su * 100) + s)
    des <- build_multi_design(sim)

    rp <- fit_proreg(des)
    rt <- fit_tmb(des)
    ep <- extract_est(rp$fit)
    et <- extract_est(rt$fit)

    P$time[s] <- rp$time
    Tm$time[s] <- rt$time

    if (!is.null(ep)) {
      P$ok[s] <- TRUE
      P$g0[s, ] <- ep$g0
      P$g1[s, ] <- ep$g1
      P$phi[s, ] <- ep$phi
      P$sigma[s] <- ep$sigma
      P$se_g0[s, ] <- ep$se_g0
      P$se_g1[s, ] <- ep$se_g1
    }
    if (!is.null(et)) {
      Tm$ok[s] <- TRUE
      Tm$g0[s, ] <- et$g0
      Tm$g1[s, ] <- et$g1
      Tm$phi[s, ] <- et$phi
      Tm$sigma[s] <- et$sigma
      Tm$se_g0[s, ] <- et$se_g0
      Tm$se_g1[s, ] <- et$se_g1
    }

    # replication row (slopes + sigma + times)
    rep_rows[[length(rep_rows) + 1]] <- data.frame(
      sigma_u = su, sim = s, n = N,
      ok_proreg = P$ok[s], ok_tmb = Tm$ok[s],
      g1_1_P = P$g1[s, 1], g1_2_P = P$g1[s, 2], g1_3_P = P$g1[s, 3],
      g1_1_T = Tm$g1[s, 1], g1_2_T = Tm$g1[s, 2], g1_3_T = Tm$g1[s, 3],
      sigma_P = P$sigma[s], sigma_T = Tm$sigma[s],
      time_P = P$time[s], time_T = Tm$time[s],
      stringsAsFactors = FALSE
    )
    utils::write.csv(do.call(rbind, rep_rows), rep_path, row.names = FALSE)

    if (s %% 5 == 0 || s == NSIM) {
      cat(sprintf(
        "  sim %d/%d | conv P:%d T:%d | tP=%.2fs tT=%.2fs\n",
        s, NSIM, sum(P$ok), sum(Tm$ok),
        rp$time, rt$time
      ))
      flush.console()
    }
  }

  # summary for each method / parameter
  add_sum <- function(engine, st, which_ok) {
    for (l in seq_len(L)) {
      sg0 <- summarize_param(st$g0[, l], st$se_g0[, l], truth$gamma0[l], NSIM)
      sg1 <- summarize_param(st$g1[, l], st$se_g1[, l], truth$gamma1[l], NSIM)
      sph <- summarize_param(st$phi[, l], rep(NA_real_, NSIM), truth$phi[l], NSIM)
      sum_rows[[length(sum_rows) + 1]] <<- data.frame(
        sigma_u = su, engine = engine, param = paste0("g0_", l),
        true = truth$gamma0[l], t(sg0),
        time_mean = mean(st$time, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      sum_rows[[length(sum_rows) + 1]] <<- data.frame(
        sigma_u = su, engine = engine, param = paste0("g1_", l),
        true = truth$gamma1[l], t(sg1),
        time_mean = mean(st$time, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      sum_rows[[length(sum_rows) + 1]] <<- data.frame(
        sigma_u = su, engine = engine, param = paste0("phi_", l),
        true = truth$phi[l], t(sph),
        time_mean = mean(st$time, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
    ss <- summarize_param(st$sigma, rep(NA_real_, NSIM), su, NSIM)
    sum_rows[[length(sum_rows) + 1]] <<- data.frame(
      sigma_u = su, engine = engine, param = "sigma",
      true = su, t(ss),
      time_mean = mean(st$time, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  add_sum("PROreg", P)
  add_sum("PROregTMB", Tm)

  utils::write.csv(do.call(rbind, sum_rows), sum_path, row.names = FALSE)

  cat(sprintf(
    "  mean time PROreg=%.2fs TMB=%.2fs speedup=%.1fx | conv P=%.0f%% T=%.0f%%\n",
    mean(P$time, na.rm = TRUE), mean(Tm$time, na.rm = TRUE),
    mean(P$time, na.rm = TRUE) / max(mean(Tm$time, na.rm = TRUE), 1e-8),
    100 * mean(P$ok), 100 * mean(Tm$ok)
  ))
  # slope bias snapshot
  for (l in 1:L) {
    cat(sprintf(
      "  g1_%d true=%.2f | bias P=%+.3f T=%+.3f | CP P=%.0f T=%.0f\n",
      l, truth$gamma1[l],
      mean(P$g1[P$ok, l], na.rm = TRUE) - truth$gamma1[l],
      mean(Tm$g1[Tm$ok, l], na.rm = TRUE) - truth$gamma1[l],
      100 * mean(P$g1[P$ok, l] - 1.96 * P$se_g1[P$ok, l] <= truth$gamma1[l] &
                   P$g1[P$ok, l] + 1.96 * P$se_g1[P$ok, l] >= truth$gamma1[l],
                 na.rm = TRUE),
      100 * mean(Tm$g1[Tm$ok, l] - 1.96 * Tm$se_g1[Tm$ok, l] <= truth$gamma1[l] &
                   Tm$g1[Tm$ok, l] + 1.96 * Tm$se_g1[Tm$ok, l] >= truth$gamma1[l],
                 na.rm = TRUE)
    ))
  }
}

sum_df <- do.call(rbind, sum_rows)
cat("\n=== Summary (slopes + sigma) ===\n")
print(sum_df[sum_df$param %in% c("g1_1", "g1_2", "g1_3", "sigma"),
             c("sigma_u", "engine", "param", "true", "bias", "RMSE",
               "ESD", "CP95", "time_mean")], digits = 3)
cat("\nWrote:\n ", rep_path, "\n ", sum_path, "\n")
