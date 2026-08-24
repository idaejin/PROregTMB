#!/usr/bin/env Rscript
# Shared helpers for SMMR WW-BW multivariate BB simulation.
# Design: 01_PROJECTS/BB/ref/WWBW-SGRQ/SIMULATION-SMMR-DESIGN.md

`%||%` <- function(a, b) if (!is.null(a)) a else b

.sim_root <- function() {
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in seq_len(8L)) {
    if (file.exists(file.path(d, "DESCRIPTION")) &&
        file.exists(file.path(d, "R", "BBmm.R"))) return(d)
    nd <- dirname(d)
    if (identical(nd, d)) break
    d <- nd
  }
  d
}

load_proregtmb <- function() {
  root <- .sim_root()
  if (!requireNamespace("PROregTMB", quietly = TRUE)) {
    suppressMessages(devtools::load_all(root, quiet = TRUE))
  } else {
    suppressPackageStartupMessages(library(PROregTMB))
  }
  invisible(root)
}

default_truth <- function() {
  list(
    domains = c("Impacts", "Symptoms", "Activity"),
    alpha = c(-0.85, -0.40, -0.20),
    beta_B_fev = c(-0.20, -0.10, -0.30),
    beta_B_wt  = c(-0.18, -0.10, -0.35),
    # within slopes on UNIT within-SD scale (X^{w,z})
    beta_W_fev = c(-0.15, -0.10, -0.20),
    beta_W_wt  = c(-0.08, -0.05, -0.18),
    gamma_fev = c(0.04, 0.03, 0.05),
    gamma_wt  = c(0.03, 0.02, 0.04),
    phi = c(0.022, 0.074, 0.041),
    m = 24L,
    G = matrix(c(
      1.00, 0.80, 0.90,
      0.80, 0.80, 0.70,
      0.90, 0.70, 1.00
    ), 3L, 3L, byrow = TRUE),
    # tau on unit within-SD scale (comparable across FEV/WT)
    tau_fev = 0.12,
    tau_wt  = 0.06,
    sigma_w = 0.45
  )
}

rmvnorm_chol <- function(n, Sigma) {
  R <- chol(Sigma)
  matrix(rnorm(n * ncol(Sigma)), n, ncol(Sigma)) %*% R
}

#' Marginal AIC parameter count (fixed + phi + G + optional smooth_sd + RS).
count_q <- function(model_code, fit_obj = NULL) {
  # Unstructured 3x3 G: 3 SD + 3 corr
  q_G <- 6L
  q_phi <- 3L
  q_fix <- switch(model_code,
    M1 = 3L * 3L,           # intercept + fev + wt per domain
    M2 = 3L * 5L,           # intercept + 2 between + 2 within
    M3 = 3L * 5L,           # intercept + 2 between + 2 within nulls per domain
    M4 = 3L * 5L,
    M5 = 3L * 5L,           # M3 mean + RS (within nulls domain-specific)
    3L * 5L
  )
  # Prefer empirical fixed length when available
  if (!is.null(fit_obj) && isTRUE(fit_obj$ok)) {
    f <- fit_obj$fit
    if (!is.null(f$fixed.coef)) q_fix <- length(f$fixed.coef)
  }
  q_smooth <- if (model_code %in% c("M3", "M5")) 6L else 0L
  q_rs <- if (model_code %in% c("M4", "M5")) 3L else 0L  # 2 SD + 1 corr
  as.integer(q_fix + q_phi + q_G + q_smooth + q_rs)
}

aic_marg <- function(nll, q) 2 * as.numeric(nll) + 2 * as.numeric(q)

generate_panel <- function(scenario = c("A", "B", "C", "D"),
                           N = 500L, T = 4L,
                           truth = default_truth(),
                           seed = NULL) {
  scenario <- match.arg(scenario)
  if (!is.null(seed)) set.seed(seed)
  tr <- truth
  L <- length(tr$domains)
  stopifnot(all(eigen(tr$G, symmetric = TRUE, only.values = TRUE)$values > 1e-8))

  use_rs <- scenario %in% c("A", "C")
  use_quad <- scenario %in% c("B", "C")
  if (identical(scenario, "D")) {
    tr$beta_B_fev <- c(-0.30, -0.20, -0.40)
    tr$beta_B_wt  <- c(-0.28, -0.18, -0.38)
    tr$beta_W_fev <- c(-0.10, -0.08, -0.12)
    tr$beta_W_wt  <- c(-0.06, -0.04, -0.10)
    use_rs <- FALSE
    use_quad <- FALSE
  }

  id <- rep(seq_len(N), each = T)
  time <- rep(seq_len(T), times = N)

  fev_bar <- rnorm(N)
  wt_bar <- rnorm(N)
  fev_w_raw <- rnorm(N * T, 0, tr$sigma_w)
  wt_w_raw <- rnorm(N * T, 0, tr$sigma_w)
  # unit within-SD scale used in the mean and in fitted competitors
  fev_w <- fev_w_raw / tr$sigma_w
  wt_w <- wt_w_raw / tr$sigma_w
  fev <- fev_bar[id] + fev_w_raw
  wt <- wt_bar[id] + wt_w_raw

  a <- rmvnorm_chol(N, tr$G)
  colnames(a) <- tr$domains
  if (use_rs) {
    b_fev <- rnorm(N, 0, tr$tau_fev)
    b_wt <- rnorm(N, 0, tr$tau_wt)
  } else {
    b_fev <- rep(0, N)
    b_wt <- rep(0, N)
  }

  rows <- vector("list", L)
  for (ell in seq_len(L)) {
    eta <- tr$alpha[ell] +
      tr$beta_B_fev[ell] * fev_bar[id] +
      tr$beta_B_wt[ell] * wt_bar[id] +
      tr$beta_W_fev[ell] * fev_w +
      tr$beta_W_wt[ell] * wt_w +
      b_fev[id] * fev_w +
      b_wt[id] * wt_w +
      a[id, ell]
    if (use_quad) {
      eta <- eta + tr$gamma_fev[ell] * fev_w^2 + tr$gamma_wt[ell] * wt_w^2
    }
    p <- plogis(eta)
    y <- PROregTMB::rBB(length(p), tr$m, p, tr$phi[ell])
    rows[[ell]] <- data.frame(
      id = factor(id), time = time,
      dim = factor(tr$domains[ell], levels = tr$domains),
      y = as.integer(y), m = as.integer(tr$m),
      fev = fev, wt = wt,
      fev_bar = fev_bar[id], wt_bar = wt_bar[id],
      fev_w = fev_w, wt_w = wt_w,
      stringsAsFactors = FALSE
    )
  }
  long <- do.call(rbind, rows)
  rownames(long) <- NULL

  list(
    long = long, scenario = scenario,
    N = as.integer(N), T = as.integer(T), L = L,
    truth = tr, use_rs = use_rs, use_quad = use_quad,
    a = a, b_fev = b_fev, b_wt = b_wt
  )
}

fit_competitors <- function(panel,
                            models = c("M1", "M2", "M3", "M4"),
                            maxiter = 200L,
                            k_smooth = 8L) {
  long <- panel$long
  out <- list()

  run_one <- function(code, label, form, random) {
    t0 <- proc.time()[["elapsed"]]
    fit <- tryCatch(
      BBmm(form, random = random, dim = "dim", corr = "unstructured",
           m = "m", data = long, maxiter = maxiter),
      error = function(e) e
    )
    elapsed <- proc.time()[["elapsed"]] - t0
    if (inherits(fit, "error")) {
      return(list(code = code, label = label, ok = FALSE,
                  error = conditionMessage(fit), elapsed = elapsed,
                  q = count_q(code), nll = NA_real_, aic = NA_real_))
    }
    nll <- as.numeric(fit$nll %||% NA_real_)
    q <- count_q(code, list(ok = TRUE, fit = fit))
    list(
      code = code, label = label, ok = TRUE,
      conv = as.character(fit$conv %||% NA_character_),
      nll = nll, q = q, aic = aic_marg(nll, q),
      elapsed = elapsed,
      beta = fit$fixed.coef, fixed.vcov = fit$fixed.vcov,
      phi = as.numeric(fit$phi.coef %||% fit$phi),
      sigma = as.numeric(fit$sigma.coef %||% fit$sigma),
      Corr = fit$Corr,
      smooth_sd = fit$smooth_sd, lambda = fit$lambda,
      fit = fit
    )
  }

  if ("M1" %in% models)
    out$M1 <- run_one("M1", "M1_contemporaneous", y ~ fev + wt,
                      ~ (0 + dim | id))
  if ("M2" %in% models)
    out$M2 <- run_one("M2", "M2_WWBW_linear",
                      y ~ fev_bar + wt_bar + fev_w + wt_w, ~ (0 + dim | id))
  if ("M3" %in% models)
    out$M3 <- run_one("M3", "M3_WWBW_spline",
                      as.formula(sprintf(
                        "y ~ fev_bar + wt_bar + s(fev_w, by = dim, k = %d) + s(wt_w, by = dim, k = %d)",
                        k_smooth, k_smooth)),
                      ~ (0 + dim | id))
  if ("M4" %in% models)
    out$M4 <- run_one("M4", "M4_WWBW_RS",
                      y ~ fev_bar + wt_bar + fev_w + wt_w,
                      ~ (0 + dim | id) + (0 + fev_w + wt_w | id))
  if ("M5" %in% models)
    out$M5 <- run_one("M5", "M5_WWBW_spline_RS",
                      as.formula(sprintf(
                        "y ~ fev_bar + wt_bar + s(fev_w, by = dim, k = %d) + s(wt_w, by = dim, k = %d)",
                        k_smooth, k_smooth)),
                      ~ (0 + dim | id) + (0 + fev_w + wt_w | id))
  out
}

coef_by_domain <- function(beta, domains, term) {
  nm <- names(beta)
  sapply(domains, function(d) {
    hit <- which(nm == paste0(d, ".", term))
    if (!length(hit)) NA_real_ else as.numeric(beta[hit[1]])
  })
}

se_by_domain <- function(vcov, beta, domains, term) {
  if (is.null(vcov) || is.null(beta)) return(rep(NA_real_, length(domains)))
  nm <- names(beta)
  sapply(domains, function(d) {
    j <- which(nm == paste0(d, ".", term))
    if (!length(j)) return(NA_real_)
    sqrt(max(as.numeric(vcov[j[1], j[1]]), 0))
  })
}

#' Discrete curvature proxy: mean abs second difference of predicted within effect.
curvature_proxy_linear_fit <- function(res, channel = c("fev_w", "wt_w"),
                                       grid = seq(-2, 2, length.out = 41)) {
  channel <- match.arg(channel)
  if (!isTRUE(res$ok) || is.null(res$beta)) return(NA_real_)
  # For M2/M4: linear slope -> second diff ~ 0
  # For M3: use domain-average of linear null coef only as crude; better use smooth
  # Prefer evaluating from fitted object if smooth present
  f <- res$fit
  if (!is.null(f$smooth) && !is.null(f$s)) {
    # fall back: use max |s| / smooth_sd as roughness score if predict unavailable
    return(mean(abs(as.numeric(f$s)), na.rm = TRUE))
  }
  0
}

extract_G_hat <- function(res) {
  if (!isTRUE(res$ok) || is.null(res$Corr) || is.null(res$sigma)) {
    return(list(G = NULL, corr_rmse = NA_real_, sd_rmse = NA_real_))
  }
  # first Corr block is domain intercepts
  C <- res$Corr
  if (is.list(C)) C <- C[[1]]
  sd <- as.numeric(res$sigma)[1:3]
  if (length(sd) < 3 || !is.matrix(C)) {
    return(list(G = NULL, corr_rmse = NA_real_, sd_rmse = NA_real_))
  }
  G <- diag(sd) %*% C[1:3, 1:3, drop = FALSE] %*% diag(sd)
  list(G = G, C = C[1:3, 1:3, drop = FALSE], sd = sd)
}

G_metrics <- function(res, G_true) {
  eg <- extract_G_hat(res)
  if (is.null(eg$G)) {
    return(c(corr_rmse = NA_real_, sd_rmse = NA_real_, frob = NA_real_))
  }
  C_true <- cov2cor(G_true)
  sd_true <- sqrt(diag(G_true))
  c(
    corr_rmse = sqrt(mean((eg$C[upper.tri(eg$C)] - C_true[upper.tri(C_true)])^2)),
    sd_rmse = sqrt(mean((eg$sd - sd_true)^2)),
    frob = sqrt(sum((eg$G - G_true)^2))
  )
}

extract_tau <- function(res) {
  if (!isTRUE(res$ok) || is.null(res$sigma)) {
    return(c(tau_fev = NA_real_, tau_wt = NA_real_))
  }
  # M4: sigma = (3 domain SD, tau_fev, tau_wt)
  s <- as.numeric(res$sigma)
  if (length(s) >= 5) c(tau_fev = s[4], tau_wt = s[5])
  else c(tau_fev = NA_real_, tau_wt = NA_real_)
}

#' Full replicate summary for simulation metrics.
summarise_replicate <- function(panel, fits) {
  tr <- panel$truth
  domains <- tr$domains
  coef_rows <- list()
  model_rows <- list()

  for (nm in names(fits)) {
    res <- fits[[nm]]
    q <- res$q %||% count_q(nm, res)
    nll <- res$nll %||% NA_real_
    aic <- res$aic %||% aic_marg(nll, q)
    gm <- G_metrics(res, tr$G)
    tau <- extract_tau(res)
    phi <- res$phi %||% rep(NA_real_, 3)

    model_rows[[nm]] <- data.frame(
      scenario = panel$scenario, model = nm,
      label = res$label %||% nm,
      ok = isTRUE(res$ok),
      conv = res$conv %||% NA_character_,
      nll = nll, q = q, aic = aic,
      elapsed = res$elapsed %||% NA_real_,
      corr_rmse_G = unname(gm["corr_rmse"]),
      sd_rmse_G = unname(gm["sd_rmse"]),
      frob_G = unname(gm["frob"]),
      tau_fev = unname(tau["tau_fev"]),
      tau_wt = unname(tau["tau_wt"]),
      bias_tau_fev = unname(tau["tau_fev"]) - tr$tau_fev,
      bias_tau_wt = unname(tau["tau_wt"]) - tr$tau_wt,
      phi_I = phi[1], phi_S = phi[2], phi_A = phi[3],
      bias_phi_I = phi[1] - tr$phi[1],
      bias_phi_S = phi[2] - tr$phi[2],
      bias_phi_A = phi[3] - tr$phi[3],
      roughness_s = curvature_proxy_linear_fit(res),
      stringsAsFactors = FALSE
    )

    if (!isTRUE(res$ok) || is.null(res$beta)) next
    for (term_pair in list(
      list(est = "fev_bar", truth = tr$beta_B_fev, kind = "B_fev"),
      list(est = "wt_bar",  truth = tr$beta_B_wt,  kind = "B_wt"),
      list(est = "fev_w",   truth = tr$beta_W_fev, kind = "W_fev"),
      list(est = "wt_w",    truth = tr$beta_W_wt,  kind = "W_wt"),
      list(est = "fev",     truth = NA,            kind = "cont_fev"),
      list(est = "wt",      truth = NA,            kind = "cont_wt")
    )) {
      est <- coef_by_domain(res$beta, domains, term_pair$est)
      se <- se_by_domain(res$fixed.vcov, res$beta, domains, term_pair$est)
      for (i in seq_along(domains)) {
        if (!is.finite(est[i])) next
        truth_i <- if (length(term_pair$truth) == length(domains))
          term_pair$truth[i] else NA_real_
        lo <- est[i] - 1.96 * se[i]
        hi <- est[i] + 1.96 * se[i]
        cover <- if (is.finite(truth_i) && is.finite(se[i]))
          as.integer(lo <= truth_i && truth_i <= hi) else NA_integer_
        coef_rows[[length(coef_rows) + 1L]] <- data.frame(
          scenario = panel$scenario, model = nm, domain = domains[i],
          effect = term_pair$kind,
          truth = truth_i, est = est[i], se = se[i],
          bias = est[i] - truth_i, cover95 = cover,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  list(
    models = do.call(rbind, model_rows),
    coefs = if (length(coef_rows)) do.call(rbind, coef_rows) else NULL
  )
}
