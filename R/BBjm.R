#' One-stage joint model: Beta-Binomial mixed + Weibull survival (TMB)
#'
#' Fits a shared-parameter joint model with Laplace approximation:
#' \deqn{y_{ij}\mid u_i \sim \mathrm{BB}(m,p_{ij},\phi)}
#' \deqn{\mathrm{logit}(p_{ij})=(\beta_0+u_{i0})+(\beta_1+u_{i1})t_{ij}}
#' \deqn{h_i(t)=w_1 w_2 t^{w_2-1}\exp\{\alpha\, m_i(t)\}}
#' with \eqn{m_i(t)=m\,\mathrm{logit}^{-1}((\beta_0+u_{i0})+(\beta_1+u_{i1})t)}
#' and \eqn{u_{i0}\sim N(0,\sigma_0^2)}, \eqn{u_{i1}\sim N(0,\sigma_1^2)}.
#'
#' @param long Long-format data frame with columns `id`, `time`, `y`.
#' @param surv Survival data frame with one row per `id`, columns
#'   `id`, `time`, `status` (1 = event, 0 = censored). IDs must match `long`.
#' @param m Maximum BB score (scalar).
#' @param n_quad Quadrature points for cumulative hazard (default 40).
#' @param maxiter `nlminb` iterations.
#' @param silent Suppress TMB tracing.
#' @param control Extra `nlminb` control.
#' @param start Optional named list of starting values.
#' @return Object of class `BBjm`.
#' @export
BBjm <- function(long, surv, m, n_quad = 40L, maxiter = 200,
                 silent = TRUE, control = list(), start = list()) {
  need <- c("id", "time", "y")
  if (!all(need %in% names(long))) {
    stop("long must have columns: id, time, y", call. = FALSE)
  }
  if (!all(c("id", "time", "status") %in% names(surv))) {
    stop("surv must have columns: id, time, status", call. = FALSE)
  }
  m <- as.numeric(m)
  if (length(m) != 1L || m <= 0 || m != as.integer(m)) {
    stop("m must be a positive integer scalar", call. = FALSE)
  }

  long <- as.data.frame(long)
  surv <- as.data.frame(surv)
  long$id <- factor(long$id)
  surv$id <- factor(surv$id, levels = levels(long$id))
  if (any(is.na(surv$id))) {
    stop("surv$id must be a subset of long$id levels", call. = FALSE)
  }
  # Keep subjects that appear in both; order by factor levels present in surv
  ids <- levels(droplevels(surv$id))
  long <- long[long$id %in% ids, , drop = FALSE]
  long$id <- factor(long$id, levels = ids)
  surv <- surv[match(ids, as.character(surv$id)), , drop = FALSE]
  n_subj <- length(ids)

  if (any(long$y < 0 | long$y > m) || any(long$y != as.integer(long$y))) {
    stop("long$y must be integers in 0..m", call. = FALSE)
  }
  if (any(surv$time <= 0)) stop("surv$time must be positive", call. = FALSE)
  if (!all(surv$status %in% c(0L, 1L))) {
    stop("surv$status must be 0/1", call. = FALSE)
  }

  ensure_tmb_dll("bb_jm")

  # Starting values from two-stage-ish warm start: BBmm then rough alpha
  zz0 <- model.matrix(~ long$id - 1)
  zz1 <- zz0 * as.numeric(long$time)
  Z <- cbind(zz0, zz1)
  X <- model.matrix(~ time, data = long)
  bb0 <- tryCatch(
    BBmm(X = X, y = as.numeric(long$y), Z = Z,
         nRandComp = c(n_subj, n_subj), m = m, silent = TRUE),
    error = function(e) NULL
  )

  if (!is.null(bb0) && identical(bb0$conv, "yes")) {
    beta0 <- as.numeric(bb0$fixed.coef)
    if (length(beta0) < 2L) beta0 <- c(beta0, 0)
    phi0 <- as.numeric(bb0$phi.coef)[1]
    sig0 <- as.numeric(bb0$sigma.coef)
    if (length(sig0) < 2L) sig0 <- c(sig0, 0.1)
    u_all <- as.numeric(bb0$random.coef)
    u0_s <- u_all[seq_len(n_subj)]
    u1_s <- u_all[n_subj + seq_len(n_subj)]
  } else {
    beta0 <- c(qlogis(mean(long$y) / m), 0)
    phi0 <- 0.5
    sig0 <- c(0.5, 0.2)
    u0_s <- rep(0, n_subj)
    u1_s <- rep(0, n_subj)
  }

  parameters <- list(
    beta = beta0[1:2],
    log_phi = log(max(phi0, 1e-3)),
    log_sigma = log(pmax(sig0[1:2], 1e-3)),
    alpha = 0,
    log_w1 = log(0.1),
    log_w2 = log(1.6),
    u0 = u0_s,
    u1 = u1_s
  )
  # user overrides
  for (nm in names(start)) {
    if (nm %in% names(parameters)) parameters[[nm]] <- start[[nm]]
  }

  data_tmb <- list(
    y = as.numeric(long$y),
    time_long = as.numeric(long$time),
    m_long = rep(m, nrow(long)),
    subj_long = as.integer(long$id) - 1L,
    surv_time = as.numeric(surv$time),
    status = as.integer(surv$status),
    m = m,
    n_quad = as.integer(n_quad)
  )

  obj <- TMB::MakeADFun(
    data = data_tmb,
    parameters = parameters,
    random = c("u0", "u1"),
    DLL = "bb_jm",
    silent = silent
  )

  ctrl <- modifyList(list(iter.max = maxiter, eval.max = maxiter * 2L), control)
  opt <- tryCatch(
    stats::nlminb(obj$par, obj$fn, obj$gr, control = ctrl),
    error = function(e) e
  )
  if (inherits(opt, "error")) {
    return(list(conv = "no", message = conditionMessage(opt)))
  }

  conv <- if (opt$convergence == 0) "yes" else "no"
  sdr <- tryCatch(TMB::sdreport(obj), error = function(e) NULL)

  pf <- opt$par
  nm <- names(pf)
  beta <- unname(pf[grep("^beta", nm)])
  names(beta) <- c("(Intercept)", "time")[seq_along(beta)]
  log_phi <- unname(pf["log_phi"])
  log_sigma <- unname(pf[grep("^log_sigma", nm)])
  alpha <- unname(pf["alpha"])
  w1 <- exp(unname(pf["log_w1"]))
  w2 <- exp(unname(pf["log_w2"]))
  phi <- exp(log_phi)
  sigma <- exp(log_sigma)
  names(sigma) <- c("RI", "RS")

  # Random effects at mode
  lp <- obj$env$last.par.best
  u0_hat <- as.numeric(lp[names(lp) == "u0"])
  u1_hat <- as.numeric(lp[names(lp) == "u1"])
  names(u0_hat) <- names(u1_hat) <- ids

  # Fixed-effect SE block
  fixed.vcov <- NULL
  se <- setNames(rep(NA_real_, length(pf)), nm)
  if (!is.null(sdr)) {
    fixed.vcov <- as.matrix(sdr$cov.fixed)
    se_fixed <- sqrt(diag(fixed.vcov))
    names(se_fixed) <- names(sdr$par.fixed)
    se[names(se_fixed)] <- se_fixed
  }

  alpha_se <- unname(se["alpha"])
  out <- list(
    fixed.coef = beta,
    phi = phi,
    sigma = sigma,
    alpha = alpha,
    alpha.se = alpha_se,
    w1 = w1,
    w2 = w2,
    random.intercept = u0_hat,
    random.slope = u1_hat,
    fixed.vcov = fixed.vcov,
    sdreport = sdr,
    conv = conv,
    nll = opt$objective,
    iter = as.integer(opt$iterations),
    nObs = nrow(long),
    nSubj = n_subj,
    nEvent = sum(surv$status),
    m = m,
    n_quad = as.integer(n_quad),
    opt = opt,
    obj = obj,
    call = match.call()
  )
  class(out) <- "BBjm"
  out
}

#' @export
print.BBjm <- function(x, ...) {
  cat("Call:\t")
  print(x$call)
  cat("\nOne-stage BB–Weibull joint model (TMB Laplace)\n")
  cat("Longitudinal fixed effects:\n")
  print(x$fixed.coef)
  cat("\nphi:", x$phi, "\n")
  cat("sigma (RI, RS):", paste(round(x$sigma, 4), collapse = ", "), "\n")
  cat("\nAssociation alpha:", x$alpha)
  if (is.finite(x$alpha.se)) cat(" (SE", round(x$alpha.se, 4), ")")
  cat("\nWeibull (w1, w2):", round(x$w1, 4), ",", round(x$w2, 4), "\n")
  cat("n subjects:", x$nSubj, "| n obs:", x$nObs, "| events:", x$nEvent, "\n")
  cat("Laplace nll:", round(x$nll, 3), "| conv:", x$conv,
      "| iter:", x$iter, "\n")
  invisible(x)
}

#' @export
summary.BBjm <- function(object, ...) {
  beta.se <- if (!is.null(object$fixed.vcov)) {
    nm <- names(object$sdreport$par.fixed)
    b_idx <- grep("^beta", nm)
    sqrt(diag(object$fixed.vcov)[b_idx])
  } else {
    rep(NA_real_, length(object$fixed.coef))
  }
  beta.tab <- cbind(
    Estimate = object$fixed.coef,
    StdErr = beta.se,
    z = object$fixed.coef / beta.se,
    p.value = 2 * pnorm(-abs(object$fixed.coef / beta.se))
  )

  alpha.tab <- cbind(
    Estimate = object$alpha,
    StdErr = object$alpha.se,
    z = object$alpha / object$alpha.se,
    p.value = 2 * pnorm(-abs(object$alpha / object$alpha.se))
  )
  rownames(alpha.tab) <- "alpha"

  out <- list(
    call = object$call,
    coefficients = beta.tab,
    alpha = alpha.tab,
    phi = object$phi,
    sigma = object$sigma,
    w1 = object$w1,
    w2 = object$w2,
    nll = object$nll,
    conv = object$conv,
    nSubj = object$nSubj,
    nObs = object$nObs,
    nEvent = object$nEvent
  )
  class(out) <- "summary.BBjm"
  out
}

#' @export
print.summary.BBjm <- function(x, ...) {
  cat("Call:\t")
  print(x$call)
  cat("\nLongitudinal fixed effects:\n\n")
  printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
  cat("\nAssociation (current-value) parameter:\n\n")
  printCoefmat(x$alpha, P.values = TRUE, has.Pvalue = TRUE)
  cat("\nphi:", x$phi, " | sigma RI,RS:", paste(round(x$sigma, 4), collapse = ", "))
  cat("\nWeibull w1,w2:", round(x$w1, 4), ",", round(x$w2, 4))
  cat("\nnll:", round(x$nll, 3), "| conv:", x$conv, "\n")
  invisible(x)
}

#' @export
coef.BBjm <- function(object, ...) {
  c(object$fixed.coef, alpha = object$alpha, phi = object$phi,
    sigma_RI = object$sigma[1], sigma_RS = object$sigma[2],
    w1 = object$w1, w2 = object$w2)
}
