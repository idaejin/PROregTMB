.validate_bb_inputs <- function(y, m, X, context = "model") {
  if (missing(y) || missing(m) || missing(X)) {
    stop(sprintf("%s: y, m and X are required.", context), call. = FALSE)
  }

  X <- as.matrix(X)
  y <- as.numeric(y)
  m <- as.numeric(m)

  if (nrow(X) != length(y)) {
    stop(sprintf("%s: nrow(X) must equal length(y).", context), call. = FALSE)
  }

  if (length(m) == 1L) {
    m <- rep(m, length(y))
  } else if (length(m) != length(y)) {
    stop(sprintf("%s: m must be scalar or length(y).", context), call. = FALSE)
  }

  if (any(!is.finite(y))) {
    stop(sprintf("%s: y must contain only finite values.", context), call. = FALSE)
  }
  if (any(!is.finite(m)) || any(m <= 0) || any(m != as.integer(m))) {
    stop(sprintf("%s: m must be positive integer(s).", context), call. = FALSE)
  }
  if (any(y != as.integer(y))) {
    stop(sprintf("%s: y must be integer-valued counts.", context), call. = FALSE)
  }
  if (any(y < 0 | y > m)) {
    stop(sprintf("%s: y must be bounded between 0 and m.", context), call. = FALSE)
  }

  list(y = y, m = m)
}

.bbreg_diagnostics <- function(opt, obj, par = NULL) {
  if (is.null(par)) {
    par <- opt$par
  }
  grad <- tryCatch(obj$gr(par), error = function(e) rep(NA_real_, length(par)))
  list(
    convergence_code = opt$convergence,
    converged = isTRUE(opt$convergence == 0L),
    iterations = as.integer(opt$iterations),
    objective = opt$objective,
    gradient_norm = if (length(grad) > 0L) sqrt(sum(as.numeric(grad)^2, na.rm = TRUE)) else NA_real_,
    max_abs_gradient = if (length(grad) > 0L) max(abs(as.numeric(grad)), na.rm = TRUE) else NA_real_
  )
}

#' Fit a beta-binomial logistic regression via TMB
#'
#' Model
#' \deqn{y_i \sim \mathrm{BB}(m_i, p_i, \phi),\quad
#' \mathrm{logit}(p_i) = x_i^\top\beta + \sum_j f_j(x_{ij})}
#' with optional additive P-splines via Eilers (1999) mixed reparameterisation
#' \eqn{f_j = X_{\mathrm{null}}\beta_{\mathrm{null}} + Z_j s_j},
#' \eqn{s_j \sim N(0,\sigma_{s_j}^2 I)}, \eqn{Z_j = B D'(DD')^{-1}}
#' residualised against the fixed design.
#'
#' @param formula Model formula for the mean (logit link). Use
#'   \code{s(x, ndx, pord)} for smooths.
#' @param m Maximum score (scalar or vector of length n).
#' @param data Optional data frame.
#' @param method `"mle"` (default): maximize the TMB likelihood with
#'   `nlminb`. `"bayes"`: same point fit, then NUTS via
#'   [tmbstan::tmbstan()] on a template with weak priors
#'   (parametric models only; no `s()` yet).
#' @param maxiter Maximum `nlminb` iterations.
#' @param control Passed to [stats::nlminb()].
#' @param silent Suppress TMB tracing.
#' @param chains,iter,warmup,seed Passed to [tmbstan::tmbstan()] when
#'   `method = "bayes"`.
#' @return Object of class `BBreg`. With `method = "bayes"`, also contains
#'   `posterior` (summary table), `stanfit`, and `time_bayes`.
#' @export
BBreg <- function(formula, m, data = list(),
                  method = c("mle", "bayes"),
                  maxiter = 100, control = list(), silent = TRUE,
                  chains = 2L, iter = 1000L, warmup = 400L, seed = 1L) {
  method <- match.arg(method)
  empty_data <- is.null(data) || (is.list(data) && !is.data.frame(data) && !length(data))
  if (!empty_data) data <- as.data.frame(data)

  sm <- .build_smooth_design(formula, data = if (empty_data) NULL else data,
                             envir = parent.frame())
  mf <- if (empty_data) {
    model.frame(formula = sm$fixed, data = parent.frame())
  } else {
    model.frame(formula = sm$fixed, data = data)
  }
  X <- model.matrix(attr(mf, "terms"), data = mf)
  y <- as.numeric(model.response(mf))
  validated <- .validate_bb_inputs(y = y, m = m, X = X, context = "BBreg")
  y <- validated$y
  n <- length(y)
  X <- .merge_X_null(X, sm$X_null)
  if (isTRUE(sm$n_smooth > 0L)) {
    sm <- .reorth_smooth_Zs(sm, X)
  }
  if (sm$n_smooth > 0L && nrow(sm$Zs) != n) {
    stop("smooth design nrow does not match response length", call. = FALSE)
  }

  if (length(m) == 1L) {
    balanced <- "yes"
    m. <- rep(as.numeric(m), n)
  } else {
    m. <- as.numeric(m)
    if (length(m.) != n) stop("m must be scalar or length(y)", call. = FALSE)
    balanced <- if (length(unique(m.)) == 1L) "yes" else "no"
  }

  validated <- .validate_bb_inputs(y = y, m = m., X = X, context = "BBreg")
  y <- validated$y
  m. <- validated$m

  if (identical(method, "bayes") && sm$n_smooth > 0L) {
    stop('method = "bayes" does not support s() smooths yet', call. = FALSE)
  }

  ensure_tmb_dll("bb_reg")

  # Starting values: binomial GLM on parametric+null part + moment phi
  glm0 <- stats::glm.fit(X, y / m., family = stats::binomial(), weights = m.)
  beta0 <- as.numeric(glm0$coefficients)
  beta0[!is.finite(beta0)] <- 0
  if (length(beta0) < ncol(X)) beta0 <- c(beta0, rep(0, ncol(X) - length(beta0)))
  if (length(beta0) > ncol(X)) beta0 <- beta0[seq_len(ncol(X))]
  mu <- mean(y)
  v <- stats::var(y)
  mbar <- mean(m.)
  phi_mm2 <- {
    ph <- (v - mu * (1 - mu / mbar)) /
      (mu * (1 - mu / mbar) * (mbar - 1) - (v - mu * (1 - mu / mbar)))
    if (!is.finite(ph) || ph <= 0) 0.5 else ph
  }

  spar <- .tmb_smooth_parameters(sm)
  parameters <- c(
    list(beta = beta0, log_phi = log(max(phi_mm2, 1e-3))),
    spar
  )
  data_tmb <- c(
    list(y = y, m = m., X = X),
    .tmb_smooth_data(sm, n = n)
  )

  random <- if (sm$n_smooth > 0L) "s" else NULL
  obj <- TMB::MakeADFun(
    data = data_tmb,
    parameters = parameters,
    random = random,
    DLL = "bb_reg",
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
  diagnostics <- .bbreg_diagnostics(opt, obj)
  sdr <- tryCatch(
    TMB::sdreport(obj, getJointPrecision = sm$n_smooth > 0L),
    error = function(e) NULL
  )

  beta <- opt$par[grep("^beta", names(opt$par))]
  if (length(beta) == 0L) beta <- opt$par[seq_len(ncol(X))]
  names(beta) <- colnames(X)
  log_phi <- unname(opt$par["log_phi"])
  phi <- exp(log_phi)

  lambda <- numeric(0)
  smooth_sd <- numeric(0)
  s_hat <- numeric(0)
  s.vcov <- NULL
  fhat <- list()
  if (sm$n_smooth > 0L) {
    sd_hat <- opt$par[grep("^log_sds", names(opt$par))]
    smooth_sd <- setNames(exp(unname(sd_hat)), sm$labels)
    lambda <- setNames(1 / (smooth_sd^2), sm$labels)
    s_hat <- tryCatch(
      as.numeric(obj$env$last.par.best[obj$env$random]),
      error = function(e) rep(NA_real_, ncol(sm$Zs))
    )
    names(s_hat) <- colnames(sm$Zs)
    for (j in seq_len(sm$n_smooth)) {
      idx <- sm$blocks_idx[[j]]
      fhat[[j]] <- as.numeric(sm$Zs[, idx, drop = FALSE] %*% s_hat[idx])
    }
    names(fhat) <- sm$labels
    s.vcov <- .s_vcov_from_joint(
      if (!is.null(sdr)) sdr$jointPrecision else NULL,
      n_s = length(s_hat),
      n_u = 0L
    )
  }

  if (!is.null(sdr)) {
    vcov.b <- as.matrix(sdr$cov.fixed)
    nm <- names(sdr$par.fixed)
    b_idx <- grep("^beta", nm)
    phi_idx <- which(nm == "log_phi")
    if (length(b_idx)) {
      vcov.b <- vcov.b[b_idx, b_idx, drop = FALSE]
      rownames(vcov.b) <- colnames(vcov.b) <- colnames(X)
    }
    psi.var <- if (length(phi_idx)) as.numeric(sdr$cov.fixed[phi_idx, phi_idx]) else NA_real_
  } else {
    H <- optimHess(opt$par, obj$fn, obj$gr)
    invH <- tryCatch(solve(H), error = function(e) matrix(NA_real_, length(opt$par), length(opt$par)))
    b_idx <- seq_len(length(beta))
    vcov.b <- invH[b_idx, b_idx, drop = FALSE]
    rownames(vcov.b) <- colnames(vcov.b) <- colnames(X)
    psi.var <- invH[length(opt$par), length(opt$par)]
  }

  eta <- as.numeric(X %*% beta)
  if (sm$n_smooth > 0L) eta <- eta + as.numeric(sm$Zs %*% s_hat)
  fitted.values <- 1 / (1 + exp(-eta))

  e <- sum(y) / sum(m.)
  loglik <- function(p_hat, phi_hat) {
    a <- p_hat / phi_hat
    b <- (1 - p_hat) / phi_hat
    sum(
      lgamma(m. + 1) - lgamma(y + 1) - lgamma(m. - y + 1) +
        lgamma(a + y) - lgamma(a) +
        lgamma(b + m. - y) - lgamma(b) +
        lgamma(a + b) - lgamma(a + b + m.)
    )
  }
  p_sat <- pmin(pmax(y / m., 1e-8), 1 - 1e-8)
  deviance <- as.numeric(2 * (loglik(p_sat, phi) - loglik(fitted.values, phi)))
  null.p <- rep(e, n)
  null.p <- pmin(pmax(null.p, 1e-8), 1 - 1e-8)
  null.deviance <- as.numeric(2 * (loglik(p_sat, phi) - loglik(null.p, phi)))

  df <- n - length(beta) - 1L - sm$n_smooth
  null.df <- n - 1L

  coef.b <- matrix(beta, ncol = 1L, dimnames = list(names(beta), NULL))

  out <- list(
    coefficients = coef.b,
    vcov = vcov.b,
    phi = phi,
    psi = log_phi,
    psi.var = psi.var,
    conv = conv,
    diagnostics = diagnostics,
    fitted.values = fitted.values,
    deviance = deviance,
    df = df,
    null.deviance = null.deviance,
    null.df = null.df,
    iter = as.integer(opt$iterations),
    X = X,
    y = y,
    m = if (balanced == "yes") m.[1] else m.,
    balanced = balanced,
    nObs = n,
    method = method,
    opt = opt,
    obj = obj,
    sdreport = sdr,
    nll = opt$objective,
    smooth = if (sm$n_smooth > 0L) sm else NULL,
    s = if (length(s_hat)) s_hat else NULL,
    s.vcov = s.vcov,
    smooth_sd = if (length(smooth_sd)) smooth_sd else NULL,
    lambda = if (length(lambda)) lambda else NULL,
    fhat = if (length(fhat)) fhat else NULL,
    posterior = NULL,
    stanfit = NULL,
    time_bayes = NA_real_
  )
  class(out) <- "BBreg"
  out$call <- match.call()
  out <- .bbreg_attach_formulation(out, formula = formula)

  if (identical(method, "bayes")) {
    # parametric-only payload for bb_reg_prior
    data_bayes <- list(y = y, m = m., X = X)
    out <- .BBreg_add_bayes(
      out, data_tmb = data_bayes,
      chains = chains, iter = iter, warmup = warmup,
      seed = seed, silent = silent
    )
  }
  out
}

# NUTS via tmbstan on bb_reg_prior (weak priors); attaches to a BBreg object
.BBreg_add_bayes <- function(fit, data_tmb, chains, iter, warmup, seed, silent) {
  if (!requireNamespace("tmbstan", quietly = TRUE) ||
      !requireNamespace("rstan", quietly = TRUE)) {
    stop('method = "bayes" requires packages tmbstan and rstan', call. = FALSE)
  }
  ensure_tmb_dll("bb_reg_prior")
  beta0 <- as.numeric(fit$coefficients)
  phi0 <- as.numeric(fit$phi)[1]
  obj <- TMB::MakeADFun(
    data = data_tmb,
    parameters = list(beta = beta0, log_phi = log(max(phi0, 1e-6))),
    DLL = "bb_reg_prior",
    silent = silent
  )
  init_fn <- function() {
    list(
      beta = beta0 + stats::rnorm(length(beta0), 0, 0.05),
      log_phi = log(max(phi0, 1e-6)) + stats::rnorm(1, 0, 0.05)
    )
  }
  t1 <- proc.time()[["elapsed"]]
  stanfit <- tmbstan::tmbstan(
    obj,
    chains = as.integer(chains),
    iter = as.integer(iter),
    warmup = as.integer(warmup),
    seed = as.integer(seed),
    refresh = 0,
    init = init_fn,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
  )
  fit$time_bayes <- proc.time()[["elapsed"]] - t1
  fit$stanfit <- stanfit

  sm <- as.data.frame(rstan::summary(stanfit)$summary)
  b_rows <- grep("^beta", rownames(sm))
  phi_draws <- exp(as.matrix(stanfit)[, "log_phi"])
  post_beta <- data.frame(
    param = rownames(fit$coefficients),
    mean = sm[b_rows, "mean"],
    sd = sm[b_rows, "sd"],
    q025 = sm[b_rows, "2.5%"],
    q975 = sm[b_rows, "97.5%"],
    n_eff = sm[b_rows, "n_eff"],
    Rhat = sm[b_rows, "Rhat"],
    stringsAsFactors = FALSE
  )
  post_phi <- data.frame(
    param = "phi",
    mean = mean(phi_draws),
    sd = stats::sd(phi_draws),
    q025 = unname(stats::quantile(phi_draws, 0.025)),
    q975 = unname(stats::quantile(phi_draws, 0.975)),
    n_eff = sm["log_phi", "n_eff"],
    Rhat = sm["log_phi", "Rhat"],
    stringsAsFactors = FALSE
  )
  fit$posterior <- rbind(post_beta, post_phi)
  fit
}

#' @export
print.BBreg <- function(x, ...) {
  cat("Call:\t")
  print(x$call)
  .cat_model_block("Model (BBreg):", if (!is.null(x$model)) x$model else .bb_model_lines_reg())
  cat("Method:", if (!is.null(x$method)) x$method else "mle")
  if (identical(x$method, "bayes")) cat(" (point estimates below are MLE)")
  cat("\n\n")
  cat("beta (mean, logit scale):\n")
  print(if (!is.null(x$beta)) x$beta else setNames(as.numeric(x$coefficients),
                                                   rownames(x$coefficients)))
  cat("\nphi (BB dispersion):", as.numeric(x$phi)[1L], "\n")
  if (!is.null(x$log_phi)) {
    cat("log(phi):", as.numeric(x$log_phi)[1L], "\n")
  }
  if (!is.null(x$smooth_sd)) {
    cat("\nsmooth_sd (Eilers P-spline RE SD):\n")
    print(x$smooth_sd)
    if (!is.null(x$lambda)) {
      cat("lambda (= 1/smooth_sd^2):\n")
      print(x$lambda)
    }
  } else if (!is.null(x$lambda)) {
    cat("\nlambda (P-spline smoothing):\n")
    print(x$lambda)
  }
  cat("\nDeviance:", x$deviance, " on ", x$df, " degrees of freedom\n")
  cat("Null deviance:", x$null.deviance, "on", x$null.df, " degrees of freedom\n")
  if (identical(x$method, "bayes") && !is.null(x$posterior)) {
    cat("\nPosterior summary (tmbstan):\n")
    print(x$posterior, row.names = FALSE, digits = 4)
    if (is.finite(x$time_bayes)) {
      cat(sprintf("Bayes sampling time: %.2fs\n", x$time_bayes))
    }
  }
  if (x$balanced == "yes") {
    cat("\nBalanced data, m =", x$m, "\n")
  } else {
    cat("\nUnbalanced m (vector).\n")
  }
  invisible(x)
}

#' @export
summary.BBreg <- function(object, ...) {
  beta <- if (!is.null(object$beta)) object$beta else {
    setNames(as.numeric(object$coefficients), rownames(object$coefficients))
  }
  beta.se <- sqrt(diag(object$vcov))
  beta.tval <- as.numeric(beta) / beta.se
  beta.TAB <- cbind(
    Estimate = as.numeric(beta),
    StdErr = beta.se,
    t.value = beta.tval,
    p.value = 2 * pt(-abs(beta.tval), df = max(object$df, 1L))
  )
  rownames(beta.TAB) <- names(beta)

  log_phi <- if (!is.null(object$log_phi)) {
    as.numeric(object$log_phi)[1L]
  } else {
    as.numeric(object$psi)[1L]
  }
  phi_tab <- cbind(
    Estimate = c(phi = as.numeric(object$phi)[1L], log_phi = log_phi),
    StdErr = c(NA_real_, sqrt(as.numeric(object$psi.var)[1L]))
  )
  # SE for phi via delta method if log_phi SE available
  se_log <- sqrt(as.numeric(object$psi.var)[1L])
  if (is.finite(se_log)) {
    phi_tab[1L, "StdErr"] <- as.numeric(object$phi)[1L] * se_log
  }

  coefficients.psi <- cbind(
    Estimate = log_phi,
    StdErr = se_log
  )
  rownames(coefficients.psi) <- "log(phi)"

  Chi <- object$null.deviance - object$deviance
  Chi.p.value <- 1 - pchisq(Chi, object$null.df - object$df)

  out <- list(
    call = object$call,
    model = if (!is.null(object$model)) object$model else .bb_model_lines_reg(),
    coefficients = beta.TAB,
    beta.table = beta.TAB,
    phi.table = phi_tab,
    psi.table = coefficients.psi,
    deviance = object$deviance,
    df = object$df,
    null.deviance = object$null.deviance,
    null.df = object$null.df,
    Goodness.of.fit = Chi.p.value,
    iter = object$iter,
    nObs = object$nObs,
    m = object$m,
    balanced = object$balanced,
    conv = object$conv,
    method = object$method
  )
  class(out) <- "summary.BBreg"
  out
}

#' @export
print.summary.BBreg <- function(x, ...) {
  cat("Call:\t")
  print(x$call)
  .cat_model_block("Model (BBreg):", x$model)
  cat("beta (logit mean):\n\n")
  printCoefmat(x$beta.table, P.values = TRUE, has.Pvalue = TRUE)
  cat("\n---------------------------------------------------------------\n")
  cat("phi (BB dispersion) and log(phi):\n\n")
  print(x$phi.table)
  cat("\n---------------------------------------------------------------\n")
  cat("Deviance:", x$deviance, " on ", x$df, " degrees of freedom\n")
  cat("Null deviance:", x$null.deviance, " on ", x$null.df, " degrees of freedom\n")
  cat("Deviance test p-value:", x$Goodness.of.fit, "\n")
  if (x$balanced == "yes") {
    cat("\nBalanced data, m =", x$m)
  } else {
    cat("\nUnbalanced m.")
  }
  cat("\nNumber of iterations:", x$iter, "\n")
  invisible(x)
}

#' @export
coef.BBreg <- function(object, ...) {
  if (!is.null(object$beta)) return(object$beta)
  setNames(as.numeric(object$coefficients), rownames(object$coefficients))
}

#' @export
vcov.BBreg <- function(object, ...) object$vcov

#' @export
fitted.BBreg <- function(object, ...) {
  if (!is.null(object$p)) return(object$p)
  object$fitted.values
}