#' Fit a beta-binomial logistic regression via TMB
#'
#' Model
#' \deqn{y_i \sim \mathrm{BB}(m_i, p_i, \phi),\quad
#' \mathrm{logit}(p_i) = x_i^\top\beta}
#' estimated by Template Model Builder.
#'
#' @param formula Model formula for the mean (logit link).
#' @param m Maximum score (scalar or vector of length n).
#' @param data Optional data frame.
#' @param method `"mle"` (default): maximize the TMB likelihood with
#'   `nlminb`. `"bayes"`: same point fit, then NUTS via
#'   [tmbstan::tmbstan()] on a template with weak priors
#'   \eqn{\beta_j\sim N(0,5^2)}, \eqn{\log\phi\sim N(\log 0.3,1)}
#'   (requires suggested packages **tmbstan** and **rstan**).
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
  if (any(m != as.integer(m)) || min(m) <= 0) {
    stop("m must be positive integer(s)", call. = FALSE)
  }

  mf <- model.frame(formula = formula, data = data)
  X <- model.matrix(attr(mf, "terms"), data = mf)
  y <- as.numeric(model.response(mf))
  n <- length(y)

  if (length(m) == 1L) {
    balanced <- "yes"
    m. <- rep(as.numeric(m), n)
  } else {
    m. <- as.numeric(m)
    if (length(m.) != n) stop("m must be scalar or length(y)", call. = FALSE)
    balanced <- if (length(unique(m.)) == 1L) "yes" else "no"
  }

  if (any(y != as.integer(y))) stop("y must be integer", call. = FALSE)
  if (any(y < 0 | y > m.)) stop("y must be bounded between 0 and m", call. = FALSE)

  ensure_tmb_dll("bb_reg")

  # Starting values: binomial GLM + moment phi
  glm0 <- stats::glm.fit(X, y / m., family = stats::binomial(), weights = m.)
  beta0 <- as.numeric(glm0$coefficients)
  beta0[!is.finite(beta0)] <- 0
  p0 <- as.numeric(glm0$fitted.values)
  p0 <- pmin(pmax(p0, 1e-4), 1 - 1e-4)
  mu <- mean(y)
  v <- stats::var(y)
  mbar <- mean(m.)
  # Simpler MM matching PROreg spirit
  phi_mm2 <- {
    ph <- (v - mu * (1 - mu / mbar)) / (mu * (1 - mu / mbar) * (mbar - 1) - (v - mu * (1 - mu / mbar)))
    if (!is.finite(ph) || ph <= 0) 0.5 else ph
  }

  parameters <- list(beta = beta0, log_phi = log(max(phi_mm2, 1e-3)))
  data_tmb <- list(y = y, m = m., X = X)

  obj <- TMB::MakeADFun(
    data = data_tmb,
    parameters = parameters,
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
  sdr <- tryCatch(TMB::sdreport(obj), error = function(e) NULL)

  beta <- opt$par[grep("^beta", names(opt$par))]
  if (length(beta) == 0L) beta <- opt$par[seq_len(ncol(X))]
  names(beta) <- colnames(X)
  log_phi <- unname(opt$par["log_phi"])
  phi <- exp(log_phi)

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

  df <- n - length(beta) - 1L
  null.df <- n - 1L

  coef.b <- matrix(beta, ncol = 1L, dimnames = list(names(beta), NULL))

  out <- list(
    coefficients = coef.b,
    vcov = vcov.b,
    phi = phi,
    psi = log_phi,
    psi.var = psi.var,
    conv = conv,
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
    posterior = NULL,
    stanfit = NULL,
    time_bayes = NA_real_
  )
  class(out) <- "BBreg"
  out$call <- match.call()
  out <- .bbreg_attach_formulation(out, formula = formula)

  if (identical(method, "bayes")) {
    out <- .BBreg_add_bayes(
      out, data_tmb = data_tmb,
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