#' Fit a beta-binomial logistic regression via TMB
#'
#' Marginal MLE of
#' \deqn{y_i \sim \mathrm{BB}(m_i, p_i, \phi),\quad
#' \mathrm{logit}(p_i) = x_i^\top\beta}
#' using Template Model Builder automatic differentiation.
#'
#' @param formula Model formula for the mean (logit link).
#' @param m Maximum score (scalar or vector of length n).
#' @param data Optional data frame.
#' @param maxiter Maximum `nlminb` iterations.
#' @param control Passed to [stats::nlminb()].
#' @param silent Suppress TMB tracing.
#' @return Object of class `BBreg` (API aligned with PROreg::BBreg).
#' @export
BBreg <- function(formula, m, data = list(), maxiter = 100,
                  control = list(), silent = TRUE) {
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
  phi_mm <- (v - mbar * mean(p0) * (1 - mean(p0))) /
    (mbar * mean(p0) * (1 - mean(p0)) * (mbar - 1) - (v - mbar * mean(p0) * (1 - mean(p0))) + 1e-8)
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
    # Drop log_phi row/col for beta block
    nm <- names(sdr$par.fixed)
    b_idx <- grep("^beta", nm)
    phi_idx <- which(nm == "log_phi")
    if (length(b_idx)) {
      vcov.b <- vcov.b[b_idx, b_idx, drop = FALSE]
      rownames(vcov.b) <- colnames(vcov.b) <- colnames(X)
    }
    psi.var <- if (length(phi_idx)) as.numeric(sdr$cov.fixed[phi_idx, phi_idx]) else NA_real_
    # ADREPORT phi SE if available
    if (!is.null(sdr$sd) && "phi" %in% names(sdr$value)) {
      # keep psi.var from Hessian of log_phi
    }
  } else {
    # Numerical Hessian fallback
    H <- optimHess(opt$par, obj$fn, obj$gr)
    invH <- tryCatch(solve(H), error = function(e) matrix(NA_real_, length(opt$par), length(opt$par)))
    b_idx <- seq_len(length(beta))
    vcov.b <- invH[b_idx, b_idx, drop = FALSE]
    rownames(vcov.b) <- colnames(vcov.b) <- colnames(X)
    psi.var <- invH[length(opt$par), length(opt$par)]
  }

  eta <- as.numeric(X %*% beta)
  fitted.values <- 1 / (1 + exp(-eta))

  # Deviance (same definition as PROreg)
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
  # Saturated: p = y/m (with care at 0/m)
  p_sat <- pmin(pmax(y / m., 1e-8), 1 - 1e-8)
  deviance <- as.numeric(2 * (loglik(p_sat, phi) - loglik(fitted.values, phi)))
  null.p <- rep(e, n)
  null.p <- pmin(pmax(null.p, 1e-8), 1 - 1e-8)
  # Null model: re-estimate phi? PROreg uses same phi; keep same phi for parity of df
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
    opt = opt,
    obj = obj,
    sdreport = sdr,
    nll = opt$objective
  )
  class(out) <- "BBreg"
  out$call <- match.call()
  out$formula <- formula
  out
}

#' @export
print.BBreg <- function(x, ...) {
  cat("Call:\t")
  print(x$call)
  cat("\nBeta coefficients:\n")
  print(t(x$coefficients))
  cat("\nDispersion parameter:", x$phi, "\n")
  cat("\nDeviance:", x$deviance, " on ", x$df, " degrees of freedom\n")
  cat("Null deviance:", x$null.deviance, "on", x$null.df, " degrees of freedom\n")
  if (x$balanced == "yes") {
    cat("\nBalanced data, maximum score in the trials:", x$m, "\n")
  } else {
    cat("\nNo balanced data.\n")
  }
  invisible(x)
}

#' @export
summary.BBreg <- function(object, ...) {
  beta.se <- sqrt(diag(object$vcov))
  beta.tval <- as.numeric(object$coefficients) / beta.se
  beta.TAB <- cbind(
    Estimate = as.numeric(object$coefficients),
    StdErr = beta.se,
    t.value = beta.tval,
    p.value = 2 * pt(-abs(beta.tval), df = max(object$df, 1L))
  )
  rownames(beta.TAB) <- rownames(object$coefficients)

  coefficients.psi <- cbind(
    Estimate = object$psi,
    StdErr = sqrt(object$psi.var)
  )
  rownames(coefficients.psi) <- "log(phi)"

  Chi <- object$null.deviance - object$deviance
  Chi.p.value <- 1 - pchisq(Chi, object$null.df - object$df)

  out <- list(
    call = object$call,
    coefficients = beta.TAB,
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
    conv = object$conv
  )
  class(out) <- "summary.BBreg"
  out
}

#' @export
print.summary.BBreg <- function(x, ...) {
  cat("Call:\t")
  print(x$call)
  cat("\nBeta coefficients:\n\n")
  printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
  cat("\n---------------------------------------------------------------\n")
  cat("Dispersion parameter coefficients:\n\n")
  print(x$psi.table)
  cat("\n---------------------------------------------------------------\n")
  cat("Deviance:", x$deviance, " on ", x$df, " degrees of freedom\n")
  cat("Null deviance:", x$null.deviance, " on ", x$null.df, " degrees of freedom\n")
  cat("Deviance test p-value:", x$Goodness.of.fit, "\n")
  if (x$balanced == "yes") {
    cat("\nBalanced data, maximum score in the trials:", x$m)
  } else {
    cat("\nNo balanced data.")
  }
  cat("\nNumber of iterations:", x$iter, "\n")
  invisible(x)
}

#' @export
coef.BBreg <- function(object, ...) {
  setNames(as.numeric(object$coefficients), rownames(object$coefficients))
}

#' @export
vcov.BBreg <- function(object, ...) object$vcov

#' @export
fitted.BBreg <- function(object, ...) object$fitted.values
