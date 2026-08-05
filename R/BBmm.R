#' Beta-binomial mixed-effects model via TMB
#'
#' Fits
#' \deqn{y \mid u \sim \mathrm{BB}(m, p, \phi),\quad
#' \mathrm{logit}(p) = X\beta + Zu,\quad u \sim N(0, D)}
#' with Laplace approximation of the marginal likelihood (TMB).
#'
#' Specify the random structure either with `random.formula` or with
#' design matrix `Z` plus `nRandComp`.
#'
#' @param fixed.formula Fixed-effects formula (response ~ covariates).
#' @param X Fixed-effects design matrix (alternative to `fixed.formula`).
#' @param y Response vector (required if `X` is supplied). For
#'   multidimensional models, stack outcomes one after another.
#' @param random.formula Random-effects formula (e.g. `~ group`).
#' @param Z Random-effects design matrix.
#' @param nRandComp Number of random effects per variance component
#'   (required with `Z`).
#' @param m Maximum score (scalar or vector).
#' @param data Optional data frame.
#' @param maxiter Maximum `nlminb` iterations.
#' @param show Logical; print progress.
#' @param nDim Number of response dimensions (default 1).
#' @param silent Suppress TMB tracing.
#' @param control Extra [stats::nlminb()] control.
#' @return Object of class `BBmm`.
#' @export
BBmm <- function(fixed.formula, X, y, random.formula = NULL, Z = NULL,
                 nRandComp = NULL, m, data = list(),
                 maxiter = 100, show = FALSE, nDim = 1L,
                 silent = TRUE, control = list()) {

  nDim <- as.integer(nDim)
  if (nDim < 1L) stop("nDim must be >= 1", call. = FALSE)

  # ----- Fixed part -----
  if (!missing(fixed.formula)) {
    if (!missing(X) || !missing(y)) {
      stop("Specify either fixed.formula or (X, y), not both", call. = FALSE)
    }
    mf <- model.frame(formula = fixed.formula, data = data)
    X <- model.matrix(attr(mf, "terms"), data = mf)
    y <- as.numeric(model.response(mf))
    formula_out <- fixed.formula
  } else {
    if (missing(X) || missing(y)) {
      stop("Provide fixed.formula or both X and y", call. = FALSE)
    }
    X <- as.matrix(X)
    y <- as.numeric(y)
    formula_out <- NULL
  }

  nObs <- length(y)
  if (nrow(X) != nObs) stop("nrow(X) must equal length(y)", call. = FALSE)

  if (any(m != as.integer(m)) || min(m) <= 0) {
    stop("m must be positive integer(s)", call. = FALSE)
  }
  if (length(m) == 1L) {
    balanced <- "yes"
    m. <- rep(as.numeric(m), nObs)
  } else {
    m. <- as.numeric(m)
    if (length(m.) != nObs) stop("m must be scalar or length(y)", call. = FALSE)
    balanced <- if (length(unique(m.)) == 1L) "yes" else "no"
  }
  if (any(y != as.integer(y))) stop("y must be integer", call. = FALSE)
  if (any(y < 0 | y > m.)) stop("y must be bounded between 0 and m", call. = FALSE)
  if (nObs %% nDim != 0L) {
    stop("length(y) must be divisible by nDim", call. = FALSE)
  }

  # ----- Random part -----
  if (is.null(random.formula) && is.null(Z)) {
    stop("Random part must be specified (random.formula or Z)", call. = FALSE)
  }
  if (!is.null(random.formula) && !is.null(Z)) {
    stop("Random part specified twice", call. = FALSE)
  }
  if (!is.null(Z) && is.null(nRandComp)) {
    stop("nRandComp must be specified when Z is given", call. = FALSE)
  }
  if (is.null(Z) && !is.null(nRandComp)) {
    stop("nRandComp only used when Z is given", call. = FALSE)
  }

  if (is.null(Z)) {
    # ~ factor1 + factor2  -> one variance component per factor
    rf <- stats::update(random.formula, ~ . - 1)
    random.mf <- model.frame(formula = rf, data = data)
    nComp <- ncol(random.mf)
    nRandComp <- integer(nComp)
    Z <- NULL
    namesRand <- names(random.mf)
    for (i in seq_len(nComp)) {
      z <- model.matrix(~ random.mf[[i]] - 1)
      Z <- cbind(Z, z)
      nRandComp[i] <- ncol(z)
    }
    Z <- as.matrix(Z)
  } else {
    Z <- as.matrix(Z)
    nComp <- length(nRandComp)
    nRandComp <- as.integer(nRandComp)
    if (ncol(Z) != sum(nRandComp)) {
      stop("sum(nRandComp) must equal ncol(Z)", call. = FALSE)
    }
    namesRand <- as.character(seq_len(nComp))
  }

  if (nrow(Z) != nObs) stop("nrow(Z) must equal length(y)", call. = FALSE)
  nRand <- ncol(Z)
  re_comp <- rep(seq_len(nComp) - 1L, times = nRandComp)

  # Dimension id for stacked multi-response
  n_per_dim <- nObs / nDim
  dim_id <- rep(seq_len(nDim) - 1L, each = n_per_dim)

  ensure_tmb_dll("bb_mm")

  # ----- Starting values -----
  # Intercept-only / fixed-only BB as warm start (ignore RE)
  if (!is.null(formula_out) && nDim == 1L) {
    bb0 <- tryCatch(
      BBreg(formula_out, m = m., data = data, silent = TRUE),
      error = function(e) NULL
    )
  } else {
    bb0 <- NULL
  }

  if (!is.null(bb0) && identical(bb0$conv, "yes") && length(bb0$coefficients) == ncol(X)) {
    beta0 <- as.numeric(bb0$coefficients)
    phi0 <- rep(bb0$phi, nDim)
  } else {
    glm0 <- tryCatch(
      stats::glm.fit(X, pmin(pmax(y / m., 1e-4), 1 - 1e-4),
                     family = stats::binomial(), weights = m.),
      error = function(e) NULL
    )
    if (!is.null(glm0)) {
      beta0 <- as.numeric(glm0$coefficients)
      beta0[!is.finite(beta0)] <- 0
    } else {
      beta0 <- rep(0, ncol(X))
      beta0[1] <- stats::qlogis(mean(y / m.))
    }
    phi0 <- rep(0.5, nDim)
  }

  parameters <- list(
    beta = beta0,
    log_phi = log(pmax(phi0, 1e-3)),
    log_sigma = rep(log(0.5), nComp),
    u = rep(0, nRand)
  )

  data_tmb <- list(
    y = y,
    m = m.,
    X = X,
    Z = Z,
    re_comp = as.integer(re_comp),
    dim_id = as.integer(dim_id),
    nDim = nDim
  )

  obj <- TMB::MakeADFun(
    data = data_tmb,
    parameters = parameters,
    random = "u",
    DLL = "bb_mm",
    silent = silent
  )

  ctrl <- modifyList(list(iter.max = maxiter, eval.max = maxiter * 2L), control)
  if (show) cat("Optimizing Laplace approximate marginal likelihood...\n")

  opt <- tryCatch(
    stats::nlminb(obj$par, obj$fn, obj$gr, control = ctrl),
    error = function(e) e
  )
  if (inherits(opt, "error")) {
    return(list(conv = "no", message = conditionMessage(opt)))
  }

  conv <- if (opt$convergence == 0) "yes" else "no"
  sdr <- tryCatch(TMB::sdreport(obj, getJointPrecision = FALSE), error = function(e) NULL)

  # Extract estimates
  par_fixed <- opt$par
  nm <- names(par_fixed)
  beta <- unname(par_fixed[grep("^beta", nm)])
  if (length(beta) == 0L) beta <- unname(par_fixed[seq_len(ncol(X))])
  names(beta) <- colnames(X)

  log_phi <- unname(par_fixed[grep("^log_phi", nm)])
  log_sigma <- unname(par_fixed[grep("^log_sigma", nm)])
  phi <- exp(log_phi)
  all.sigma <- exp(log_sigma)
  names(all.sigma) <- namesRand

  # Random effects: last mode
  u_hat <- tryCatch(as.numeric(obj$env$last.par.best[obj$env$random]),
                    error = function(e) rep(NA_real_, nRand))
  names(u_hat) <- seq_len(nRand)

  # Variances from sdreport
  fixed.vcov <- matrix(NA_real_, length(beta), length(beta),
                       dimnames = list(names(beta), names(beta)))
  psi.var <- rep(NA_real_, nDim)
  all.sigma.var <- rep(NA_real_, nComp)

  if (!is.null(sdr)) {
    covf <- as.matrix(sdr$cov.fixed)
    fnm <- names(sdr$par.fixed)
    b_idx <- grep("^beta", fnm)
    if (length(b_idx) == length(beta)) {
      fixed.vcov <- covf[b_idx, b_idx, drop = FALSE]
      rownames(fixed.vcov) <- colnames(fixed.vcov) <- names(beta)
    }
    phi_idx <- grep("^log_phi", fnm)
    if (length(phi_idx)) psi.var <- diag(covf)[phi_idx]
    sig_idx <- grep("^log_sigma", fnm)
    if (length(sig_idx)) {
      # delta method for sigma = exp(log_sigma): Var(sigma) = sigma^2 Var(log_sigma)
      all.sigma.var <- (all.sigma^2) * diag(covf)[sig_idx]
    }
  }

  eta <- as.numeric(X %*% beta + Z %*% u_hat)
  fitted <- 1 / (1 + exp(-eta))

  # Simple deviance using conditional BB (like PROreg conditional fit)
  loglik_bb <- function(p_hat, phi_hat, y_, m_) {
    # phi_hat recycled by dim
    a <- p_hat / phi_hat
    b <- (1 - p_hat) / phi_hat
    sum(
      lgamma(m_ + 1) - lgamma(y_ + 1) - lgamma(m_ - y_ + 1) +
        lgamma(a + y_) - lgamma(a) +
        lgamma(b + m_ - y_) - lgamma(b) +
        lgamma(a + b) - lgamma(a + b + m_)
    )
  }
  phi_i <- phi[dim_id + 1L]
  p_sat <- pmin(pmax(y / m., 1e-8), 1 - 1e-8)
  deviance <- as.numeric(2 * (loglik_bb(p_sat, phi_i, y, m.) -
                                loglik_bb(fitted, phi_i, y, m.)))
  e <- sum(y) / sum(m.)
  null.p <- pmin(pmax(rep(e, nObs), 1e-8), 1 - 1e-8)
  null.deviance <- as.numeric(2 * (loglik_bb(p_sat, phi_i, y, m.) -
                                     loglik_bb(null.p, phi_i, y, m.)))
  df <- nObs - length(beta) - length(all.sigma) - nDim
  null.df <- nObs - 1L - length(all.sigma) - nDim

  # D matrix
  d <- rep(all.sigma^2, times = nRandComp)
  D <- diag(d, nrow = nRand)

  psi <- log(phi)
  names(psi) <- if (nDim == 1L) "log(phi)" else paste0("log(phi)", seq_len(nDim))

  out <- list(
    fixed.coef = beta,
    fixed.vcov = fixed.vcov,
    random.coef = u_hat,
    sigma.coef = all.sigma,
    sigma.var = all.sigma.var,
    phi.coef = if (nDim == 1L) phi[1] else phi,
    psi.coef = if (nDim == 1L) psi[1] else psi,
    psi.var = if (nDim == 1L) psi.var[1] else psi.var,
    fitted.values = fitted,
    conv = conv,
    deviance = deviance,
    df = df,
    null.deviance = null.deviance,
    null.df = null.df,
    nRand = nRand,
    nComp = nComp,
    nRandComp = nRandComp,
    namesRand = namesRand,
    iter = as.integer(opt$iterations),
    nObs = nObs,
    nDim = nDim,
    y = y,
    X = X,
    Z = Z,
    D = D,
    balanced = balanced,
    m = if (balanced == "yes") m.[1] else m.,
    opt = opt,
    obj = obj,
    sdreport = sdr,
    nll = opt$objective
  )
  class(out) <- "BBmm"
  out$call <- match.call()
  out$formula <- formula_out
  out
}

#' @export
print.BBmm <- function(x, ...) {
  cat("Call:\t")
  print(x$call)
  cat("\nFixed effects estimation:\n")
  print(x$fixed.coef)
  cat("\nStandard deviation of normal random effects:\n")
  for (i in seq_along(x$namesRand)) {
    cat(x$namesRand[i], x$sigma.coef[i], "\n")
  }
  cat("\nBeta-binomial dispersion parameter:",
      paste(x$phi.coef, collapse = ", "), "\n")
  cat("\nDeviance of the model:", x$deviance)
  cat("\nNumber of iterations:", x$iter)
  if (x$balanced == "yes") {
    cat("\nBalanced data, maximum score number:", x$m, "\n")
  } else {
    cat("\nNo balanced data.\n")
  }
  invisible(x)
}

#' @export
summary.BBmm <- function(object, ...) {
  fixed.se <- sqrt(diag(object$fixed.vcov))
  fixed.tval <- object$fixed.coef / fixed.se
  fixed.TAB <- cbind(
    Estimate = object$fixed.coef,
    StdErr = fixed.se,
    t.value = fixed.tval,
    p.value = 2 * pnorm(-abs(fixed.tval))
  )

  psi <- object$psi.coef
  psi.se <- sqrt(object$psi.var)
  if (length(psi) == 1L) {
    psi.table <- cbind(Estimate = psi, StdErr = psi.se)
    rownames(psi.table) <- "log(phi)"
  } else {
    psi.table <- cbind(Estimate = psi, StdErr = psi.se)
  }

  sigma.table <- cbind(
    Estimate = object$sigma.coef,
    StdErr = sqrt(object$sigma.var)
  )
  rownames(sigma.table) <- object$namesRand

  Chi <- object$null.deviance - object$deviance
  Chi.p.value <- 1 - pchisq(Chi, object$null.df - object$df)

  res <- list(
    call = object$call,
    fixed.coefficients = fixed.TAB,
    sigma.table = sigma.table,
    psi.table = psi.table,
    random.coef = object$random.coef,
    iter = object$iter,
    nObs = object$nObs,
    nRand = object$nRand,
    nComp = object$nComp,
    nRandComp = object$nRandComp,
    deviance = object$deviance,
    df = object$df,
    null.deviance = object$null.deviance,
    null.df = object$null.df,
    Goodness.of.fit = Chi.p.value,
    balanced = object$balanced,
    m = object$m,
    conv = object$conv,
    nll = object$nll
  )
  class(res) <- "summary.BBmm"
  res
}

#' @export
print.summary.BBmm <- function(x, ...) {
  cat("Call:\t")
  print(x$call)
  cat("\nFixed effects coefficients:\n\n")
  printCoefmat(x$fixed.coefficients, P.values = TRUE, has.Pvalue = TRUE)
  cat("\n---------------------------------------------------------------\n")
  cat("Random effects dispersion parameter(s):\n\n")
  print(x$sigma.table)
  cat("\n---------------------------------------------------------------\n")
  cat("Logarithm of beta-binomial dispersion parameter log(phi):\n\n")
  print(x$psi.table)
  cat("\n---------------------------------------------------------------\n")
  cat("Deviance of the model:", x$deviance, "; with", x$df, "degrees of freedom.\n")
  cat("Deviance of the null model", x$null.deviance, "; with", x$null.df, "degrees of freedom.\n")
  cat("Deviance goodness-of-fit test p-value:", x$Goodness.of.fit, "\n")
  if (!is.null(x$nll)) cat("Laplace approx. nll:", x$nll, "\n")
  cat("\nNumber of observations:", x$nObs)
  cat("\nNumber of iterations:", x$iter)
  if (x$balanced == "yes") {
    cat("\nBalanced data, maximum score number:", x$m)
  } else {
    cat("\nNo balanced data.")
  }
  cat("\nNumber of random effects in each random component:", x$nRandComp, "\n\n")
  invisible(x)
}

#' @export
coef.BBmm <- function(object, ...) object$fixed.coef

#' @export
vcov.BBmm <- function(object, ...) object$fixed.vcov

#' @export
fitted.BBmm <- function(object, ...) object$fitted.values
