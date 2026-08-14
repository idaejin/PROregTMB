#' Beta-binomial mixed-effects model via TMB
#'
#' Unified shared-latent formulation
#' \deqn{\mathrm{logit}(p_{ij}^{(\ell)})
#'   = x_{ij}^{(\ell)\top}\beta^{(\ell)} + z_{ij}^{\top} a_i,
#'   \quad a_i \sim N(0,G),}
#' with Laplace approximation of the marginal likelihood. Here \(\ell\) is the
#' outcome dimension (\(L=1\) univariate), \(j\) indexes visits
#' (\(n_i=1\) cross-sectional), and \(z_{ij}=(1)\) or \((1,t_{ij})^\top\).
#'
#' **Random effects**
#' \itemize{
#'   \item \code{random = ~ (1 | id)} — random intercept
#'   \item \code{random = ~ (1 + time | id)} — RI + RS;
#'     \code{corr = "unstructured"} (default if \(q>1\)) or \code{"diag"}
#'   \item Legacy: \code{random.formula = ~ id}
#'   \item Advanced: \code{Z} + \code{nRandComp}
#' }
#'
#' **Multivariate (shared \(a_i\))**
#' \itemize{
#'   \item Wide (CS): \code{cbind(y1,y2,y3) ~ x} with \code{random = ~ (1|id)}
#'   \item Long (CS or longitudinal): \code{y ~ x + time} with
#'     \code{dim = "domain"} (column naming the dimension \(\ell\))
#' }
#' Low-level stacking: [multi_bb_stack()].
#'
#' @param fixed.formula Fixed-effects formula. Multivariate wide form:
#'   \code{cbind(y1,y2) ~ x}.
#' @param X Fixed-effects design matrix (alternative to `fixed.formula`).
#' @param y Response vector (required if `X` is supplied).
#' @param random Random-effects formula with `|` bars, or a named list.
#' @param random.formula Legacy PROreg-style grouping factors only.
#' @param Z Random-effects design matrix.
#' @param nRandComp Number of RE per variance component (with `Z`).
#' @param corr Within-subject \(G\) for multi-term RE: `"unstructured"` /
#'   `"us"` or `"diag"`. Aliases: `"cor"`, `"correlated"`.
#' @param dim Character: name of the dimension/domain column for **long**
#'   multivariate data (shared latent \(a_i\)). Mutually exclusive with
#'   \code{cbind()} responses.
#' @param m Maximum score (scalar, vector, length-\(L\), or column name).
#' @param data Data frame.
#' @param method `"mle"` (default) or `"bayes"`.
#' @param kappa Weight for sum-to-zero centering of each `s()` (default `1e8`).
#' @param maxiter Maximum `nlminb` iterations.
#' @param show Logical; print progress.
#' @param nDim Number of dimensions (usually auto-set from `cbind` / `dim`).
#' @param silent Suppress TMB tracing.
#' @param control Extra [stats::nlminb()] control.
#' @param chains,iter,warmup,seed Stan controls when `method = "bayes"`.
#' @param laplace_bayes If `TRUE`, tmbstan uses Laplace for `u`.
#' @return Object of class `BBmm`.
#' @export
BBmm <- function(fixed.formula, X, y, random = NULL,
                 random.formula = NULL, Z = NULL,
                 nRandComp = NULL,
                 corr = c("unstructured", "us", "diag"),
                 dim = NULL,
                 m, data = list(),
                 method = c("mle", "bayes"),
                 kappa = 1e6,
                 maxiter = 100, show = FALSE, nDim = 1L,
                 silent = TRUE, control = list(),
                 chains = 2L, iter = 1000L, warmup = 400L, seed = 1L,
                 laplace_bayes = TRUE) {
  method <- match.arg(method)
  corr <- .normalize_corr(if (length(corr)) corr[1] else "unstructured")
  nDim <- as.integer(nDim)
  if (nDim < 1L) stop("nDim must be >= 1", call. = FALSE)

  formula_out <- NULL
  multi <- NULL
  sm <- NULL

  # ----- Multivariate clean API (cbind / dim=) -----
  if (!missing(fixed.formula) &&
      ( !is.null(dim) || .is_cbind_response(fixed.formula) )) {
    if (length(grep("^s\\s*\\(",
                    attr(stats::terms(fixed.formula, specials = "s"),
                         "term.labels"))) ) {
      stop("s() smooths are not yet supported with cbind()/dim= multivariate API",
           call. = FALSE)
    }
    if (is.null(data) || !(is.data.frame(data) || is.list(data))) {
      stop("multivariate BBmm requires data = ...", call. = FALSE)
    }
    if (!is.null(Z) || !is.null(nRandComp)) {
      stop("Do not pass Z/nRandComp with cbind()/dim=; use random=",
           call. = FALSE)
    }
    multi <- .prepare_multivariate_bbmm(
      fixed.formula = fixed.formula,
      data = data,
      m = m,
      dim = dim,
      random = random,
      random.formula = random.formula,
      corr = corr
    )
    y <- multi$y
    X <- multi$X
    m. <- multi$m
    nDim <- multi$nDim
    data <- multi$data_re
    re <- multi$re
    formula_out <- multi$formula_out
    balanced <- if (length(unique(m.)) == 1L) "yes" else "no"
    nObs <- length(y)
  } else {
    # ----- Univariate / manual X,y -----
    if (!missing(fixed.formula)) {
      if (!missing(X) || !missing(y)) {
        stop("Specify either fixed.formula or (X, y), not both", call. = FALSE)
      }
      if (is.null(data) || (is.list(data) && !is.data.frame(data) && !length(data))) {
        stop("fixed.formula requires data", call. = FALSE)
      }
      data <- as.data.frame(data)
      m <- .resolve_m_from_data(m, data, nrow(data))
      sm <- .build_smooth_design(fixed.formula, data = data)
      mf <- model.frame(formula = sm$fixed, data = data)
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
      sm <- list(n_smooth = 0L)
    }

    nObs <- length(y)
    if (nrow(X) != nObs) stop("nrow(X) must equal length(y)", call. = FALSE)

    if (is.character(m)) m <- .resolve_m_from_data(m, data, nObs)
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

    if (is.null(data) || (is.list(data) && !is.data.frame(data) && !length(data))) {
      data <- data.frame(y = y)
    } else {
      data <- as.data.frame(data)
    }

    re <- .build_re_structure(
      random = random,
      random.formula = random.formula,
      Z = Z,
      nRandComp = nRandComp,
      data = data,
      corr = corr,
      nObs = nObs
    )
  }

  if (any(y != as.integer(y))) stop("y must be integer", call. = FALSE)
  if (any(m. != as.integer(m.)) || min(m.) <= 0) {
    stop("m must be positive integer(s)", call. = FALSE)
  }
  if (any(y < 0 | y > m.)) stop("y must be bounded between 0 and m", call. = FALSE)
  if (is.null(multi) && (nObs %% nDim != 0L)) {
    stop("length(y) must be divisible by nDim", call. = FALSE)
  }
  if (is.null(sm)) sm <- list(n_smooth = 0L)
  if (identical(method, "bayes") && isTRUE(sm$n_smooth > 0L)) {
    stop('method = "bayes" does not support s() smooths yet', call. = FALSE)
  }

  Z <- re$Z
  nRand <- re$nRand
  namesRand <- re$namesRand
  nRandComp <- re$nRandComp
  nComp <- re$nComp

  # Dimension id for phi_ell
  if (!is.null(multi) && identical(multi$mode, "long")) {
    # blocks may have unequal size: build dim_id from ordered factor
    dcol <- multi$dim_names
    # y stacked by dimension in .prepare; use lengths of unique runs
    # Reconstruct from nDim equal blocks if balanced; else from data
    if (!is.null(dim) && dim %in% names(data)) {
      dim_id <- as.integer(droplevels(as.factor(data[[dim]]))) - 1L
    } else {
      n_per_dim <- nObs / nDim
      dim_id <- rep(seq_len(nDim) - 1L, each = n_per_dim)
    }
  } else {
    n_per_dim <- nObs / nDim
    dim_id <- rep(seq_len(nDim) - 1L, each = n_per_dim)
  }
  if (length(dim_id) != nObs) {
    stop("internal dim_id length mismatch", call. = FALSE)
  }

  ensure_tmb_dll("bb_mm")

  # ----- Starting values -----
  if (!is.null(formula_out) && nDim == 1L) {
    bb0 <- tryCatch(
      BBreg(formula_out, m = m., data = data, kappa = kappa, silent = TRUE),
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

  spar <- .tmb_smooth_parameters(sm)
  parameters <- c(
    list(
      beta = beta0,
      log_phi = log(pmax(phi0, 1e-3)),
      theta_re = re$theta0,
      u = rep(0, nRand)
    ),
    spar
  )

  data_tmb <- c(
    list(
      y = y,
      m = m.,
      X = X,
      Z = .as_dense_Z(Z),
      dim_id = as.integer(dim_id),
      nDim = nDim,
      n_blocks = as.integer(re$n_blocks),
      block_G = as.integer(re$block_G),
      block_q = as.integer(re$block_q),
      block_corr = as.integer(re$block_corr),
      block_theta0 = as.integer(re$block_theta0)
    ),
    .tmb_smooth_data(sm, n = nObs, kappa = kappa)
  )

  random_tmb <- if (isTRUE(sm$n_smooth > 0L)) c("u", "alpha") else "u"
  obj <- TMB::MakeADFun(
    data = data_tmb,
    parameters = parameters,
    random = random_tmb,
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
  sdr <- tryCatch(
    TMB::sdreport(obj, getJointPrecision = isTRUE(sm$n_smooth > 0L)),
    error = function(e) NULL
  )

  # Extract estimates
  par_fixed <- opt$par
  nm <- names(par_fixed)
  beta <- unname(par_fixed[grep("^beta", nm)])
  if (length(beta) == 0L) beta <- unname(par_fixed[seq_len(ncol(X))])
  names(beta) <- colnames(X)

  log_phi <- unname(par_fixed[grep("^log_phi", nm)])
  theta_hat <- unname(par_fixed[grep("^theta_re", nm)])
  if (length(theta_hat) == 0L) theta_hat <- re$theta0
  phi <- exp(log_phi)

  # Per-block Sigma / sd / Corr
  Sigma_list <- vector("list", re$n_blocks)
  names(Sigma_list) <- vapply(re$blocks, `[[`, character(1), "name")
  Corr_list <- Sigma_list
  sd_all <- numeric(0)
  sd_names <- character(0)
  for (b in seq_len(re$n_blocks)) {
    q <- re$block_q[b]
    nt <- .n_theta_sigma(q, re$blocks[[b]]$corr)
    i0 <- re$block_theta0[b] + 1L
    th <- theta_hat[i0:(i0 + nt - 1L)]
    unpacked <- .sigma_from_theta(th, q, re$blocks[[b]]$corr)
    tn <- re$blocks[[b]]$term_names
    rownames(unpacked$Sigma) <- colnames(unpacked$Sigma) <- tn
    rownames(unpacked$Corr) <- colnames(unpacked$Corr) <- tn
    names(unpacked$sd) <- tn
    Sigma_list[[b]] <- unpacked$Sigma
    Corr_list[[b]] <- unpacked$Corr
    sd_all <- c(sd_all, unpacked$sd)
    sd_names <- c(sd_names, paste0(re$blocks[[b]]$name, ".", tn))
  }
  names(sd_all) <- sd_names
  all.sigma <- sd_all

  # Random effects: subject u and optional smooth alpha
  rand_hat <- tryCatch(as.numeric(obj$env$last.par.best[obj$env$random]),
                       error = function(e) rep(NA_real_, nRand))
  rand_names <- names(obj$env$last.par.best[obj$env$random])
  if (is.null(rand_names)) rand_names <- rep("u", length(rand_hat))
  # TMB often names all "u" / "alpha" without index; split by length
  u_hat <- rand_hat[seq_len(nRand)]
  names(u_hat) <- colnames(Z)
  alpha <- NULL
  alpha.vcov <- NULL
  lambda <- NULL
  fhat <- NULL
  if (isTRUE(sm$n_smooth > 0L)) {
    n_alpha <- sum(sm$smooth_K)
    alpha <- rand_hat[nRand + seq_len(n_alpha)]
    names(alpha) <- colnames(sm$B)
    lam_hat <- par_fixed[grep("^log_lambda", nm)]
    lambda <- setNames(exp(unname(lam_hat)), sm$labels)
    fhat <- lapply(seq_len(sm$n_smooth), function(j) {
      idx <- sm$blocks_idx[[j]]
      as.numeric(sm$specs[[j]]$B %*% alpha[idx])
    })
    names(fhat) <- sm$labels
    alpha.vcov <- .alpha_vcov_from_joint(
      if (!is.null(sdr)) sdr$jointPrecision else NULL,
      n_alpha = n_alpha,
      n_u = nRand
    )
  }

  fixed.vcov <- matrix(NA_real_, length(beta), length(beta),
                       dimnames = list(names(beta), names(beta)))
  psi.var <- rep(NA_real_, nDim)
  all.sigma.var <- rep(NA_real_, length(all.sigma))

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
  }

  eta <- as.numeric(X %*% beta + as.numeric(Z %*% u_hat))
  if (isTRUE(sm$n_smooth > 0L)) eta <- eta + as.numeric(sm$B %*% alpha)
  fitted <- 1 / (1 + exp(-eta))

  loglik_bb <- function(p_hat, phi_hat, y_, m_) {
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
  n_re_par <- length(theta_hat)
  df <- nObs - length(beta) - n_re_par - nDim
  null.df <- nObs - 1L - n_re_par - nDim

  # Full D (block-diagonal copies of each Sigma)
  D <- matrix(0, nRand, nRand)
  u_off <- 0L
  for (b in seq_len(re$n_blocks)) {
    G <- re$block_G[b]
    q <- re$block_q[b]
    Sb <- Sigma_list[[b]]
    for (g in seq_len(G)) {
      idx <- u_off + (g - 1L) * q + seq_len(q)
      D[idx, idx] <- Sb
    }
    u_off <- u_off + G * q
  }

  psi <- log(phi)
  names(psi) <- if (nDim == 1L) "log(phi)" else paste0("log(phi)", seq_len(nDim))

  # Legacy sigma.coef: for pure RI blocks use block SDs named by factor;
  # for RI+RS keep term-level names
  sigma.coef <- all.sigma

  out <- list(
    fixed.coef = beta,
    fixed.vcov = fixed.vcov,
    random.coef = u_hat,
    sigma.coef = sigma.coef,
    sigma.var = all.sigma.var,
    Sigma = Sigma_list,
    Corr = Corr_list,
    theta_re = theta_hat,
    re_blocks = re$blocks,
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
    method = method,
    corr = corr,
    dim = dim,
    structure = if (!is.null(multi)) "shared" else NULL,
    dim_names = if (!is.null(multi)) multi$dim_names else NULL,
    random = random,
    smooth = if (isTRUE(sm$n_smooth > 0L)) sm else NULL,
    alpha = alpha,
    alpha.vcov = alpha.vcov,
    lambda = lambda,
    fhat = fhat,
    kappa = kappa,
    opt = opt,
    obj = obj,
    sdreport = sdr,
    nll = opt$objective,
    posterior = NULL,
    stanfit = NULL,
    time_bayes = NA_real_
  )
  class(out) <- "BBmm"
  out$call <- match.call()
  out$formula <- formula_out
  out <- .bbmm_attach_formulation(
    out,
    fixed.formula = formula_out,
    random.formula = if (!is.null(random)) random else random.formula
  )

  if (identical(method, "bayes")) {
    out <- .BBmm_add_bayes(
      out, data_tmb = data_tmb, parameters = parameters,
      chains = chains, iter = iter, warmup = warmup, seed = seed,
      silent = silent, laplace_bayes = laplace_bayes
    )
  }
  out
}

# Bayes via tmbstan on bb_mm_prior
.BBmm_add_bayes <- function(fit, data_tmb, parameters, chains, iter,
                            warmup, seed, silent, laplace_bayes = TRUE) {
  if (!requireNamespace("tmbstan", quietly = TRUE) ||
      !requireNamespace("rstan", quietly = TRUE)) {
    stop('method = "bayes" requires packages tmbstan and rstan', call. = FALSE)
  }
  # bb_mm_prior is parametric + subject RE only
  drop <- c("n_smooth", "B", "S", "C", "smooth_K", "smooth_off", "kappa")
  data_tmb <- data_tmb[setdiff(names(data_tmb), drop)]
  parameters <- parameters[setdiff(names(parameters), c("log_lambda", "alpha"))]
  ensure_tmb_dll("bb_mm_prior")
  beta0 <- as.numeric(fit$fixed.coef)
  phi0 <- as.numeric(fit$phi.coef)
  if (length(phi0) == 1L) phi0 <- rep(phi0, fit$nDim)
  th0 <- as.numeric(fit$theta_re)
  u0 <- as.numeric(fit$random.coef)
  pars <- list(
    beta = beta0,
    log_phi = log(pmax(phi0, 1e-6)),
    theta_re = th0,
    u = u0
  )
  obj <- TMB::MakeADFun(
    data = data_tmb,
    parameters = pars,
    random = "u",
    DLL = "bb_mm_prior",
    silent = silent
  )
  init_fn <- function() {
    list(
      beta = beta0 + stats::rnorm(length(beta0), 0, 0.05),
      log_phi = log(pmax(phi0, 1e-6)) + stats::rnorm(length(phi0), 0, 0.05),
      theta_re = th0 + stats::rnorm(length(th0), 0, 0.05),
      u = u0
    )
  }
  if (isTRUE(laplace_bayes)) {
    init_fn <- function() {
      list(
        beta = beta0 + stats::rnorm(length(beta0), 0, 0.05),
        log_phi = log(pmax(phi0, 1e-6)) + stats::rnorm(length(phi0), 0, 0.05),
        theta_re = th0 + stats::rnorm(length(th0), 0, 0.05)
      )
    }
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
    laplace = isTRUE(laplace_bayes),
    control = list(adapt_delta = 0.95, max_treedepth = 12)
  )
  fit$time_bayes <- proc.time()[["elapsed"]] - t1
  fit$stanfit <- stanfit

  phi_names <- if (fit$nDim == 1L) "phi" else paste0("phi", seq_len(fit$nDim))
  fit$posterior <- .posterior_from_stanfit(
    stanfit,
    beta_names = names(fit$fixed.coef),
    phi_names = phi_names,
    sigma_names = names(fit$sigma.coef)
  )
  fit
}

#' @export
print.BBmm <- function(x, ...) {
  cat("Call:\t")
  print(x$call)
  .cat_model_block(
    "Model (BBmm):",
    if (!is.null(x$model)) x$model else .bb_model_lines_mm(x$nDim, x$namesRand)
  )
  if (!is.null(x$fixed.formula) || !is.null(x$random.formula)) {
    cat("Formulas: fixed = ")
    print(if (!is.null(x$fixed.formula)) x$fixed.formula else NA)
    cat("          random = ")
    print(if (!is.null(x$random.formula)) x$random.formula else NA)
    cat("\n")
  }
  cat("Method:", if (!is.null(x$method)) x$method else "mle")
  if (identical(x$method, "bayes")) cat(" (point estimates below are MLE / Laplace)")
  cat("\n\n")
  cat("beta (fixed effects, logit scale):\n")
  print(if (!is.null(x$beta)) x$beta else x$fixed.coef)

  cat("\nRandom-effect SDs (sigma):\n")
  sig <- if (!is.null(x$sigma)) x$sigma else x$sigma.coef
  for (i in seq_along(sig)) {
    nm <- if (!is.null(names(sig))) names(sig)[i] else as.character(i)
    cat("  ", nm, ": ", sig[i], "\n", sep = "")
  }
  if (!is.null(x$Corr) && length(x$Corr)) {
    for (nm in names(x$Corr)) {
      C <- x$Corr[[nm]]
      if (is.matrix(C) && nrow(C) > 1L) {
        cat("\nCorr (", nm, "):\n", sep = "")
        print(round(C, 3))
      }
    }
  }
  cat("\nphi (BB dispersion):", paste(x$phi, collapse = ", "), "\n")
  if (identical(x$method, "bayes") && !is.null(x$posterior)) {
    cat("\nPosterior summary (tmbstan):\n")
    print(x$posterior, row.names = FALSE, digits = 4)
    if (is.finite(x$time_bayes)) {
      cat(sprintf("Bayes sampling time: %.2fs\n", x$time_bayes))
    }
  }
  cat("\nDeviance:", x$deviance)
  cat("\nNumber of iterations:", x$iter)
  if (x$balanced == "yes") {
    cat("\nBalanced data, m =", x$m, "\n")
  } else {
    cat("\nUnbalanced m.\n")
  }
  invisible(x)
}

#' @export
summary.BBmm <- function(object, ...) {
  beta <- if (!is.null(object$beta)) object$beta else object$fixed.coef
  fixed.se <- sqrt(diag(object$fixed.vcov))
  fixed.tval <- as.numeric(beta) / fixed.se
  fixed.TAB <- cbind(
    Estimate = as.numeric(beta),
    StdErr = fixed.se,
    t.value = fixed.tval,
    p.value = 2 * pnorm(-abs(fixed.tval))
  )
  rownames(fixed.TAB) <- names(beta)

  log_phi <- object$log_phi
  if (is.null(log_phi)) log_phi <- object$psi.coef
  phi <- object$phi
  if (is.null(phi)) phi <- object$phi.coef
  psi.se <- sqrt(object$psi.var)
  if (length(as.numeric(log_phi)) == 1L) {
    psi.table <- cbind(Estimate = as.numeric(log_phi), StdErr = as.numeric(psi.se)[1L])
    rownames(psi.table) <- "log(phi)"
    phi.table <- cbind(
      Estimate = c(phi = as.numeric(phi)[1L], log_phi = as.numeric(log_phi)[1L]),
      StdErr = c(as.numeric(phi)[1L] * as.numeric(psi.se)[1L], as.numeric(psi.se)[1L])
    )
  } else {
    psi.table <- cbind(Estimate = as.numeric(log_phi), StdErr = as.numeric(psi.se))
    rownames(psi.table) <- names(log_phi)
    phi.table <- cbind(
      Estimate = as.numeric(phi),
      StdErr = as.numeric(phi) * as.numeric(psi.se)
    )
    rownames(phi.table) <- paste0("phi", seq_along(phi))
  }

  sigma <- if (!is.null(object$sigma)) object$sigma else object$sigma.coef
  sigma.table <- cbind(
    Estimate = as.numeric(sigma),
    StdErr = sqrt(object$sigma.var)
  )
  rownames(sigma.table) <- names(sigma)

  Chi <- object$null.deviance - object$deviance
  Chi.p.value <- 1 - pchisq(Chi, object$null.df - object$df)

  res <- list(
    call = object$call,
    model = if (!is.null(object$model)) object$model else .bb_model_lines_mm(object$nDim, object$namesRand),
    fixed.coefficients = fixed.TAB,
    beta.table = fixed.TAB,
    sigma.table = sigma.table,
    Sigma = object$Sigma,
    Corr = object$Corr,
    phi.table = phi.table,
    psi.table = psi.table,
    random.coef = if (!is.null(object$u)) object$u else object$random.coef,
    u = if (!is.null(object$u)) object$u else object$random.coef,
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
    nll = object$nll,
    method = object$method,
    fixed.formula = object$fixed.formula,
    random.formula = object$random.formula
  )
  class(res) <- "summary.BBmm"
  res
}

#' @export
print.summary.BBmm <- function(x, ...) {
  cat("Call:\t")
  print(x$call)
  .cat_model_block("Model (BBmm):", x$model)
  cat("beta (fixed effects):\n\n")
  printCoefmat(x$beta.table, P.values = TRUE, has.Pvalue = TRUE)
  cat("\n---------------------------------------------------------------\n")
  cat("sigma (RE SD):\n\n")
  print(x$sigma.table)
  if (!is.null(x$Corr) && length(x$Corr)) {
    for (nm in names(x$Corr)) {
      C <- x$Corr[[nm]]
      if (is.matrix(C) && nrow(C) > 1L) {
        cat("\nCorr (", nm, "):\n", sep = "")
        print(round(C, 3))
      }
    }
  }
  cat("\n---------------------------------------------------------------\n")
  cat("phi (BB dispersion) and log(phi):\n\n")
  print(x$phi.table)
  cat("\n---------------------------------------------------------------\n")
  cat("Deviance:", x$deviance, "; df =", x$df, "\n")
  cat("Null deviance:", x$null.deviance, "; df =", x$null.df, "\n")
  cat("Deviance goodness-of-fit p-value:", x$Goodness.of.fit, "\n")
  if (!is.null(x$nll)) cat("Laplace approx. nll:", x$nll, "\n")
  cat("\nNumber of observations:", x$nObs)
  cat("\nNumber of iterations:", x$iter)
  if (x$balanced == "yes") {
    cat("\nBalanced data, m =", x$m)
  } else {
    cat("\nUnbalanced m.")
  }
  cat("\nRE counts per component (nRandComp):", paste(x$nRandComp, collapse = ", "), "\n\n")
  invisible(x)
}

#' @export
coef.BBmm <- function(object, ...) {
  if (!is.null(object$beta)) return(object$beta)
  object$fixed.coef
}

#' @export
vcov.BBmm <- function(object, ...) object$fixed.vcov

#' @export
fitted.BBmm <- function(object, ...) {
  if (!is.null(object$p)) return(object$p)
  object$fitted.values
}
