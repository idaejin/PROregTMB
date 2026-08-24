## Eilers–Marx P-splines via mixed-model reparameterisation
## (Eilers 1999; same engine as BB-GAM / Report 4):
##   f(x) = X_null beta_null + Z s,   s ~ N(0, sigma_s^2 I)
##   Z = B D' (D D')^{-1}, residualised against the fixed design.
## Difference from the former GMRF-on-B path: the polynomial null space
## lives in X (once), not inside a penalised GMRF.

#' Truncated power (Eilers & Marx)
#' @keywords internal
.tpower <- function(x, t, p) (x - t)^p * (x > t)

#' B-spline via truncated powers (BB-GAM notes)
#' @keywords internal
.bspline_eilers <- function(x, xl, xr, ndx, bdeg) {
  dx <- (xr - xl) / ndx
  knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
  P <- outer(x, knots, .tpower, bdeg)
  n <- ncol(P)
  Dd <- diff(diag(n), diff = bdeg + 1) / (gamma(bdeg + 1) * dx^bdeg)
  B <- (-1)^(bdeg + 1) * P %*% t(Dd)
  list(B = B, knots = knots, dx = dx, xl = xl, xr = xr, ndx = ndx, bdeg = bdeg)
}

#' @keywords internal
.diff_penalty_eilers <- function(n_basis, pord = 2L) {
  diff(diag(as.integer(n_basis)), differences = as.integer(pord))
}

#' Eilers Z = B D'(DD')^{-1}
#' @keywords internal
.eilers_Z <- function(B, D) {
  Z <- B %*% t(D) %*% solve(tcrossprod(D))
  colnames(Z) <- paste0("z", seq_len(ncol(Z)))
  Z
}

#' Residualise Z against X (X'Z ≈ 0)
#' @keywords internal
.orth_to_X <- function(Z, X) {
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  if (!ncol(Z) || !nrow(X)) return(Z)
  Q <- qr.Q(qr(X))
  Z - Q %*% crossprod(Q, Z)
}

#' Drop near-zero / dependent columns
#' @keywords internal
.reduce_Z_rank <- function(Z, tol = 1e-8) {
  Z <- as.matrix(Z)
  if (!ncol(Z)) return(Z)
  nrm <- sqrt(colSums(Z * Z))
  keep <- nrm > tol * max(nrm, 1)
  Z <- Z[, keep, drop = FALSE]
  if (!ncol(Z)) return(Z)
  qrz <- qr(Z, tol = tol)
  Z[, qrz$pivot[seq_len(qrz$rank)], drop = FALSE]
}

#' B-spline basis (Eilers & Marx, equally spaced knots) — public helper
#'
#' Uses \pkg{splines} `splineDesign` (compatible with earlier PROregTMB API).
#' The model engine uses the truncated-power construction internally.
#'
#' @inheritParams bbase
#' @export
bbase <- function(x, xl = min(x, na.rm = TRUE), xr = max(x, na.rm = TRUE),
                  ndx = 10L, bdeg = 3L, nseg = NULL) {
  x <- as.numeric(x)
  if (!is.null(nseg)) ndx <- nseg
  ndx <- as.integer(ndx)
  bdeg <- as.integer(bdeg)
  if (ndx < 1L) stop("ndx must be >= 1", call. = FALSE)
  if (!is.finite(xl) || !is.finite(xr) || xr <= xl) {
    stop("need finite xl < xr for B-spline basis", call. = FALSE)
  }
  dx <- (xr - xl) / ndx
  knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
  B <- splines::splineDesign(knots, x, ord = bdeg + 1L, outer.ok = TRUE)
  attr(B, "ndx") <- ndx
  attr(B, "nseg") <- ndx
  attr(B, "bdeg") <- bdeg
  attr(B, "xl") <- xl
  attr(B, "xr") <- xr
  attr(B, "knots") <- knots
  B
}

#' Difference penalty matrix S = D'D (Eilers & Marx)
#' @export
penalty_diff <- function(K, pord = 2L, d = NULL) {
  K <- as.integer(K)
  if (!is.null(d)) pord <- d
  pord <- as.integer(pord)
  if (K <= pord) stop("K must be > pord", call. = FALSE)
  D <- diff(diag(K), differences = pord)
  crossprod(D)
}

#' Rank-1 sum-to-zero centering penalty (legacy helper; unused by Eilers engine)
#' @export
penalty_center <- function(B) {
  B <- as.matrix(B)
  Bt1 <- crossprod(B, rep(1, nrow(B)))
  tcrossprod(Bt1)
}

#' P-spline smooth term (Eilers notation)
#'
#' Engine: mixed-model map \eqn{Z = B D'(DD')^{-1}} with polynomial null
#' space in the fixed design (Eilers 1999), matching BB-GAM / Report 4.
#'
#' @param x Numeric covariate.
#' @param ndx Number of equal intervals (or set `k` = basis dimension).
#' @param pord Difference penalty order (default 2).
#' @param bdeg B-spline degree (default 3).
#' @param xl,xr Optional knot boundaries.
#' @param by Optional grouping factor (domain-specific smooths).
#' @param k Basis dimension; if `ndx` missing, `ndx = k - bdeg`.
#' @param ... Ignored.
#' @export
s <- function(x, ndx = 10L, pord = 2L, bdeg = 3L,
              xl = NULL, xr = NULL, by = NULL, k = NULL, ...) {
  as.numeric(x)
}

#' @keywords internal
.parse_smooth_formula <- function(formula) {
  if (!inherits(formula, "formula")) {
    stop("formula must be a formula", call. = FALSE)
  }
  tt <- stats::terms(formula, specials = "s")
  sp <- attr(tt, "specials")$s
  resp <- if (length(formula) == 3L) formula[[2L]] else NULL
  labs <- attr(tt, "term.labels")
  if (is.null(sp)) {
    return(list(fixed = formula, s_calls = list(), s_labels = character()))
  }
  vars <- as.list(attr(tt, "variables"))[-1L]
  s_calls <- lapply(sp, function(i) vars[[i]])
  is_s <- vapply(labs, function(lab) grepl("^s\\s*\\(", lab), logical(1))
  keep <- labs[!is_s]
  if (length(keep)) {
    fixed <- stats::reformulate(keep, response = resp)
  } else if (!is.null(resp)) {
    fixed <- stats::reformulate("1", response = resp)
  } else {
    fixed <- ~ 1
  }
  environment(fixed) <- environment(formula)
  list(fixed = fixed, s_calls = s_calls, s_labels = labs[is_s])
}

#' @keywords internal
.parametric_rhs <- function(formula, drop = character()) {
  parts <- .parse_smooth_formula(formula)
  tl <- attr(stats::terms(parts$fixed), "term.labels")
  tl <- setdiff(tl, drop)
  if (!length(tl)) return(~ 1)
  stats::reformulate(tl)
}

#' Parse one s() call (covariate + args); no basis yet
#' @keywords internal
.parse_s_call <- function(cl, data, envir = parent.frame()) {
  if (!is.call(cl) || !identical(cl[[1L]], as.name("s"))) {
    stop("expected a call to s()", call. = FALSE)
  }
  args <- as.list(cl)[-1L]
  if (!length(args)) stop("s() needs a covariate", call. = FALSE)
  x_expr <- args[[1L]]
  xname <- paste(deparse(x_expr, width.cutoff = 500L), collapse = "")
  x <- as.numeric(eval(x_expr, envir = data, enclos = envir))
  if (anyNA(x)) stop("NA in smooth covariate: ", xname, call. = FALSE)

  get_arg <- function(name, default) {
    if (name %in% names(args)) eval(args[[name]], envir = data, enclos = envir)
    else default
  }
  bdeg <- as.integer(get_arg("bdeg", 3L))
  ndx <- get_arg("ndx", NULL)
  if (is.null(ndx)) ndx <- get_arg("nseg", NULL)
  if (is.null(ndx)) {
    k_basis <- get_arg("k", NULL)
    if (!is.null(k_basis)) {
      ndx <- as.integer(k_basis) - bdeg
      if (ndx < 1L) stop("k must be > bdeg", call. = FALSE)
    } else ndx <- 10L
  }
  pord <- get_arg("pord", NULL)
  if (is.null(pord)) pord <- get_arg("m", 2L)
  if (is.null(pord)) pord <- get_arg("d", 2L)
  xl <- get_arg("xl", NULL)
  xr <- get_arg("xr", NULL)
  eps <- as.numeric(get_arg("eps", 0.01))
  if (is.null(xl)) xl <- min(x) - eps
  if (is.null(xr)) xr <- max(x) + eps
  if (xr <= xl) xr <- xl + 1e-6

  by_expr <- if ("by" %in% names(args)) args[["by"]] else NULL
  by_fac <- NULL
  by_name <- NA_character_
  if (!is.null(by_expr)) {
    by_fac <- droplevels(as.factor(eval(by_expr, envir = data, enclos = envir)))
    if (length(by_fac) != length(x)) {
      stop("by= length must match covariate", call. = FALSE)
    }
    by_name <- paste(deparse(by_expr, width.cutoff = 500L), collapse = "")
  }

  list(
    x = x, xname = xname, x_expr = x_expr,
    ndx = as.integer(ndx), pord = as.integer(pord), bdeg = bdeg,
    xl = xl, xr = xr, by_fac = by_fac, by_name = by_name
  )
}

#' Build Eilers mixed smooth design from a formula
#'
#' @return list with `fixed`, `X_null` (linear null pieces not already named
#'   in a later merge), `Zs`, `s_comp`, `n_smooth`, `labels`, `specs`,
#'   `engine = "eilers"`.
#' @keywords internal
.build_smooth_design <- function(formula, data, envir = parent.frame(),
                                 orth_Z = TRUE) {
  parts <- .parse_smooth_formula(formula)
  if (is.null(data) || (is.list(data) && !is.data.frame(data) && !length(data))) {
    av <- all.vars(formula)
    data <- as.data.frame(mget(av, envir = envir, inherits = TRUE),
                          optional = TRUE)
  } else {
    data <- as.data.frame(data)
  }
  n <- nrow(data)
  empty <- list(
    fixed = parts$fixed,
    n_smooth = 0L,
    engine = "eilers",
    X_null = matrix(0, n, 0L),
    Zs = matrix(0, n, 0L),
    s_comp = integer(0),
    smooth_K = integer(0),
    labels = character(),
    specs = list(),
    blocks_idx = list(),
    # legacy aliases so old checks don't explode
    B = matrix(0, n, 0L),
    S = matrix(0, 0L, 0L),
    C = matrix(0, 0L, 0L)
  )
  if (!length(parts$s_calls)) return(empty)

  parsed <- lapply(parts$s_calls, .parse_s_call, data = data, envir = envir)

  # Null-space linear terms (pord = 2): one column per covariate for shared
  # s(x); for s(x, by = g), one column per level, named "level.x" with
  # values x * I(g == level). Domain-specific by= therefore matches a
  # fully domain-specific mixed P-spline (linear null + wiggly).
  X_null_cols <- list()
  for (j in seq_along(parsed)) {
    pj <- parsed[[j]]
    if (is.null(pj$by_fac)) {
      nm <- pj$xname
      if (!nm %in% names(X_null_cols)) {
        X_null_cols[[nm]] <- matrix(pj$x, ncol = 1L, dimnames = list(NULL, nm))
      }
    } else {
      for (lev in levels(pj$by_fac)) {
        nm <- paste0(lev, ".", pj$xname)
        w <- as.numeric(pj$by_fac == lev)
        X_null_cols[[nm]] <- matrix(
          pj$x * w, ncol = 1L, dimnames = list(NULL, nm)
        )
      }
    }
  }
  X_null <- if (length(X_null_cols)) {
    do.call(cbind, X_null_cols)
  } else {
    matrix(0, n, 0L)
  }

  # Temporary X for orthogonalisation: intercept + null + will merge later
  # Callers should re-orth against full model X; we orth against null+1 here
  # and BBreg/BBmm re-orth against cbind(X_param, X_null) after merge.
  X_tmp <- cbind("(Intercept)" = rep(1, n), X_null)

  specs <- list()
  Z_list <- list()
  labels <- character()
  s_comp <- integer()
  comp_id <- 0L

  for (j in seq_along(parsed)) {
    pj <- parsed[[j]]
    bb <- .bspline_eilers(pj$x, pj$xl, pj$xr, pj$ndx, pj$bdeg)
    D <- .diff_penalty_eilers(ncol(bb$B), pj$pord)
    Z_raw <- .eilers_Z(bb$B, D)

    if (is.null(pj$by_fac)) {
      lab <- paste0("s(", pj$xname, ")")
      Zj <- .orth_to_X(Z_raw, X_tmp)
      Zj <- .reduce_Z_rank(Zj)
      if (!ncol(Zj)) next
      colnames(Zj) <- paste0(lab, ".", seq_len(ncol(Zj)))
      Z_list[[length(Z_list) + 1L]] <- Zj
      labels <- c(labels, lab)
      s_comp <- c(s_comp, rep(comp_id, ncol(Zj)))
      specs[[length(specs) + 1L]] <- list(
        label = lab, xname = pj$xname, x = pj$x,
        by_level = NA_character_, by_var = NA_character_,
        ndx = pj$ndx, pord = pj$pord, bdeg = pj$bdeg,
        xl = pj$xl, xr = pj$xr, D = D, B = bb$B,
        K = ncol(Zj), null_name = pj$xname
      )
      comp_id <- comp_id + 1L
    } else {
      for (lev in levels(pj$by_fac)) {
        lab <- paste0("s(", pj$xname, "):", lev)
        w <- as.numeric(pj$by_fac == lev)
        Zj <- Z_raw * w
        Zj <- .orth_to_X(Zj, X_tmp)
        Zj <- .reduce_Z_rank(Zj)
        if (!ncol(Zj)) next
        colnames(Zj) <- paste0(lab, ".", seq_len(ncol(Zj)))
        Z_list[[length(Z_list) + 1L]] <- Zj
        labels <- c(labels, lab)
        s_comp <- c(s_comp, rep(comp_id, ncol(Zj)))
        specs[[length(specs) + 1L]] <- list(
          label = lab, xname = pj$xname, x = pj$x[w > 0],
          by_level = lev, by_var = pj$by_name,
          ndx = pj$ndx, pord = pj$pord, bdeg = pj$bdeg,
          xl = pj$xl, xr = pj$xr, D = D, B = bb$B,
          K = ncol(Zj),
          null_name = paste0(lev, ".", pj$xname),
          by_mask = w
        )
        comp_id <- comp_id + 1L
      }
    }
  }

  if (!length(Z_list)) {
    empty$fixed <- parts$fixed
    empty$X_null <- X_null
    return(empty)
  }

  Zs <- do.call(cbind, Z_list)
  # Joint rank reduction may drop columns; rebuild s_comp carefully
  # (skip joint drop to keep s_comp aligned — orth already done)

  Ks <- vapply(specs, `[[`, integer(1), "K")
  off <- if (length(Ks) <= 1L) {
    if (length(Ks) == 1L) 0L else integer(0)
  } else {
    c(0L, cumsum(Ks[-length(Ks)]))
  }
  blocks_idx <- lapply(seq_along(Ks), function(j) off[j] + seq_len(Ks[j]))

  list(
    fixed = parts$fixed,
    n_smooth = length(specs),
    engine = "eilers",
    X_null = X_null,
    Zs = Zs,
    s_comp = as.integer(s_comp),
    smooth_K = as.integer(Ks),
    smooth_off = as.integer(off),
    labels = labels,
    specs = specs,
    blocks_idx = blocks_idx,
    B = Zs, # alias used by some nrow checks
    S = matrix(0, 0L, 0L),
    C = matrix(0, 0L, 0L)
  )
}

#' Re-orthogonalise Zs against the full fixed design X
#' @keywords internal
.reorth_smooth_Zs <- function(sm, X) {
  if (is.null(sm) || !isTRUE(sm$n_smooth > 0L)) return(sm)
  Zs <- .orth_to_X(sm$Zs, X)
  # rebuild by block to keep labels, dropping empty blocks
  new_Z <- list()
  new_specs <- list()
  new_labels <- character()
  new_comp <- integer()
  comp <- 0L
  for (j in seq_len(sm$n_smooth)) {
    idx <- sm$blocks_idx[[j]]
    Zj <- .reduce_Z_rank(Zs[, idx, drop = FALSE])
    if (!ncol(Zj)) next
    colnames(Zj) <- paste0(sm$labels[j], ".", seq_len(ncol(Zj)))
    new_Z[[length(new_Z) + 1L]] <- Zj
    sp <- sm$specs[[j]]
    sp$K <- ncol(Zj)
    new_specs[[length(new_specs) + 1L]] <- sp
    new_labels <- c(new_labels, sm$labels[j])
    new_comp <- c(new_comp, rep(comp, ncol(Zj)))
    comp <- comp + 1L
  }
  if (!length(new_Z)) {
    sm$n_smooth <- 0L
    sm$Zs <- matrix(0, nrow(X), 0L)
    sm$B <- sm$Zs
    sm$s_comp <- integer(0)
    sm$labels <- character()
    sm$specs <- list()
    sm$blocks_idx <- list()
    sm$smooth_K <- integer(0)
    return(sm)
  }
  Zs2 <- do.call(cbind, new_Z)
  Ks <- vapply(new_specs, `[[`, integer(1), "K")
  off <- if (length(Ks) <= 1L) 0L else c(0L, cumsum(Ks[-length(Ks)]))
  if (length(Ks) == 0L) off <- integer(0)
  sm$Zs <- Zs2
  sm$B <- Zs2
  sm$n_smooth <- length(new_specs)
  sm$specs <- new_specs
  sm$labels <- new_labels
  sm$s_comp <- as.integer(new_comp)
  sm$smooth_K <- as.integer(Ks)
  sm$smooth_off <- as.integer(off)
  sm$blocks_idx <- lapply(seq_along(Ks), function(j) off[j] + seq_len(Ks[j]))
  sm
}

#' Merge X_null columns into parametric X (skip name collisions)
#' @keywords internal
.merge_X_null <- function(X, X_null) {
  if (is.null(X_null) || !ncol(X_null)) return(X)
  cn <- colnames(X)
  add <- setdiff(colnames(X_null), cn)
  if (!length(add)) return(X)
  # Skip a shared null name (e.g. "fev_w") when domain-prefixed versions
  # already exist in X (e.g. "Impacts.fev_w"); do not drop domain-prefixed
  # X_null columns that are the intended by= null space.
  if (length(cn)) {
    bare <- sub("^.*\\.", "", cn)
    shared_add <- add[!grepl("\\.", add)]
    drop_shared <- shared_add[shared_add %in% bare]
    add <- setdiff(add, drop_shared)
  }
  if (!length(add)) return(X)
  cbind(X, X_null[, add, drop = FALSE])
}

#' Stack Eilers smooth rows for wide multivariate (shared smooth)
#' @keywords internal
.rep_smooth_rows <- function(sm, L) {
  L <- as.integer(L)
  if (is.null(sm) || !isTRUE(sm$n_smooth > 0L) || L <= 1L) return(sm)
  sm$Zs <- do.call(rbind, replicate(L, sm$Zs, simplify = FALSE))
  sm$B <- sm$Zs
  if (!is.null(sm$X_null) && ncol(sm$X_null)) {
    sm$X_null <- do.call(rbind, replicate(L, sm$X_null, simplify = FALSE))
  }
  sm
}

#' TMB data block for Eilers smooths
#' @keywords internal
.tmb_smooth_data <- function(sm, n) {
  if (is.null(sm) || sm$n_smooth < 1L) {
    return(list(
      n_smooth = 0L,
      Zs = matrix(0, n, 0L),
      s_comp = integer(0)
    ))
  }
  list(
    n_smooth = as.integer(sm$n_smooth),
    Zs = as.matrix(sm$Zs),
    s_comp = as.integer(sm$s_comp)
  )
}

#' TMB parameters for Eilers smooth block
#' @keywords internal
.tmb_smooth_parameters <- function(sm) {
  if (is.null(sm) || sm$n_smooth < 1L) {
    return(list(log_sds = numeric(0), s = numeric(0)))
  }
  list(
    log_sds = rep(0, sm$n_smooth),
    s = rep(0, ncol(sm$Zs))
  )
}

#' Predict an additive smooth (Eilers mixed; wiggly + null if available)
#'
#' Pointwise bands use the marginal Laplace covariance from TMB
#' \code{jointPrecision}. With \code{centered = TRUE} (recommended for
#' within WW--BW displays), the target is the contrast
#' \eqn{f(x)-f(0)} and the SE uses the joint covariance of the null-space
#' linear coefficient and the wiggly coefficients \eqn{s} for that smooth
#' (fixed--random cross terms included; hyperparameters held at their MLE).
#'
#' @param object A `BBreg` or `BBmm` fit with smooths.
#' @param which Integer index or label of the smooth.
#' @param x Optional grid; default uses observed covariate.
#' @param level Confidence level (default 0.95).
#' @param se If `TRUE` (default), attach `se`, `lwr`, `upr`.
#' @param centered If `TRUE`, return \eqn{f(x)-f(0)} with joint SE; if `FALSE`,
#'   return the uncentered smooth (wiggly SE only, legacy).
#' @param method Kept for API compatibility (`"marrawood"`).
#' @return Data frame with `x`, `fit`, and usually `se`, `lwr`, `upr`.
#' @references
#' Eilers, P. H. C. (1999). Discussion of Currie & Durban. JRSS-C.
#' Marra, G. and Wood, S. N. (2012). Scand J Statist 39, 53–74.
#' @export
predict_smooth <- function(object, which = 1L, x = NULL, level = 0.95,
                           se = TRUE, centered = FALSE,
                           method = c("marrawood")) {
  method <- match.arg(method)
  sm <- object$smooth
  if (is.null(sm) || is.null(sm$n_smooth) || sm$n_smooth < 1L) {
    stop("object has no smooth terms", call. = FALSE)
  }
  if (is.character(which)) {
    j <- match(which, sm$labels)
    if (is.na(j)) stop("unknown smooth: ", which, call. = FALSE)
  } else {
    j <- as.integer(which)
  }
  if (j < 1L || j > sm$n_smooth) stop("which out of range", call. = FALSE)

  spec <- sm$specs[[j]]
  idx <- sm$blocks_idx[[j]]
  s_hat <- as.numeric(object$s)[idx]
  if (is.null(x)) {
    x <- spec$x
  } else {
    x <- as.numeric(x)
  }

  bb <- .bspline_eilers(x, spec$xl, spec$xr, spec$ndx, spec$bdeg)
  Zg <- .eilers_Z(bb$B, spec$D)
  if (ncol(Zg) > length(s_hat)) {
    Zg <- Zg[, seq_len(length(s_hat)), drop = FALSE]
  } else if (ncol(Zg) < length(s_hat)) {
    Zg <- cbind(Zg, matrix(0, nrow(Zg), length(s_hat) - ncol(Zg)))
  }
  f_w <- as.numeric(Zg %*% s_hat)

  # Null linear coefficient (shared or domain-prefixed)
  b_null <- 0
  bn <- spec$null_name
  beta <- object$beta
  null_nm <- NA_character_
  if (!is.null(beta) && !is.null(bn)) {
    cand <- c(
      bn,
      if (!is.null(spec$by_level) && !is.na(spec$by_level))
        paste0(spec$by_level, ".", bn)
    )
    hit <- cand[cand %in% names(beta)]
    if (length(hit)) {
      null_nm <- hit[1]
      b_null <- unname(beta[[null_nm]])
    }
  }
  f_null <- b_null * x
  f <- f_null + f_w

  if (isTRUE(centered)) {
    Z0 <- .eilers_Z(
      .bspline_eilers(0, spec$xl, spec$xr, spec$ndx, spec$bdeg)$B,
      spec$D
    )
    if (ncol(Z0) > length(s_hat)) {
      Z0 <- Z0[, seq_len(length(s_hat)), drop = FALSE]
    } else if (ncol(Z0) < length(s_hat)) {
      Z0 <- cbind(Z0, matrix(0, 1L, length(s_hat) - ncol(Z0)))
    }
    f0 <- b_null * 0 + as.numeric(Z0 %*% s_hat)
    f <- f - f0
    Zd <- Zg - matrix(as.numeric(Z0), nrow(Zg), ncol(Zg), byrow = TRUE)
  } else {
    Zd <- Zg
  }

  out <- data.frame(x = x, fit = f)
  if (!isTRUE(se)) return(out)

  se_hat <- rep(NA_real_, length(x))
  se_method <- "none"

  if (isTRUE(centered)) {
    meta <- .smooth_contrast_jp_index(object, j)
    V <- NULL
    if (!is.null(meta)) {
      V <- .marginal_cov_from_joint(object$sdreport$jointPrecision, meta$idx)
    }
    if (!is.null(V) && all(is.finite(V))) {
      if (isTRUE(meta$has_null)) {
        C <- cbind(x, Zd)
      } else {
        C <- Zd
      }
      if (ncol(C) == ncol(V)) {
        se_hat <- sqrt(pmax(0, rowSums((C %*% V) * C)))
        se_method <- "joint_null_s"
      }
    }
  }

  if (se_method == "none") {
    # Legacy / fallback: wiggly block only
    Vfull <- object$s.vcov
    if (is.null(Vfull) && !is.null(object$sdreport)) {
      n_u <- if (inherits(object, "BBmm")) object$nRand else 0L
      Vfull <- .s_vcov_from_joint(
        object$sdreport$jointPrecision,
        n_s = length(object$s),
        n_u = n_u
      )
    }
    if (!is.null(Vfull)) {
      Vj <- Vfull[idx, idx, drop = FALSE]
      if (all(is.finite(Vj))) {
        se_hat <- as.numeric(smooth_se_cpp(Zd, Vj))
        se_method <- if (isTRUE(centered)) "s_only_centered" else "marrawood"
        # Partial upgrade: add null variance without cross term if available
        if (isTRUE(centered) && !is.na(null_nm) && !is.null(object$fixed.vcov) &&
            null_nm %in% rownames(object$fixed.vcov)) {
          se_hat <- sqrt(pmax(0, se_hat^2 + (x^2) * object$fixed.vcov[null_nm, null_nm]))
          se_method <- "s_plus_null_var"
        }
      }
    }
  }

  z <- stats::qnorm(1 - (1 - level) / 2)
  out$se <- se_hat
  out$lwr <- f - z * se_hat
  out$upr <- f + z * se_hat
  attr(out, "level") <- level
  attr(out, "which") <- j
  attr(out, "label") <- sm$labels[j]
  attr(out, "method") <- se_method
  attr(out, "centered") <- isTRUE(centered)
  out
}

#' Plot an additive P-spline with pointwise bands
#' @inheritParams predict_smooth
#' @param col,shade,add,xlab,ylab,main,ylim,... Plot options.
#' @export
plot_smooth <- function(object, which = "all", x = NULL, level = 0.95,
                        col = "darkorange2", shade = grDevices::adjustcolor(col, 0.30),
                        add = FALSE, xlab = NULL, ylab = NULL, main = NULL,
                        ylim = NULL, ...) {
  sm <- object$smooth
  if (is.null(sm) || sm$n_smooth < 1L) {
    stop("object has no smooth terms", call. = FALSE)
  }
  if (identical(which, "all")) which <- seq_len(sm$n_smooth)
  if (length(which) > 1L) {
    if (isTRUE(add)) stop("add = TRUE requires a single smooth", call. = FALSE)
    n <- length(which)
    nc <- if (n > 3L) ceiling(sqrt(n)) else n
    nr <- ceiling(n / nc)
    op <- graphics::par(mfrow = c(nr, nc), mar = c(4, 4, 2.5, 1))
    on.exit(graphics::par(op), add = TRUE)
    out <- lapply(which, function(w) {
      plot_smooth(object, which = w, x = x, level = level, col = col,
                  shade = shade, add = FALSE, xlab = xlab, ylab = ylab,
                  main = main, ylim = ylim, ...)
    })
    names(out) <- if (is.numeric(which)) sm$labels[which] else which
    return(invisible(out))
  }
  pr <- predict_smooth(object, which = which, x = x, level = level, se = TRUE)
  o <- order(pr$x)
  pr <- pr[o, , drop = FALSE]
  lab <- attr(pr, "label")
  if (is.null(lab)) lab <- paste0("smooth ", which)
  if (is.null(xlab)) xlab <- "x"
  if (is.null(ylab)) ylab <- "f(x)"
  if (is.null(main)) {
    main <- paste0(lab, " (", round(100 * level), "% CI)")
  }
  if (is.null(ylim)) ylim <- range(pr$fit, pr$lwr, pr$upr, na.rm = TRUE)
  if (!isTRUE(add)) {
    graphics::plot(pr$x, pr$fit, type = "n", xlab = xlab, ylab = ylab,
                   main = main, ylim = ylim, ...)
    graphics::abline(h = 0, col = "grey70")
  }
  if (!is.null(pr$lwr) && !is.null(pr$upr)) {
    graphics::polygon(c(pr$x, rev(pr$x)), c(pr$lwr, rev(pr$upr)),
                      col = shade, border = NA)
  }
  graphics::lines(pr$x, pr$fit, col = col, lwd = 2, ...)
  invisible(pr)
}
