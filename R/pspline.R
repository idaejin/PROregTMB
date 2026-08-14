## Additive P-splines for BBreg / BBmm (Eilers & Marx + sum-to-zero)
##
## Model (additive on the logit scale):
##   logit(p) = X beta + sum_j f_j(x_j),   f_j = B_j gamma_j
##
## Penalty per smooth (as in pspline_additive.R / _ci.R):
##   lambda_j * || D_{pord} gamma_j ||^2
##   + kappa * || 1' B_j gamma_j ||^2     # sum_i f_j(x_ij) = 0
##
## In TMB: gamma_j ~ GMRF(Q_j),  Q_j = lambda_j S_j + kappa C_j,
## with S_j = D'D and C_j = (B_j'1)(1'B_j); lambda_j via Laplace.

#' B-spline basis (Eilers & Marx, equally spaced knots)
#'
#' @param x Covariate values.
#' @param xl,xr Left/right boundaries (default range of `x`).
#' @param ndx Number of equal-width intervals on \eqn{[xl, xr]} (Eilers).
#' @param bdeg B-spline degree (default 3 = cubic).
#' @param nseg Deprecated alias of `ndx`.
#' @return Basis matrix with attributes `ndx`, `bdeg`, `xl`, `xr`, `knots`.
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
#'
#' @param K Basis dimension.
#' @param pord Difference / penalty order (default 2).
#' @param d Deprecated alias of `pord`.
#' @return Symmetric `K x K` penalty matrix.
#' @export
penalty_diff <- function(K, pord = 2L, d = NULL) {
  K <- as.integer(K)
  if (!is.null(d)) pord <- d
  pord <- as.integer(pord)
  if (K <= pord) stop("K must be > pord", call. = FALSE)
  D <- diff(diag(K), differences = pord)
  crossprod(D)
}

#' Rank-1 sum-to-zero centering penalty C = (B'1)(1'B)
#'
#' Enforces \eqn{1^\top B\gamma \approx 0} (mean of the smooth at the
#' observed design is zero), as in the additive P-spline notes.
#'
#' @param B Basis matrix (n x K).
#' @return Symmetric `K x K` matrix.
#' @export
penalty_center <- function(B) {
  B <- as.matrix(B)
  Bt1 <- crossprod(B, rep(1, nrow(B)))
  tcrossprod(Bt1)
}

#' P-spline smooth term (Eilers notation)
#'
#' @param x Numeric covariate (univariate).
#' @param ndx Number of equal intervals on the covariate domain.
#' @param pord Order of the difference penalty (default 2).
#' @param bdeg B-spline degree (default 3).
#' @param xl,xr Optional knot boundaries.
#' @param ... Ignored.
#' @export
s <- function(x, ndx = 10L, pord = 2L, bdeg = 3L,
              xl = NULL, xr = NULL, ...) {
  as.numeric(x)
}

#' Split formula into parametric part and s() smooth calls
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

#' Evaluate one s() call into basis + penalties
#' @keywords internal
.eval_s_call <- function(cl, data, envir = parent.frame()) {
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
  ndx <- get_arg("ndx", NULL)
  if (is.null(ndx)) ndx <- get_arg("nseg", NULL)
  if (is.null(ndx)) ndx <- get_arg("k", 10L)
  pord <- get_arg("pord", NULL)
  if (is.null(pord)) pord <- get_arg("m", 2L)
  if (is.null(pord)) pord <- get_arg("d", 2L)
  bdeg <- as.integer(get_arg("bdeg", 3L))
  ndx <- as.integer(ndx)
  pord <- as.integer(pord)
  xl <- get_arg("xl", NULL)
  xr <- get_arg("xr", NULL)
  if (is.null(xl)) xl <- min(x)
  if (is.null(xr)) xr <- max(x)
  if (xr <= xl) xr <- xl + 1e-6

  B <- bbase(x, xl = xl, xr = xr, ndx = ndx, bdeg = bdeg)
  K <- ncol(B)
  S <- penalty_diff(K, pord = pord)
  C <- penalty_center(B)
  list(
    label = paste0("s(", xname, ")"),
    xname = xname,
    x = x,
    ndx = ndx,
    pord = pord,
    bdeg = bdeg,
    xl = xl,
    xr = xr,
    knots = attr(B, "knots"),
    B = B,
    S = S,
    C = C,
    K = K
  )
}

#' Build additive smooth design: stacked B, block-diag S and C
#' @keywords internal
.build_smooth_design <- function(formula, data, envir = parent.frame()) {
  parts <- .parse_smooth_formula(formula)
  if (is.null(data) || (is.list(data) && !is.data.frame(data) && !length(data))) {
    av <- all.vars(formula)
    data <- as.data.frame(mget(av, envir = envir, inherits = TRUE),
                          optional = TRUE)
  } else {
    data <- as.data.frame(data)
  }
  n <- nrow(data)
  if (!length(parts$s_calls)) {
    return(list(
      fixed = parts$fixed,
      n_smooth = 0L,
      B = matrix(0, n, 0L),
      S = matrix(0, 0L, 0L),
      C = matrix(0, 0L, 0L),
      smooth_K = integer(0),
      smooth_off = integer(0),
      specs = list(),
      labels = character(),
      blocks_idx = list()
    ))
  }
  specs <- lapply(parts$s_calls, .eval_s_call, data = data, envir = envir)
  labels <- vapply(specs, `[[`, character(1), "label")
  B <- do.call(cbind, lapply(specs, `[[`, "B"))
  S <- as.matrix(Matrix::bdiag(lapply(specs, `[[`, "S")))
  C <- as.matrix(Matrix::bdiag(lapply(specs, `[[`, "C")))
  Ks <- vapply(specs, `[[`, integer(1), "K")
  off <- if (length(Ks) <= 1L) {
    if (length(Ks) == 1L) 0L else integer(0)
  } else {
    c(0L, cumsum(Ks[-length(Ks)]))
  }
  blocks_idx <- lapply(seq_along(Ks), function(j) {
    off[j] + seq_len(Ks[j])
  })
  colnames(B) <- unlist(lapply(seq_along(specs), function(j) {
    paste0(labels[j], ".", seq_len(Ks[j]))
  }))
  list(
    fixed = parts$fixed,
    n_smooth = length(specs),
    B = B,
    S = S,
    C = C,
    smooth_K = as.integer(Ks),
    smooth_off = as.integer(off),
    specs = specs,
    labels = labels,
    blocks_idx = blocks_idx
  )
}

#' TMB data block for smooths
#' @keywords internal
.tmb_smooth_data <- function(sm, n, kappa = 1e6) {
  if (is.null(sm) || sm$n_smooth < 1L) {
    return(list(
      n_smooth = 0L,
      B = matrix(0, n, 0L),
      S = matrix(0, 0L, 0L),
      C = matrix(0, 0L, 0L),
      smooth_K = integer(0),
      smooth_off = integer(0),
      kappa = as.numeric(kappa)
    ))
  }
  list(
    n_smooth = as.integer(sm$n_smooth),
    B = sm$B,
    S = sm$S,
    C = sm$C,
    smooth_K = as.integer(sm$smooth_K),
    smooth_off = as.integer(sm$smooth_off),
    kappa = as.numeric(kappa)
  )
}

#' TMB parameters for smooth block
#' @keywords internal
.tmb_smooth_parameters <- function(sm) {
  if (is.null(sm) || sm$n_smooth < 1L) {
    return(list(log_lambda = numeric(0), alpha = numeric(0)))
  }
  list(
    log_lambda = rep(log(10), sm$n_smooth),
    alpha = rep(0, sum(sm$smooth_K))
  )
}

#' Predict an additive smooth with Marra–Wood (2012) pointwise bands
#'
#' Pointwise Bayesian confidence intervals for a centered smooth
#' \eqn{f_j(x)=b(x)^\top\gamma_j} as in Marra & Wood (2012), Biometrika:
#' \deqn{\widehat f_j(x)\pm z_{1-\alpha/2}\sqrt{b(x)^\top V_{\gamma_j} b(x)},}
#' where \eqn{V_{\gamma}} is the **marginal** Laplace / Bayesian posterior
#' covariance of the spline coefficients from TMB's joint precision
#' (Wahba–Silverman prior induced by the P-spline penalty). These bands
#' target good *across-the-function* frequentist coverage.
#'
#' @param object A `BBreg` or `BBmm` fit with smooths.
#' @param which Integer index or label of the smooth.
#' @param x Optional grid; default uses observed covariate.
#' @param level Confidence level (default 0.95).
#' @param se If `TRUE` (default), attach Marra–Wood `se`, `lwr`, `upr`.
#' @param method Only `"marrawood"` (default) is implemented; accepted for
#'   API clarity.
#' @return Data frame with `x`, `fit`, and usually `se`, `lwr`, `upr`.
#' @references
#' Marra, G. and Wood, S. N. (2012). Coverage properties of confidence
#' intervals for generalized additive model components. *Scandinavian
#' Journal of Statistics*, 39, 53–74.
#' @export
#' @seealso [plot_smooth()]
predict_smooth <- function(object, which = 1L, x = NULL, level = 0.95,
                           se = TRUE, method = c("marrawood")) {
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
  alpha <- as.numeric(object$alpha)[idx]
  if (is.null(x)) {
    x <- spec$x
  } else {
    x <- as.numeric(x)
  }
  Bg <- bbase(x, xl = spec$xl, xr = spec$xr, ndx = spec$ndx, bdeg = spec$bdeg)
  f <- as.numeric(Bg %*% alpha)

  out <- data.frame(x = x, fit = f)
  if (!isTRUE(se)) return(out)

  Vfull <- object$alpha.vcov
  if (is.null(Vfull) && !is.null(object$sdreport)) {
    n_u <- if (inherits(object, "BBmm")) object$nRand else 0L
    Vfull <- .alpha_vcov_from_joint(
      object$sdreport$jointPrecision,
      n_alpha = length(object$alpha),
      n_u = n_u
    )
  }
  if (is.null(Vfull) && !is.null(object$obj)) {
    sdr <- tryCatch(
      TMB::sdreport(object$obj, getJointPrecision = TRUE),
      error = function(e) NULL
    )
    if (!is.null(sdr)) {
      n_u <- if (inherits(object, "BBmm")) object$nRand else 0L
      Vfull <- .alpha_vcov_from_joint(
        sdr$jointPrecision,
        n_alpha = length(object$alpha),
        n_u = n_u
      )
    }
  }
  if (is.null(Vfull)) return(out)

  Vj <- Vfull[idx, idx, drop = FALSE]
  if (!all(is.finite(Vj))) return(out)

  se_hat <- as.numeric(smooth_se_cpp(Bg, Vj))
  z <- stats::qnorm(1 - (1 - level) / 2)
  out$se <- se_hat
  out$lwr <- f - z * se_hat
  out$upr <- f + z * se_hat
  attr(out, "level") <- level
  attr(out, "which") <- j
  attr(out, "label") <- sm$labels[j]
  attr(out, "method") <- "marrawood"
  out
}

#' Plot an additive P-spline with Marra–Wood (2012) confidence bands
#'
#' Draws \eqn{\hat f_j} and a shaded pointwise Bayesian band
#' (Marra & Wood, 2012). See [predict_smooth()].
#'
#' @param object A `BBreg` or `BBmm` fit with smooths.
#' @param which Integer index or label of the smooth (or `"all"`).
#' @param x Optional evaluation grid.
#' @param level Confidence level.
#' @param col Line color for the estimate.
#' @param shade Band fill color.
#' @param add If `TRUE`, add to the current plot (`which` must be length 1).
#' @param xlab,ylab,main Axis labels / title (`NULL` = defaults).
#' @param ylim y-limits; default spans band and fit.
#' @param ... Passed to [graphics::plot()] / [graphics::lines()].
#' @return Invisibly, the data frame from [predict_smooth()] (or a list if
#'   `which = "all"`).
#' @export
plot_smooth <- function(object, which = "all", x = NULL, level = 0.95,
                        col = "darkorange2", shade = grDevices::adjustcolor(col, 0.30),
                        add = FALSE, xlab = NULL, ylab = NULL, main = NULL,
                        ylim = NULL, ...) {
  sm <- object$smooth
  if (is.null(sm) || sm$n_smooth < 1L) {
    stop("object has no smooth terms", call. = FALSE)
  }
  if (identical(which, "all")) {
    which <- seq_len(sm$n_smooth)
  }
  if (length(which) > 1L) {
    if (isTRUE(add)) stop("add = TRUE requires a single smooth", call. = FALSE)
    n <- length(which)
    nr <- 1L
    nc <- n
    if (n > 3L) {
      nc <- ceiling(sqrt(n))
      nr <- ceiling(n / nc)
    }
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

  pr <- predict_smooth(object, which = which, x = x, level = level, se = TRUE,
                       method = "marrawood")
  o <- order(pr$x)
  pr <- pr[o, , drop = FALSE]
  lab <- attr(pr, "label")
  if (is.null(lab)) lab <- paste0("smooth ", which)
  if (is.null(xlab)) xlab <- "x"
  if (is.null(ylab)) ylab <- "f(x)"
  if (is.null(main)) {
    main <- paste0(lab, " (Marra-Wood ", round(100 * level), "% CI)")
  }
  if (is.null(ylim)) {
    ylim <- range(pr$fit, pr$lwr, pr$upr, na.rm = TRUE)
  }

  if (!isTRUE(add)) {
    graphics::plot(pr$x, pr$fit, type = "n", xlab = xlab, ylab = ylab,
                   main = main, ylim = ylim, ...)
    graphics::abline(h = 0, col = "grey70")
  }
  if (!is.null(pr$lwr) && !is.null(pr$upr)) {
    graphics::polygon(
      c(pr$x, rev(pr$x)),
      c(pr$lwr, rev(pr$upr)),
      col = shade, border = NA
    )
  }
  graphics::lines(pr$x, pr$fit, col = col, lwd = 2, ...)
  invisible(pr)
}

