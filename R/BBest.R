#' Estimate beta-binomial parameters (no covariates)
#'
#' @param y Response counts.
#' @param m Maximum score (scalar or vector).
#' @param method `"MM"` (moments) or `"MLE"` (via TMB / BBreg intercept-only).
#' @return Object of class `BBest`.
#' @export
BBest <- function(y, m, method = c("MM", "MLE")) {
  method <- match.arg(method)
  y <- as.numeric(y)
  n <- length(y)
  m. <- if (length(m) == 1L) rep(m, n) else as.numeric(m)
  if (length(m.) != n) stop("m must be scalar or length(y)", call. = FALSE)
  if (any(y < 0 | y > m.)) stop("y must be in 0..m", call. = FALSE)

  balanced <- if (length(unique(m.)) == 1L) "yes" else "no"

  if (method == "MM") {
    E <- mean(y)
    V <- stats::var(y)
    mbar <- mean(m.)
    p <- E / mbar
    p <- min(max(p, 1e-8), 1 - 1e-8)
    denom <- mbar * p * (1 - p) * mbar - V
    phi <- if (abs(denom) < .Machine$double.eps) NA_real_ else
      (V - mbar * p * (1 - p)) / denom
    if (is.na(phi) || phi <= 0) phi <- 1e-4
    out <- list(
      p = p, phi = phi, pVar = NA_real_,
      psi = log(phi), psiVar = NA_real_,
      m = if (balanced == "yes") m.[1] else m.,
      balanced = balanced, method = method
    )
    class(out) <- "BBest"
    return(out)
  }

  # MLE via intercept-only BBreg
  dat <- data.frame(y = y)
  fit <- BBreg(y ~ 1, m = m., data = dat)
  p <- as.numeric(fit$fitted.values[1])
  out <- list(
    p = p, phi = fit$phi, pVar = as.numeric(fit$vcov[1, 1]) * (p * (1 - p))^2,
    psi = fit$psi, psiVar = fit$psi.var,
    m = if (balanced == "yes") m.[1] else m.,
    balanced = balanced, method = method,
    tmb = fit
  )
  class(out) <- "BBest"
  out
}

#' @export
print.BBest <- function(x, ...) {
  cat("Beta-binomial parameter estimates (", x$method, ")\n", sep = "")
  cat("  p   :", x$p, "\n")
  cat("  phi :", x$phi, "\n")
  invisible(x)
}

#' @export
summary.BBest <- function(object, ...) {
  tab <- cbind(
    Estimate = c(p = object$p, phi = object$phi, `log(phi)` = object$psi)
  )
  out <- list(coefficients = tab, method = object$method, object = object)
  class(out) <- "summary.BBest"
  out
}

#' @export
print.summary.BBest <- function(x, ...) {
  cat("Beta-binomial estimates (", x$method, ")\n\n", sep = "")
  print(x$coefficients)
  invisible(x)
}
