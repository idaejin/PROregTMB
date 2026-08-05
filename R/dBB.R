#' Density of the beta-binomial distribution
#'
#' @param m Maximum score (number of trials), positive integer.
#' @param p Probability parameter in (0, 1).
#' @param phi Dispersion parameter (> 0).
#' @return Numeric vector of length `m + 1` with P(Y = 0), ..., P(Y = m).
#' @export
dBB <- function(m, p, phi) {
  if (length(m) > 1L || length(p) > 1L || length(phi) > 1L) {
    stop("m, p and phi must be scalars", call. = FALSE)
  }
  if (m != as.integer(m) || m <= 0) stop("m must be a positive integer", call. = FALSE)
  if (p < 0 || p > 1) stop("p must be in [0, 1]", call. = FALSE)
  if (phi <= 0) stop("phi must be positive", call. = FALSE)

  tt <- seq.int(0L, m)
  # Log-space for numerical stability
  a <- p / phi
  b <- (1 - p) / phi
  log_pmf <- lgamma(m + 1) - lgamma(tt + 1) - lgamma(m - tt + 1) +
    lgamma(a + tt) - lgamma(a) +
    lgamma(b + m - tt) - lgamma(b) +
    lgamma(a + b) - lgamma(a + b + m)
  exp(log_pmf)
}

#' Random generation from the beta-binomial distribution
#'
#' @param k Number of draws.
#' @param m Maximum score (scalar or length-k).
#' @param p Probability parameter (scalar or length-k).
#' @param phi Dispersion parameter (> 0).
#' @return Integer vector of length `k`.
#' @export
rBB <- function(k, m, p, phi) {
  if (k != as.integer(k) || k <= 0) stop("k must be a positive integer", call. = FALSE)
  if (phi <= 0) stop("phi must be positive", call. = FALSE)
  if (min(p) < 0 || max(p) > 1) stop("p must be in [0, 1]", call. = FALSE)
  if (min(m) <= 0 || any(m != as.integer(m))) stop("m must be positive integer(s)", call. = FALSE)
  if (length(m) > 1L && length(m) < k) stop("m must be scalar or length k", call. = FALSE)
  if (length(p) > 1L && length(p) < k) stop("p must be scalar or length k", call. = FALSE)

  alpha <- p / phi
  beta <- (1 - p) / phi
  u <- rbeta(k, alpha, beta)
  rbinom(k, m, u)
}
