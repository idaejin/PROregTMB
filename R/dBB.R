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
  as.numeric(dBB_cpp(as.integer(m), as.numeric(p), as.numeric(phi)))
}

#' Random generation from the beta-binomial distribution
#'
#' Implemented in Rcpp for large `k` (simulations / benchmarks).
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
  as.integer(rBB_cpp(as.integer(k), as.numeric(m), as.numeric(p), as.numeric(phi)[1]))
}

#' Pointwise beta-binomial log-PMF (Rcpp)
#'
#' @param y Counts.
#' @param m Maximum scores (scalar or length `y`).
#' @param p Probabilities (scalar or length `y`).
#' @param phi Dispersion (scalar or length `y`).
#' @return Numeric vector of log-probabilities.
#' @export
ldBB <- function(y, m, p, phi) {
  y <- as.numeric(y)
  m <- as.numeric(m)
  p <- as.numeric(p)
  phi <- as.numeric(phi)
  as.numeric(ldbb_cpp(y, m, p, phi))
}
