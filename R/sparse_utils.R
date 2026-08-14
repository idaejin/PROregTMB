## Sparse helpers: Z RE design and Marra–Wood Bayesian Cov(alpha)

#' Marginal Cov(alpha) from TMB jointPrecision (Marra–Wood Bayesian CI)
#'
#' For random effects stacked as \code{(u, alpha)}, the Marra–Wood (2012)
#' Bayesian covariance of the spline coefficients is the **marginal**
#' block of \eqn{Q^{-1}}, i.e. \code{solve(Q)[alpha, alpha]}, obtained via a
#' sparse solve against the alpha identity columns (not
#' \code{solve(Q[alpha,alpha])}, which ignores dependence on subject RE).
#'
#' @param jointPrecision Sparse/dense joint precision from `sdreport`.
#' @param n_alpha Number of spline coefficients.
#' @param n_u Number of subject random effects preceding alpha (BBmm).
#' @return Dense `n_alpha x n_alpha` covariance, or `NULL`.
#' @keywords internal
.alpha_vcov_from_joint <- function(jointPrecision, n_alpha, n_u = 0L) {
  if (is.null(jointPrecision) || n_alpha < 1L) return(NULL)
  Jp <- jointPrecision
  if (!inherits(Jp, "Matrix")) {
    Jp <- Matrix::Matrix(as.matrix(Jp), sparse = TRUE)
  }
  rn <- colnames(Jp)
  if (is.null(rn)) rn <- rownames(Jp)
  a_idx <- integer(0)
  if (!is.null(rn)) {
    a_idx <- which(rn == "alpha" | grepl("^alpha($|\\[)", rn))
    if (!length(a_idx)) a_idx <- grep("alpha", rn)
  }
  if (length(a_idx) != n_alpha) {
    nr <- nrow(Jp)
    if (as.integer(n_u) + n_alpha == nr) {
      a_idx <- seq.int(as.integer(n_u) + 1L, as.integer(n_u) + n_alpha)
    } else if (n_alpha <= nr) {
      a_idx <- seq.int(nr - n_alpha + 1L, nr)
    } else {
      return(NULL)
    }
  }
  # Marginal Cov(alpha) = (Q^{-1})_{aa} via sparse solve Q W = I[, a]
  tryCatch({
    nq <- nrow(Jp)
    E <- Matrix::sparseMatrix(
      i = a_idx,
      j = seq_along(a_idx),
      x = rep(1, length(a_idx)),
      dims = c(nq, length(a_idx))
    )
    W <- Matrix::solve(Jp, E)
    as.matrix(W[a_idx, , drop = FALSE])
  }, error = function(e) NULL)
}

#' Densify Z for TMB DATA_MATRIX (keeps attributes)
#' @keywords internal
.as_dense_Z <- function(Z) {
  if (inherits(Z, "sparseMatrix")) {
    Zn <- as.matrix(Z)
    colnames(Zn) <- colnames(Z)
    return(Zn)
  }
  as.matrix(Z)
}
