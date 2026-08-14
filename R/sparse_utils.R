## Sparse helpers: Z RE design and alpha covariance from jointPrecision

#' Extract Cov(alpha) from TMB jointPrecision without densifying all RE
#'
#' Uses only the alpha block of the sparse precision (plus optional subject
#' RE block size to locate indices when names are ambiguous).
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
    Jp <- Matrix::Matrix(Jp, sparse = TRUE)
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
      # last n_alpha (common when random = c(u, alpha))
      a_idx <- seq.int(nr - n_alpha + 1L, nr)
    } else {
      return(NULL)
    }
  }
  Qa <- Jp[a_idx, a_idx, drop = FALSE]
  V <- tryCatch(as.matrix(Matrix::solve(Qa)), error = function(e) NULL)
  V
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
