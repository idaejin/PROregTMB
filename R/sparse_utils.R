## Sparse helpers: Z RE design and Marra–Wood Bayesian Cov(alpha)

#' Marginal covariance block from TMB jointPrecision
#'
#' Returns \eqn{(Q^{-1})[I,I]} for index set \code{idx} via a sparse solve
#' \code{Q W = E}, not by inverting the submatrix \code{Q[I,I]} alone.
#'
#' @param jointPrecision Sparse/dense joint precision from `sdreport`.
#' @param idx Integer indices into the joint parameter order.
#' @return Dense \code{length(idx) x length(idx)} covariance, or `NULL`.
#' @keywords internal
.marginal_cov_from_joint <- function(jointPrecision, idx) {
  if (is.null(jointPrecision) || !length(idx)) return(NULL)
  idx <- as.integer(idx)
  if (anyNA(idx) || any(idx < 1L)) return(NULL)
  Jp <- jointPrecision
  if (!inherits(Jp, "Matrix")) {
    Jp <- Matrix::Matrix(as.matrix(Jp), sparse = TRUE)
  }
  nq <- nrow(Jp)
  if (any(idx > nq)) return(NULL)
  tryCatch({
    E <- Matrix::sparseMatrix(
      i = idx,
      j = seq_along(idx),
      x = rep(1, length(idx)),
      dims = c(nq, length(idx))
    )
    W <- Matrix::solve(Jp, E)
    as.matrix(W[idx, , drop = FALSE])
  }, error = function(e) NULL)
}

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
    if (!length(a_idx)) a_idx <- which(rn == "s")
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
  .marginal_cov_from_joint(Jp, a_idx)
}

#' Joint indices of fixed beta + wiggly s for one smooth (for contrast SE)
#' @keywords internal
.smooth_contrast_jp_index <- function(object, j) {
  sdr <- object$sdreport
  if (is.null(sdr) || is.null(sdr$jointPrecision)) return(NULL)
  Jp <- sdr$jointPrecision
  rn <- rownames(Jp)
  if (is.null(rn)) rn <- colnames(Jp)
  if (is.null(rn)) return(NULL)
  sm <- object$smooth
  spec <- sm$specs[[j]]
  s_block <- sm$blocks_idx[[j]]
  s_jp <- which(rn == "s")
  if (length(s_jp) < max(s_block)) return(NULL)
  idx_s <- s_jp[s_block]
  beta <- object$beta
  bn <- spec$null_name
  b_i <- integer(0)
  null_name <- NA_character_
  if (!is.null(beta) && !is.null(bn) && nzchar(bn)) {
    cand <- c(
      bn,
      if (!is.null(spec$by_level) && !is.na(spec$by_level))
        paste0(spec$by_level, ".", bn)
    )
    hit <- cand[cand %in% names(beta)]
    if (length(hit)) {
      b_pos <- match(hit[1], names(beta))
      beta_jp <- which(rn == "beta")
      if (length(beta_jp) >= b_pos) {
        b_i <- beta_jp[b_pos]
        null_name <- hit[1]
      }
    }
  }
  list(idx = c(b_i, idx_s), has_null = length(b_i) == 1L, null_name = null_name)
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
