# Formulation-aligned field aliases and print helpers for BBreg / BBmm.
# Primary names match the math (beta, phi, u, sigma, p); PROreg-style
# fields are kept as aliases for compatibility.

#' @keywords internal
.bb_model_lines_reg <- function(n_smooth = 0L) {
  if (n_smooth > 0L) {
    c(
      "Y_i ~ BB(m_i, p_i, phi)",
      "logit(p_i) = x_i' beta + sum_j f_j(x_ij)",
      "f_j = X_null beta_null + Z_j s_j;  s_j ~ N(0, sigma_sj^2 I)",
      "Z_j = B D'(DD')^{-1} (Eilers mixed; residualised vs X)"
    )
  } else {
    c(
      "Y_i ~ BB(m_i, p_i, phi)",
      "logit(p_i) = x_i' beta"
    )
  }
}

#' @keywords internal
.bb_model_lines_mm <- function(nDim = 1L, namesRand = NULL, re_blocks = NULL) {
  phi_line <- if (as.integer(nDim) > 1L) {
    "Y_i | u ~ BB(m_i, p_i, phi_{d(i)})"
  } else {
    "Y_i | u ~ BB(m_i, p_i, phi)"
  }
  if (length(re_blocks)) {
    bits <- vapply(re_blocks, function(b) {
      tn <- paste(b$term_names, collapse = " + ")
      corr <- if (!.is_corr_diag(b$corr) && b$q > 1L) "unstructured" else "diag"
      sprintf("(%s | %s) [%s]", tn, b$name, corr)
    }, character(1))
    re <- paste0("u ~ N(0, block-diag Sigma); ", paste(bits, collapse = " + "))
  } else if (length(namesRand)) {
    re <- paste0(
      "u ~ N(0, D),  D = block-diag(sigma_c^2);  c in {",
      paste(namesRand, collapse = ", "), "}"
    )
  } else {
    re <- "u ~ N(0, D)"
  }
  c(
    phi_line,
    "logit(p_i) = x_i' beta + z_i' u",
    re
  )
}

#' @keywords internal
.cat_model_block <- function(title, lines) {
  cat(title, "\n", sep = "")
  for (ln in lines) cat("  ", ln, "\n", sep = "")
  cat("\n")
}

#' Attach formulation-primary names on a BBreg fit (keeps PROreg aliases).
#' @keywords internal
.bbreg_attach_formulation <- function(out, formula) {
  beta <- setNames(as.numeric(out$coefficients), rownames(out$coefficients))
  out$beta <- beta
  out$phi <- as.numeric(out$phi)
  out$log_phi <- as.numeric(out$psi)
  out$p <- as.numeric(out$fitted.values)
  ns <- if (!is.null(out$smooth)) out$smooth$n_smooth else 0L
  out$model <- .bb_model_lines_reg(n_smooth = ns)
  out$formula <- formula
  out
}

#' Attach formulation-primary names on a BBmm fit (keeps PROreg aliases).
#' @keywords internal
.bbmm_attach_formulation <- function(out, fixed.formula = NULL,
                                     random.formula = NULL) {
  out$beta <- setNames(as.numeric(out$fixed.coef), names(out$fixed.coef))
  out$u <- as.numeric(out$random.coef)
  names(out$u) <- names(out$random.coef)
  out$sigma <- setNames(as.numeric(out$sigma.coef), names(out$sigma.coef))
  out$phi <- out$phi.coef
  out$log_phi <- out$psi.coef
  out$p <- as.numeric(out$fitted.values)
  out$fixed.formula <- fixed.formula
  out$random.formula <- random.formula
  out$model <- .bb_model_lines_mm(out$nDim, out$namesRand, out$re_blocks)
  out
}
