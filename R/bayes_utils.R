#' @keywords internal
.match_fit_method <- function(method) {
  method <- tolower(as.character(method)[1])
  if (!method %in% c("mle", "bayes")) {
    stop("method must be 'mle' or 'bayes'", call. = FALSE)
  }
  method
}

#' Run tmbstan on a TMB object (optional Suggests).
#' @keywords internal
.run_tmbstan <- function(obj, init_fn, chains = 2L, iter = 1000L,
                         warmup = 400L, seed = 1L, laplace = FALSE,
                         adapt_delta = 0.95, max_treedepth = 12L) {
  if (!requireNamespace("tmbstan", quietly = TRUE) ||
      !requireNamespace("rstan", quietly = TRUE)) {
    stop("method = 'bayes' requires packages tmbstan and rstan", call. = FALSE)
  }
  tmbstan::tmbstan(
    obj,
    chains = as.integer(chains),
    iter = as.integer(iter),
    warmup = as.integer(warmup),
    seed = as.integer(seed),
    refresh = 0,
    init = init_fn,
    laplace = isTRUE(laplace),
    control = list(
      adapt_delta = adapt_delta,
      max_treedepth = as.integer(max_treedepth)
    )
  )
}

#' Posterior summary for fixed-effect style parameters from a stanfit.
#' @keywords internal
.posterior_from_stanfit <- function(fit, beta_names, phi_names = "phi",
                                    sigma_names = NULL) {
  sm <- as.data.frame(rstan::summary(fit)$summary)
  draws <- as.matrix(fit)
  rows <- list()

  b_idx <- grep("^beta(\\[|$)", rownames(sm))
  if (!length(b_idx) && any(grepl("^beta", colnames(draws)))) {
    b_idx <- grep("^beta", rownames(sm))
  }
  if (length(b_idx)) {
    nm <- beta_names
    if (length(nm) != length(b_idx)) nm <- rownames(sm)[b_idx]
    rows[[length(rows) + 1L]] <- data.frame(
      param = nm,
      mean = sm[b_idx, "mean"],
      sd = sm[b_idx, "sd"],
      q025 = sm[b_idx, "2.5%"],
      q975 = sm[b_idx, "97.5%"],
      n_eff = sm[b_idx, "n_eff"],
      Rhat = sm[b_idx, "Rhat"],
      stringsAsFactors = FALSE
    )
  }

  # phi via exp(log_phi*)
  lp_cols <- grep("^log_phi", colnames(draws), value = TRUE)
  if (length(lp_cols)) {
    for (j in seq_along(lp_cols)) {
      ph <- exp(draws[, lp_cols[j]])
      nm <- if (length(phi_names) >= j) phi_names[j] else paste0("phi", j)
      rname <- lp_cols[j]
      rows[[length(rows) + 1L]] <- data.frame(
        param = nm,
        mean = mean(ph),
        sd = stats::sd(ph),
        q025 = unname(stats::quantile(ph, 0.025)),
        q975 = unname(stats::quantile(ph, 0.975)),
        n_eff = if (rname %in% rownames(sm)) sm[rname, "n_eff"] else NA_real_,
        Rhat = if (rname %in% rownames(sm)) sm[rname, "Rhat"] else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }

  ls_cols <- grep("^log_sigma", colnames(draws), value = TRUE)
  if (length(ls_cols) && !is.null(sigma_names)) {
    for (j in seq_along(ls_cols)) {
      sg <- exp(draws[, ls_cols[j]])
      nm <- if (length(sigma_names) >= j) sigma_names[j] else paste0("sigma", j)
      rname <- ls_cols[j]
      rows[[length(rows) + 1L]] <- data.frame(
        param = nm,
        mean = mean(sg),
        sd = stats::sd(sg),
        q025 = unname(stats::quantile(sg, 0.025)),
        q975 = unname(stats::quantile(sg, 0.975)),
        n_eff = if (rname %in% rownames(sm)) sm[rname, "n_eff"] else NA_real_,
        Rhat = if (rname %in% rownames(sm)) sm[rname, "Rhat"] else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}
