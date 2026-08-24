#!/usr/bin/env Rscript
# Tests for fitted WW-BW within P-splines (Smooth / M3_Hps) vs Linear within.
#
# 1) Global Laplace nℓ comparison: Smooth vs Linear (and vs RS if available)
# 2) Per-smooth Wald test H0: s_g = 0 using Laplace s.vcov
#    (naive chi^2; under null sigma_s -> 0 this is on the boundary — interpret
#     as a descriptive diagnostic, not an exact classical p-value)
#
# Usage (from PROregTMB root):
#   Rscript scripts/test_wwbw_splines.R

suppressPackageStartupMessages({
  .root <- local({
    d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
    for (i in seq_len(6L)) {
      if (file.exists(file.path(d, "DESCRIPTION")) &&
          file.exists(file.path(d, "R", "BBmm.R"))) return(d)
      nd <- dirname(d)
      if (identical(nd, d)) break
      d <- nd
    }
    d
  })
  if (!requireNamespace("PROregTMB", quietly = TRUE)) {
    # try local load
    suppressMessages(devtools::load_all(.root, quiet = TRUE))
  } else {
    library(PROregTMB)
  }
})

out_dir <- file.path(.root, "scripts", "out_spline_tests")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rea_path <- file.path(.root, "scripts", "out_report4_reanalysis", "reanalysis_report4.rds")
rs_path <- file.path(.root, "scripts", "out_wwbw_rs_within", "rs_within_full.rds")
stopifnot(file.exists(rea_path))

rea <- readRDS(rea_path)
hps <- rea$results$M3_Hps
hlin <- rea$results$M3_Hlin
stopifnot(isTRUE(hps$ok), isTRUE(hlin$ok))

fit <- hps$fit
stopifnot(!is.null(fit$s), !is.null(fit$s.vcov), !is.null(fit$smooth_sd))

# ---- 1) Global nℓ comparisons ----
nll_smooth <- as.numeric(hps$nll)
nll_linear <- as.numeric(hlin$nll)
delta_lin <- nll_linear - nll_smooth
lr_lin <- 2 * delta_lin

global_rows <- data.frame(
  comparison = "Smooth_vs_Linear",
  nll_alt = nll_smooth,
  nll_null = nll_linear,
  delta_nll = delta_lin,
  twice_delta = lr_lin,
  note = "Linear within is nested in the P-spline null space (pord=2); 6 smoothing SDs on boundary under H0",
  stringsAsFactors = FALSE
)

if (file.exists(rs_path)) {
  rs <- readRDS(rs_path)
  # tolerate list layouts
  nll_rs <- NA_real_
  if (!is.null(rs$nll)) nll_rs <- as.numeric(rs$nll)
  if (is.list(rs) && !is.null(rs$results)) {
    nm <- names(rs$results)
    hit <- grep("RS$|Hlin\\+RS|shared", nm, ignore.case = TRUE, value = TRUE)
    if (length(hit)) nll_rs <- as.numeric(rs$results[[hit[1]]]$nll)
  }
  if (is.list(rs) && !is.null(rs$compare)) {
    cmp <- rs$compare
    if (is.data.frame(cmp) && "nll" %in% names(cmp)) {
      i <- grep("^RS$|shared", cmp$model, ignore.case = TRUE)
      if (length(i)) nll_rs <- as.numeric(cmp$nll[i[1]])
    }
  }
  # fallback: compare_full.csv
  cmp_csv <- file.path(dirname(rs_path), "compare_full.csv")
  if (!is.finite(nll_rs) && file.exists(cmp_csv)) {
    cmp <- read.csv(cmp_csv, stringsAsFactors = FALSE)
    i <- which(grepl("^RS$|shared", cmp$model, ignore.case = TRUE) &
                 !grepl("domain|by.?dim", cmp$model, ignore.case = TRUE))
    if (!length(i)) i <- which(grepl("RS", cmp$model, ignore.case = TRUE))
    if (length(i)) nll_rs <- as.numeric(cmp$nll[i[1]])
  }
  if (is.finite(nll_rs)) {
    global_rows <- rbind(
      global_rows,
      data.frame(
        comparison = "Smooth_vs_RS",
        nll_alt = nll_rs,
        nll_null = nll_smooth,
        delta_nll = nll_smooth - nll_rs,
        twice_delta = 2 * (nll_smooth - nll_rs),
        note = "Not nested; descriptive ranking only (RS preferred on panel)",
        stringsAsFactors = FALSE
      )
    )
  }
}

# ---- 2) Per-smooth Wald: H0 s_g = 0 ----
s <- as.numeric(fit$s)
names(s) <- names(fit$s)
V <- as.matrix(fit$s.vcov)
stopifnot(nrow(V) == length(s), ncol(V) == length(s))

group <- sub("\\.[0-9]+$", "", names(s))
ug <- unique(group)

wald_one <- function(idx) {
  sj <- s[idx]
  Vj <- V[idx, idx, drop = FALSE]
  # ridge tiny eigenvalues for numerical stability
  ev <- eigen(Vj, symmetric = TRUE)
  lam <- pmax(ev$values, 1e-10)
  Vinv <- ev$vectors %*% diag(1 / lam, length(lam)) %*% t(ev$vectors)
  W <- as.numeric(t(sj) %*% Vinv %*% sj)
  df <- length(sj)
  p <- pchisq(W, df = df, lower.tail = FALSE)
  list(W = W, df = df, p = p, max_abs_s = max(abs(sj)))
}

wald_rows <- lapply(ug, function(g) {
  idx <- which(group == g)
  w <- wald_one(idx)
  sd_g <- if (!is.null(fit$smooth_sd) && g %in% names(fit$smooth_sd)) {
    as.numeric(fit$smooth_sd[g])
  } else {
    NA_real_
  }
  lam_g <- if (!is.null(fit$lambda) && g %in% names(fit$lambda)) {
    as.numeric(fit$lambda[g])
  } else if (is.finite(sd_g) && sd_g > 0) {
    1 / sd_g^2
  } else {
    NA_real_
  }
  data.frame(
    smooth = g,
    n_coef = length(idx),
    smooth_sd = sd_g,
    lambda = lam_g,
    Wald = w$W,
    df = w$df,
    p_chisq = w$p,
    max_abs_s = w$max_abs_s,
    stringsAsFactors = FALSE
  )
})
wald_tab <- do.call(rbind, wald_rows)
rownames(wald_tab) <- NULL

# joint Wald all wiggly coefficients
w_all <- wald_one(seq_along(s))
joint <- data.frame(
  smooth = "ALL_within_smooths",
  n_coef = length(s),
  smooth_sd = NA_real_,
  lambda = NA_real_,
  Wald = w_all$W,
  df = w_all$df,
  p_chisq = w_all$p,
  max_abs_s = w_all$max_abs_s,
  stringsAsFactors = FALSE
)
wald_tab <- rbind(wald_tab, joint)

# ---- write ----
write.csv(global_rows, file.path(out_dir, "spline_global_nll.csv"), row.names = FALSE)
write.csv(wald_tab, file.path(out_dir, "spline_wald_by_smooth.csv"), row.names = FALSE)

cat("=== Global Laplace nll comparisons ===\n")
print(global_rows, row.names = FALSE)
cat("\n=== Per-smooth Wald (H0: wiggly s_g = 0; Laplace s.vcov) ===\n")
print(wald_tab, row.names = FALSE, digits = 4)
cat("\nCaveat: p_chisq is a Laplace / Wald diagnostic. Under H0 the smoothing\n")
cat("variance sits on the boundary, so classical chi^2 p-values are approximate.\n")
cat("Prefer ranking by nll (Smooth vs Linear vs RS) for model choice.\n")
cat("\nWrote:\n  ", file.path(out_dir, "spline_global_nll.csv"), "\n  ",
    file.path(out_dir, "spline_wald_by_smooth.csv"), "\n", sep = "")
