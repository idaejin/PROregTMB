#!/usr/bin/env Rscript
# WW-BW sensitivity plots: linear within + random slopes vs Hlin / Hps.
#
# Outputs under scripts/out_wwbw_rs_within/:
#   rs_nll_compare.{png,pdf}
#   rs_sigma.{png,pdf}
#   rs_within_fan.{png,pdf}
#
# Usage (from PROregTMB root):
#   Rscript scripts/plot_wwbw_rs_sensitivity.R
#   Rscript scripts/plot_wwbw_rs_sensitivity.R --refit   # force RS refit for betas

args <- commandArgs(trailingOnly = TRUE)
force_refit <- "--refit" %in% args

suppressPackageStartupMessages({
  .pkg <- normalizePath(".", winslash = "/", mustWork = FALSE)
  if (!file.exists(file.path(.pkg, "DESCRIPTION"))) {
    .pkg <- normalizePath("..", winslash = "/", mustWork = FALSE)
  }
  if (file.exists(file.path(.pkg, "DESCRIPTION")) &&
      requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(.pkg, quiet = TRUE)
  } else {
    library(PROregTMB)
  }
})

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

out_dir <- file.path(.root, "scripts", "out_wwbw_rs_within")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
sum_rds <- file.path(out_dir, "rs_within_full.rds")
beta_rds <- file.path(out_dir, "rs_fit_betas.rds")
stopifnot(file.exists(sum_rds))
summ <- readRDS(sum_rds)

# Reference Hps / Hlin from Report-4 reanalysis
rea <- readRDS(file.path(.root, "scripts", "out_report4_reanalysis",
                         "reanalysis_report4.rds"))
hlin_beta <- rea$results$M3_Hlin$beta
hps_nll <- as.numeric(rea$results$M3_Hps$nll)

dims <- c("Impacts", "Symptoms", "Activity")
col_pop <- "#2166ac"
col_hlin <- "grey35"
col_band1 <- adjustcolor("#2166ac", 0.18)
col_band2 <- adjustcolor("#2166ac", 0.08)

# ---- 1) nll comparison ----
nll_df <- data.frame(
  model = c("Linear", "Smooth", "RS", "RS-domain"),
  nll = c(summ$Hlin$nll, hps_nll, summ$RS$nll, summ$RS_bydim$nll),
  stringsAsFactors = FALSE
)
nll_df$delta <- nll_df$nll - nll_df$nll[1]
# Display top = best (lowest nll)
ord <- order(nll_df$nll)
nll_show <- nll_df[ord, , drop = FALSE]
cols_nll <- c("Linear" = "grey70", "Smooth" = "#1b9e77",
              "RS" = "#2166ac", "RS-domain" = "#7570b3")

plot_nll <- function() {
  op <- par(mar = c(5, 11, 3, 2))
  on.exit(par(op), add = TRUE)
  xlim <- c(min(nll_show$nll) - 80, max(nll_show$nll) + 30)
  bp <- barplot(
    nll_show$nll, names.arg = nll_show$model, horiz = TRUE,
    las = 1, col = unname(cols_nll[nll_show$model]),
    xlim = xlim, xlab = "Approximate marginal nll (lower better)",
    main = "WW-BW mean structures (M3 intercepts)"
  )
  abline(v = nll_df$nll[nll_df$model == "Linear"], lty = 3, col = "grey40")
  text(nll_show$nll, bp,
       sprintf("  %.1f (D %+.1f)", nll_show$nll, nll_show$delta),
       adj = c(0, 0.5), cex = 0.8, xpd = TRUE)
}

png(file.path(out_dir, "rs_nll_compare.png"), width = 900, height = 420, res = 120)
plot_nll()
dev.off()
pdf(file.path(out_dir, "rs_nll_compare.pdf"), width = 8, height = 3.6)
plot_nll()
dev.off()
message("Wrote rs_nll_compare.{png,pdf}")

# ---- 2) sigma comparison (slope panel on comparable within-SD scale) ----
sig_h <- summ$Hlin$sigma
sig_r <- summ$RS$sigma
# Align names for intercepts + slopes
nm_int <- c("Impacts", "Symptoms", "Activity")

# Within-covariate SDs on the visit panel (same prep as RS fits)
data_path <- normalizePath(
  file.path(.root, "..", "ref", "BB-GAM", "longitudinal", "data", "COPD_All.txt"),
  mustWork = TRUE
)
raw0 <- read.table(data_path, header = TRUE, stringsAsFactors = FALSE)
raw0$id <- factor(raw0$id)
raw0$fev <- as.numeric(raw0$fev)
raw0$wt <- as.numeric(raw0$wt)
ok0 <- complete.cases(raw0[, c("fev", "wt", "age", "sex", "cluster_INMA",
                               "IMPACTS", "SYMPTOMS", "ACTIVITY")])
raw0 <- raw0[ok0, , drop = FALSE]
raw0$fev_w <- raw0$fev - ave(raw0$fev, raw0$id, FUN = mean)
raw0$wt_w <- raw0$wt - ave(raw0$wt, raw0$id, FUN = mean)
sd_fev_w <- sd(raw0$fev_w)
sd_wt_w <- sd(raw0$wt_w)
tau_fev <- as.numeric(sig_r[4])
tau_wt <- as.numeric(sig_r[5])
# Comparable scales: logit SD per 1 within-SD of the covariate
ys_scaled <- c(tau_fev * sd_fev_w, tau_wt * sd_wt_w)
nm_sl_scaled <- c("FEV (per 1 within-SD)", "6MWT (per 1 within-SD)")

plot_sigma <- function() {
  op <- par(mfrow = c(1, 2), mar = c(6.5, 4.5, 3, 1))
  on.exit(par(op), add = TRUE)

  # Intercept SDs
  yh <- as.numeric(sig_h[1:3])
  yr <- as.numeric(sig_r[1:3])
  mat <- rbind(Linear = yh, RS = yr)
  barplot(mat, beside = TRUE, names.arg = nm_int, las = 1,
          col = c("grey70", col_pop), ylab = "SD",
          main = "Domain random intercepts",
          ylim = c(0, max(mat) * 1.15), cex.names = 0.9)
  legend("topright", c("Linear", "RS"), fill = c("grey70", col_pop), bty = "n")

  # Slope SDs on within-SD scale (comparable across covariates)
  bp2 <- barplot(ys_scaled, names.arg = nm_sl_scaled, las = 2, col = col_pop,
                 ylab = "SD on logit scale",
                 main = "Within random slopes (RS)",
                 ylim = c(0, max(ys_scaled) * 1.45), cex.names = 0.85)
  text(bp2, ys_scaled, sprintf("%.3f", ys_scaled), pos = 3, cex = 0.9)
  mtext(
    sprintf(
      "Raw: tau_FEV=%.3f /pp, tau_6MWT=%.4f /m; within SD = %.1f pp, %.0f m",
      tau_fev, tau_wt, sd_fev_w, sd_wt_w
    ),
    side = 1, line = 5.2, cex = 0.7, col = "grey30"
  )
}

png(file.path(out_dir, "rs_sigma.png"), width = 1000, height = 520, res = 120)
plot_sigma()
dev.off()
pdf(file.path(out_dir, "rs_sigma.pdf"), width = 9, height = 4.6)
plot_sigma()
dev.off()
message("Wrote rs_sigma.{png,pdf} (slope panel = per within-SD)")
message(sprintf("  scaled tau: FEV=%.3f  6MWT=%.3f", ys_scaled[1], ys_scaled[2]))

# ---- 3) Within fan: need RS betas ----
get_or_fit_rs_beta <- function() {
  if (!force_refit && file.exists(beta_rds)) {
    message("Loading cached RS betas: ", beta_rds)
    return(readRDS(beta_rds))
  }
  message("Fitting M3_Hlin+RS to recover betas (approx 3-4 min)...")
  data_path <- normalizePath(
    file.path(.root, "..", "ref", "BB-GAM", "longitudinal", "data", "COPD_All.txt"),
    mustWork = TRUE
  )
  zscore <- function(x) {
    s <- sd(x, na.rm = TRUE)
    if (!is.finite(s) || s < 1e-8) return(x - mean(x, na.rm = TRUE))
    (x - mean(x, na.rm = TRUE)) / s
  }
  sgrq_bin <- function(x) pmin(24L, as.integer(floor(as.numeric(x) / 4)))

  raw <- read.table(data_path, header = TRUE, stringsAsFactors = FALSE)
  raw$id <- factor(raw$id)
  raw$time <- as.integer(raw$time)
  raw$sex <- factor(raw$sex, levels = c("1", "2"), labels = c("Male", "Female"))
  raw$fev <- as.numeric(raw$fev)
  raw$wt <- as.numeric(raw$wt)
  raw$age <- as.numeric(raw$age)
  raw$time_c <- raw$time - 1L
  raw$age_z <- zscore(raw$age)
  inna_map <- c("1" = "B", "2" = "C", "3" = "A", "4" = "D")
  raw$hclus4 <- factor(inna_map[as.character(raw$cluster_INMA)],
                       levels = c("B", "A", "C", "D"))
  raw$fev_bar <- ave(raw$fev, raw$id, FUN = mean)
  raw$wt_bar <- ave(raw$wt, raw$id, FUN = mean)
  raw$fev_w <- raw$fev - raw$fev_bar
  raw$wt_w <- raw$wt - raw$wt_bar
  raw$fev_bar_z <- zscore(raw$fev_bar)
  raw$wt_bar_z <- zscore(raw$wt_bar)
  ok <- complete.cases(raw[, c("fev", "wt", "age", "sex", "cluster_INMA",
                               "IMPACTS", "SYMPTOMS", "ACTIVITY")])
  raw <- raw[ok, , drop = FALSE]
  raw$id <- droplevels(raw$id)

  y_cols <- c("IMPACTS", "SYMPTOMS", "ACTIVITY")
  long <- do.call(rbind, lapply(seq_along(dims), function(d) {
    data.frame(
      id = raw$id, time = raw$time, time_c = raw$time_c,
      dim = factor(dims[d], levels = dims),
      y = sgrq_bin(raw[[y_cols[d]]]), m = 24L,
      sex = raw$sex, hclus4 = raw$hclus4, age_z = raw$age_z,
      fev_bar_z = raw$fev_bar_z, wt_bar_z = raw$wt_bar_z,
      fev_w = raw$fev_w, wt_w = raw$wt_w,
      stringsAsFactors = FALSE
    )
  }))
  rownames(long) <- NULL

  fit <- BBmm(
    y ~ hclus4 + time_c + sex + age_z + fev_bar_z + wt_bar_z + fev_w + wt_w,
    random = ~ (0 + dim | id) + (0 + fev_w + wt_w | id),
    dim = "dim", corr = "unstructured", m = "m", data = long,
    silent = TRUE, show = TRUE, maxiter = 600L
  )
  out <- list(
    beta = fit$beta, sigma = fit$sigma, nll = fit$nll, conv = fit$conv,
    fev_w = raw$fev_w, wt_w = raw$wt_w
  )
  saveRDS(out, beta_rds)
  message("Cached ", beta_rds)
  out
}

rsb <- get_or_fit_rs_beta()
beta_rs <- rsb$beta
sig_fev <- as.numeric(summ$RS$sigma["id.fev_w"])
sig_wt  <- as.numeric(summ$RS$sigma["id.wt_w"])

slope <- function(beta, dm, v) {
  nm <- paste0(dm, ".", v)
  if (nm %in% names(beta)) unname(beta[[nm]]) else NA_real_
}
a0 <- function(beta, dm) {
  nm <- paste0(dm, ".(Intercept)")
  if (nm %in% names(beta)) unname(beta[[nm]]) else 0
}

# grids from data
if (!is.null(rsb$fev_w)) {
  fev_w_all <- rsb$fev_w
  wt_w_all <- rsb$wt_w
} else {
  # fallback: rebuild quickly from COPD
  data_path <- normalizePath(
    file.path(.root, "..", "ref", "BB-GAM", "longitudinal", "data", "COPD_All.txt"),
    mustWork = TRUE
  )
  raw0 <- read.table(data_path, header = TRUE, stringsAsFactors = FALSE)
  raw0$fev <- as.numeric(raw0$fev); raw0$wt <- as.numeric(raw0$wt)
  raw0$id <- factor(raw0$id)
  fev_bar <- ave(raw0$fev, raw0$id, FUN = mean)
  wt_bar <- ave(raw0$wt, raw0$id, FUN = mean)
  fev_w_all <- raw0$fev - fev_bar
  wt_w_all <- raw0$wt - wt_bar
}
fev_grid <- seq(quantile(fev_w_all, 0.02, na.rm = TRUE),
                quantile(fev_w_all, 0.98, na.rm = TRUE), length.out = 80)
wt_grid <- seq(quantile(wt_w_all, 0.02, na.rm = TRUE),
               quantile(wt_w_all, 0.98, na.rm = TRUE), length.out = 80)
sd_fev_w <- sd(fev_w_all, na.rm = TRUE)
sd_wt_w <- sd(wt_w_all, na.rm = TRUE)

fan_curve <- function(dm, v, x, sig) {
  b_rs <- slope(beta_rs, dm, v)
  b_hl <- slope(hlin_beta, dm, v)
  a <- a0(beta_rs, dm)
  # population + heterogeneity bands on logit, then plogis
  eta0 <- a + b_rs * x
  list(
    x = x,
    mu = plogis(eta0),
    mu_p1 = plogis(a + (b_rs + sig) * x),
    mu_m1 = plogis(a + (b_rs - sig) * x),
    mu_p2 = plogis(a + (b_rs + 2 * sig) * x),
    mu_m2 = plogis(a + (b_rs - 2 * sig) * x),
    mu_hlin = plogis(a0(hlin_beta, dm) + b_hl * x),
    b_rs = b_rs, b_hl = b_hl, sig = sig
  )
}

plot_fan <- function() {
  op <- par(mfrow = c(2, 3), mar = c(4.2, 4, 2.8, 0.8))
  on.exit(par(op), add = TRUE)
  for (v in c("fev_w", "wt_w")) {
    # Plot on within-SD x-scale so FEV and 6MWT fan widths are comparable
    sd_w <- if (v == "fev_w") sd_fev_w else sd_wt_w
    x_raw <- if (v == "fev_w") fev_grid else wt_grid
    x_z <- x_raw / sd_w
    sig_raw <- if (v == "fev_w") sig_fev else sig_wt
    sig_z <- sig_raw * sd_w
    xlab <- if (v == "fev_w") {
      "FEV within (SD units from person mean)"
    } else {
      "6MWT within (SD units from person mean)"
    }
    for (dm in dims) {
      # Use slopes on the z-scale: b_z = b_raw * sd_w
      b_rs <- slope(beta_rs, dm, v) * sd_w
      b_hl <- slope(hlin_beta, dm, v) * sd_w
      a <- a0(beta_rs, dm)
      eta0 <- a + b_rs * x_z
      sp <- list(
        x = x_z,
        mu = plogis(eta0),
        mu_p1 = plogis(a + (b_rs + sig_z) * x_z),
        mu_m1 = plogis(a + (b_rs - sig_z) * x_z),
        mu_p2 = plogis(a + (b_rs + 2 * sig_z) * x_z),
        mu_m2 = plogis(a + (b_rs - 2 * sig_z) * x_z),
        mu_hlin = plogis(a0(hlin_beta, dm) + b_hl * x_z),
        b_rs = b_rs, sig = sig_z
      )
      ylim <- range(c(sp$mu_m2, sp$mu_p2, sp$mu_hlin), na.rm = TRUE)
      # Extra headroom so the legend does not sit on the fan
      pad <- 0.06 * diff(ylim)
      ylim <- c(max(0, ylim[1] - 0.02), min(1, ylim[2] + pad))
      plot(sp$x, sp$mu, type = "n", xlab = xlab, ylab = expression(mu),
           main = paste0(dm, if (v == "fev_w") " - FEV within" else " - 6MWT within"),
           ylim = ylim)
      polygon(c(sp$x, rev(sp$x)), c(sp$mu_m2, rev(sp$mu_p2)),
              col = col_band2, border = NA)
      polygon(c(sp$x, rev(sp$x)), c(sp$mu_m1, rev(sp$mu_p1)),
              col = col_band1, border = NA)
      lines(sp$x, sp$mu, lwd = 2.2, col = col_pop)
      lines(sp$x, sp$mu_hlin, lwd = 1.8, lty = 2, col = col_hlin)
      abline(v = 0, lty = 3, col = "grey50")
      grid()
      if (dm == dims[1] && v == "fev_w") {
        legend("topright",
               c("RS population", "+/-1 SD", "+/-2 SD", "Linear (no RS)"),
               col = c(col_pop, adjustcolor(col_pop, 0.45),
                       adjustcolor(col_pop, 0.2), col_hlin),
               lty = c(1, NA, NA, 2), lwd = c(2.2, NA, NA, 1.8),
               pch = c(NA, 15, 15, NA), pt.cex = 1.4,
               bg = "white", box.col = "grey70", box.lwd = 0.6,
               cex = 0.72, y.intersp = 1.25, x.intersp = 0.85,
               inset = c(0.01, 0.01))
      }
    }
  }
}

png(file.path(out_dir, "rs_within_fan.png"), width = 1100, height = 780, res = 120)
plot_fan()
dev.off()
pdf(file.path(out_dir, "rs_within_fan.pdf"), width = 9.5, height = 6.5)
plot_fan()
dev.off()
message("Wrote rs_within_fan.{png,pdf}")

# coef snapshot
coef_tab <- do.call(rbind, lapply(dims, function(dm) {
  data.frame(
    dim = dm,
    fev_w_Hlin = slope(hlin_beta, dm, "fev_w"),
    fev_w_RS = slope(beta_rs, dm, "fev_w"),
    wt_w_Hlin = slope(hlin_beta, dm, "wt_w"),
    wt_w_RS = slope(beta_rs, dm, "wt_w"),
    stringsAsFactors = FALSE
  )
}))
write.csv(coef_tab, file.path(out_dir, "rs_vs_hlin_within_slopes.csv"),
          row.names = FALSE)
print(coef_tab)
message("Done. Figures in ", out_dir)
