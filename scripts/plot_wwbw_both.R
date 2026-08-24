## Plot WW (P-spline) and BW (linear) partial effects from M3_Hps
## Preferred model: between = fev_bar_z, wt_bar_z (parametric);
##                  within  = s(fev_w), s(wt_w) (Eilers).
##
## Usage: Rscript scripts/plot_wwbw_both.R

suppressPackageStartupMessages({
  .pkg <- normalizePath(".", winslash = "/", mustWork = FALSE)
  if (!file.exists(file.path(.pkg, "DESCRIPTION")))
    .pkg <- normalizePath("..", winslash = "/", mustWork = FALSE)
  if (file.exists(file.path(.pkg, "DESCRIPTION")) &&
      requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(.pkg, quiet = TRUE)
  } else library(PROregTMB)
})

.root <- local({
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in seq_len(6L)) {
    if (file.exists(file.path(d, "R", "BBmm.R"))) return(d)
    nd <- dirname(d); if (identical(nd, d)) break; d <- nd
  }
  d
})

out_dir <- file.path(.root, "scripts", "out_report4_reanalysis")
fit <- readRDS(file.path(out_dir, "reanalysis_report4.rds"))$results$M3_Hps$fit
stopifnot(!is.null(fit))

data_path <- normalizePath(
  file.path(.root, "..", "ref", "BB-GAM", "longitudinal", "data", "COPD_All.txt"),
  mustWork = TRUE
)
raw <- read.table(data_path, header = TRUE, stringsAsFactors = FALSE)
raw$id <- factor(raw$id)
raw$fev <- as.numeric(raw$fev); raw$wt <- as.numeric(raw$wt)
raw$fev_bar <- ave(raw$fev, raw$id, FUN = mean)
raw$wt_bar <- ave(raw$wt, raw$id, FUN = mean)
raw$fev_w <- raw$fev - raw$fev_bar
raw$wt_w <- raw$wt - raw$wt_bar
zscore <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s < 1e-8) x - mean(x, na.rm = TRUE) else (x - mean(x, na.rm = TRUE)) / s
}
# person-level between (one row per id)
uid <- !duplicated(raw$id)
fev_bar_z <- zscore(raw$fev_bar[uid])
wt_bar_z <- zscore(raw$wt_bar[uid])
fev_w_all <- raw$fev_w
wt_w_all <- raw$wt_w

dims <- c("Impacts", "Symptoms", "Activity")
beta <- fit$beta
Vbeta <- fit$fixed.vcov
zcrit <- qnorm(0.975)

get_b <- function(dm, v) {
  nm <- paste0(dm, ".", v)
  if (nm %in% names(beta)) unname(beta[[nm]]) else NA_real_
}
get_se_b <- function(dm, v) {
  nm <- paste0(dm, ".", v)
  if (!is.null(Vbeta) && nm %in% rownames(Vbeta))
    sqrt(max(0, Vbeta[nm, nm])) else NA_real_
}
get_a0 <- function(dm) {
  nm <- paste0(dm, ".(Intercept)")
  if (nm %in% names(beta)) unname(beta[[nm]]) else 0
}

# ---- Within contrast: joint Cov(beta_null, s) via predict_smooth ----
within_contrast <- function(lab, x, dm) {
  pr <- predict_smooth(fit, which = lab, x = x, se = TRUE, centered = TRUE)
  a0 <- get_a0(dm)
  eta <- a0 + pr$fit
  data.frame(
    x = pr$x,
    mu = plogis(eta),
    mu_lo = plogis(eta - zcrit * pr$se),
    mu_hi = plogis(eta + zcrit * pr$se),
    se_method = attr(pr, "method")
  )
}

# ---- Between linear partial on mu ----
between_line <- function(v, dm, x) {
  b <- get_b(dm, v); se_b <- get_se_b(dm, v); a0 <- get_a0(dm)
  # contrast vs mean 0 of the z-scored between: delta = b * x
  eta <- a0 + b * x
  se <- abs(x) * se_b
  data.frame(x = x, mu = plogis(eta),
             mu_lo = plogis(eta - zcrit * se),
             mu_hi = plogis(eta + zcrit * se),
             b = b, se_b = se_b)
}

fev_w_grid <- seq(quantile(fev_w_all, 0.02), quantile(fev_w_all, 0.98), length.out = 80)
wt_w_grid  <- seq(quantile(wt_w_all, 0.02),  quantile(wt_w_all, 0.98),  length.out = 80)
fev_b_grid <- seq(quantile(fev_bar_z, 0.02), quantile(fev_bar_z, 0.98), length.out = 80)
wt_b_grid  <- seq(quantile(wt_bar_z, 0.02),  quantile(wt_bar_z, 0.98),  length.out = 80)

out_png <- file.path(out_dir, "wwbw_both_effects.png")
out_pdf <- file.path(out_dir, "wwbw_both_effects.pdf")

plot_all <- function() {
  op <- par(mfrow = c(4, 3), mar = c(4, 4, 2.2, 0.8))
  on.exit(par(op), add = TRUE)

  # Row 1-2: BETWEEN (linear)
  for (v in c("fev_bar_z", "wt_bar_z")) {
    xlab <- if (v == "fev_bar_z") "FEV person-mean (z)" else "6MWT person-mean (z)"
    xseq <- if (v == "fev_bar_z") fev_b_grid else wt_b_grid
    for (dm in dims) {
      sp <- between_line(v, dm, xseq)
      ylim <- range(c(sp$mu_lo, sp$mu_hi), na.rm = TRUE)
      pad <- if (dm == dims[1] && v == "fev_bar_z") 0.08 * diff(ylim) else 0.02
      ylim <- c(max(0, ylim[1] - 0.02), min(1, ylim[2] + pad))
      plot(sp$x, sp$mu, type = "n", xlab = xlab, ylab = expression(mu),
           main = paste0(dm, " - BW linear"), ylim = ylim)
      polygon(c(sp$x, rev(sp$x)), c(sp$mu_lo, rev(sp$mu_hi)),
              col = adjustcolor("#1b9e77", 0.25), border = NA)
      lines(sp$x, sp$mu, lwd = 2, col = "#1b9e77")
      abline(v = 0, lty = 3, col = "grey50")
      grid()
      if (dm == dims[1] && v == "fev_bar_z") {
        legend("topright", c("between (linear)", "95% band"),
               col = c("#1b9e77", adjustcolor("#1b9e77", 0.25)),
               lty = c(1, NA), lwd = 2, pch = c(NA, 15), pt.cex = 1.6,
               bg = "white", box.col = "grey70", box.lwd = 0.6,
               cex = 0.72, y.intersp = 1.25, x.intersp = 0.85,
               inset = c(0.01, 0.01))
      }
    }
  }

  # Row 3-4: WITHIN (P-spline)
  for (v in c("fev_w", "wt_w")) {
    xlab <- if (v == "fev_w") "FEV within (pp from person mean)" else
      "6MWT within (m from person mean)"
    xseq <- if (v == "fev_w") fev_w_grid else wt_w_grid
    for (dm in dims) {
      lab <- paste0("s(", v, "):", dm)
      sp <- within_contrast(lab, xseq, dm)
      ylim <- range(c(sp$mu_lo, sp$mu_hi), na.rm = TRUE)
      pad <- if (dm == dims[1] && v == "fev_w") 0.08 * diff(ylim) else 0.02
      ylim <- c(max(0, ylim[1] - 0.02), min(1, ylim[2] + pad))
      plot(sp$x, sp$mu, type = "n", xlab = xlab, ylab = expression(mu),
           main = paste0(dm, " - WW P-spline"), ylim = ylim)
      polygon(c(sp$x, rev(sp$x)), c(sp$mu_lo, rev(sp$mu_hi)),
              col = adjustcolor("#d95f02", 0.25), border = NA)
      lines(sp$x, sp$mu, lwd = 2, col = "#d95f02")
      abline(v = 0, lty = 3, col = "grey50")
      grid()
      if (dm == dims[1] && v == "fev_w") {
        legend("topright", c("within P-spline", "95% band"),
               col = c("#d95f02", adjustcolor("#d95f02", 0.25)),
               lty = c(1, NA), lwd = 2, pch = c(NA, 15), pt.cex = 1.6,
               bg = "white", box.col = "grey70", box.lwd = 0.6,
               cex = 0.72, y.intersp = 1.25, x.intersp = 0.85,
               inset = c(0.01, 0.01))
      }
    }
  }
}

png(out_png, width = 1100, height = 1400, res = 120)
plot_all()
dev.off()
message("Wrote ", out_png)

pdf(out_pdf, width = 9, height = 11)
plot_all()
dev.off()
message("Wrote ", out_pdf)

# coef table for between
bw <- do.call(rbind, lapply(dims, function(dm) {
  data.frame(
    dim = dm,
    fev_bar_z = get_b(dm, "fev_bar_z"),
    se_fev = get_se_b(dm, "fev_bar_z"),
    wt_bar_z = get_b(dm, "wt_bar_z"),
    se_wt = get_se_b(dm, "wt_bar_z")
  )
}))
write.csv(bw, file.path(out_dir, "between_coefs.csv"), row.names = FALSE)
message("Wrote ", file.path(out_dir, "between_coefs.csv"))
print(bw, row.names = FALSE, digits = 3)
message(
  "Note: preferred model has BW *linear* and WW *P-spline* ",
  "(not s() on person means)."
)
