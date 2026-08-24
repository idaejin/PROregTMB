#!/usr/bin/env Rscript
# WW-BW both-effects figure for the best-by-nll model:
#   M3 Hlin + within random slopes (shared).
# Layout mirrors plot_wwbw_both.R:
#   rows 1-2: between linear (green)
#   rows 3-4: within linear + +/-1, +/-2 SD random-slope fan (blue)
#
# Usage (from PROregTMB root):
#   Rscript scripts/plot_wwbw_rs_both.R

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

out_dir <- file.path(.root, "scripts", "out_wwbw_rs_within")
beta_rds <- file.path(out_dir, "rs_fit_betas.rds")
sum_rds <- file.path(out_dir, "rs_within_full.rds")
stopifnot(file.exists(beta_rds), file.exists(sum_rds))

rsb <- readRDS(beta_rds)
summ <- readRDS(sum_rds)

# Need fixed.vcov for between 95% bands; refit once if missing
if (is.null(rsb$fixed.vcov)) {
  message("fixed.vcov missing: refitting Hlin+RS to recover SEs (~3-4 min)...")
  data_path0 <- normalizePath(
    file.path(.root, "..", "ref", "BB-GAM", "longitudinal", "data", "COPD_All.txt"),
    mustWork = TRUE
  )
  zscore0 <- function(x) {
    s <- sd(x, na.rm = TRUE)
    if (!is.finite(s) || s < 1e-8) return(x - mean(x, na.rm = TRUE))
    (x - mean(x, na.rm = TRUE)) / s
  }
  sgrq_bin0 <- function(x) pmin(24L, as.integer(floor(as.numeric(x) / 4)))
  raw0 <- read.table(data_path0, header = TRUE, stringsAsFactors = FALSE)
  raw0$id <- factor(raw0$id)
  raw0$time <- as.integer(raw0$time)
  raw0$sex <- factor(raw0$sex, levels = c("1", "2"), labels = c("Male", "Female"))
  raw0$fev <- as.numeric(raw0$fev)
  raw0$wt <- as.numeric(raw0$wt)
  raw0$age <- as.numeric(raw0$age)
  raw0$time_c <- raw0$time - 1L
  raw0$age_z <- zscore0(raw0$age)
  inna_map <- c("1" = "B", "2" = "C", "3" = "A", "4" = "D")
  raw0$hclus4 <- factor(inna_map[as.character(raw0$cluster_INMA)],
                        levels = c("B", "A", "C", "D"))
  raw0$fev_bar <- ave(raw0$fev, raw0$id, FUN = mean)
  raw0$wt_bar <- ave(raw0$wt, raw0$id, FUN = mean)
  raw0$fev_w <- raw0$fev - raw0$fev_bar
  raw0$wt_w <- raw0$wt - raw0$wt_bar
  raw0$fev_bar_z <- zscore0(raw0$fev_bar)
  raw0$wt_bar_z <- zscore0(raw0$wt_bar)
  ok <- complete.cases(raw0[, c("fev", "wt", "age", "sex", "cluster_INMA",
                                "IMPACTS", "SYMPTOMS", "ACTIVITY")])
  raw0 <- raw0[ok, , drop = FALSE]
  raw0$id <- droplevels(raw0$id)
  dims0 <- c("Impacts", "Symptoms", "Activity")
  y_cols <- c("IMPACTS", "SYMPTOMS", "ACTIVITY")
  long0 <- do.call(rbind, lapply(seq_along(dims0), function(d) {
    data.frame(
      id = raw0$id, time = raw0$time, time_c = raw0$time_c,
      dim = factor(dims0[d], levels = dims0),
      y = sgrq_bin0(raw0[[y_cols[d]]]), m = 24L,
      sex = raw0$sex, hclus4 = raw0$hclus4, age_z = raw0$age_z,
      fev_bar_z = raw0$fev_bar_z, wt_bar_z = raw0$wt_bar_z,
      fev_w = raw0$fev_w, wt_w = raw0$wt_w,
      stringsAsFactors = FALSE
    )
  }))
  fit0 <- BBmm(
    y ~ hclus4 + time_c + sex + age_z + fev_bar_z + wt_bar_z + fev_w + wt_w,
    random = ~ (0 + dim | id) + (0 + fev_w + wt_w | id),
    dim = "dim", corr = "unstructured", m = "m", data = long0,
    silent = TRUE, show = TRUE, maxiter = 600L
  )
  rsb <- list(
    beta = fit0$beta, fixed.vcov = fit0$fixed.vcov,
    sigma = fit0$sigma, nll = fit0$nll, conv = fit0$conv,
    fev_w = raw0$fev_w, wt_w = raw0$wt_w
  )
  saveRDS(rsb, beta_rds)
  message("Updated ", beta_rds, " with fixed.vcov")
}

beta <- rsb$beta
sig_fev <- as.numeric(summ$RS$sigma["id.fev_w"])
sig_wt  <- as.numeric(summ$RS$sigma["id.wt_w"])
Vbeta <- rsb$fixed.vcov
zcrit <- qnorm(0.975)
stopifnot(!is.null(Vbeta))

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
  if (!is.finite(s) || s < 1e-8) x - mean(x, na.rm = TRUE) else
    (x - mean(x, na.rm = TRUE)) / s
}
uid <- !duplicated(raw$id)
fev_bar_z <- zscore(raw$fev_bar[uid])
wt_bar_z <- zscore(raw$wt_bar[uid])
fev_w_all <- raw$fev_w
wt_w_all <- raw$wt_w

dims <- c("Impacts", "Symptoms", "Activity")
col_bw <- "#1b9e77"
col_ww <- "#2166ac"

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

between_line <- function(v, dm, x) {
  b <- get_b(dm, v); se_b <- get_se_b(dm, v); a0 <- get_a0(dm)
  eta <- a0 + b * x
  se <- if (is.finite(se_b)) abs(x) * se_b else rep(NA_real_, length(x))
  data.frame(
    x = x, mu = plogis(eta),
    mu_lo = if (all(is.na(se))) NA_real_ else plogis(eta - zcrit * se),
    mu_hi = if (all(is.na(se))) NA_real_ else plogis(eta + zcrit * se)
  )
}

within_fan <- function(v, dm, x, sig) {
  b <- get_b(dm, v); a0 <- get_a0(dm)
  eta0 <- a0 + b * x
  list(
    x = x,
    mu = plogis(eta0),
    mu_p1 = plogis(a0 + (b + sig) * x),
    mu_m1 = plogis(a0 + (b - sig) * x),
    mu_p2 = plogis(a0 + (b + 2 * sig) * x),
    mu_m2 = plogis(a0 + (b - 2 * sig) * x)
  )
}

fev_w_grid <- seq(quantile(fev_w_all, 0.02, na.rm = TRUE),
                  quantile(fev_w_all, 0.98, na.rm = TRUE), length.out = 80)
wt_w_grid  <- seq(quantile(wt_w_all, 0.02, na.rm = TRUE),
                  quantile(wt_w_all, 0.98, na.rm = TRUE), length.out = 80)
sd_fev_w <- sd(fev_w_all, na.rm = TRUE)
sd_wt_w <- sd(wt_w_all, na.rm = TRUE)
fev_b_grid <- seq(quantile(fev_bar_z, 0.02, na.rm = TRUE),
                  quantile(fev_bar_z, 0.98, na.rm = TRUE), length.out = 80)
wt_b_grid  <- seq(quantile(wt_bar_z, 0.02, na.rm = TRUE),
                  quantile(wt_bar_z, 0.98, na.rm = TRUE), length.out = 80)

out_png <- file.path(out_dir, "wwbw_rs_both_effects.png")
out_pdf <- file.path(out_dir, "wwbw_rs_both_effects.pdf")

plot_all <- function() {
  op <- par(mfrow = c(4, 3), mar = c(4, 4, 2.2, 0.8))
  on.exit(par(op), add = TRUE)

  for (v in c("fev_bar_z", "wt_bar_z")) {
    xlab <- if (v == "fev_bar_z") "FEV person-mean (z)" else "6MWT person-mean (z)"
    xseq <- if (v == "fev_bar_z") fev_b_grid else wt_b_grid
    for (dm in dims) {
      sp <- between_line(v, dm, xseq)
      if (all(is.na(sp$mu_lo))) {
        ylim <- range(sp$mu, na.rm = TRUE)
      } else {
        ylim <- range(c(sp$mu_lo, sp$mu_hi), na.rm = TRUE)
      }
      pad <- if (dm == dims[1] && v == "fev_bar_z") 0.08 * diff(ylim) else 0.02
      ylim <- c(max(0, ylim[1] - 0.02), min(1, ylim[2] + pad))
      plot(sp$x, sp$mu, type = "n", xlab = xlab, ylab = expression(mu),
           main = paste0(dm, " - BW linear"), ylim = ylim)
      if (!all(is.na(sp$mu_lo))) {
        polygon(c(sp$x, rev(sp$x)), c(sp$mu_lo, rev(sp$mu_hi)),
                col = adjustcolor(col_bw, 0.25), border = NA)
      }
      lines(sp$x, sp$mu, lwd = 2, col = col_bw)
      abline(v = 0, lty = 3, col = "grey50")
      grid()
      if (dm == dims[1] && v == "fev_bar_z") {
        legend("topright", c("between (linear)", "95% band"),
               col = c(col_bw, adjustcolor(col_bw, 0.25)),
               lty = c(1, NA), lwd = 2, pch = c(NA, 15), pt.cex = 1.6,
               bg = "white", box.col = "grey70", box.lwd = 0.6,
               cex = 0.72, y.intersp = 1.25, x.intersp = 0.85,
               inset = c(0.01, 0.01))
      }
    }
  }

  for (v in c("fev_w", "wt_w")) {
    sd_w <- if (v == "fev_w") sd_fev_w else sd_wt_w
    x_raw <- if (v == "fev_w") fev_w_grid else wt_w_grid
    x_z <- x_raw / sd_w
    sig_z <- (if (v == "fev_w") sig_fev else sig_wt) * sd_w
    xlab <- if (v == "fev_w") {
      "FEV within (SD units from person mean)"
    } else {
      "6MWT within (SD units from person mean)"
    }
    for (dm in dims) {
      # slopes on within-SD scale
      b_z <- get_b(dm, v) * sd_w
      a0 <- get_a0(dm)
      eta0 <- a0 + b_z * x_z
      sp <- list(
        x = x_z,
        mu = plogis(eta0),
        mu_p1 = plogis(a0 + (b_z + sig_z) * x_z),
        mu_m1 = plogis(a0 + (b_z - sig_z) * x_z),
        mu_p2 = plogis(a0 + (b_z + 2 * sig_z) * x_z),
        mu_m2 = plogis(a0 + (b_z - 2 * sig_z) * x_z)
      )
      ylim <- range(c(sp$mu_m2, sp$mu_p2), na.rm = TRUE)
      pad <- 0.06 * diff(ylim)
      ylim <- c(max(0, ylim[1] - 0.02), min(1, ylim[2] + pad))
      plot(sp$x, sp$mu, type = "n", xlab = xlab, ylab = expression(mu),
           main = paste0(dm, " - WW linear + RS"), ylim = ylim)
      polygon(c(sp$x, rev(sp$x)), c(sp$mu_m2, rev(sp$mu_p2)),
              col = adjustcolor(col_ww, 0.10), border = NA)
      polygon(c(sp$x, rev(sp$x)), c(sp$mu_m1, rev(sp$mu_p1)),
              col = adjustcolor(col_ww, 0.22), border = NA)
      lines(sp$x, sp$mu, lwd = 2.2, col = col_ww)
      abline(v = 0, lty = 3, col = "grey50")
      grid()
      if (dm == dims[1] && v == "fev_w") {
        legend("topright",
               c("pop. slope", "+/-1 SD", "+/-2 SD"),
               col = c(col_ww, adjustcolor(col_ww, 0.45),
                       adjustcolor(col_ww, 0.2)),
               lty = c(1, NA, NA), lwd = c(2.2, NA, NA),
               pch = c(NA, 15, 15), pt.cex = 1.4,
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

bw <- do.call(rbind, lapply(dims, function(dm) {
  data.frame(
    dim = dm,
    fev_bar_z = get_b(dm, "fev_bar_z"),
    wt_bar_z = get_b(dm, "wt_bar_z"),
    fev_w = get_b(dm, "fev_w"),
    wt_w = get_b(dm, "wt_w"),
    stringsAsFactors = FALSE
  )
}))
write.csv(bw, file.path(out_dir, "rs_both_coefs.csv"), row.names = FALSE)
print(bw, row.names = FALSE, digits = 4)
message("sigma FEV^w=", round(sig_fev, 5), "  6MWT^w=", round(sig_wt, 5))
