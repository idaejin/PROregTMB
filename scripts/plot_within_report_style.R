## WW-BW within smooths — Report-style partial effects
##
## Inferential target: within contrast relative to own mean (x = 0)
##   delta(x) = f(x) - f(0)
##   SE from joint Laplace Cov(beta_null, s) via TMB jointPrecision
##   (hyperparameters held at MLE; empirical-Bayes / conditional bands)
##
## Display: mu(x) = plogis(a0 + delta(x)), x on p2–p98, linear within overlay.
##
## Usage (from PROregTMB root):
##   Rscript scripts/plot_within_report_style.R

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
        file.exists(file.path(d, "R", "BBmm.R"))) {
      return(d)
    }
    nd <- dirname(d)
    if (identical(nd, d)) break
    d <- nd
  }
  d
})

out_dir <- file.path(.root, "scripts", "out_report4_reanalysis")
rds <- file.path(out_dir, "reanalysis_report4.rds")
stopifnot(file.exists(rds))

obj <- readRDS(rds)
fit <- obj$results$M3_Hps$fit
fit_lin <- obj$results$M3_Hlin$fit
stopifnot(!is.null(fit), !is.null(fit_lin))

data_path <- normalizePath(
  file.path(.root, "..", "ref", "BB-GAM", "longitudinal", "data", "COPD_All.txt"),
  mustWork = TRUE
)
raw <- read.table(data_path, header = TRUE, stringsAsFactors = FALSE)
raw$id <- factor(raw$id)
raw$fev <- as.numeric(raw$fev)
raw$wt <- as.numeric(raw$wt)
raw$fev_bar <- ave(raw$fev, raw$id, FUN = mean)
raw$wt_bar <- ave(raw$wt, raw$id, FUN = mean)
raw$fev_w <- raw$fev - raw$fev_bar
raw$wt_w <- raw$wt - raw$wt_bar

dims <- c("Impacts", "Symptoms", "Activity")
fev_grid <- seq(stats::quantile(raw$fev_w, 0.02),
                stats::quantile(raw$fev_w, 0.98), length.out = 80)
wt_grid <- seq(stats::quantile(raw$wt_w, 0.02),
               stats::quantile(raw$wt_w, 0.98), length.out = 80)

beta <- fit$beta
get_intercept <- function(dm) {
  nm <- paste0(dm, ".(Intercept)")
  if (nm %in% names(beta)) unname(beta[[nm]]) else 0
}

lin_slope <- function(v, dm) {
  nm <- paste0(dm, ".", v)
  b <- fit_lin$beta
  if (nm %in% names(b)) unname(b[[nm]]) else NA_real_
}

within_contrast <- function(lab, x, dm, level = 0.95) {
  pr <- predict_smooth(fit, which = lab, x = x, level = level,
                       se = TRUE, centered = TRUE)
  z <- stats::qnorm(1 - (1 - level) / 2)
  a0 <- get_intercept(dm)
  eta <- a0 + pr$fit
  data.frame(
    x = pr$x,
    delta = pr$fit,
    se_delta = pr$se,
    eta = eta,
    mu = stats::plogis(eta),
    mu_lo = stats::plogis(eta - z * pr$se),
    mu_hi = stats::plogis(eta + z * pr$se),
    se_method = attr(pr, "method")
  )
}

out_png <- file.path(out_dir, "within_smooths_report_style.png")
out_pdf <- file.path(out_dir, "within_smooths_report_style.pdf")

plot_panels <- function() {
  op <- par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1))
  on.exit(par(op), add = TRUE)
  for (v in c("fev_w", "wt_w")) {
    xlab <- if (identical(v, "fev_w")) {
      "FEV within (pp from person mean)"
    } else {
      "6MWT within (m from person mean)"
    }
    xseq <- if (identical(v, "fev_w")) fev_grid else wt_grid
    for (dm in dims) {
      lab <- paste0("s(", v, "):", dm)
      sp <- within_contrast(lab, xseq, dm)
      b1 <- lin_slope(v, dm)
      a0 <- get_intercept(dm)
      mu_lin <- if (is.finite(b1)) {
        stats::plogis(a0 + b1 * xseq)
      } else {
        rep(NA_real_, length(xseq))
      }

      ylim <- range(c(sp$mu_lo, sp$mu_hi, mu_lin), na.rm = TRUE)
      pad <- if (identical(dm, dims[1]) && identical(v, "fev_w")) {
        0.10 * diff(ylim)
      } else {
        0.02
      }
      ylim <- c(max(0, ylim[1] - 0.02), min(1, ylim[2] + pad))
      plot(sp$x, sp$mu, type = "n", xlab = xlab, ylab = expression(mu),
           main = dm, ylim = ylim)
      if (all(is.finite(c(sp$mu_lo, sp$mu_hi)))) {
        polygon(c(sp$x, rev(sp$x)), c(sp$mu_lo, rev(sp$mu_hi)),
                col = adjustcolor("#d95f02", alpha.f = 0.25), border = NA)
      }
      lines(sp$x, sp$mu, lwd = 2, col = "#d95f02")
      if (all(is.finite(mu_lin))) {
        lines(sp$x, mu_lin, lwd = 2, lty = 2, col = "#666666")
      }
      abline(v = 0, lty = 3, col = "grey50")
      if (identical(dm, dims[1]) && identical(v, "fev_w")) {
        legend("topright",
               c("within P-spline", "95% band", "within linear"),
               col = c("#d95f02", adjustcolor("#d95f02", 0.25), "#666666"),
               lty = c(1, NA, 2), lwd = c(2, NA, 2), pch = c(NA, 15, NA),
               pt.cex = 1.5, bg = "white", box.col = "grey70", box.lwd = 0.6,
               cex = 0.75, y.intersp = 1.25, x.intersp = 0.85,
               inset = c(0.01, 0.01))
      }
      grid()
    }
  }
}

png(out_png, width = 1100, height = 720, res = 120)
plot_panels()
dev.off()
message("Wrote ", out_png)

pdf(out_pdf, width = 9, height = 5.5)
plot_panels()
dev.off()
message("Wrote ", out_pdf)

rows <- list()
se_methods <- character(0)
for (v in c("fev_w", "wt_w")) {
  xseq <- if (identical(v, "fev_w")) fev_grid else wt_grid
  for (dm in dims) {
    lab <- paste0("s(", v, "):", dm)
    sp <- within_contrast(lab, xseq, dm)
    se_methods <- c(se_methods, unique(as.character(sp$se_method)))
    sp$cov <- v
    sp$dim <- dm
    rows[[length(rows) + 1L]] <- sp
  }
}
curves <- do.call(rbind, rows)
csv <- file.path(out_dir, "within_smooths_report_style.csv")
write.csv(curves, csv, row.names = FALSE)
message("Wrote ", csv)
message(
  "Bands = 95% for within contrast f(x)-f(0); SE method(s): ",
  paste(unique(se_methods), collapse = ", ")
)
