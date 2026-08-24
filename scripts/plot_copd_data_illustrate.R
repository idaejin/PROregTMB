#!/usr/bin/env Rscript
# Illustrative COPD panel plots (data description; ASCII labels).
#
# Usage (from PROregTMB root):
#   Rscript scripts/plot_copd_data_illustrate.R

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
})

data_path <- normalizePath(
  file.path(.root, "..", "ref", "BB-GAM", "longitudinal", "data", "COPD_All.txt"),
  mustWork = TRUE
)
out_dir <- file.path(.root, "scripts", "out_data_illustrate")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

raw <- read.table(data_path, header = TRUE, stringsAsFactors = FALSE)
raw$id <- factor(raw$id)
raw$time <- as.integer(raw$time)
raw$fev <- as.numeric(raw$fev)
raw$wt <- as.numeric(raw$wt)
raw$age <- as.numeric(raw$age)
raw$sex <- factor(raw$sex, levels = c("1", "2"), labels = c("Male", "Female"))
inna_map <- c("1" = "B", "2" = "C", "3" = "A", "4" = "D")
raw$hclus4 <- factor(inna_map[as.character(raw$cluster_INMA)],
                     levels = c("B", "A", "C", "D"))

# Panel structure uses the full visit file (matches manuscript N=1772).
# Physiology / SGRQ plots use complete cases only.
panel <- raw
panel$id <- droplevels(factor(panel$id))
panel$n_vis <- ave(rep(1L, nrow(panel)), panel$id, FUN = length)

ok <- is.finite(raw$fev) & is.finite(raw$wt) & is.finite(raw$age) &
  raw$fev > 0 & raw$wt > 0 &
  is.finite(raw$SYMPTOMS) & is.finite(raw$IMPACTS) & is.finite(raw$ACTIVITY) &
  !is.na(raw$sex) & !is.na(raw$hclus4)
dat <- raw[ok, , drop = FALSE]
dat$id <- droplevels(dat$id)

dat$fev_bar <- ave(dat$fev, dat$id, FUN = mean)
dat$wt_bar <- ave(dat$wt, dat$id, FUN = mean)
dat$fev_w <- dat$fev - dat$fev_bar
dat$wt_w <- dat$wt - dat$wt_bar
dat$n_vis <- ave(rep(1L, nrow(dat)), dat$id, FUN = length)

uid <- !duplicated(dat$id)
pers <- dat[uid, , drop = FALSE]
uid_panel <- !duplicated(panel$id)
pers_panel <- panel[uid_panel, , drop = FALSE]

col_dom <- c(Impacts = "#1b9e77", Symptoms = "#d95f02", Activity = "#7570b3")
col_phys <- c(FEV = "#2166ac", Walk = "#b2182b")

save_fig <- function(name, width_in, height_in, plot_fun, png_w = NULL, png_h = NULL) {
  stopifnot(is.function(plot_fun))
  pdf(file.path(out_dir, paste0(name, ".pdf")), width = width_in, height = height_in)
  plot_fun()
  dev.off()
  png(file.path(out_dir, paste0(name, ".png")),
      width = if (is.null(png_w)) as.integer(width_in * 120) else png_w,
      height = if (is.null(png_h)) as.integer(height_in * 120) else png_h,
      res = 120)
  plot_fun()
  dev.off()
  message("Wrote ", name)
}

# ---- 1) Panel structure (full visit file; avoid label overlap) ----
save_fig("data_panel_structure", 9, 4.5, function() {
  op <- par(mfrow = c(1, 2), mar = c(4.2, 4.5, 3.6, 1))
  on.exit(par(op), add = TRUE)

  n_id <- nlevels(panel$id)
  n_vis_tot <- nrow(panel)
  tab_n <- table(pers_panel$n_vis)
  bp <- barplot(as.numeric(tab_n), names.arg = names(tab_n),
                col = "#4d4d4d", border = NA,
                xlab = "Visits per subject", ylab = "Number of subjects",
                main = "",
                ylim = c(0, max(tab_n) * 1.18))
  text(bp, as.numeric(tab_n), labels = as.numeric(tab_n), pos = 3, cex = 0.85)
  title(main = "Follow-up length", line = 2.2, cex.main = 1.1)
  title(main = sprintf("N subjects = %d", n_id), line = 0.7, cex.main = 0.95,
        font.main = 1, col.main = "grey25")

  tab_t <- table(panel$time)
  bp2 <- barplot(as.numeric(tab_t), names.arg = paste0("t=", names(tab_t)),
                 col = "#4393c3", border = NA,
                 xlab = "Visit index", ylab = "Number of records",
                 main = "",
                 ylim = c(0, max(tab_t) * 1.18))
  text(bp2, as.numeric(tab_t), labels = as.numeric(tab_t), pos = 3, cex = 0.85)
  title(main = "Sample size by visit", line = 2.2, cex.main = 1.1)
  title(main = sprintf("N visits = %d", n_vis_tot), line = 0.7, cex.main = 0.95,
        font.main = 1, col.main = "grey25")
})

# ---- 2) SGRQ domain distributions ----
save_fig("data_sgrq_domains", 9.5, 4.5, function() {
  op <- par(mfrow = c(1, 3), mar = c(4.5, 4.2, 2.5, 0.8))
  on.exit(par(op), add = TRUE)
  for (nm in c("IMPACTS", "SYMPTOMS", "ACTIVITY")) {
    lab <- c(IMPACTS = "Impacts", SYMPTOMS = "Symptoms", ACTIVITY = "Activity")[nm]
    hist(dat[[nm]], breaks = 20, col = adjustcolor(col_dom[lab], 0.7),
         border = "white", main = paste0("SGRQ ", lab),
         xlab = "Score (0-100; higher = worse)", ylab = "Count",
         xlim = c(0, 100))
    abline(v = median(dat[[nm]]), lty = 2, lwd = 2, col = "grey20")
  }
})

# ---- 3) Physiology: raw, between, within ----
save_fig("data_wwbw_decomp", 9.5, 6.8, function() {
  op <- par(mfrow = c(2, 3), mar = c(4.2, 4.2, 2.4, 0.8))
  on.exit(par(op), add = TRUE)

  hist(dat$fev, breaks = 25, col = adjustcolor(col_phys["FEV"], 0.65),
       border = "white", main = "FEV1 % pred (raw)",
       xlab = "FEV1 % predicted", ylab = "Count")
  hist(pers$fev_bar, breaks = 20, col = adjustcolor(col_phys["FEV"], 0.65),
       border = "white", main = "FEV between (person mean)",
       xlab = "Person-mean FEV1 % pred", ylab = "Count")
  hist(dat$fev_w, breaks = 25, col = adjustcolor(col_phys["FEV"], 0.65),
       border = "white", main = "FEV within (visit - mean)",
       xlab = "Deviation (percentage points)", ylab = "Count")
  abline(v = 0, lty = 3, col = "grey40")

  hist(dat$wt, breaks = 25, col = adjustcolor(col_phys["Walk"], 0.55),
       border = "white", main = "6MWT (raw)",
       xlab = "Walk distance (m)", ylab = "Count")
  hist(pers$wt_bar, breaks = 20, col = adjustcolor(col_phys["Walk"], 0.55),
       border = "white", main = "6MWT between (person mean)",
       xlab = "Person-mean 6MWT (m)", ylab = "Count")
  hist(dat$wt_w, breaks = 25, col = adjustcolor(col_phys["Walk"], 0.55),
       border = "white", main = "6MWT within (visit - mean)",
       xlab = "Deviation (m)", ylab = "Count")
  abline(v = 0, lty = 3, col = "grey40")
})

# ---- 4a) Spaghetti trajectories (sample; mean line, no SE band) ----
set.seed(2)
ids4 <- as.character(pers$id[pers$n_vis == 4L])
ids_sp <- sample(ids4, size = min(40L, length(ids4)))
sp <- dat[as.character(dat$id) %in% ids_sp, , drop = FALSE]

save_fig("data_trajectories", 9.5, 7.2, function() {
  op <- par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.4, 0.8))
  on.exit(par(op), add = TRUE)

  plot_spag <- function(ycol, ylab, main, col) {
    yr <- range(sp[[ycol]], na.rm = TRUE)
    plot(NA, xlim = c(1, 4), ylim = yr,
         xlab = "Visit", ylab = ylab, main = main, xaxt = "n")
    axis(1, at = 1:4)
    for (id in ids_sp) {
      w <- sp[as.character(sp$id) == id, , drop = FALSE]
      w <- w[order(w$time), , drop = FALSE]
      lines(w$time, w[[ycol]], col = adjustcolor(col, 0.35), lwd = 1)
      points(w$time, w[[ycol]], col = adjustcolor(col, 0.55), pch = 16, cex = 0.45)
    }
    mu <- tapply(sp[[ycol]], sp$time, mean, na.rm = TRUE)
    lines(as.numeric(names(mu)), as.numeric(mu), lwd = 2.5, col = "black")
  }

  plot_spag("fev", "FEV1 % predicted", "FEV trajectories (n=40 with 4 visits)",
            col_phys["FEV"])
  plot_spag("wt", "6MWT (m)", "6MWT trajectories", col_phys["Walk"])
  plot_spag("ACTIVITY", "SGRQ Activity (0-100)", "Activity trajectories",
            col_dom["Activity"])
  plot_spag("IMPACTS", "SGRQ Impacts (0-100)", "Impacts trajectories",
            col_dom["Impacts"])
})

# ---- 4b) Boxplots by visit (full sample; kept as separate figure) ----
save_fig("data_by_visit_boxplots", 9.5, 7.2, function() {
  op <- par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.4, 0.8))
  on.exit(par(op), add = TRUE)

  plot_box_visit <- function(ycol, ylab, main, col) {
    dat$time_f <- factor(dat$time)
    boxplot(dat[[ycol]] ~ dat$time_f,
            col = adjustcolor(col, 0.55), border = "grey25",
            xlab = "Visit", ylab = ylab, main = main,
            outline = TRUE, outpch = 16, outcex = 0.35,
            outcol = adjustcolor(col, 0.35),
            medcol = "black", medlwd = 2)
    mu <- tapply(dat[[ycol]], dat$time, mean, na.rm = TRUE)
    points(seq_along(mu), as.numeric(mu), pch = 18, cex = 1.2, col = "black")
    if (identical(ycol, "fev")) {
      legend("bottomleft", c("median (box)", "mean"),
             pch = c(NA, 18), lty = c(1, NA), lwd = c(2, NA),
             col = "black", bty = "n", cex = 0.75)
    }
  }

  plot_box_visit("fev", "FEV1 % predicted", "FEV by visit",
                 col_phys["FEV"])
  plot_box_visit("wt", "6MWT (m)", "6MWT by visit", col_phys["Walk"])
  plot_box_visit("ACTIVITY", "SGRQ Activity (0-100)", "Activity by visit",
                 col_dom["Activity"])
  plot_box_visit("IMPACTS", "SGRQ Impacts (0-100)", "Impacts by visit",
                 col_dom["Impacts"])
})

# ---- 5) Between vs within association (scatter) ----
save_fig("data_between_within_scatter", 9.5, 6.2, function() {
  op <- par(mfrow = c(2, 3), mar = c(4.2, 4.2, 2.4, 0.8))
  on.exit(par(op), add = TRUE)
  act_bar <- ave(dat$ACTIVITY, dat$id, FUN = mean)

  plot(pers$fev_bar, act_bar[uid], pch = 16, cex = 0.5,
       col = adjustcolor(col_dom["Activity"], 0.4),
       xlab = "Person-mean FEV1 % pred", ylab = "Person-mean Activity",
       main = "Between: FEV vs Activity")
  plot(pers$wt_bar, act_bar[uid], pch = 16, cex = 0.5,
       col = adjustcolor(col_dom["Activity"], 0.4),
       xlab = "Person-mean 6MWT (m)", ylab = "Person-mean Activity",
       main = "Between: 6MWT vs Activity")
  plot(pers$fev_bar, pers$wt_bar, pch = 16, cex = 0.5,
       col = adjustcolor("grey20", 0.35),
       xlab = "Person-mean FEV1 % pred", ylab = "Person-mean 6MWT (m)",
       main = "Between: FEV vs 6MWT")

  plot(dat$fev_w, dat$ACTIVITY, pch = 16, cex = 0.35,
       col = adjustcolor(col_dom["Activity"], 0.25),
       xlab = "FEV within (pp)", ylab = "Activity score",
       main = "Within: FEV deviation vs Activity")
  abline(v = 0, lty = 3, col = "grey50")
  plot(dat$wt_w, dat$ACTIVITY, pch = 16, cex = 0.35,
       col = adjustcolor(col_dom["Activity"], 0.25),
       xlab = "6MWT within (m)", ylab = "Activity score",
       main = "Within: 6MWT deviation vs Activity")
  abline(v = 0, lty = 3, col = "grey50")
  plot(dat$fev_w, dat$wt_w, pch = 16, cex = 0.35,
       col = adjustcolor("grey20", 0.2),
       xlab = "FEV within (pp)", ylab = "6MWT within (m)",
       main = "Within: FEV vs 6MWT deviations")
  abline(h = 0, v = 0, lty = 3, col = "grey50")
})

# ---- 6) Domain joint view at baseline (time==1) ----
b1 <- dat[dat$time == 1L, , drop = FALSE]
save_fig("data_domains_joint", 9, 4.5, function() {
  op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1))
  on.exit(par(op), add = TRUE)
  plot(b1$IMPACTS, b1$ACTIVITY, pch = 16, cex = 0.5,
       col = adjustcolor("#4d4d4d", 0.35),
       xlab = "Impacts", ylab = "Activity",
       main = "Baseline visit: Impacts vs Activity",
       xlim = c(0, 100), ylim = c(0, 100))
  abline(0, 1, lty = 3, col = "grey50")
  plot(b1$SYMPTOMS, b1$ACTIVITY, pch = 16, cex = 0.5,
       col = adjustcolor("#4d4d4d", 0.35),
       xlab = "Symptoms", ylab = "Activity",
       main = "Baseline visit: Symptoms vs Activity",
       xlim = c(0, 100), ylim = c(0, 100))
  abline(0, 1, lty = 3, col = "grey50")
})

# ---- 7) Covariate snapshot (age, sex, INMA) ----
save_fig("data_covariates", 9.5, 4.2, function() {
  op <- par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2.5, 1))
  on.exit(par(op), add = TRUE)
  hist(pers$age, breaks = 15, col = "grey70", border = "white",
       main = "Age (person level)", xlab = "Age (years)", ylab = "Subjects")
  tab_s <- table(pers$sex)
  barplot(tab_s, col = c("#4393c3", "#d6604d"), border = NA,
          main = "Sex", ylab = "Subjects", xlab = "")
  tab_c <- table(pers$hclus4)
  barplot(tab_c, col = "#74add1", border = NA,
          main = "INMA phenotype", ylab = "Subjects", xlab = "Cluster")
})

message("All figures in ", out_dir)
writeLines(
  c(
    sprintf("n_subjects=%d", nlevels(dat$id)),
    sprintf("n_visits=%d", nrow(dat)),
    sprintf("visits_per_id: %s",
            paste(sprintf("%s=%s", names(table(pers$n_vis)),
                          as.integer(table(pers$n_vis))), collapse = ", "))
  ),
  file.path(out_dir, "README.txt")
)
