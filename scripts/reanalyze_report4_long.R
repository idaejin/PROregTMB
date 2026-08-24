## Reanalysis of Report 4 (report_long_mv_tmb.pdf) with PROregTMB API
## Models (all long + dim=, M3 unless noted):
##   L1  linear contemporaneous: hclus4 + time_c + sex + age_z + fev_z + wt_z
##   M1  same mean, shared RI (1|id) — ranking check vs M3
##   P   pooled P-spline: + s(fev_z, by=dim) + s(wt_z, by=dim)
##   Hlin WW–BW linear: fev_bar_z + wt_bar_z + fev_w + wt_w
##   Hps  WW–BW semi:   fev_bar_z + wt_bar_z + s(fev_w,by=dim) + s(wt_w,by=dim)
##
## Note: Report 4 uses *shared* phi; PROregTMB uses domain-specific phi_ell.
##       BB-GAM phi is brms-style concentration; PROregTMB phi is 1/concentration scale.
##
## Usage (from PROregTMB root):
##   Rscript scripts/reanalyze_report4_long.R
##   Rscript scripts/reanalyze_report4_long.R --quick

args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args

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

data_path <- normalizePath(
  file.path(.root, "..", "ref", "BB-GAM", "longitudinal", "data", "COPD_All.txt"),
  mustWork = FALSE
)
stopifnot(file.exists(data_path))

out_dir <- file.path(.root, "scripts", "out_report4_reanalysis")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sgrq_bin <- function(x) pmin(24L, as.integer(floor(as.numeric(x) / 4)))
zscore <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s < 1e-8) return(x - mean(x, na.rm = TRUE))
  (x - mean(x, na.rm = TRUE)) / s
}

raw <- read.table(data_path, header = TRUE, stringsAsFactors = FALSE)
raw$id <- factor(raw$id)
raw$time <- as.integer(raw$time)
raw$sex <- factor(raw$sex, levels = c("1", "2"), labels = c("Male", "Female"))
raw$fev <- as.numeric(raw$fev)
raw$wt <- as.numeric(raw$wt)
raw$age <- as.numeric(raw$age)
raw$time_c <- raw$time - 1L
raw$age_z <- zscore(raw$age)
raw$fev_z <- zscore(raw$fev)
raw$wt_z <- zscore(raw$wt)
inna_map <- c("1" = "B", "2" = "C", "3" = "A", "4" = "D")
raw$hclus4 <- factor(inna_map[as.character(raw$cluster_INMA)],
                     levels = c("B", "A", "C", "D"))
raw$fev_bar <- ave(raw$fev, raw$id, FUN = mean)
raw$wt_bar <- ave(raw$wt, raw$id, FUN = mean)
raw$fev_w <- raw$fev - raw$fev_bar
raw$wt_w <- raw$wt - raw$wt_bar
raw$fev_bar_z <- zscore(raw$fev_bar)
raw$wt_bar_z <- zscore(raw$wt_bar)

if (quick) {
  set.seed(1)
  keep_id <- sample(levels(raw$id), size = min(150L, nlevels(raw$id)))
  raw <- raw[raw$id %in% keep_id, , drop = FALSE]
  raw$id <- droplevels(raw$id)
  message("Quick subset: ", nlevels(raw$id), " subjects, ", nrow(raw), " visits")
}

dims <- c("Impacts", "Symptoms", "Activity")
y_cols <- c("IMPACTS", "SYMPTOMS", "ACTIVITY")
long <- do.call(rbind, lapply(seq_along(dims), function(d) {
  data.frame(
    id = raw$id, time = raw$time, time_c = raw$time_c,
    dim = factor(dims[d], levels = dims),
    y = sgrq_bin(raw[[y_cols[d]]]), m = 24L,
    sex = raw$sex, hclus4 = raw$hclus4, age_z = raw$age_z,
    fev_z = raw$fev_z, wt_z = raw$wt_z,
    fev_bar_z = raw$fev_bar_z, wt_bar_z = raw$wt_bar_z,
    fev_w = raw$fev_w, wt_w = raw$wt_w,
    stringsAsFactors = FALSE
  )
}))
rownames(long) <- NULL
message(sprintf("Long: N=%d  subjects=%d  visits(raw)=%d",
                nrow(long), nlevels(long$id), nrow(raw)))

fit_one <- function(label, form, random, corr = "unstructured", maxiter = 400L) {
  message("\n===== ", label, " =====")
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    BBmm(
      form, random = random, dim = "dim", corr = corr,
      m = "m", data = long, silent = TRUE, show = TRUE, maxiter = maxiter
    ),
    error = function(e) e
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  if (inherits(fit, "error")) {
    message("FAILED: ", conditionMessage(fit))
    return(list(label = label, ok = FALSE, message = conditionMessage(fit),
                elapsed = elapsed))
  }
  message(sprintf("conv=%s  nll=%.2f  time=%.1fs", fit$conv, fit$nll, elapsed))
  list(
    label = label, ok = TRUE, conv = fit$conv, nll = fit$nll,
    elapsed = elapsed, beta = fit$beta, phi = fit$phi, sigma = fit$sigma,
    Corr = fit$Corr$id, Sigma = fit$Sigma$id, lambda = fit$lambda,
    fit = fit
  )
}

re_m3 <- ~ (0 + dim | id)
re_m1 <- ~ (1 | id)

results <- list()
results$M1_lin <- fit_one(
  "M1 linear contemporaneous",
  y ~ hclus4 + time_c + sex + age_z + fev_z + wt_z,
  re_m1, corr = "diag", maxiter = 300L
)
results$M3_lin <- fit_one(
  "M3 linear contemporaneous (Report Table 4)",
  y ~ hclus4 + time_c + sex + age_z + fev_z + wt_z,
  re_m3, maxiter = 400L
)
results$M3_pool <- fit_one(
  "M3 pooled P-spline raw fev/wt (Report Sec 4.4)",
  y ~ hclus4 + time_c + sex + age_z +
    s(fev_z, by = dim, k = 8) + s(wt_z, by = dim, k = 8),
  re_m3, maxiter = 500L
)
results$M3_Hlin <- fit_one(
  "M3 WW-BW linear within (Report hybrid linear)",
  y ~ hclus4 + time_c + sex + age_z +
    fev_bar_z + wt_bar_z + fev_w + wt_w,
  re_m3, maxiter = 400L
)
results$M3_Hps <- fit_one(
  "M3 WW-BW P-spline within (Report preferred)",
  y ~ hclus4 + time_c + sex + age_z +
    fev_bar_z + wt_bar_z +
    s(fev_w, by = dim, k = 8) + s(wt_w, by = dim, k = 8),
  re_m3, maxiter = 500L
)

# ---- comparison table vs Report 4 Laplace nll (shared-phi BB-GAM) ----
report_nll <- c(
  M1_lin = 14342.2,
  M3_lin = 14139.4,
  M3_pool = 14091.2,
  M3_Hlin = 14120.9,
  M3_Hps = 14089.4
)

cmp <- do.call(rbind, lapply(names(results), function(nm) {
  r <- results[[nm]]
  data.frame(
    model = nm,
    label = r$label,
    ok = isTRUE(r$ok),
    conv = if (isTRUE(r$ok)) r$conv else NA_character_,
    nll_PROregTMB = if (isTRUE(r$ok)) r$nll else NA_real_,
    nll_Report4 = unname(report_nll[nm]),
    elapsed_s = r$elapsed,
    stringsAsFactors = FALSE
  )
}))
cmp$delta_nll <- cmp$nll_PROregTMB - cmp$nll_Report4

message("\n========== COMPARISON ==========\n")
print(cmp, row.names = FALSE)

# Fixed effects excerpt for M3 linear (vs Table 4)
if (isTRUE(results$M3_lin$ok)) {
  message("\n--- M3 linear beta (compare Table 4) ---\n")
  print(round(results$M3_lin$beta, 3))
  message("\nphi: ", paste(round(results$M3_lin$phi, 4), collapse = ", "))
  message("\nCorr:\n")
  print(round(results$M3_lin$Corr, 3))
  message("\nsigma:\n")
  print(round(results$M3_lin$sigma, 3))
}

if (isTRUE(results$M3_Hps$ok)) {
  message("\n--- M3 WW-BW preferred: between/within betas + lambda ---\n")
  b <- results$M3_Hps$beta
  keep <- grepl("fev_bar|wt_bar|hclus4|time_c|Intercept|sex|age", names(b))
  print(round(b[keep], 3))
  if (!is.null(results$M3_Hps$lambda)) {
    message("\nlambda (within smooths):\n")
    print(round(results$M3_Hps$lambda, 2))
  }
  message("\nCorr:\n")
  print(round(results$M3_Hps$Corr, 3))
}

# Within smooth plots for preferred model
if (isTRUE(results$M3_Hps$ok)) {
  fit <- results$M3_Hps$fit
  png(file.path(out_dir, "within_smooths.png"), width = 1100, height = 720, res = 120)
  op <- par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1))
  for (v in c("fev_w", "wt_w")) {
    for (dm in dims) {
      lab <- paste0("s(", v, "):", dm)
      if (!lab %in% fit$smooth$labels) next
      spec <- fit$smooth$specs[[match(lab, fit$smooth$labels)]]
      g <- seq(spec$xl, spec$xr, length.out = 80)
      pr <- tryCatch(predict_smooth(fit, which = lab, x = g), error = function(e) NULL)
      if (is.null(pr)) next
      ylim <- range(pr$fit, pr$lwr, pr$upr, 0, finite = TRUE)
      if (!all(is.finite(ylim))) ylim <- range(pr$fit, 0, finite = TRUE)
      plot(pr$x, pr$fit, type = "n", ylim = ylim, xlab = v, ylab = "f", main = lab)
      if (!is.null(pr$lwr) && all(is.finite(c(pr$lwr, pr$upr)))) {
        polygon(c(pr$x, rev(pr$x)), c(pr$lwr, rev(pr$upr)),
                col = adjustcolor("darkorange2", 0.3), border = NA)
      }
      lines(pr$x, pr$fit, col = "darkorange2", lwd = 2)
      abline(h = 0, col = "grey70")
    }
  }
  par(op)
  dev.off()
  message("Wrote ", file.path(out_dir, "within_smooths.png"))
}

saveRDS(
  list(cmp = cmp, results = results, quick = quick, n_id = nlevels(long$id),
       n_obs = nrow(long), report_nll = report_nll),
  file.path(out_dir, "reanalysis_report4.rds")
)
write.csv(cmp, file.path(out_dir, "compare_nll.csv"), row.names = FALSE)
message("Wrote ", file.path(out_dir, "reanalysis_report4.rds"))
