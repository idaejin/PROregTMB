#!/usr/bin/env Rscript
# Fit M2 (domain RI, diagonal G) vs existing M3 numbers.
# Means: Visit-level and Linear WW-BW.
#
# Usage (from PROregTMB root):
#   Rscript scripts/fit_m2_compare.R

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

data_path <- normalizePath(
  file.path(.root, "..", "ref", "BB-GAM", "longitudinal", "data", "COPD_All.txt"),
  mustWork = TRUE
)
out_dir <- file.path(.root, "scripts", "out_m2_compare")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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
message(sprintf("Long: N=%d  subjects=%d", nrow(long), nlevels(long$id)))

re_dom <- ~ (0 + dim | id)
form_visit <- y ~ hclus4 + time_c + sex + age_z + fev_z + wt_z
form_linear <- y ~ hclus4 + time_c + sex + age_z +
  fev_bar_z + wt_bar_z + fev_w + wt_w

fit_one <- function(label, form, corr, maxiter = 400L) {
  message("\n==== ", label, "  corr=", corr, " ====")
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    BBmm(form, random = re_dom, dim = "dim", corr = corr,
         m = "m", data = long, silent = TRUE, show = TRUE, maxiter = maxiter),
    error = function(e) e
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  if (inherits(fit, "error")) {
    message("FAILED: ", conditionMessage(fit))
    return(list(label = label, ok = FALSE, error = conditionMessage(fit),
                elapsed = elapsed))
  }
  message(sprintf("conv=%s  nll=%.4f  elapsed=%.1fs", fit$conv, fit$nll, elapsed))
  message("sigma: ", paste(round(fit$sigma, 4), collapse = ", "))
  if (!is.null(fit$Corr)) print(fit$Corr)
  list(label = label, ok = TRUE, conv = fit$conv, nll = as.numeric(fit$nll),
       sigma = fit$sigma, Corr = fit$Corr, elapsed = elapsed)
}

# Reference M3 from prior reanalysis CSV
ref <- read.csv(file.path(.root, "scripts", "out_report4_reanalysis", "compare_nll.csv"))
nll_visit_m3 <- ref$nll_PROregTMB[ref$model == "M3_lin"]
nll_linear_m3 <- ref$nll_PROregTMB[ref$model == "M3_Hlin"]

res_visit_m2 <- fit_one("Visit-M2", form_visit, corr = "diag", maxiter = 400L)
res_linear_m2 <- fit_one("Linear-M2", form_linear, corr = "diag", maxiter = 400L)

tab <- data.frame(
  model = c("Visit-shared (M1)", "Visit-M2", "Visit (M3)",
            "Linear-M2", "Linear (M3)"),
  nll = c(
    ref$nll_PROregTMB[ref$model == "M1_lin"],
    if (isTRUE(res_visit_m2$ok)) res_visit_m2$nll else NA_real_,
    nll_visit_m3,
    if (isTRUE(res_linear_m2$ok)) res_linear_m2$nll else NA_real_,
    nll_linear_m3
  ),
  stringsAsFactors = FALSE
)
tab$delta_vs_Visit_M3 <- tab$nll - nll_visit_m3
tab$delta_vs_Linear_M3 <- tab$nll - nll_linear_m3
print(tab)

write.csv(tab, file.path(out_dir, "m2_compare.csv"), row.names = FALSE)
saveRDS(
  list(table = tab, Visit_M2 = res_visit_m2, Linear_M2 = res_linear_m2,
       Visit_M3_nll = nll_visit_m3, Linear_M3_nll = nll_linear_m3),
  file.path(out_dir, "m2_compare.rds")
)
message("Wrote ", out_dir)
