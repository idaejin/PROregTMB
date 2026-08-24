#!/usr/bin/env Rscript
# Sensitivity: BW linear + WW linear + random slopes on within (fev_w, wt_w),
# keeping M3 domain-specific random intercepts.
#
# Data prep matches scripts/reanalyze_report4_long.R
#
# Usage (from package root):
#   Rscript scripts/fit_wwbw_rs_within.R
#   Rscript scripts/fit_wwbw_rs_within.R --quick

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

out_dir <- file.path(.root, "scripts", "out_wwbw_rs_within")
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
inna_map <- c("1" = "B", "2" = "C", "3" = "A", "4" = "D")
raw$hclus4 <- factor(inna_map[as.character(raw$cluster_INMA)],
                     levels = c("B", "A", "C", "D"))
raw$fev_bar <- ave(raw$fev, raw$id, FUN = mean)
raw$wt_bar <- ave(raw$wt, raw$id, FUN = mean)
raw$fev_w <- raw$fev - raw$fev_bar
raw$wt_w <- raw$wt - raw$wt_bar
raw$fev_bar_z <- zscore(raw$fev_bar)
raw$wt_bar_z <- zscore(raw$wt_bar)

# Drop incomplete rows needed for WW–BW
ok <- complete.cases(raw[, c("fev", "wt", "age", "sex", "cluster_INMA",
                             "IMPACTS", "SYMPTOMS", "ACTIVITY")])
raw <- raw[ok, , drop = FALSE]
raw$id <- droplevels(raw$id)

if (quick) {
  set.seed(1)
  keep_id <- sample(levels(raw$id), size = min(80L, nlevels(raw$id)))
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
    fev_bar_z = raw$fev_bar_z, wt_bar_z = raw$wt_bar_z,
    fev_w = raw$fev_w, wt_w = raw$wt_w,
    stringsAsFactors = FALSE
  )
}))
rownames(long) <- NULL
message(sprintf("Long: N=%d  subjects=%d  visits(raw)=%d",
                nrow(long), nlevels(long$id), nrow(raw)))

form_lin <- y ~ hclus4 + time_c + sex + age_z +
  fev_bar_z + wt_bar_z + fev_w + wt_w

fit_one <- function(label, random, maxiter = 500L) {
  message("\n==== ", label, " ====")
  message("random: ", paste(deparse(random), collapse = " "))
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    BBmm(
      form_lin, random = random, dim = "dim", corr = "unstructured",
      m = "m", data = long, silent = TRUE, show = TRUE, maxiter = maxiter
    ),
    error = function(e) e
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  if (inherits(fit, "error")) {
    message("FAILED: ", conditionMessage(fit))
    return(list(label = label, ok = FALSE, error = conditionMessage(fit),
                elapsed = elapsed))
  }
  nll <- as.numeric(fit$nll)
  phi <- as.numeric(fit$phi)
  names(phi) <- fit$dim_names
  message(sprintf("conv=%s  nll=%.4f  elapsed=%.1fs", fit$conv, nll, elapsed))
  message("phi: ", paste(sprintf("%s=%.4f", names(phi), phi), collapse = "  "))
  if (!is.null(fit$sigma) && length(fit$sigma)) {
    message("sigma: ", paste(sprintf("%.5f", fit$sigma), collapse = ", "))
  }
  if (!is.null(fit$Corr)) {
    message("Corr:")
    print(fit$Corr)
  }
  list(
    label = label, ok = TRUE, conv = fit$conv, nll = nll, phi = phi,
    sigma = fit$sigma, Corr = fit$Corr, Sigma = fit$Sigma,
    elapsed = elapsed, fit = fit
  )
}

res_hlin <- fit_one("M3_Hlin", random = ~ (0 + dim | id), maxiter = 400L)

res_rs <- fit_one(
  "M3_Hlin_RS_within",
  random = ~ (0 + dim | id) + (0 + fev_w + wt_w | id),
  maxiter = 600L
)

# Heavier: domain-specific within slopes
res_rs_dom <- fit_one(
  "M3_Hlin_RS_within_bydim",
  random = ~ (0 + dim | id) + (0 + dim:fev_w + dim:wt_w | id),
  maxiter = 700L
)

summarize <- function(r) {
  if (!isTRUE(r$ok)) {
    return(data.frame(model = r$label, ok = FALSE, conv = NA_character_,
                      nll = NA_real_, elapsed = r$elapsed, note = r$error,
                      stringsAsFactors = FALSE))
  }
  data.frame(
    model = r$label, ok = TRUE, conv = as.character(r$conv),
    nll = r$nll, elapsed = r$elapsed, note = "",
    stringsAsFactors = FALSE
  )
}

tab <- rbind(summarize(res_hlin), summarize(res_rs), summarize(res_rs_dom))
if (isTRUE(res_hlin$ok)) {
  tab$dnll_vs_Hlin <- tab$nll - res_hlin$nll
}
message("\n========== COMPARISON ==========\n")
print(tab, row.names = FALSE)

# Also print reference Hps nll from prior reanalysis if present
cmp_path <- file.path(.root, "scripts", "out_report4_reanalysis", "compare_nll.csv")
if (file.exists(cmp_path)) {
  message("\nPrior reanalysis nll (domain-specific phi):\n")
  print(read.csv(cmp_path))
}

saveRDS(
  list(
    quick = quick, n_id = nlevels(long$id), n_visit = nrow(raw),
    table = tab,
    Hlin = if (isTRUE(res_hlin$ok)) {
      list(nll = res_hlin$nll, phi = res_hlin$phi, sigma = res_hlin$sigma,
           Corr = res_hlin$Corr, Sigma = res_hlin$Sigma, conv = res_hlin$conv)
    } else NULL,
    RS = if (isTRUE(res_rs$ok)) {
      list(nll = res_rs$nll, phi = res_rs$phi, sigma = res_rs$sigma,
           Corr = res_rs$Corr, Sigma = res_rs$Sigma, conv = res_rs$conv)
    } else NULL,
    RS_bydim = if (isTRUE(res_rs_dom$ok)) {
      list(nll = res_rs_dom$nll, phi = res_rs_dom$phi, sigma = res_rs_dom$sigma,
           Corr = res_rs_dom$Corr, Sigma = res_rs_dom$Sigma,
           conv = res_rs_dom$conv)
    } else NULL
  ),
  file.path(out_dir, if (quick) "rs_within_quick.rds" else "rs_within_full.rds")
)
write.csv(tab, file.path(out_dir, if (quick) "compare_quick.csv" else "compare_full.csv"),
          row.names = FALSE)
message("Wrote ", out_dir)
