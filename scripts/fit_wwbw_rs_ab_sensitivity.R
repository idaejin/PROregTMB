#!/usr/bin/env Rscript
# Sensitivity: allow Corr(a, b) via one unstructured 5x5 RE block
#   random = ~ (0 + dim + fev_w + wt_w | id)
# Compare to preferred RS (block-diagonal a ⊥ b) from out_wwbw_rs_within/.
#
# Usage (from PROregTMB root):
#   Rscript scripts/fit_wwbw_rs_ab_sensitivity.R
#   Rscript scripts/fit_wwbw_rs_ab_sensitivity.R --quick

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
  mustWork = TRUE
)
out_dir <- file.path(.root, "scripts", "out_wwbw_rs_ab_sens")
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

# Sanity: single block should have q = 5
re_chk <- PROregTMB:::.build_re_structure(
  random = ~ (0 + dim + fev_w + wt_w | id),
  data = long, corr = "unstructured", nObs = nrow(long)
)
message("RE check: n_blocks=", re_chk$n_blocks,
        " q=", re_chk$block_q[1],
        " terms=", paste(re_chk$blocks[[1]]$term_names, collapse = ","),
        " n_theta=", length(re_chk$theta0))
stopifnot(re_chk$n_blocks == 1L, re_chk$block_q[1] == 5L)

message("\n==== RS_ab (5x5 unstructured) ====")
t0 <- proc.time()[["elapsed"]]
fit_ab <- tryCatch(
  BBmm(
    form_lin,
    random = ~ (0 + dim + fev_w + wt_w | id),
    dim = "dim", corr = "unstructured",
    m = "m", data = long, silent = TRUE, show = TRUE,
    maxiter = if (quick) 200L else 800L
  ),
  error = function(e) e
)
elapsed <- proc.time()[["elapsed"]] - t0
if (inherits(fit_ab, "error")) {
  message("FAILED: ", conditionMessage(fit_ab))
  quit(status = 1)
}

nll_ab <- as.numeric(fit_ab$nll)
q_ab <- length(fit_ab$opt$par)
aic_ab <- 2 * nll_ab + 2 * q_ab
message(sprintf("conv=%s  nll=%.4f  q=%d  AIC=%.1f  elapsed=%.1fs",
                fit_ab$conv, nll_ab, q_ab, aic_ab, elapsed))
message("sigma: ", paste(sprintf("%s=%.5f", names(fit_ab$sigma), fit_ab$sigma),
                         collapse = ", "))
Corr_ab <- fit_ab$Corr[[1]]
message("Corr (5x5):")
print(round(Corr_ab, 3))

# Cross-block correlations a vs b (rows/cols 1:3 vs 4:5 if term order is dim*, fev_w, wt_w)
tn <- colnames(Corr_ab)
a_idx <- grep("dim|Impacts|Symptoms|Activity", tn, ignore.case = TRUE)
# term names from model.matrix(~0+dim+fev_w+wt_w) typically dimImpacts,...
if (!length(a_idx)) a_idx <- 1:3
b_idx <- which(grepl("fev_w|wt_w", tn))
if (!length(b_idx)) b_idx <- 4:5
C_block <- Corr_ab[a_idx, b_idx, drop = FALSE]
message("Cross Corr(a,b):")
print(round(C_block, 3))

# Preferred RS reference
rs_ref <- file.path(.root, "scripts", "out_wwbw_rs_within", "rs_fit_betas.rds")
rs_nll <- 13935.327387811
rs_q <- 45L
rs_beta <- NULL
if (file.exists(rs_ref) && !quick) {
  rs <- readRDS(rs_ref)
  rs_nll <- as.numeric(rs$nll)
  rs_beta <- rs$beta
  message(sprintf("RS reference nll=%.4f (from rs_fit_betas.rds)", rs_nll))
} else if (quick) {
  message("Quick mode: fitting RS block-diag for local comparison...")
  fit_rs <- BBmm(
    form_lin,
    random = ~ (0 + dim | id) + (0 + fev_w + wt_w | id),
    dim = "dim", corr = "unstructured",
    m = "m", data = long, silent = TRUE, show = TRUE, maxiter = 200L
  )
  rs_nll <- as.numeric(fit_rs$nll)
  rs_q <- length(fit_rs$opt$par)
  rs_beta <- fit_rs$beta
}

aic_rs <- 2 * rs_nll + 2 * rs_q
beta_ab <- fit_ab$beta

# Physiology slopes comparison
phys <- c(
  "Impacts.fev_bar_z", "Impacts.wt_bar_z", "Impacts.fev_w", "Impacts.wt_w",
  "Symptoms.fev_bar_z", "Symptoms.wt_bar_z", "Symptoms.fev_w", "Symptoms.wt_w",
  "Activity.fev_bar_z", "Activity.wt_bar_z", "Activity.fev_w", "Activity.wt_w"
)
phys <- phys[phys %in% names(beta_ab)]
cmp <- data.frame(
  coef = phys,
  RS = if (!is.null(rs_beta)) unname(rs_beta[phys]) else NA_real_,
  RS_ab = unname(beta_ab[phys]),
  stringsAsFactors = FALSE
)
if (!is.null(rs_beta)) {
  cmp$diff <- cmp$RS_ab - cmp$RS
  cmp$rel <- cmp$diff / pmax(abs(cmp$RS), 1e-8)
}

tab <- data.frame(
  model = c("RS (a ⊥ b)", "RS_ab (5x5)"),
  nll = c(rs_nll, nll_ab),
  q = c(rs_q, q_ab),
  AIC = c(aic_rs, aic_ab),
  conv = c("yes", fit_ab$conv),
  stringsAsFactors = FALSE
)
tab$dnll_vs_RS <- tab$nll - rs_nll
tab$dAIC_vs_RS <- tab$AIC - aic_rs

write.csv(tab, file.path(out_dir, "compare_ab_sens.csv"), row.names = FALSE)
write.csv(cmp, file.path(out_dir, "beta_phys_compare.csv"), row.names = FALSE)
write.csv(as.data.frame(as.table(Corr_ab)),
          file.path(out_dir, "corr_5x5_long.csv"), row.names = FALSE)
write.csv(as.data.frame(as.table(C_block)),
          file.path(out_dir, "corr_cross_ab.csv"), row.names = FALSE)

saveRDS(
  list(
    quick = quick, nll = nll_ab, q = q_ab, AIC = aic_ab, conv = fit_ab$conv,
    elapsed = elapsed, sigma = fit_ab$sigma, Corr = Corr_ab, C_cross = C_block,
    beta = beta_ab, compare = tab, beta_phys = cmp,
    # keep opt$par length; drop heavy obj if huge
    opt_par = fit_ab$opt$par,
    maxgrad = if (!is.null(fit_ab$sdreport$gradient.fixed))
      max(abs(fit_ab$sdreport$gradient.fixed)) else NA_real_
  ),
  file.path(out_dir, "rs_ab_sens.rds")
)

message("\n========== SUMMARY ==========")
print(tab, row.names = FALSE, digits = 4)
message("\nPhysiology fixed effects:")
print(cmp, row.names = FALSE, digits = 4)
message("\nCross |Corr(a,b)| max = ",
        sprintf("%.3f", max(abs(C_block), na.rm = TRUE)))
message("Wrote outputs in ", out_dir)
