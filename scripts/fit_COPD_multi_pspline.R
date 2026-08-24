## Port of BB-GAM COPD multivariate P-spline (CS visit-1 / multi2)
## to PROregTMB BBmm API: domain-specific s(..., by = dim) + M3 RE.
##
## BB-GAM reference: ref/BB-GAM/R/fit_pspline_M123.R (model_flag = 3)
## Mean (per domain via dim=):
##   logit(p) ~ hclus4 + wt + fev1p + s(wt, by=dim, k=8) + s(fev1p, by=dim, k=8)
##             + u_i^(dim),   u_i ~ N(0, Sigma)   # M3
## Linear wt/fev stay in X (P-spline null-space role in BB-GAM).
##
## Usage (from PROregTMB root):
##   Rscript scripts/fit_COPD_multi_pspline.R
##   Rscript scripts/fit_COPD_multi_pspline.R --quick   # ~120 subjects
##   Rscript scripts/fit_COPD_multi_pspline.R --m1      # shared RI instead of M3

args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args
use_m1 <- "--m1" %in% args

suppressPackageStartupMessages({
  .pkg <- normalizePath(".", winslash = "/", mustWork = FALSE)
  # allow running from PROregTMB/ or PROregTMB/scripts/
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
  normalizePath(".", winslash = "/", mustWork = FALSE)
})

bb_root <- normalizePath(file.path(.root, ".."), winslash = "/", mustWork = FALSE)
data_path <- file.path(bb_root, "ref", "BB-GAM", "multi2.RData")
if (!file.exists(data_path)) {
  data_path <- file.path(.root, "..", "ref", "BB-GAM", "multi2.RData")
  data_path <- normalizePath(data_path, mustWork = FALSE)
}
stopifnot(file.exists(data_path))

out_dir <- file.path(.root, "scripts", "out_COPD_multi_pspline")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

load(data_path) # multi2
dat0 <- multi2
dat0$wt <- as.numeric(dat0$wt)
dat0$fev1p <- as.numeric(dat0$fev1p)
if (all(is.na(dat0$fev1p))) dat0$fev1p <- as.numeric(dat0$fev1crudo)
dat0$hclus4 <- factor(dat0$hclus4.aldatuta, levels = c("B", "A", "C", "D"))
dat0$id <- factor(seq_len(nrow(dat0)))

if (quick) {
  set.seed(1)
  keep <- sort(sample.int(nrow(dat0), size = min(120L, nrow(dat0))))
  dat0 <- dat0[keep, , drop = FALSE]
  dat0$id <- factor(dat0$id)
  message("Quick subset: n = ", nrow(dat0))
}

dims <- c("Impacts", "Symptoms", "Activity")
y_cols <- c("STGIMPN", "STGSYMN", "STGACTN")
m_trials <- 24L

long <- do.call(rbind, lapply(seq_along(dims), function(d) {
  data.frame(
    id = dat0$id,
    dim = factor(dims[d], levels = dims),
    y = as.integer(dat0[[y_cols[d]]]),
    m = m_trials,
    wt = dat0$wt,
    fev1p = dat0$fev1p,
    hclus4 = dat0$hclus4,
    stringsAsFactors = FALSE
  )
}))
rownames(long) <- NULL

# Standardize smooth covariates (same knot logic, stabler Laplace)
long$wt_z <- as.numeric(scale(long$wt))
long$fev1p_z <- as.numeric(scale(long$fev1p))

message(sprintf(
  "COPD long: N=%d  subjects=%d  domains=%d",
  nrow(long), nlevels(long$id), nlevels(long$dim)
))

t0 <- proc.time()[["elapsed"]]
# Linear wt/fev per domain stay in X (BB-GAM null space in X);
# s() = nonlinear deviation. Default RE = M3 (0+dim|id); --m1 for shared RI.
re_form <- if (use_m1) ~ (1 | id) else ~ (0 + dim | id)
message("RE: ", if (use_m1) "M1 shared (1|id)" else "M3 (0+dim|id) unstructured")
fit <- BBmm(
  y ~ hclus4 + wt_z + fev1p_z +
    s(wt_z, by = dim, k = 8) + s(fev1p_z, by = dim, k = 8),
  random = re_form,
  dim = "dim",
  corr = if (use_m1) "diag" else "unstructured",
  m = "m",
  data = long,
  silent = TRUE,
  show = TRUE,
  maxiter = 500
)
elapsed <- proc.time()[["elapsed"]] - t0

message("conv: ", fit$conv)
if (!is.null(fit$opt$message)) message("opt: ", fit$opt$message)
message(sprintf("time: %.1fs", elapsed))
message("nDim: ", fit$nDim)
message("lambda:\n")
print(round(fit$lambda, 4))
message("phi:\n")
print(round(fit$phi, 4))
message("sigma:\n")
print(round(fit$sigma, 4))
if (!is.null(fit$Corr$id)) {
  message("Corr (domain RE):\n")
  print(round(fit$Corr$id, 3))
}
if (!is.null(fit$Sigma$id)) {
  message("Sigma (domain RE):\n")
  print(round(fit$Sigma$id, 3))
}
message("beta (head):\n")
print(round(head(fit$beta, 12), 4))

sum_f <- vapply(fit$fhat, sum, numeric(1))
message("sum fhat (should be ~0):\n")
print(round(sum_f, 5))

# Smooth curves + bands (on standardized scale)
pdf(file.path(out_dir, "COPD_smooths.pdf"), width = 9, height = 6)
op <- par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1))
for (v in c("wt_z", "fev1p_z")) {
  for (dm in dims) {
    lab <- paste0("s(", v, "):", dm)
    if (!lab %in% fit$smooth$labels) next
    xr <- range(long[[v]], na.rm = TRUE)
    g <- seq(xr[1], xr[2], length.out = 80)
    pr <- tryCatch(
      predict_smooth(fit, which = lab, x = g, method = "marrawood"),
      error = function(e) NULL
    )
    if (is.null(pr)) next
    ylim <- range(pr$fit, pr$lwr, pr$upr, finite = TRUE)
    if (!all(is.finite(ylim))) ylim <- range(pr$fit, finite = TRUE)
    plot(pr$x, pr$fit, type = "n", ylim = ylim,
         xlab = v, ylab = "f", main = lab)
    if (all(c("lwr", "upr") %in% names(pr)) &&
        all(is.finite(c(pr$lwr, pr$upr)))) {
      polygon(c(pr$x, rev(pr$x)), c(pr$lwr, rev(pr$upr)),
              col = adjustcolor("steelblue", 0.30), border = NA)
    }
    lines(pr$x, pr$fit, col = "steelblue", lwd = 2)
    abline(h = 0, col = "grey70")
  }
}
par(op)
dev.off()
message("Wrote ", file.path(out_dir, "COPD_smooths.pdf"))

saveRDS(
  list(fit = fit, elapsed = elapsed, quick = quick, n = nlevels(long$id)),
  file.path(out_dir, "fit_COPD_multi_pspline.rds")
)
message("Wrote ", file.path(out_dir, "fit_COPD_multi_pspline.rds"))
