# Small comparison: PROregTMB versus PROreg
# Run from the package root with:
#   Rscript scripts/compare_PROreg.R

if (!requireNamespace("PROreg", quietly = TRUE)) {
  stop("Install the PROreg package before running this script.", call. = FALSE)
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Install devtools to load the local PROregTMB package.", call. = FALSE)
}

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
pkg_root <- if (length(file_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
} else {
  normalizePath(".")
}

devtools::load_all(pkg_root, quiet = TRUE)
options(PROregTMB.verbose_compile = FALSE)

compare <- function(label, fit_tmb, fit_proreg, tmb_time, proreg_time,
                    tmb_coef, proreg_coef) {
  cmp <- data.frame(
    model = label,
    parameter = names(tmb_coef),
    PROregTMB = unname(tmb_coef),
    PROreg = unname(proreg_coef),
    difference = unname(tmb_coef - proreg_coef),
    row.names = NULL,
    check.names = FALSE
  )
  print(cmp, digits = 4)
  cat(
    "Time (PROregTMB / PROreg): ", round(tmb_time, 3), " / ",
    round(proreg_time, 3), " seconds; ratio = ",
    round(tmb_time / max(proreg_time, 1e-8), 2), "\n",
    sep = ""
  )
  cat(
    "Convergence (PROregTMB / PROreg): ", fit_tmb$conv, " / ",
    fit_proreg$conv, "\n\n", sep = ""
  )
  invisible(cmp)
}

set.seed(2026)
m <- 10L

# Marginal beta-binomial regression.
n <- 500L
x <- rnorm(n)
beta <- c("(Intercept)" = -0.5, x = 0.8)
phi <- 0.25
p <- plogis(beta[1] + beta[2] * x)
y <- PROregTMB::rBB(n, m, p, phi)
dat_reg <- data.frame(y = y, x = x)

# Warm up TMB so compilation is not included in the comparison.
invisible(PROregTMB::BBreg(y ~ x, m = m, data = dat_reg, silent = TRUE))

tmb_start <- proc.time()[["elapsed"]]
fit_tmb_reg <- PROregTMB::BBreg(y ~ x, m = m, data = dat_reg, silent = TRUE)
tmb_time_reg <- proc.time()[["elapsed"]] - tmb_start

proreg_start <- proc.time()[["elapsed"]]
fit_proreg_reg <- PROreg::BBreg(y ~ x, m = m, data = dat_reg)
proreg_time_reg <- proc.time()[["elapsed"]] - proreg_start

coef_tmb_reg <- c(fit_tmb_reg$beta, phi = fit_tmb_reg$phi)
coef_proreg_reg <- c(fit_proreg_reg$coefficients, phi = fit_proreg_reg$phi)
names(coef_proreg_reg) <- names(coef_tmb_reg)

cat("=== BBreg: marginal beta-binomial regression ===\n")
compare("BBreg", fit_tmb_reg, fit_proreg_reg, tmb_time_reg,
        proreg_time_reg, coef_tmb_reg, coef_proreg_reg)

# Mixed beta-binomial regression with a random intercept.
n_group <- 30L
n_per_group <- 8L
group <- factor(rep(seq_len(n_group), each = n_per_group))
x <- rnorm(n_group * n_per_group)
u <- rnorm(n_group, sd = 0.8)
p <- plogis(0.4 - 0.6 * x + u[group])
y <- PROregTMB::rBB(length(x), m, p, phi = 0.15)
dat_mm <- data.frame(y = y, x = x, group = group)

# Warm up the mixed-model DLL as well.
invisible(PROregTMB::BBmm(
  y ~ x, random.formula = ~ group, m = m, data = dat_mm, silent = TRUE
))

tmb_start <- proc.time()[["elapsed"]]
fit_tmb_mm <- PROregTMB::BBmm(
  y ~ x, random.formula = ~ group, m = m, data = dat_mm, silent = TRUE
)
tmb_time_mm <- proc.time()[["elapsed"]] - tmb_start

proreg_start <- proc.time()[["elapsed"]]
fit_proreg_mm <- PROreg::BBmm(
  fixed.formula = y ~ x, random.formula = ~ group,
  m = m, data = dat_mm
)
proreg_time_mm <- proc.time()[["elapsed"]] - proreg_start

coef_tmb_mm <- c(fit_tmb_mm$beta, phi = fit_tmb_mm$phi,
                 sigma = fit_tmb_mm$sigma[1])
coef_proreg_mm <- c(fit_proreg_mm$fixed.coef, phi = fit_proreg_mm$phi.coef,
                    sigma = fit_proreg_mm$sigma.coef[1])
names(coef_proreg_mm) <- names(coef_tmb_mm)

cat("=== BBmm: random-intercept beta-binomial model ===\n")
compare("BBmm", fit_tmb_mm, fit_proreg_mm, tmb_time_mm,
        proreg_time_mm, coef_tmb_mm, coef_proreg_mm)

cat("Note: estimates may differ when the random-effect signal is weak because\n")
cat("PROreg and PROregTMB use different likelihood approximations.\n")