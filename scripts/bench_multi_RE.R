# Benchmark: PROreg vs PROregTMB with several random-effect components
# Rscript scripts/bench_multi_RE.R

suppressPackageStartupMessages({
  library(PROreg)
})

args <- commandArgs(trailingOnly = FALSE)
f <- grep("^--file=", args, value = TRUE)
pkg_root <- if (length(f)) {
  normalizePath(file.path(dirname(sub("^--file=", "", f)), ".."))
} else {
  normalizePath(".")
}
devtools::load_all(pkg_root, quiet = TRUE)

sim_two_re <- function(n_site, n_doc_per_site, n_per_doc,
                       beta = c(0.3, -0.25),
                       phi = 0.12,
                       sigma_site = 0.55,
                       sigma_doc = 0.35,
                       m = 10,
                       seed = 1) {
  set.seed(seed)
  n_doc <- n_site * n_doc_per_site
  k <- n_doc * n_per_doc
  site <- factor(rep(seq_len(n_site), each = n_doc_per_site * n_per_doc))
  doc <- factor(rep(seq_len(n_doc), each = n_per_doc))
  x <- rnorm(k)
  u_site <- rnorm(n_site, 0, sigma_site)
  u_doc <- rnorm(n_doc, 0, sigma_doc)
  eta <- beta[1] + beta[2] * x +
    u_site[as.integer(site)] + u_doc[as.integer(doc)]
  p <- 1 / (1 + exp(-eta))
  y <- PROregTMB::rBB(k, m, p, phi)
  list(
    data = data.frame(y = y, x = x, site = site, doc = doc),
    m = m,
    truth = list(
      beta = beta, phi = phi,
      sigma_site = sigma_site, sigma_doc = sigma_doc,
      n = k, n_site = n_site, n_doc = n_doc,
      n_re = n_site + n_doc
    )
  )
}

fit_both <- function(dat, m, reps = 3L) {
  # PROreg
  t_old <- numeric(reps)
  fit_old <- NULL
  for (r in seq_len(reps)) {
    t0 <- proc.time()[["elapsed"]]
    fit_old <- PROreg::BBmm(
      fixed.formula = y ~ x,
      random.formula = ~ site + doc,
      m = m,
      data = dat
    )
    t_old[r] <- proc.time()[["elapsed"]] - t0
  }

  # TMB
  t_tmb <- numeric(reps)
  fit_tmb <- NULL
  for (r in seq_len(reps)) {
    t0 <- proc.time()[["elapsed"]]
    fit_tmb <- PROregTMB::BBmm(
      fixed.formula = y ~ x,
      random.formula = ~ site + doc,
      m = m,
      data = dat
    )
    t_tmb[r] <- proc.time()[["elapsed"]] - t0
  }

  list(
    old = fit_old, tmb = fit_tmb,
    time_old = t_old, time_tmb = t_tmb
  )
}

summarize_scenario <- function(label, sim, fits) {
  old <- fits$old
  tmb <- fits$tmb
  tr <- sim$truth
  data.frame(
    scenario = label,
    n = tr$n,
    n_site = tr$n_site,
    n_doc = tr$n_doc,
    n_re = tr$n_re,
    true_beta0 = tr$beta[1],
    true_beta1 = tr$beta[2],
    true_phi = tr$phi,
    true_sigma_site = tr$sigma_site,
    true_sigma_doc = tr$sigma_doc,
    proreg_beta0 = as.numeric(old$fixed.coef[1]),
    proreg_beta1 = as.numeric(old$fixed.coef[2]),
    proreg_phi = as.numeric(old$phi.coef),
    # PROreg leaves sigma.coef unnamed; order follows namesRand
    proreg_sigma_site = as.numeric(old$sigma.coef[match("site", old$namesRand)]),
    proreg_sigma_doc = as.numeric(old$sigma.coef[match("doc", old$namesRand)]),
    proreg_conv = old$conv,
    proreg_iter = old$iter,
    tmb_beta0 = as.numeric(tmb$fixed.coef[1]),
    tmb_beta1 = as.numeric(tmb$fixed.coef[2]),
    tmb_phi = as.numeric(tmb$phi.coef),
    tmb_sigma_site = as.numeric(tmb$sigma.coef[match("site", tmb$namesRand)]),
    tmb_sigma_doc = as.numeric(tmb$sigma.coef[match("doc", tmb$namesRand)]),
    tmb_conv = tmb$conv,
    tmb_iter = tmb$iter,
    tmb_nll = tmb$nll,
    time_proreg_mean = mean(fits$time_old),
    time_proreg_sd = stats::sd(fits$time_old),
    time_tmb_mean = mean(fits$time_tmb),
    time_tmb_sd = stats::sd(fits$time_tmb),
    speedup = mean(fits$time_old) / mean(fits$time_tmb),
    cor_u = {
      # Align RE predictions: both use same Z order (site then doc)
      suppressWarnings(cor(tmb$random.coef, old$random.coef))
    },
    stringsAsFactors = FALSE
  )
}

# Warm-up compile (exclude from timings)
cat("Warm-up compile...\n")
wu <- sim_two_re(4, 2, 4, seed = 0)
invisible(PROregTMB::BBmm(
  y ~ x, random.formula = ~ site + doc, m = wu$m, data = wu$data
))

scenarios <- list(
  list(label = "S1_small",  n_site = 10, n_doc_per_site = 3, n_per_doc = 5, seed = 11),
  list(label = "S2_medium", n_site = 20, n_doc_per_site = 4, n_per_doc = 6, seed = 22),
  list(label = "S3_large",  n_site = 30, n_doc_per_site = 5, n_per_doc = 8, seed = 33),
  list(label = "S4_xlarge", n_site = 40, n_doc_per_site = 6, n_per_doc = 8, seed = 44)
)

rows <- list()
for (sc in scenarios) {
  cat("\n===", sc$label, "===\n")
  sim <- sim_two_re(
    sc$n_site, sc$n_doc_per_site, sc$n_per_doc,
    seed = sc$seed
  )
  cat("n =", sim$truth$n, "| sites =", sim$truth$n_site,
      "| docs =", sim$truth$n_doc, "| #RE =", sim$truth$n_re, "\n")
  fits <- fit_both(sim$data, sim$m, reps = 3L)
  rows[[sc$label]] <- summarize_scenario(sc$label, sim, fits)
  print(rows[[sc$label]][, c(
    "n", "n_re",
    "proreg_sigma_site", "tmb_sigma_site",
    "proreg_sigma_doc", "tmb_sigma_doc",
    "time_proreg_mean", "time_tmb_mean", "speedup", "cor_u"
  )], digits = 4)
}

res <- do.call(rbind, rows)
rownames(res) <- NULL

out_dir <- file.path(pkg_root, "scripts", "bench_out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
csv_path <- file.path(out_dir, "bench_multi_RE.csv")
json_path <- file.path(out_dir, "bench_multi_RE.json")
utils::write.csv(res, csv_path, row.names = FALSE)

# JSON without needing jsonlite
to_json <- function(df) {
  lines <- vapply(seq_len(nrow(df)), function(i) {
    row <- df[i, ]
    parts <- vapply(names(row), function(nm) {
      v <- row[[nm]]
      if (is.character(v) || is.factor(v)) {
        sprintf('"%s":"%s"', nm, gsub('"', '\\"', as.character(v)))
      } else if (is.logical(v)) {
        sprintf('"%s":%s', nm, if (isTRUE(v)) "true" else "false")
      } else if (is.na(v)) {
        sprintf('"%s":null', nm)
      } else {
        sprintf('"%s":%s', nm, format(v, scientific = FALSE, digits = 8, trim = TRUE))
      }
    }, character(1))
    paste0("{", paste(parts, collapse = ","), "}")
  }, character(1))
  paste0("[\n", paste(lines, collapse = ",\n"), "\n]\n")
}
writeLines(to_json(res), json_path)

cat("\n=== Full results ===\n")
print(res, digits = 4)
cat("\nWrote:\n ", csv_path, "\n ", json_path, "\n")
