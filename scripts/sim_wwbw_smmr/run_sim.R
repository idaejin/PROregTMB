#!/usr/bin/env Rscript
# SMMR WW-BW simulation runner (AIC + recovery metrics).
#
# Usage (from PROregTMB root):
#   Rscript scripts/sim_wwbw_smmr/run_sim.R --quick
#   Rscript scripts/sim_wwbw_smmr/run_sim.R --scenario=A --N=200 --nsim=50 --models=M1,M2,M3,M4

suppressPackageStartupMessages({
  args_all <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args_all, value = TRUE)
  this_dir <- if (length(f)) {
    dirname(normalizePath(sub("^--file=", "", f[1])))
  } else normalizePath("scripts/sim_wwbw_smmr")
  source(file.path(this_dir, "helpers.R"))
})

user_args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), user_args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}
has_flag <- function(name) any(user_args == paste0("--", name))

QUICK <- has_flag("quick")
SCENARIO <- get_arg("scenario", "A")
NSIM <- as.integer(get_arg("nsim", if (QUICK) "3" else "50"))
N <- as.integer(get_arg("N", if (QUICK) "80" else "500"))
Tvis <- as.integer(get_arg("T", "4"))
SEED0 <- as.integer(get_arg("seed", "1"))
MODELS <- strsplit(get_arg("models", if (QUICK) "M1,M2,M4" else "M1,M2,M3,M4"),
                   ",", fixed = TRUE)[[1]]
MAXITER <- as.integer(get_arg("maxiter", if (QUICK) "80" else "250"))

root <- load_proregtmb()
out_dir <- file.path(root, "scripts", "sim_wwbw_smmr", "out",
                     sprintf("sc%s_N%d_T%d", SCENARIO, N, Tvis))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message(sprintf(
  "SMMR sim scenario=%s  N=%d T=%d nsim=%d models=%s\n out=%s",
  SCENARIO, N, Tvis, NSIM, paste(MODELS, collapse = ","), out_dir
))

mod_rows <- list()
coef_rows <- list()

for (r in seq_len(NSIM)) {
  seed <- SEED0 + r - 1L
  message(sprintf("--- replicate %d/%d (seed=%d) ---", r, NSIM, seed))
  panel <- generate_panel(SCENARIO, N = N, T = Tvis, seed = seed)
  fits <- fit_competitors(panel, models = MODELS, maxiter = MAXITER)
  sumr <- summarise_replicate(panel, fits)

  mr <- sumr$models
  mr$rep <- r
  mr$seed <- seed
  mod_rows[[r]] <- mr

  if (!is.null(sumr$coefs)) {
    cr <- sumr$coefs
    cr$rep <- r
    cr$seed <- seed
    coef_rows[[r]] <- cr
  }

  write.csv(do.call(rbind, mod_rows),
            file.path(out_dir, "models_by_rep.csv"), row.names = FALSE)
  if (length(coef_rows)) {
    write.csv(do.call(rbind, coef_rows),
              file.path(out_dir, "coefs_by_rep.csv"), row.names = FALSE)
  }
}

models_df <- do.call(rbind, mod_rows)
coefs_df <- if (length(coef_rows)) do.call(rbind, coef_rows) else NULL

# --- Selection: AIC primary, nll secondary ---
okm <- models_df[models_df$ok %in% TRUE, , drop = FALSE]
best_aic <- by(okm, okm$rep, function(d) d$model[which.min(d$aic)])
best_nll <- by(okm, okm$rep, function(d) d$model[which.min(d$nll)])
sel <- data.frame(
  criterion = c(rep("AIC", length(unlist(best_aic))),
                rep("nll", length(unlist(best_nll)))),
  model = c(as.character(unlist(best_aic)), as.character(unlist(best_nll)))
)
sel_tab <- as.data.frame(table(sel$criterion, sel$model), stringsAsFactors = FALSE)
names(sel_tab) <- c("criterion", "model", "n_best")

# --- Aggregate coef bias/RMSE/coverage ---
coef_sum <- NULL
if (!is.null(coefs_df)) {
  okc <- coefs_df[is.finite(coefs_df$bias), , drop = FALSE]
  if (nrow(okc)) {
    coef_sum <- do.call(rbind, lapply(
      split(okc, list(okc$model, okc$domain, okc$effect), drop = TRUE),
      function(d) {
        data.frame(
          model = d$model[1], domain = d$domain[1], effect = d$effect[1],
          n = nrow(d),
          bias = mean(d$bias, na.rm = TRUE),
          rmse = sqrt(mean(d$bias^2, na.rm = TRUE)),
          cover95 = mean(d$cover95, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }
    ))
    rownames(coef_sum) <- NULL
  }
}

# --- Aggregate recovery of G, tau, phi ---
rec_sum <- do.call(rbind, lapply(split(okm, okm$model), function(d) {
  data.frame(
    model = d$model[1], n = nrow(d),
    mean_aic = mean(d$aic, na.rm = TRUE),
    mean_nll = mean(d$nll, na.rm = TRUE),
    mean_q = mean(d$q, na.rm = TRUE),
    corr_rmse_G = mean(d$corr_rmse_G, na.rm = TRUE),
    sd_rmse_G = mean(d$sd_rmse_G, na.rm = TRUE),
    bias_tau_fev = mean(d$bias_tau_fev, na.rm = TRUE),
    rmse_tau_fev = sqrt(mean(d$bias_tau_fev^2, na.rm = TRUE)),
    bias_tau_wt = mean(d$bias_tau_wt, na.rm = TRUE),
    rmse_tau_wt = sqrt(mean(d$bias_tau_wt^2, na.rm = TRUE)),
    bias_phi_I = mean(d$bias_phi_I, na.rm = TRUE),
    bias_phi_S = mean(d$bias_phi_S, na.rm = TRUE),
    bias_phi_A = mean(d$bias_phi_A, na.rm = TRUE),
    mean_roughness_s = mean(d$roughness_s, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
rownames(rec_sum) <- NULL

write.csv(sel_tab, file.path(out_dir, "selection_counts.csv"), row.names = FALSE)
write.csv(rec_sum, file.path(out_dir, "recovery_summary.csv"), row.names = FALSE)
if (!is.null(coef_sum))
  write.csv(coef_sum, file.path(out_dir, "coef_summary.csv"), row.names = FALSE)

saveRDS(list(
  scenario = SCENARIO, N = N, T = Tvis, nsim = NSIM, models = MODELS,
  models_df = models_df, coefs_df = coefs_df,
  selection = sel_tab, recovery = rec_sum, coef_summary = coef_sum
), file.path(out_dir, "sim_result.rds"))

message("\n=== Selection (AIC primary) ===")
print(sel_tab)
message("\n=== Recovery summary ===")
print(rec_sum)
if (!is.null(coef_sum)) {
  message("\n=== Coef summary (head) ===")
  print(utils::head(coef_sum, 20))
}
message("Wrote ", out_dir)
