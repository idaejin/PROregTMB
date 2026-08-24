#!/usr/bin/env Rscript
# Recode COPD SF-36 (0-100) to binomial form for beta-binomial models.
# Method: Arostegui et al. (2013), same cuts as PROreg::SF36rec
# (PROreg validates name "MT" but implements "MH"; we accept "MH".)
#
# Usage (from PROregTMB root or anywhere):
#   Rscript scripts/prep_copd_sf36_bb.R
#
# Outputs:
#   scripts/out_copd_sf36_bb/sf36_bb_wide.rds
#   scripts/out_copd_sf36_bb/sf36_bb_long.rds
#   scripts/out_copd_sf36_bb/sf36_bb_summary.csv

suppressPackageStartupMessages({
  if (!requireNamespace("car", quietly = TRUE)) {
    stop("Package 'car' required for recode()", call. = FALSE)
  }
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

.root <- local({
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in seq_len(8L)) {
    if (file.exists(file.path(d, "DESCRIPTION")) &&
        file.exists(file.path(d, "R", "BBmm.R"))) return(d)
    nd <- dirname(d)
    if (identical(nd, d)) break
    d <- nd
  }
  normalizePath(
    "/Users/daejin/Dropbox/IE/research/01_PROJECTS/BB/PROregTMB",
    mustWork = FALSE
  )
})

data_path <- normalizePath(
  file.path(.root, "..", "ref", "BB-GAM", "longitudinal", "data", "COPD_All.txt"),
  mustWork = TRUE
)
out_dir <- file.path(.root, "scripts", "out_copd_sf36_bb")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Domain order + binomial maxima after SF36rec
sf36_m <- c(
  PF = 20L, RP = 4L, BP = 9L, GH = 20L,
  VT = 20L, SF = 8L, RE = 3L, MH = 13L
)

#' Arostegui / PROreg SF-36 recode (fixed MH name check)
SF36rec_bb <- function(x, name) {
  name <- as.character(name)[1]
  if (!is.numeric(x)) stop("x must be numeric", call. = FALSE)
  if (any(x < 0 | x > 100, na.rm = TRUE)) {
    stop("x must be bounded between 0 and 100", call. = FALSE)
  }
  if (!name %in% names(sf36_m)) {
    stop("name must be one of ", paste(names(sf36_m), collapse = ", "),
         call. = FALSE)
  }
  out <- switch(
    name,
    PF = car::recode(x, "0:2.5=0; 2.5:7.5=1; 7.5:12.5=2; 12.5:17.5=3; 17.5:22.5=4;
      22.5:27.5=5; 27.5:32.5=6; 32.5:37.5=7; 37.5:42.5=8; 42.5:47.5=9;
      47.5:52.5=10; 52.5:57.5=11; 57.5:62.5=12; 62.5:67.5=13; 67.5:72.5=14;
      72.5:77.5=15; 77.5:82.5=16; 82.5:87.5=17; 87.5:92.5=18; 92.5:97.5=19;
      97.5:100=20"),
    RP = car::recode(x, "0:12.5=0 ; 12.5:37.5=1 ; 37.5:62.5=2 ; 62.5:87.5=3 ;
      87.5:100=4"),
    BP = car::recode(x, "0:5=0; 5:16=1; 16:26=2; 26:36=3; 36:47=4; 47:57=5; 57:67=6;
      67:77=7; 77:92=8; 92:100=9"),
    GH = car::recode(x, "0:2.5=0; 2.5:7.5=1; 7.5:13.5=2; 13.5:18.5=3; 18.5:23.5=4;
      23.5:28.5=5; 28.5:33.5=6; 33.5:38.5=7; 38.5:43.5=8; 43.5:48.5=9;
      48.5:53.5=10; 53.5:58.5=11; 58.5:63.5=12; 63.5:68.5=13; 68.5:73.5=14;
      73.5:78.5=15; 78.5:83.5=16; 83.5:88.5=17; 88.5:93.5=18; 93.5:98.5=19;
      98.5:100=20"),
    VT = car::recode(x, "0:2.5=0; 2.5:7.5=1; 7.5:12.5=2; 12.5:17.5=3; 17.5:22.5=4;
      22.5:27.5=5; 27.5:32.5=6; 32.5:37.5=7; 37.5:42.5=8; 42.5:47.5=9;
      47.5:52.5=10; 52.5:57.5=11; 57.5:62.5=12; 62.5:67.5=13; 67.5:72.5=14;
      72.5:77.5=15; 77.5:82.5=16; 82.5:87.5=17; 87.5:92.5=18; 92.5:97.5=19;
      97.5:100=20"),
    SF = car::recode(x, "0:6.25=0; 6.25:18.75=1; 18.75:31.25=2; 31.25:43.75=3;
      43.75:56.25=4; 56.25:68.75=5; 68.75:81.25=6; 81.25:93.75=7;
      93.75:100=8"),
    RE = car::recode(x, "0:16.67=0 ; 16.67:50=1 ; 50:83.33=2 ; 83.33:100=3"),
    MH = car::recode(x, "0:2=0;2:10=1;10:18=2;18:26=3;26:34=4;34:42=5;42:50=6;
      50:58=7;58:66=8;66:74=9;74:82=10;82:90=11;90:98=12;98:100=13")
  )
  as.integer(out)
}

raw <- read.table(data_path, header = TRUE, stringsAsFactors = FALSE)
message("COPD rows: ", nrow(raw), "  subjects: ", length(unique(raw$id)))

# Recode each domain
for (nm in names(sf36_m)) {
  raw[[paste0(nm, "_bb")]] <- SF36rec_bb(as.numeric(raw[[nm]]), nm)
}

# Quick integrity: max <= m, min >= 0
sum_rows <- lapply(names(sf36_m), function(nm) {
  y <- raw[[paste0(nm, "_bb")]]
  data.frame(
    domain = nm,
    m = sf36_m[[nm]],
    n = sum(is.finite(y)),
    n_na = sum(!is.finite(y)),
    ymin = min(y, na.rm = TRUE),
    ymax = max(y, na.rm = TRUE),
    n_unique = length(unique(y[is.finite(y)])),
    mean = mean(y, na.rm = TRUE),
    ok = isTRUE(min(y, na.rm = TRUE) >= 0 &&
                  max(y, na.rm = TRUE) <= sf36_m[[nm]]),
    stringsAsFactors = FALSE
  )
})
sum_df <- do.call(rbind, sum_rows)
print(sum_df)
if (!all(sum_df$ok)) stop("Recoded values outside 0..m for some domain", call. = FALSE)

# Covariates for WW-BW (same as SGRQ pipeline)
zscore <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s < 1e-8) return(x - mean(x, na.rm = TRUE))
  (x - mean(x, na.rm = TRUE)) / s
}
raw$id <- factor(raw$id)
raw$time <- as.integer(raw$time)
raw$time_c <- raw$time - 1L
raw$sex <- factor(raw$sex, levels = c("1", "2"), labels = c("Male", "Female"))
raw$age_z <- zscore(as.numeric(raw$age))
inna_map <- c("1" = "B", "2" = "C", "3" = "A", "4" = "D")
raw$hclus4 <- factor(inna_map[as.character(raw$cluster_INMA)],
                     levels = c("B", "A", "C", "D"))
raw$fev <- as.numeric(raw$fev)
raw$wt <- as.numeric(raw$wt)
raw$fev_bar <- ave(raw$fev, raw$id, FUN = mean)
raw$wt_bar <- ave(raw$wt, raw$id, FUN = mean)
raw$fev_w <- raw$fev - raw$fev_bar
raw$wt_w <- raw$wt - raw$wt_bar
raw$fev_bar_z <- zscore(raw$fev_bar)
raw$wt_bar_z <- zscore(raw$wt_bar)
raw$fev_z <- zscore(raw$fev)
raw$wt_z <- zscore(raw$wt)

ok <- complete.cases(raw[, c("fev", "wt", "age", "sex", "cluster_INMA",
                             paste0(names(sf36_m), "_bb"))])
raw <- raw[ok, , drop = FALSE]
raw$id <- droplevels(raw$id)
message("After complete cases: ", nrow(raw), " visits, ",
        nlevels(raw$id), " subjects")

# Long stacked for BBmm (domain-specific m)
long <- do.call(rbind, lapply(names(sf36_m), function(nm) {
  data.frame(
    id = raw$id,
    time = raw$time,
    time_c = raw$time_c,
    dim = factor(nm, levels = names(sf36_m)),
    y = raw[[paste0(nm, "_bb")]],
    m = as.integer(sf36_m[[nm]]),
    sex = raw$sex,
    hclus4 = raw$hclus4,
    age_z = raw$age_z,
    fev_z = raw$fev_z,
    wt_z = raw$wt_z,
    fev_bar_z = raw$fev_bar_z,
    wt_bar_z = raw$wt_bar_z,
    fev_w = raw$fev_w,
    wt_w = raw$wt_w,
    # keep 0-100 originals for reference
    score_100 = raw[[nm]],
    stringsAsFactors = FALSE
  )
}))
rownames(long) <- NULL

saveRDS(raw, file.path(out_dir, "sf36_bb_wide.rds"))
saveRDS(long, file.path(out_dir, "sf36_bb_long.rds"))
utils::write.csv(sum_df, file.path(out_dir, "sf36_bb_summary.csv"),
                 row.names = FALSE)
utils::write.csv(
  data.frame(domain = names(sf36_m), m = as.integer(sf36_m)),
  file.path(out_dir, "sf36_bb_m.csv"), row.names = FALSE
)

message("Wrote ", out_dir)
message("Long: ", nrow(long), " rows (= visits x 8 domains)")
message("Example BBmm call:")
message('  BBmm(y ~ hclus4 + time_c + sex + age_z + fev_bar_z + wt_bar_z + fev_w + wt_w,')
message('       random = ~ (0 + dim | id), dim = "dim", corr = "unstructured",')
message('       m = "m", data = long)')
