.tmb_env <- new.env(parent = emptyenv())

.tmb_src_dir <- function() {
  # Installed package: inst/TMB -> system.file("TMB", ...)
  cand <- system.file("TMB", package = "PROregTMB")
  if (nzchar(cand) && file.exists(file.path(cand, "bb_reg.cpp"))) {
    return(cand)
  }
  # devtools::load_all / source checkout
  pkg <- tryCatch(find.package("PROregTMB", quiet = TRUE), error = function(e) "")
  if (nzchar(pkg)) {
    # Sometimes load_all keeps path to source
    for (d in c(
      file.path(pkg, "inst", "TMB"),
      file.path(pkg, "TMB"),
      file.path(dirname(pkg), "inst", "TMB")
    )) {
      if (file.exists(file.path(d, "bb_reg.cpp"))) return(d)
    }
  }
  # Working directory = package root
  for (d in c(
    file.path(getwd(), "inst", "TMB"),
    file.path(getwd(), "TMB")
  )) {
    if (file.exists(file.path(d, "bb_reg.cpp"))) return(d)
  }
  stop("Cannot locate TMB templates for PROregTMB.", call. = FALSE)
}

# Eigen + -Wall often warns on unused-but-set locals under GCC/Clang.
.tmb_compile_flags <- function() {
  cxx <- tryCatch(
    system2(file.path(R.home("bin"), "R"),
            c("CMD", "config", "CXX17"),
            stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  )
  if (length(cxx) && any(grepl("clang|g\\+\\+|gcc", cxx, ignore.case = TRUE))) {
    return("-Wno-unused-but-set-variable")
  }
  ""
}

#' Ensure a TMB DLL is compiled and loaded
#' @param name Template basename without extension (e.g. "bb_reg")
#' @param force Recompile even if already loaded
#' @return Invisibly, the DLL name
#' @keywords internal
ensure_tmb_dll <- function(name, force = FALSE) {
  loaded <- vapply(getLoadedDLLs(), function(x) x[["name"]], character(1))
  if (!force && (name %in% loaded || isTRUE(.tmb_env[[name]]))) {
    return(invisible(name))
  }

  src_dir <- .tmb_src_dir()
  cpp <- file.path(src_dir, paste0(name, ".cpp"))
  if (!file.exists(cpp)) {
    stop("TMB template not found: ", cpp, call. = FALSE)
  }

  work <- file.path(tempdir(), "PROregTMB_tmb")
  dir.create(work, showWarnings = FALSE, recursive = TRUE)
  file.copy(cpp, file.path(work, basename(cpp)), overwrite = TRUE)

  wd <- getwd()
  on.exit(setwd(wd), add = TRUE)
  setwd(work)

  TMB::compile(basename(cpp), flags = .tmb_compile_flags())
  dyn.load(TMB::dynlib(name))
  .tmb_env[[name]] <- TRUE
  invisible(name)
}
