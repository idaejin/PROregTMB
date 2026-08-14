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
#'
#' On first use (or when the `.cpp` template changes), compiles via
#' [TMB::compile()] in a quiet subprocess so `clang++` / `make` lines do
#' not flood the console. A one-line status message is shown unless
#' `getOption("PROregTMB.verbose_compile")` is `FALSE` or `quiet = TRUE`.
#'
#' @param name Template basename without extension (e.g. "bb_reg")
#' @param force Recompile even if already loaded
#' @param quiet If `TRUE`, no status message (compiler output still suppressed).
#' @return Invisibly, the DLL name
#' @keywords internal
ensure_tmb_dll <- function(name, force = FALSE, quiet = FALSE) {
  src_dir <- .tmb_src_dir()
  cpp <- file.path(src_dir, paste0(name, ".cpp"))
  if (!file.exists(cpp)) {
    stop("TMB template not found: ", cpp, call. = FALSE)
  }
  mtime <- file.info(cpp)$mtime
  cached_mtime <- .tmb_env[[paste0(name, "_mtime")]]
  loaded <- vapply(getLoadedDLLs(), function(x) x[["name"]], character(1))
  stale <- is.null(cached_mtime) || !identical(cached_mtime, mtime)
  if (!force && !stale && (name %in% loaded || isTRUE(.tmb_env[[name]]))) {
    return(invisible(name))
  }
  # Drop previous DLL if we must rebuild
  if (name %in% loaded) {
    dll <- getLoadedDLLs()[[name]]
    try(dyn.unload(dll[["path"]]), silent = TRUE)
  }
  .tmb_env[[name]] <- NULL

  work <- file.path(tempdir(), "PROregTMB_tmb")
  dir.create(work, showWarnings = FALSE, recursive = TRUE)
  file.copy(cpp, file.path(work, basename(cpp)), overwrite = TRUE)

  show_msg <- !isTRUE(quiet) &&
    !isFALSE(getOption("PROregTMB.verbose_compile", TRUE))
  if (show_msg) {
    message(
      "PROregTMB: compiling TMB model '", name,
      "' (first use or template changed; may take a minute)..."
    )
  }

  # Subprocess keeps SHLIB / clang noise off the console
  flags <- .tmb_compile_flags()
  script <- tempfile("PROregTMB_compile_", fileext = ".R")
  writeLines(
    c(
      "Sys.setenv(MAKEFLAGS = paste(Sys.getenv('MAKEFLAGS'), '-s'))",
      paste0("setwd(", deparse(work), ")"),
      paste0(
        "TMB::compile(", deparse(basename(cpp)),
        ", flags = ", deparse(flags), ")"
      )
    ),
    script
  )
  on.exit(unlink(script), add = TRUE)
  out <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", script),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(out, "status")
  if (!is.null(status) && !identical(as.integer(status), 0L)) {
    stop(
      "TMB compile failed for '", name, "':\n",
      paste(out, collapse = "\n"),
      call. = FALSE
    )
  }
  so <- file.path(work, paste0(name, .Platform$dynlib.ext))
  if (!file.exists(so)) {
    stop(
      "TMB compile produced no DLL for '", name, "':\n",
      paste(out, collapse = "\n"),
      call. = FALSE
    )
  }

  if (show_msg) message("PROregTMB: '", name, "' ready.")

  dyn.load(so)
  .tmb_env[[name]] <- TRUE
  .tmb_env[[paste0(name, "_mtime")]] <- mtime
  invisible(name)
}
