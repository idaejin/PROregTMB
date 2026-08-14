.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "PROregTMB ", utils::packageVersion("PROregTMB"),
    ": fit objects expose formulation fields beta/phi/p (BBreg) ",
    "and beta/u/sigma/phi/p (BBmm); PROreg-style aliases remain. ",
    "Avoid loading PROreg in the same session (S3 print/summary clash)."
  )
}
