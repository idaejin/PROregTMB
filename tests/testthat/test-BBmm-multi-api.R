test_that("BBmm cbind multivariate shared RI", {
  skip_on_cran()
  set.seed(11)
  n <- 50
  x <- runif(n)
  id <- factor(seq_len(n))
  u <- rnorm(n, 0, 0.6)
  y1 <- rBB(n, 20, 1 / (1 + exp(-(1 - 1.2 * x + u))), 0.4)
  y2 <- rBB(n, 10, 1 / (1 + exp(-(-1.5 + 2 * x + u))), 0.9)
  dat <- data.frame(y1, y2, x, id)
  fit <- BBmm(
    cbind(y1, y2) ~ x, random = ~ (1 | id),
    m = c(20, 10), data = dat, silent = TRUE
  )
  expect_equal(fit$conv, "yes")
  expect_equal(fit$nDim, 2L)
  expect_equal(fit$structure, "shared")
  expect_length(fit$phi, 2L)
  expect_true(any(grepl("y1\\.", names(fit$beta))))
  expect_true(any(grepl("y2\\.", names(fit$beta))))
})

test_that("BBmm long dim= shared RI+RS", {
  skip_on_cran()
  set.seed(12)
  n <- 40
  n_per <- 5
  id <- factor(rep(seq_len(n), each = n_per))
  time <- as.numeric(rep(seq_len(n_per), times = n))
  x <- rnorm(n * n_per)
  # shared RI+RS
  Sigma <- matrix(c(0.5, -0.1, -0.1, 0.16), 2, 2)
  U <- matrix(rnorm(n * 2), n, 2) %*% chol(Sigma)
  u0 <- U[id, 1]
  u1 <- U[id, 2]
  mk_long <- function(domain, m, g0, g1) {
    eta <- g0 + g1 * x + u0 + u1 * scale(time)[, 1]
    data.frame(
      y = rBB(length(x), m, 1 / (1 + exp(-eta)), 0.2),
      x = x, time = time, id = id,
      domain = domain, m = m
    )
  }
  long <- rbind(
    mk_long("SF36", 20, 0.5, -0.3),
    mk_long("ACT", 10, -0.5, 0.4)
  )
  fit <- BBmm(
    y ~ x + time,
    random = ~ (1 + time | id),
    dim = "domain",
    corr = "unstructured",
    m = "m",
    data = long,
    silent = TRUE
  )
  expect_equal(fit$conv, "yes")
  expect_equal(fit$nDim, 2L)
  expect_equal(fit$structure, "shared")
  expect_equal(fit$re_blocks[[1]]$q, 2L)
  expect_length(fit$phi, 2L)
  expect_true(is.matrix(fit$Corr$id))
})
