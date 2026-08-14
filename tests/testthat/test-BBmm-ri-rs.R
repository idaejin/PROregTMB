test_that("BBmm random intercept via random = ~ (1|id)", {
  skip_on_cran()
  set.seed(42)
  n_g <- 25
  n_per <- 6
  m <- 10
  x <- rnorm(n_g * n_per)
  id <- factor(rep(seq_len(n_g), each = n_per))
  u <- rnorm(n_g, 0, 0.7)
  eta <- 0.4 - 0.3 * x + u[id]
  y <- rBB(length(x), m, 1 / (1 + exp(-eta)), phi = 0.2)
  dat <- data.frame(y, x, id)

  fit <- BBmm(y ~ x, random = ~ (1 | id), m = m, data = dat)
  expect_equal(fit$conv, "yes")
  expect_true(all(c("Sigma", "Corr", "sigma") %in% names(fit)))
  expect_equal(fit$re_blocks[[1]]$q, 1L)
  expect_gt(as.numeric(fit$sigma)[1], 0.2)
})

test_that("BBmm RI+RS with correlation recovers signal", {
  skip_on_cran()
  set.seed(99)
  n_g <- 40
  n_per <- 8
  m <- 10
  id <- factor(rep(seq_len(n_g), each = n_per))
  time <- as.numeric(rep(seq_len(n_per), times = n_g))
  # Correlated RI+RS
  Sigma <- matrix(c(0.64, -0.2, -0.2, 0.25), 2, 2)
  L <- chol(Sigma)
  U <- matrix(rnorm(n_g * 2), n_g, 2) %*% L
  x <- rnorm(n_g * n_per)
  eta <- 0.3 - 0.2 * x + U[id, 1] + U[id, 2] * scale(time)[, 1]
  y <- rBB(length(x), m, 1 / (1 + exp(-eta)), phi = 0.12)
  dat <- data.frame(y, x, id, time)

  fit_us <- BBmm(
    y ~ x + time, random = ~ (1 + time | id),
    corr = "unstructured", m = m, data = dat
  )
  expect_equal(fit_us$conv, "yes")
  expect_equal(fit_us$re_blocks[[1]]$q, 2L)
  expect_equal(fit_us$re_blocks[[1]]$corr, "unstructured")
  expect_true(is.matrix(fit_us$Corr$id))
  expect_equal(nrow(fit_us$Corr$id), 2L)
  # correlation should be finite in (-1,1)
  rho <- fit_us$Corr$id[1, 2]
  expect_true(is.finite(rho) && abs(rho) < 1)

  # short alias "us" and informal "cor" both normalize to unstructured
  fit_us2 <- BBmm(
    y ~ x + time, random = ~ (1 + time | id),
    corr = "us", m = m, data = dat, maxiter = 5
  )
  expect_equal(fit_us2$re_blocks[[1]]$corr, "unstructured")
  fit_alias <- BBmm(
    y ~ x + time, random = ~ (1 + time | id),
    corr = "cor", m = m, data = dat, maxiter = 5
  )
  expect_equal(fit_alias$re_blocks[[1]]$corr, "unstructured")

  fit_diag <- BBmm(
    y ~ x + time, random = list(id = ~ 1 + time),
    corr = "diag", m = m, data = dat
  )
  expect_equal(fit_diag$conv, "yes")
  expect_equal(fit_diag$re_blocks[[1]]$corr, "diag")
  expect_equal(fit_diag$Corr$id[1, 2], 0)
})

test_that("legacy random.formula still works", {
  skip_on_cran()
  set.seed(42)
  n_g <- 30
  n_per <- 6
  m <- 10
  x <- rnorm(n_g * n_per)
  z <- factor(rep(seq_len(n_g), each = n_per))
  u <- rnorm(n_g, 0, 0.8)
  eta <- 0.5 - 0.4 * x + u[z]
  y <- rBB(length(x), m, 1 / (1 + exp(-eta)), phi = 0.15)
  dat <- data.frame(y, x, z)
  fit <- BBmm(y ~ x, random.formula = ~ z, m = m, data = dat)
  expect_equal(fit$conv, "yes")
  expect_equal(unname(fit$beta), unname(as.numeric(fit$fixed.coef)))
})
