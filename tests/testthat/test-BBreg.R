test_that("dBB matches PROreg and sums to 1", {
  skip_if_not_installed("PROreg")
  m <- 10
  p <- 0.4
  phi <- 1.8
  d1 <- dBB(m, p, phi)
  d2 <- PROreg::dBB(m, p, phi)
  expect_equal(sum(d1), 1, tolerance = 1e-10)
  expect_equal(d1, d2, tolerance = 1e-8)
})

test_that("BBreg recovers known parameters", {
  skip_on_cran()
  set.seed(18)
  k <- 800
  m <- 10
  x <- rnorm(k, 5, 3)
  beta <- c(-10, 2)
  p <- 1 / (1 + exp(-(beta[1] + beta[2] * x)))
  phi <- 1.2
  y <- rBB(k, m, p, phi)
  fit <- BBreg(y ~ x, m)
  expect_equal(fit$conv, "yes")
  est <- as.numeric(fit$coefficients)
  expect_equal(est[1], beta[1], tolerance = 0.4)
  expect_equal(est[2], beta[2], tolerance = 0.15)
  expect_equal(fit$phi, phi, tolerance = 0.35)
})

test_that("BBreg close to PROreg::BBreg", {
  skip_if_not_installed("PROreg")
  skip_on_cran()
  set.seed(18)
  k <- 500
  m <- 10
  x <- rnorm(k, 5, 3)
  beta <- c(-10, 2)
  p <- 1 / (1 + exp(-(beta[1] + beta[2] * x)))
  y <- rBB(k, m, p, phi = 1.2)
  fit_tmb <- BBreg(y ~ x, m)
  fit_old <- PROreg::BBreg(y ~ x, m)
  expect_equal(as.numeric(fit_tmb$coefficients),
               as.numeric(fit_old$coefficients),
               tolerance = 0.05)
  expect_equal(fit_tmb$phi, fit_old$phi, tolerance = 0.1)
})
