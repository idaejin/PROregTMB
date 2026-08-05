test_that("BBmm recovers parameters in a simple simulation", {
  skip_on_cran()
  set.seed(42)
  n_g <- 30
  n_per <- 6
  k <- n_g * n_per
  m <- 10
  phi <- 0.15
  beta <- c(0.5, -0.4)
  sigma <- 0.8
  x <- rnorm(k)
  z <- factor(rep(seq_len(n_g), each = n_per))
  Z <- model.matrix(~ z - 1)
  u <- rnorm(n_g, 0, sigma)
  eta <- beta[1] + beta[2] * x + as.numeric(Z %*% u)
  p <- 1 / (1 + exp(-eta))
  y <- rBB(k, m, p, phi)
  dat <- data.frame(y, x, z)

  fit <- BBmm(y ~ x, random.formula = ~ z, m = m, data = dat)
  expect_equal(fit$conv, "yes")
  expect_equal(as.numeric(fit$fixed.coef[1]), beta[1], tolerance = 0.35)
  expect_equal(as.numeric(fit$fixed.coef[2]), beta[2], tolerance = 0.25)
  expect_equal(as.numeric(fit$phi.coef), phi, tolerance = 0.2)
  expect_equal(as.numeric(fit$sigma.coef), sigma, tolerance = 0.35)
})

test_that("BBmm close to PROreg::BBmm when RE is identifiable", {
  skip_if_not_installed("PROreg")
  skip_on_cran()
  set.seed(42)
  n_g <- 40
  n_per <- 8
  k <- n_g * n_per
  m <- 10
  phi <- 0.15
  beta <- c(0.5, -0.4)
  sigma <- 0.8
  x <- rnorm(k)
  z <- factor(rep(seq_len(n_g), each = n_per))
  Z <- model.matrix(~ z - 1)
  u <- rnorm(n_g, 0, sigma)
  eta <- beta[1] + beta[2] * x + as.numeric(Z %*% u)
  p <- 1 / (1 + exp(-eta))
  y <- rBB(k, m, p, phi)
  dat <- data.frame(y, x, z)

  fit_tmb <- BBmm(y ~ x, random.formula = ~ z, m = m, data = dat)
  fit_old <- PROreg::BBmm(fixed.formula = y ~ x, random.formula = ~ z,
                          m = m, data = dat)

  expect_equal(as.numeric(fit_tmb$fixed.coef),
               as.numeric(fit_old$fixed.coef),
               tolerance = 0.05)
  expect_equal(as.numeric(fit_tmb$phi.coef), fit_old$phi.coef, tolerance = 0.05)
  expect_equal(as.numeric(fit_tmb$sigma.coef),
               as.numeric(fit_old$sigma.coef),
               tolerance = 0.05)
  expect_gt(cor(fit_tmb$random.coef, fit_old$random.coef), 0.99)
})
