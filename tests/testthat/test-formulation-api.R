test_that("BBreg exposes formulation-aligned fields", {
  skip_on_cran()
  set.seed(18)
  k <- 200
  m <- 10
  x <- rnorm(k, 5, 3)
  p <- 1 / (1 + exp(-(-10 + 2 * x)))
  y <- rBB(k, m, p, phi = 1.2)
  fit <- BBreg(y ~ x, m)
  expect_true(all(c("beta", "phi", "log_phi", "p", "model") %in% names(fit)))
  expect_equal(unname(fit$beta), as.numeric(fit$coefficients), tolerance = 1e-12)
  expect_equal(fit$p, fit$fitted.values, tolerance = 1e-12)
  expect_equal(fit$log_phi, as.numeric(fit$psi), tolerance = 1e-12)
  expect_true(any(grepl("BB\\(m_i", fit$model)))
  expect_equal(coef(fit), fit$beta)
  expect_equal(fitted(fit), fit$p)
})

test_that("BBmm exposes formulation-aligned fields", {
  skip_on_cran()
  set.seed(42)
  n_g <- 20
  n_per <- 5
  m <- 10
  x <- rnorm(n_g * n_per)
  z <- factor(rep(seq_len(n_g), each = n_per))
  u <- rnorm(n_g, 0, 0.8)
  eta <- 0.5 - 0.4 * x + u[z]
  y <- rBB(length(x), m, 1 / (1 + exp(-eta)), phi = 0.15)
  dat <- data.frame(y, x, z)
  fit <- BBmm(y ~ x, random.formula = ~ z, m = m, data = dat)
  expect_true(all(c("beta", "u", "sigma", "phi", "p", "model") %in% names(fit)))
  expect_equal(unname(fit$beta), unname(as.numeric(fit$fixed.coef)), tolerance = 1e-12)
  expect_equal(unname(fit$u), unname(as.numeric(fit$random.coef)), tolerance = 1e-12)
  expect_equal(unname(fit$sigma), unname(as.numeric(fit$sigma.coef)), tolerance = 1e-12)
  expect_equal(fit$phi, fit$phi.coef)
  expect_true(!is.null(fit$fixed.formula))
  expect_true(!is.null(fit$random.formula))
  expect_equal(coef(fit), fit$beta)
})
