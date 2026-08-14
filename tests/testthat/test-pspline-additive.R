test_that("additive P-splines in BBreg with sum-to-zero", {
  skip_on_cran()
  set.seed(11)
  n <- 250
  m <- 10L
  x1 <- runif(n)
  x2 <- runif(n)
  f1 <- sin(2 * pi * x1)
  f2 <- exp(2 * x2) - (exp(2) - 1) / 2
  eta <- -0.5 + f1 + f2
  p <- 1 / (1 + exp(-eta))
  y <- rBB(n, m, p, phi = 0.2)
  dat <- data.frame(y, x1, x2)

  fit <- BBreg(y ~ s(x1, ndx = 8, pord = 2) + s(x2, ndx = 8, pord = 2),
               m = m, data = dat, silent = TRUE)
  expect_identical(fit$conv, "yes")
  expect_equal(fit$smooth$n_smooth, 2L)
  expect_true(all(is.finite(fit$lambda)))
  expect_true(all(abs(vapply(fit$fhat, sum, numeric(1))) < 1e-2))
  pr <- predict_smooth(fit, which = 1L, x = seq(0, 1, length.out = 30))
  expect_true(all(c("x", "fit") %in% names(pr)))
  expect_equal(nrow(pr), 30L)
})

test_that("BBmm with RI + one P-spline", {
  skip_on_cran()
  set.seed(12)
  n_g <- 25
  n_per <- 5
  m <- 8L
  id <- factor(rep(seq_len(n_g), each = n_per))
  x <- runif(n_g * n_per)
  u <- rnorm(n_g, 0, 0.5)
  f <- sin(2 * pi * x)
  eta <- 0.2 + f + u[id]
  y <- rBB(length(x), m, 1 / (1 + exp(-eta)), phi = 0.15)
  dat <- data.frame(y, x, id)

  fit <- BBmm(y ~ s(x, ndx = 8), random = ~ (1 | id), m = m, data = dat,
              silent = TRUE)
  expect_identical(fit$conv, "yes")
  expect_true(!is.null(fit$lambda))
  expect_true(abs(sum(fit$fhat[[1]])) < 5e-2)
})

test_that("parametric BBreg still works after smooth templates", {
  skip_on_cran()
  set.seed(13)
  n <- 200
  m <- 10L
  x <- rnorm(n)
  p <- 1 / (1 + exp(-(-1 + 0.8 * x)))
  y <- rBB(n, m, p, phi = 0.3)
  fit <- BBreg(y ~ x, m = m, silent = TRUE)
  expect_identical(fit$conv, "yes")
  expect_null(fit$smooth)
})
