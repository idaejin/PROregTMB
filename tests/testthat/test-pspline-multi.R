test_that("BBmm cbind + shared P-spline", {
  skip_on_cran()
  set.seed(21)
  n <- 60
  id <- factor(seq_len(n))
  x <- runif(n)
  a <- rnorm(n, 0, 0.5)
  f <- sin(2 * pi * x)
  y1 <- rBB(n, 20L, 1 / (1 + exp(-(0.2 + f + a))), 0.2)
  y2 <- rBB(n, 10L, 1 / (1 + exp(-(-0.3 + f + a))), 0.25)
  dat <- data.frame(y1, y2, x, id)

  fit <- BBmm(
    cbind(y1, y2) ~ s(x, ndx = 8, pord = 2),
    random = ~ (1 | id),
    m = c(20, 10),
    data = dat,
    silent = TRUE
  )
  expect_identical(fit$conv, "yes")
  expect_equal(fit$nDim, 2L)
  expect_equal(fit$smooth$n_smooth, 1L)
  expect_equal(nrow(fit$smooth$B), 2L * n)
  expect_true(abs(sum(fit$fhat[[1]])) < 0.1)
  pr <- predict_smooth(fit, which = 1L, x = seq(0, 1, length.out = 25))
  expect_equal(nrow(pr), 25L)
  expect_true(all(c("se", "lwr", "upr") %in% names(pr)))
})

test_that("BBmm long dim= + shared P-spline", {
  skip_on_cran()
  set.seed(22)
  n_g <- 40
  id <- factor(rep(seq_len(n_g), each = 2))
  domain <- factor(rep(c("SF36", "ACT"), times = n_g))
  x <- runif(n_g * 2)
  a <- rnorm(n_g, 0, 0.45)[as.integer(id)]
  f <- sin(2 * pi * x)
  m_vec <- ifelse(domain == "SF36", 20L, 10L)
  eta <- ifelse(domain == "SF36", 0.15, -0.25) + f + a
  y <- rBB(length(x), m_vec, 1 / (1 + exp(-eta)), phi = 0.2)
  long <- data.frame(y, x, id, domain, m = m_vec)

  fit <- BBmm(
    y ~ s(x, ndx = 8),
    random = ~ (1 | id),
    dim = "domain",
    m = "m",
    data = long,
    silent = TRUE
  )
  expect_identical(fit$conv, "yes")
  expect_equal(fit$nDim, 2L)
  expect_equal(fit$smooth$n_smooth, 1L)
  expect_equal(nrow(fit$smooth$B), nrow(long))
  expect_true(abs(sum(fit$fhat[[1]])) < 0.15)
})
