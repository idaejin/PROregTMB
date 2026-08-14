test_that("Rcpp rBB / dBB / ldBB", {
  set.seed(1)
  y <- rBB(5000L, m = 10L, p = 0.4, phi = 0.2)
  expect_equal(length(y), 5000L)
  expect_true(all(y >= 0L & y <= 10L))
  pmf <- dBB(10L, 0.4, 0.2)
  expect_equal(length(pmf), 11L)
  expect_equal(sum(pmf), 1, tolerance = 1e-8)
  ll <- ldBB(y[1:5], m = 10, p = 0.4, phi = 0.2)
  expect_equal(length(ll), 5L)
  expect_true(all(is.finite(ll)))
})

test_that("sparse Z from RE bars", {
  n_g <- 20L
  n_per <- 4L
  id <- factor(rep(seq_len(n_g), each = n_per))
  time <- ave(seq_along(id), id, FUN = seq_along)
  dat <- data.frame(id = id, time = time)
  re <- PROregTMB:::.build_re_structure(
    random = ~ (1 + time | id),
    data = dat,
    corr = "unstructured",
    nObs = nrow(dat)
  )
  expect_true(inherits(re$Z, "sparseMatrix"))
  expect_equal(nrow(re$Z), nrow(dat))
  expect_equal(ncol(re$Z), n_g * 2L)
  # densify works for TMB
  Zd <- PROregTMB:::.as_dense_Z(re$Z)
  expect_true(is.matrix(Zd))
  expect_equal(dim(Zd), dim(re$Z))
})

test_that("predict_smooth SE via Rcpp + sparse jointPrecision path", {
  skip_on_cran()
  set.seed(21)
  n <- 180
  m <- 8L
  x <- runif(n)
  y <- as.integer(rBB(n, m, 1 / (1 + exp(-(sin(2 * pi * x)))), 0.2))
  dat <- data.frame(y, x)
  fit <- BBreg(y ~ s(x, ndx = 8), m = m, data = dat, silent = TRUE)
  expect_identical(fit$conv, "yes")
  expect_true(!is.null(fit$alpha.vcov))
  pr <- predict_smooth(fit, which = 1L, x = seq(0, 1, length.out = 40))
  expect_true(all(c("se", "lwr", "upr") %in% names(pr)))
  expect_true(all(is.finite(pr$se)))
  expect_true(all(pr$lwr <= pr$fit & pr$fit <= pr$upr))
  # plot helper returns invisibly with bands
  pdf(NULL)
  out <- plot_smooth(fit, which = 1L, x = seq(0, 1, length.out = 30))
  dev.off()
  expect_true(all(c("lwr", "upr") %in% names(out)))
})
