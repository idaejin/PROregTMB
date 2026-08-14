test_that("multi_bb_stack + BBmm nDim recovers shared RE", {
  skip_on_cran()
  set.seed(7)
  n <- 60
  L <- 2L
  x <- runif(n)
  u <- rnorm(n, 0, 0.7)
  y1 <- rBB(n, 20, 1 / (1 + exp(-(1 - 1.2 * x + u))), 0.4)
  y2 <- rBB(n, 10, 1 / (1 + exp(-(-1.5 + 2 * x + u))), 0.9)
  des <- multi_bb_stack(list(SF36 = y1, ACT = y2), x = x, m = c(20, 10))
  fit <- BBmm(
    X = des$X, y = des$y, Z = des$Z, nRandComp = des$nRandComp,
    m = des$m, nDim = des$nDim, silent = TRUE
  )
  expect_equal(fit$conv, "yes")
  expect_equal(fit$nDim, 2L)
  expect_length(fit$phi, 2L)
  expect_true(all(grepl("SF36|ACT", names(fit$beta))))
  expect_gt(as.numeric(fit$sigma)[1], 0.2)
})
