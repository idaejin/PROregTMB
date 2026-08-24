test_that("s(x, by=) expands to one smooth per level", {
  skip_on_cran()
  set.seed(31)
  n <- 120
  id <- factor(rep(seq_len(n / 2), each = 2))
  domain <- factor(rep(c("A", "B"), length.out = n))
  x <- runif(n)
  a <- rnorm(nlevels(id), 0, 0.4)[as.integer(id)]
  f <- ifelse(domain == "A", sin(2 * pi * x), cos(2 * pi * x))
  m <- 10L
  y <- rBB(n, m, 1 / (1 + exp(-(0.1 + f + a))), 0.2)
  dat <- data.frame(y, x, id, domain, m)

  fit <- BBmm(
    y ~ s(x, by = domain, k = 8),
    random = ~ (1 | id),
    dim = "domain",
    m = "m",
    data = dat,
    silent = TRUE,
    maxiter = 300
  )
  expect_identical(fit$conv, "yes")
  expect_equal(fit$smooth$n_smooth, 2L)
  expect_true(all(grepl("^s\\(x\\):", fit$smooth$labels)))
  # Domain-specific null-space slopes (not a shared "x")
  expect_true(all(c("A.x", "B.x") %in% names(fit$beta)))
  expect_false("x" %in% names(fit$beta))
  expect_true(all(abs(vapply(fit$fhat, sum, numeric(1))) < 0.25))
  pr <- predict_smooth(fit, which = "s(x):A", x = seq(0, 1, length.out = 20))
  expect_equal(nrow(pr), 20L)
})

test_that("s(x, by=) null space is domain-specific in design", {
  set.seed(2)
  n <- 24
  domain <- factor(rep(c("A", "B", "C"), length.out = n))
  x <- seq_len(n) / n
  dat <- data.frame(y = 1L, x, domain, m = 5L)
  sm <- PROregTMB:::.build_smooth_design(
    y ~ s(x, by = domain, k = 6),
    data = dat
  )
  expect_equal(colnames(sm$X_null), c("A.x", "B.x", "C.x"))
  # Column for level A is zero outside A rows
  expect_equal(sm$X_null[domain != "A", "A.x"], rep(0, sum(domain != "A")))
  expect_equal(sm$X_null[domain == "A", "A.x"], x[domain == "A"])
  expect_equal(
    vapply(sm$specs, `[[`, character(1), "null_name"),
    c("A.x", "B.x", "C.x")
  )
})
