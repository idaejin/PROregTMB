test_that("M3: domain-specific correlated RE via (0 + dim | id)", {
  skip_on_cran()
  set.seed(42)
  n_g <- 50
  L <- 3L
  id <- factor(rep(seq_len(n_g), each = L))
  domain <- factor(rep(c("Impacts", "Symptoms", "Activity"), times = n_g),
                   levels = c("Impacts", "Symptoms", "Activity"))
  x <- rnorm(n_g * L)
  Sigma_true <- matrix(c(
    0.55, 0.22, 0.18,
    0.22, 0.45, 0.12,
    0.18, 0.12, 0.35
  ), 3, 3)
  U <- matrix(rnorm(n_g * 3), n_g, 3) %*% chol(Sigma_true)
  u <- U[cbind(as.integer(id), as.integer(domain))]
  g0 <- c(0.3, -0.2, 0.0)[as.integer(domain)]
  eta <- g0 + 0.25 * x + u
  y <- rBB(length(x), 12L, 1 / (1 + exp(-eta)), 0.18)
  dat <- data.frame(y, x, id, domain, m = 12L)

  fit <- BBmm(
    y ~ x,
    random = ~ (0 + domain | id),
    dim = "domain",
    corr = "unstructured",
    m = "m",
    data = dat,
    silent = TRUE,
    maxiter = 250
  )
  expect_equal(fit$nDim, 3L)
  expect_equal(fit$re_blocks[[1]]$q, 3L)
  expect_equal(fit$re_blocks[[1]]$corr, "unstructured")
  expect_equal(fit$nRand, n_g * 3L)
  expect_true(is.matrix(fit$Sigma$id))
  expect_equal(dim(fit$Sigma$id), c(3L, 3L))
  expect_true(is.matrix(fit$Corr$id))
  # Correlations in (-1,1), SDs positive
  expect_true(all(diag(fit$Corr$id) > 0.99))
  expect_true(all(abs(fit$Corr$id[upper.tri(fit$Corr$id)]) < 1))
  expect_true(all(fit$sigma > 0))
})

test_that("M2: domain RE independent via corr=diag", {
  skip_on_cran()
  set.seed(43)
  n_g <- 40
  id <- factor(rep(seq_len(n_g), each = 2))
  domain <- factor(rep(c("A", "B"), times = n_g))
  x <- rnorm(n_g * 2)
  uA <- rnorm(n_g, 0, 0.5)
  uB <- rnorm(n_g, 0, 0.4)
  u <- ifelse(domain == "A", uA[as.integer(id)], uB[as.integer(id)])
  y <- rBB(length(x), 10L, 1 / (1 + exp(-(0.1 * x + u))), 0.2)
  dat <- data.frame(y, x, id, domain, m = 10L)

  fit <- BBmm(
    y ~ x,
    random = ~ (0 + domain | id),
    dim = "domain",
    corr = "diag",
    m = "m",
    data = dat,
    silent = TRUE
  )
  expect_equal(fit$re_blocks[[1]]$q, 2L)
  expect_equal(fit$re_blocks[[1]]$corr, "diag")
  expect_true(max(abs(fit$Corr$id - diag(2))) < 1e-8)
})
