test_that("calc_asymptotic_variance matches manual logistic term when U absent", {
  beta_0 <- 0
  beta_a <- 0.5
  theta2 <- 0.3
  X2 <- matrix(c(1, 2), ncol = 1)
  y2 <- c(1, 0)

  result <- calc_asymptotic_variance(
    beta_0 = beta_0,
    beta_a = beta_a,
    theta2 = theta2,
    y2 = y2,
    X2 = X2,
    include_estimation_uncertainty = FALSE,
    model = "logistic"
  )

  expect_equal(result, 4.079173411852398, tolerance = 1e-10)
})

test_that("calc_asymptotic_variance handles logistic model with unpenalized block", {
  beta_0 <- 0
  beta_a <- 0.5
  theta2 <- 0.3
  X2 <- matrix(c(1, 2), ncol = 1)
  y2 <- c(1, 0)
  U2 <- matrix(1, nrow = 2, ncol = 1)
  delta <- 0.2

  result <- calc_asymptotic_variance(
    beta_0 = beta_0,
    beta_a = beta_a,
    theta2 = theta2,
    y2 = y2,
    X2 = X2,
    U2 = U2,
    delta = delta,
    include_estimation_uncertainty = FALSE,
    model = "logistic"
  )

  expect_equal(result, 4.922160276554947, tolerance = 1e-10)
})

test_that("calc_asymptotic_variance respects supplied means and variances", {
  beta_0 <- 0
  beta_a <- 0.5
  theta2 <- 0.3
  X2 <- matrix(c(1, 2), ncol = 1)
  y2 <- c(1, 0)

  beta_theta <- beta_0 + theta2 * (beta_a - beta_0)
  eta <- as.vector(X2 %*% beta_theta)
  mu <- 1 / (1 + exp(-eta))
  variance <- mu * (1 - mu)

  result <- calc_asymptotic_variance(
    beta_0 = beta_0,
    beta_a = beta_a,
    theta2 = theta2,
    y2 = y2,
    X2 = X2,
    mu2 = mu,
    variance2 = variance,
    include_estimation_uncertainty = FALSE,
    model = "logistic"
  )

  expect_equal(result, 4.079173411852398, tolerance = 1e-10)
})

test_that("calc_asymptotic_variance supports poisson responses with U", {
  beta_0 <- 0.1
  beta_a <- 0.4
  theta2 <- 0.25
  X2 <- matrix(c(1, 0.5), ncol = 1)
  y2 <- c(2, 1)
  U2 <- matrix(1, nrow = 2, ncol = 1)
  delta <- 0.3

  result <- calc_asymptotic_variance(
    beta_0 = beta_0,
    beta_a = beta_a,
    theta2 = theta2,
    y2 = y2,
    X2 = X2,
    U2 = U2,
    delta = delta,
    include_estimation_uncertainty = FALSE,
    model = "poisson"
  )

  expect_equal(result, 0.5964023253138986, tolerance = 1e-10)
})
