skip_if_not_installed("porridge")

set.seed(20240520)

test_that("theta_inf_hat_direct matches glm offset for default method", {
  n <- 400
  p <- 3
  X <- matrix(rnorm(n * p), n, p)
  beta_0 <- rnorm(p)
  beta_a <- beta_0 + rnorm(p, sd = 0.5)

  base <- drop(X %*% beta_0)
  diff <- drop(X %*% (beta_a - beta_0))

  U <- matrix(rnorm(n * 2), n, 2)
  colnames(U) <- c("z1", "z2")

  theta_true <- 0.35
  gamma_true <- c(-0.25, 0.4)
  linpred <- base + theta_true * diff + drop(U %*% gamma_true)
  prob <- plogis(linpred)
  Y <- rbinom(n, 1, prob)

  df <- data.frame(
    Y = Y,
    z1 = U[, 1],
    z2 = U[, 2],
    diff = diff,
    base = base
  )
  fit_glm <- glm(
    Y ~ z1 + z2 + diff,
    family = binomial(),
    data = df,
    offset = df$base
  )

  theta_glmoffset <- theta_inf_hat_direct(
    Y = Y,
    X = X,
    U = U,
    beta_0 = beta_0,
    beta_a = beta_a,
    model = "logistic"
  )

  expect_equal(
    theta_glmoffset,
    unname(coef(fit_glm)["diff"]),
    tolerance = 1e-6
  )
})

test_that("theta_inf_hat_direct ridge method matches glm offset", {
  n <- 400
  p <- 3
  X <- matrix(rnorm(n * p), n, p)
  beta_0 <- rnorm(p)
  beta_a <- beta_0 + rnorm(p, sd = 0.5)

  base <- drop(X %*% beta_0)
  diff <- drop(X %*% (beta_a - beta_0))

  U <- matrix(rnorm(n * 2), n, 2)
  colnames(U) <- c("z1", "z2")

  theta_true <- 0.35
  gamma_true <- c(-0.25, 0.4)
  linpred <- base + theta_true * diff + drop(U %*% gamma_true)
  prob <- plogis(linpred)
  Y <- rbinom(n, 1, prob)

  df <- data.frame(
    Y = Y,
    z1 = U[, 1],
    z2 = U[, 2],
    diff = diff,
    base = base
  )
  fit_glm <- glm(
    Y ~ z1 + z2 + diff,
    family = binomial(),
    data = df,
    offset = df$base
  )

  theta_ridge <- theta_inf_hat_direct(
    Y = Y,
    X = X,
    U = U,
    beta_0 = beta_0,
    beta_a = beta_a,
    model = "logistic",
    method = "ridge"
  )

  expect_equal(
    theta_ridge,
    unname(coef(fit_glm)["diff"]),
    tolerance = 1e-6
  )
})
