test_that("intercept_policy auto flips intercept when constants dropped", {
  Y <- rnorm(6)
  X <- cbind(const = rep(3, 6), var = seq_len(6))

  mock_solver <- function(Y, X, lambda, target, model, U, ...) {
    nU <- if (is.null(U)) 0L else ncol(U)
    nX <- if (is.null(X)) 0L else ncol(X)
    rep(0, nU + nX)
  }

  fit <- ridge_complete(
    Y = Y,
    X = X,
    lambda = 0.1,
    model = "gaussian",
    intercept = FALSE,
    intercept_policy = "auto",
    solver_fun = mock_solver
  )

  expect_true(fit$meta$intercept)
  expect_true(fit$meta$intercept_auto_flipped)
  expect_false(1L %in% fit$kept)
})

test_that("intercept_policy keep_one retains exactly one constant column", {
  Y <- rnorm(6)
  X <- cbind(const = rep(5, 6), var = seq_len(6))

  mock_solver <- function(Y, X, lambda, target, model, U, ...) {
    nU <- if (is.null(U)) 0L else ncol(U)
    nX <- if (is.null(X)) 0L else ncol(X)
    rep(0, nU + nX)
  }

  fit <- ridge_complete(
    Y = Y,
    X = X,
    lambda = 0.1,
    model = "gaussian",
    intercept = FALSE,
    intercept_policy = "keep_one",
    solver_fun = mock_solver
  )

  expect_false(fit$meta$intercept)
  expect_identical(fit$meta$kept_constant, 1L)
  expect_true(1L %in% fit$kept)
  expect_equal(fit$mu[1L], 0)
  expect_equal(fit$s[1L], 1)
})

test_that("intercept_policy drop_all warns and preserves intercept choice", {
  Y <- rnorm(6)
  X <- cbind(const = rep(4, 6), var = seq_len(6))

  mock_solver <- function(Y, X, lambda, target, model, U, ...) {
    nU <- if (is.null(U)) 0L else ncol(U)
    nX <- if (is.null(X)) 0L else ncol(X)
    rep(0, nU + nX)
  }

  fit <- expect_warning(
    ridge_complete(
      Y = Y,
      X = X,
      lambda = 0.1,
      model = "gaussian",
      intercept = FALSE,
      intercept_policy = "drop_all",
      solver_fun = mock_solver
    ),
    "Dropping"
  )

  expect_false(fit$meta$intercept)
  expect_false(fit$meta$intercept_auto_flipped)
  expect_identical(fit$meta$kept_constant, integer(0))
  expect_true(1L %in% fit$dropped)
})
