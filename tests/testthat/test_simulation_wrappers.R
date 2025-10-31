test_that("power simulations append replicate ids sequentially", {
  db_path <- tempfile(fileext = ".db")
  on.exit(unlink(db_path), add = TRUE)

  stubbed_power <- function(n, p, model, lambda, B_power, ...) {
    list(p_values = matrix(seq_len(B_power), nrow = B_power, ncol = 1))
  }

  with_mocked_bindings({
    GLMsignpostTest:::run_and_store_power_simulations(
      db_path = db_path,
      n_grid = 5,
      p = 2,
      lambda = Inf,
      test_types = "LR",
      model = "logistic",
      specifications = "well_specified",
      n_beta_values = list(NULL),
      lambda_beta = 1,
      B_power = 1,
      B_X = 1,
      gammas = 0
    )
  }, power_LR_test_extended = stubbed_power)

  get_ids <- function() {
    con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
    on.exit(DBI::dbDisconnect(con), add = TRUE)
    DBI::dbGetQuery(con, "SELECT replicate_id FROM results ORDER BY replicate_id")$replicate_id
  }

  expect_identical(get_ids(), 1L)

  with_mocked_bindings({
    GLMsignpostTest:::run_and_store_power_simulations(
      db_path = db_path,
      n_grid = 5,
      p = 2,
      lambda = Inf,
      test_types = "LR",
      model = "logistic",
      specifications = "well_specified",
      n_beta_values = list(NULL),
      lambda_beta = 1,
      B_power = 2,
      B_X = 1,
      gammas = 0
    )
  }, power_LR_test_extended = stubbed_power)

  expect_identical(get_ids(), 1:3)
})


test_that("estimation simulations append replicate ids sequentially", {
  db_path <- tempfile(fileext = ".db")
  on.exit(unlink(db_path), add = TRUE)

  stubbed_estimation <- function(n, p, lambda, model, misspecified = FALSE,
                                 n_beta = NULL, lambda_beta = 1,
                                 lambda_est = 0.5, B_power = 10,
                                 gammas = 0, scaled = TRUE) {
    ncols <- length(gammas)
    list(
      theta_hats = matrix(seq_len(B_power * ncols), nrow = B_power, ncol = ncols),
      target_loss = matrix(seq_len(B_power * ncols), nrow = B_power, ncol = ncols),
      null_loss = matrix(seq_len(B_power * ncols), nrow = B_power, ncol = ncols)
    )
  }

  with_mocked_bindings({
    GLMsignpostTest:::run_and_store_estimation_simulations(
      db_path = db_path,
      n_grid = 5,
      p = 2,
      lambda = Inf,
      model = "logistic",
      specifications = "well_specified",
      n_beta_values = list(NULL),
      lambda_beta = 1,
      lambda_est = 0.5,
      B_power = 1,
      gammas = 0
    )
  }, estimation_simulation = stubbed_estimation)

  get_ids <- function() {
    con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
    on.exit(DBI::dbDisconnect(con), add = TRUE)
    DBI::dbGetQuery(con, "SELECT replicate_id FROM results ORDER BY replicate_id")$replicate_id
  }

  expect_identical(get_ids(), 1L)

  with_mocked_bindings({
    GLMsignpostTest:::run_and_store_estimation_simulations(
      db_path = db_path,
      n_grid = 5,
      p = 2,
      lambda = Inf,
      model = "logistic",
      specifications = "well_specified",
      n_beta_values = list(NULL),
      lambda_beta = 1,
      lambda_est = 0.5,
      B_power = 2,
      gammas = 0
    )
  }, estimation_simulation = stubbed_estimation)

  expect_identical(get_ids(), 1:3)
})
