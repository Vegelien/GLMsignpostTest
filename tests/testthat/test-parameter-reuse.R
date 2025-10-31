library(DBI)
library(RSQLite)

setup_parameters_table <- function(con) {
  dbExecute(con, "CREATE TABLE parameters (
               id INTEGER PRIMARY KEY AUTOINCREMENT,
               n INTEGER NOT NULL,
               p INTEGER NOT NULL,
               lambda REAL NOT NULL,
               test_type TEXT NOT NULL,
               GLM_model TEXT NOT NULL,
               B_null_distr INTEGER NOT NULL,
               X TEXT NOT NULL,
               beta TEXT NOT NULL,
               n_beta INTEGER,
               lambda_beta REAL
             );")
}

test_that("get_or_insert_param_id reuses rows for nullable parameters", {
  con <- dbConnect(SQLite(), ":memory:")
  on.exit(dbDisconnect(con))

  setup_parameters_table(con)

  first_id <- get_or_insert_param_id(
    con = con,
    n = 100,
    p = 5,
    lambda = 0.1,
    test_type = "wald",
    GLM_model = "logistic",
    B_null_distr = 500,
    X = "design",
    beta = "coef",
    n_beta = NULL,
    lambda_beta = NULL
  )

  expect_true(is.numeric(first_id))

  second_id <- get_or_insert_param_id(
    con = con,
    n = 100,
    p = 5,
    lambda = 0.1,
    test_type = "wald",
    GLM_model = "logistic",
    B_null_distr = 500,
    X = "design",
    beta = "coef",
    n_beta = NULL,
    lambda_beta = NULL
  )

  expect_identical(second_id, first_id)

  third_id <- get_or_insert_param_id(
    con = con,
    n = 100,
    p = 5,
    lambda = 0.1,
    test_type = "wald",
    GLM_model = "logistic",
    B_null_distr = 500,
    X = "design",
    beta = "coef",
    n_beta = 20,
    lambda_beta = NULL
  )

  fourth_id <- get_or_insert_param_id(
    con = con,
    n = 100,
    p = 5,
    lambda = 0.1,
    test_type = "wald",
    GLM_model = "logistic",
    B_null_distr = 500,
    X = "design",
    beta = "coef",
    n_beta = 20,
    lambda_beta = NULL
  )

  expect_identical(fourth_id, third_id)

  ids <- dbGetQuery(con, "SELECT id FROM parameters ORDER BY id")
  expect_equal(nrow(ids), 2)
})

test_that("power parameter retrieval reuses nullable rows", {
  con <- dbConnect(SQLite(), ":memory:")
  on.exit(dbDisconnect(con))

  GLMsignpostTest:::create_power_database_schema(con)

  id_one <- GLMsignpostTest:::get_or_insert_power_param_id(
    con = con,
    n = 200,
    p = 10,
    lambda = 0.25,
    test_type = "lrt",
    GLM_model = "poisson",
    model_specification = "spec",
    n_beta = NULL,
    lambda_beta = NULL
  )

  id_two <- GLMsignpostTest:::get_or_insert_power_param_id(
    con = con,
    n = 200,
    p = 10,
    lambda = 0.25,
    test_type = "lrt",
    GLM_model = "poisson",
    model_specification = "spec",
    n_beta = NULL,
    lambda_beta = NULL
  )

  expect_identical(id_two, id_one)

  expect_equal(
    dbGetQuery(con, "SELECT COUNT(*) AS n FROM parameters")$n,
    1
  )
})

test_that("estimation parameter retrieval reuses nullable rows", {
  con <- dbConnect(SQLite(), ":memory:")
  on.exit(dbDisconnect(con))

  GLMsignpostTest:::create_estimation_database_schema(con)

  first <- GLMsignpostTest:::get_or_insert_estimation_param_id(
    con = con,
    n = 50,
    p = 4,
    lambda = 0.05,
    GLM_model = "gaussian",
    model_specification = "alt",
    n_beta = NULL,
    lambda_beta = NULL
  )

  second <- GLMsignpostTest:::get_or_insert_estimation_param_id(
    con = con,
    n = 50,
    p = 4,
    lambda = 0.05,
    GLM_model = "gaussian",
    model_specification = "alt",
    n_beta = NULL,
    lambda_beta = NULL
  )

  expect_identical(second, first)

  expect_equal(
    dbGetQuery(con, "SELECT COUNT(*) AS n FROM parameters")$n,
    1
  )
})
