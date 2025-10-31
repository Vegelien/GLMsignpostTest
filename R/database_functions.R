
# Create a table for parameter specifications
# dbExecute(con, "
# CREATE TABLE IF NOT EXISTS parameters (
#     id INTEGER PRIMARY KEY AUTOINCREMENT,
#     n INTEGER,
#     p INTEGER,
#     lambda REAL,
#     test_type TEXT,
#     GLM_model TEXT,
#     B_null_distr  INTEGER,
#     X TEXT,
#     beta TEXT,
#     UNIQUE(n, p, lambda, test_type, GLM_model, B_null_distr, X, beta)
# );")
#
# # Create a table for simulation results
# dbExecute(con, "
# CREATE TABLE IF NOT EXISTS results (
#     id INTEGER PRIMARY KEY AUTOINCREMENT,
#     param_id INTEGER,
#     setting_id INTEGER,
#     p_value00 REAL,
#     p_value01 REAL,
#     p_value02 REAL,
#     p_value03 REAL,
#     p_value04 REAL,
#     p_value05 REAL,
#     FOREIGN KEY(param_id) REFERENCES parameters(id)
# );")
#
#










#' Get or Insert Parameter ID
#'
#' Retrieves the ID of an existing parameter set in the `parameters` table.
#' If the parameters do not exist, inserts them and returns the new ID.
#'
#' @param con A database connection object.
#' @param n Integer, the number of samples.
#' @param p Integer, the number of predictors.
#' @param lambda Numeric, regularization parameter.
#' @param test_type Character, type of test used.
#' @param GLM_model Character, the type of GLM model.
#' @param B_null_distr Integer, number of bootstrap samples.
#' @param X Character, representation of predictor variables.
#' @param beta Character, representation of model coefficients.
#' @param n_beta (Optional) Integer, sample size used to estimate beta_0 and beta_a.
#' @param lambda_beta (Optional) Numeric, ridge penalty used when estimating beta_0 and beta_a.
#'
#' @return Integer, the `id` corresponding to the given parameter set.
#' @export
get_or_insert_param_id <- function(con, n, p, lambda, test_type, GLM_model, B_null_distr, X, beta, 
                                   n_beta = NULL, lambda_beta = NULL) {
  if (is.null(n_beta)) n_beta <- NA_integer_
  if (is.null(lambda_beta)) lambda_beta <- NA_real_
  # Check if parameters already exist
  query <- "SELECT id FROM parameters WHERE n = ? AND p = ? AND lambda = ? AND
            test_type = ? AND GLM_model = ? AND B_null_distr = ? AND X = ? AND beta = ? AND
            ((? IS NULL AND n_beta IS NULL) OR n_beta = ?) AND ((? IS NULL AND lambda_beta IS NULL) OR lambda_beta = ?);"

  existing_id <- dbGetQuery(con, query, params = list(n, p, lambda, test_type, GLM_model, B_null_distr, X, beta,
                                                      n_beta, n_beta, lambda_beta, lambda_beta))
  
  if (nrow(existing_id) > 0) {
    return(existing_id$id[1])  # Return existing ID
  } else {
    # Insert new parameters and return the new ID
    insert_query <- "INSERT INTO parameters (n, p, lambda, test_type, GLM_model, B_null_distr, X, beta, n_beta, lambda_beta)
                     VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?);"
    
    dbExecute(con, insert_query, params = list(n, p, lambda, test_type, GLM_model, B_null_distr, X, beta, n_beta, lambda_beta))
    
    return(dbGetQuery(con, "SELECT last_insert_rowid();")[1, 1])  # Return new ID
  }
}



#' Insert a Single Simulation Result
#'
#' Inserts a single simulation result into the `results` table.
#'
#' @param con A database connection object.
#' @param param_id Integer, the ID of the parameter set (foreign key from `parameters` table).
#' @param result_value Numeric vector of length 6, containing the result values for different test cases.
#'
#' @return None. Inserts data into the database.
#' @export
insert_simulation_result <- function(con, param_id, result_value) {
  if (length(result_value) != 6) {
    stop("Error: result_value must be a vector of exactly 6 values.")
  }

  query <- "INSERT INTO results (param_id, p_value00, p_value01, p_value02, p_value03, p_value04, p_value05)
            VALUES (?, ?, ?, ?, ?, ?, ?);"

  dbExecute(con, query, params = c(param_id, result_value))
}


#' Insert Multiple Simulation Results in Batch
#'
#' Inserts multiple rows of simulation results into the `results` table efficiently.
#'
#' @param con A database connection object.
#' @param param_id Integer, the ID of the parameter set (same for all results).
#' @param p_val_matrix A matrix where each row represents a set of result values.
#'
#' @return None. Inserts data into the database.
#' @export
insert_simulation_results_batch <- function(con, param_id, p_val_matrix) {
  if (!is.matrix(p_val_matrix)) {
    stop("Error: p_val_matrix must be a matrix.")
  }

  num_gammas <- ncol(p_val_matrix)  # Number of gamma values (columns)

  # Dynamically create column names (p_value00, p_value01, ...)
  col_names <- paste0("p_value", sprintf("%02d", seq_len(num_gammas) - 1))
  placeholders <- paste(rep("?", num_gammas + 1), collapse = ", ")  # SQL placeholders

  query <- sprintf(
    "INSERT INTO results (param_id, %s) VALUES (%s);",
    paste(col_names, collapse = ", "),
    placeholders
  )

  stmt <- dbSendStatement(con, query)

  # Loop over each row and insert values
  for (i in seq_len(nrow(p_val_matrix))) {
    row_values <- as.list(p_val_matrix[i, ])

    if (length(row_values) != num_gammas) {
      stop(
        sprintf(
          "Row %d has %d values but %d columns were expected for binding.",
          i,
          length(row_values),
          num_gammas
        )
      )
    }

    params <- c(list(param_id), row_values)
    dbBind(stmt, params)
  }

  dbClearResult(stmt)  # Free up resources
}




# Example usage:
# param_id <- get_or_insert_param_id(con, 0.5, 1.2, "test")#need more params
# print(param_id)  # Prints the ID of the parameter set
#
#
# insert_simulation_result(con, param_id = 1, result_value = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06))


#---- multiple computers ----
# Create (or connect to) a local results database
# local_con <- dbConnect(SQLite(), "local_results.db")
#
# # Create results table (if it doesn't exist)
# dbExecute(local_con, "
# CREATE TABLE IF NOT EXISTS results (
#     id INTEGER PRIMARY KEY AUTOINCREMENT,
#     param_id INTEGER,
#     result_value REAL
# );")
#









#' Retrieve Simulation Results for Given Parameters
#'
#' Queries the `results` table to retrieve results corresponding to a given parameter set.
#'
#' @param con A database connection object.
#' @param n Integer, the number of samples.
#' @param p Integer, the number of predictors.
#' @param lambda Numeric, regularization parameter.
#' @param test_type Character, type of test used.
#' @param GLM_model Character, the type of GLM model.
#' @param B_null_distr Integer, number of bootstrap samples.
#' @param X Character, representation of predictor variables.
#' @param beta Character, representation of model coefficients.
#'
#' @return A data frame containing the simulation results for the given parameters.
#' @export
get_results_by_params <- function(con, n, p, lambda, test_type, GLM_model, B_null_distr, X, beta) {
  query <- "SELECT r.* FROM results r
            JOIN parameters p ON r.param_id = p.id
            WHERE p.n = ? AND p.p = ? AND p.lambda = ?
            AND p.test_type = ? AND p.GLM_model = ?
            AND p.B_null_distr = ? AND p.X = ? AND p.beta = ?;"

  results <- dbGetQuery(con, query, params = list(n, p, lambda, test_type, GLM_model, B_null_distr, X, beta))

  if (nrow(results) == 0) {
    message("No results found for the given parameters.")
  }

  return(results)
}

# Example usage:
# results <- get_results_by_params(con, 0.5, 1.2, "test")
# print(results)  # Shows all results for the given parameter set




# handle_database <- function(db_path) {
#   con <- dbConnect(RSQLite::SQLite(), db_path)
#   on.exit(dbDisconnect(con), add = TRUE)  # Ensures disconnection when function ends
#
#   # Perform database operations here
#   dbListTables(con)
# }
#
# handle_database("simulations.db")  # Connection auto-closes after this
#

