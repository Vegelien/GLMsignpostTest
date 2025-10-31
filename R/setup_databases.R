# ============================================================================
# Database Schema Setup
# ============================================================================
# This script creates the schemas for two separate databases:
# 1. power_simulations.db - for power calculations
# 2. estimation_simulations.db - for theta_hat and loss calculations
# ============================================================================

library(DBI)
library(RSQLite)

# ============================================================================
# POWER SIMULATIONS DATABASE
# ============================================================================

create_power_database_schema <- function(con) {
  
  # Parameters table
  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS parameters (
      id INTEGER PRIMARY KEY AUTOINCREMENT,
      n INTEGER NOT NULL,
      p INTEGER NOT NULL,
      lambda REAL NOT NULL,
      test_type TEXT NOT NULL,
      GLM_model TEXT NOT NULL,
      model_specification TEXT NOT NULL,
      n_beta INTEGER,
      lambda_beta REAL,
      UNIQUE(n, p, lambda, test_type, GLM_model, model_specification, n_beta, lambda_beta)
    );
  ")
  
  # Results table - stores p-values for each gamma
  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS results (
      id INTEGER PRIMARY KEY AUTOINCREMENT,
      param_id INTEGER NOT NULL,
      replicate_id INTEGER NOT NULL,
      p_value00 REAL,
      p_value01 REAL,
      p_value02 REAL,
      p_value03 REAL,
      p_value04 REAL,
      p_value05 REAL,
      FOREIGN KEY(param_id) REFERENCES parameters(id)
    );
  ")
  
  # Create index for faster queries
  dbExecute(con, "
    CREATE INDEX IF NOT EXISTS idx_results_param_id ON results(param_id);
  ")
  
  message("✓ Power database schema created successfully")
}


# ============================================================================
# ESTIMATION SIMULATIONS DATABASE
# ============================================================================

create_estimation_database_schema <- function(con) {
  
  # Parameters table
  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS parameters (
      id INTEGER PRIMARY KEY AUTOINCREMENT,
      n INTEGER NOT NULL,
      p INTEGER NOT NULL,
      lambda REAL NOT NULL,
      GLM_model TEXT NOT NULL,
      model_specification TEXT NOT NULL,
      n_beta INTEGER,
      lambda_beta REAL,
      UNIQUE(n, p, lambda, GLM_model, model_specification, n_beta, lambda_beta)
    );
  ")
  
  # Estimation settings table
  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS estimation_settings (
      id INTEGER PRIMARY KEY AUTOINCREMENT,
      lambda_est REAL NOT NULL,
      UNIQUE(lambda_est)
    );
  ")
  
  # Results table - stores theta_hat and losses for each gamma
  # gamma values: 0, 0.1, 0.2, ..., 1.0 (11 values total)
  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS results (
      id INTEGER PRIMARY KEY AUTOINCREMENT,
      param_id INTEGER NOT NULL,
      estimation_id INTEGER NOT NULL,
      replicate_id INTEGER NOT NULL,
      theta_hat00 REAL,
      theta_hat01 REAL,
      theta_hat02 REAL,
      theta_hat03 REAL,
      theta_hat04 REAL,
      theta_hat05 REAL,
      theta_hat06 REAL,
      theta_hat07 REAL,
      theta_hat08 REAL,
      theta_hat09 REAL,
      theta_hat10 REAL,
      target_loss00 REAL,
      target_loss01 REAL,
      target_loss02 REAL,
      target_loss03 REAL,
      target_loss04 REAL,
      target_loss05 REAL,
      target_loss06 REAL,
      target_loss07 REAL,
      target_loss08 REAL,
      target_loss09 REAL,
      target_loss10 REAL,
      null_loss00 REAL,
      null_loss01 REAL,
      null_loss02 REAL,
      null_loss03 REAL,
      null_loss04 REAL,
      null_loss05 REAL,
      null_loss06 REAL,
      null_loss07 REAL,
      null_loss08 REAL,
      null_loss09 REAL,
      null_loss10 REAL,
      FOREIGN KEY(param_id) REFERENCES parameters(id),
      FOREIGN KEY(estimation_id) REFERENCES estimation_settings(id)
    );
  ")
  
  # Create indices for faster queries
  dbExecute(con, "
    CREATE INDEX IF NOT EXISTS idx_results_param_id ON results(param_id);
  ")
  
  dbExecute(con, "
    CREATE INDEX IF NOT EXISTS idx_results_estimation_id ON results(estimation_id);
  ")
  
  message("✓ Estimation database schema created successfully")
}


# ============================================================================
# HELPER FUNCTIONS FOR DATABASE OPERATIONS
# ============================================================================

#' Get or Insert Parameter ID (Power Database)
get_or_insert_power_param_id <- function(con, n, p, lambda, test_type, GLM_model, 
                                         model_specification, n_beta = NULL, lambda_beta = NULL) {
  if (is.null(n_beta)) n_beta <- NA_integer_
  if (is.null(lambda_beta)) lambda_beta <- NA_real_
  
  # Check if parameters already exist
  query <- "SELECT id FROM parameters 
            WHERE n = ? AND p = ? AND lambda = ? AND test_type = ? 
            AND GLM_model = ? AND model_specification = ? 
            AND (n_beta IS ? OR n_beta = ?) 
            AND (lambda_beta IS ? OR lambda_beta = ?);"
  
  existing_id <- dbGetQuery(con, query, params = list(
    n, p, lambda, test_type, GLM_model, model_specification, 
    n_beta, n_beta, lambda_beta, lambda_beta
  ))
  
  if (nrow(existing_id) > 0) {
    return(existing_id$id[1])
  } else {
    # Insert new parameters
    insert_query <- "INSERT INTO parameters 
                     (n, p, lambda, test_type, GLM_model, model_specification, n_beta, lambda_beta)
                     VALUES (?, ?, ?, ?, ?, ?, ?, ?);"
    
    dbExecute(con, insert_query, params = list(
      n, p, lambda, test_type, GLM_model, model_specification, n_beta, lambda_beta
    ))
    
    return(dbGetQuery(con, "SELECT last_insert_rowid();")[1, 1])
  }
}


#' Get or Insert Parameter ID (Estimation Database)
get_or_insert_estimation_param_id <- function(con, n, p, lambda, GLM_model, 
                                               model_specification, n_beta = NULL, lambda_beta = NULL) {
  if (is.null(n_beta)) n_beta <- NA_integer_
  if (is.null(lambda_beta)) lambda_beta <- NA_real_
  
  # Check if parameters already exist
  query <- "SELECT id FROM parameters 
            WHERE n = ? AND p = ? AND lambda = ? 
            AND GLM_model = ? AND model_specification = ? 
            AND (n_beta IS ? OR n_beta = ?) 
            AND (lambda_beta IS ? OR lambda_beta = ?);"
  
  existing_id <- dbGetQuery(con, query, params = list(
    n, p, lambda, GLM_model, model_specification, 
    n_beta, n_beta, lambda_beta, lambda_beta
  ))
  
  if (nrow(existing_id) > 0) {
    return(existing_id$id[1])
  } else {
    # Insert new parameters
    insert_query <- "INSERT INTO parameters 
                     (n, p, lambda, GLM_model, model_specification, n_beta, lambda_beta)
                     VALUES (?, ?, ?, ?, ?, ?, ?);"
    
    dbExecute(con, insert_query, params = list(
      n, p, lambda, GLM_model, model_specification, n_beta, lambda_beta
    ))
    
    return(dbGetQuery(con, "SELECT last_insert_rowid();")[1, 1])
  }
}


#' Get or Insert Estimation Settings ID
get_or_insert_estimation_settings_id <- function(con, lambda_est) {
  # Check if estimation settings already exist
  query <- "SELECT id FROM estimation_settings WHERE lambda_est = ?;"
  existing_id <- dbGetQuery(con, query, params = list(lambda_est))
  
  if (nrow(existing_id) > 0) {
    return(existing_id$id[1])
  } else {
    # Insert new estimation settings
    insert_query <- "INSERT INTO estimation_settings (lambda_est) VALUES (?);"
    dbExecute(con, insert_query, params = list(lambda_est))
    return(dbGetQuery(con, "SELECT last_insert_rowid();")[1, 1])
  }
}


#' Insert Power Results Batch
insert_power_results_batch <- function(con, param_id, p_val_matrix, replicate_start = 1) {
  if (!is.matrix(p_val_matrix)) {
    stop("Error: p_val_matrix must be a matrix.")
  }
  
  num_gammas <- ncol(p_val_matrix)
  col_names <- paste0("p_value", sprintf("%02d", seq_len(num_gammas) - 1))
  
  query <- sprintf(
    "INSERT INTO results (param_id, replicate_id, %s) VALUES (%s);",
    paste(col_names, collapse = ", "),
    paste(rep("?", num_gammas + 2), collapse = ", ")
  )
  
  stmt <- dbSendStatement(con, query)

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

    params <- c(list(param_id, replicate_start + i - 1), row_values)
    dbBind(stmt, params)
  }
  
  dbClearResult(stmt)
}


#' Insert Estimation Results Batch
insert_estimation_results_batch <- function(con, param_id, estimation_id, 
                                           theta_hat_matrix, target_loss_matrix, 
                                           null_loss_matrix, replicate_start = 1) {
  if (!is.matrix(theta_hat_matrix) || !is.matrix(target_loss_matrix) || !is.matrix(null_loss_matrix)) {
    stop("Error: All input matrices must be matrices.")
  }
  
  num_gammas <- ncol(theta_hat_matrix)

  if (ncol(target_loss_matrix) != num_gammas || ncol(null_loss_matrix) != num_gammas) {
    stop(
      sprintf(
        "Column count mismatch: expected %d columns for each matrix (theta_hat, target_loss, null_loss).",
        num_gammas
      )
    )
  }

  if (nrow(target_loss_matrix) != nrow(theta_hat_matrix) ||
      nrow(null_loss_matrix) != nrow(theta_hat_matrix)) {
    stop("Row count mismatch across input matrices for binding.")
  }
  theta_cols <- paste0("theta_hat", sprintf("%02d", seq_len(num_gammas) - 1))
  target_loss_cols <- paste0("target_loss", sprintf("%02d", seq_len(num_gammas) - 1))
  null_loss_cols <- paste0("null_loss", sprintf("%02d", seq_len(num_gammas) - 1))
  
  all_cols <- c(theta_cols, target_loss_cols, null_loss_cols)
  
  query <- sprintf(
    "INSERT INTO results (param_id, estimation_id, replicate_id, %s) VALUES (%s);",
    paste(all_cols, collapse = ", "),
    paste(rep("?", length(all_cols) + 3), collapse = ", ")
  )
  
  stmt <- dbSendStatement(con, query)

  for (i in seq_len(nrow(theta_hat_matrix))) {
    theta_row <- as.list(theta_hat_matrix[i, ])
    target_row <- as.list(target_loss_matrix[i, ])
    null_row <- as.list(null_loss_matrix[i, ])

    if (length(theta_row) != num_gammas ||
        length(target_row) != num_gammas ||
        length(null_row) != num_gammas) {
      stop(
        sprintf(
          "Row %d has mismatched column counts (theta: %d, target_loss: %d, null_loss: %d; expected %d).",
          i,
          length(theta_row),
          length(target_row),
          length(null_row),
          num_gammas
        )
      )
    }

    params <- c(
      list(param_id, estimation_id, replicate_start + i - 1),
      theta_row,
      target_row,
      null_row
    )

    expected_param_count <- length(all_cols) + 3
    if (length(params) != expected_param_count) {
      stop(
        sprintf(
          "Row %d binding produced %d parameters but %d columns were expected.",
          i,
          length(params),
          expected_param_count
        )
      )
    }

    dbBind(stmt, params)
  }
  
  dbClearResult(stmt)
}


# ============================================================================
# INITIALIZATION FUNCTIONS
# ============================================================================

#' Initialize Power Database
initialize_power_database <- function(db_path = "power_simulations.db") {
  con <- dbConnect(SQLite(), db_path)
  create_power_database_schema(con)
  dbDisconnect(con)
  message(paste("Power database initialized at:", db_path))
}


#' Initialize Estimation Database
initialize_estimation_database <- function(db_path = "estimation_simulations.db") {
  con <- dbConnect(SQLite(), db_path)
  create_estimation_database_schema(con)
  dbDisconnect(con)
  message(paste("Estimation database initialized at:", db_path))
}


# ============================================================================
# USAGE EXAMPLE
# ============================================================================

if (FALSE) {
  # Initialize both databases
  setwd("/Users/alexandra/Documents/Werk/Documents/Wessel - project 1/")
  initialize_power_database("power_simulations.db")
  initialize_estimation_database("estimation_simulations.db")
}
