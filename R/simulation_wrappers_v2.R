# ============================================================================
# Simulation Wrappers Using Existing Functions
# ============================================================================
# This uses the existing functions from simulationFunctions.R:
# - power_signpost_test_extended_AS()
# - power_LR_test_extended()
# Plus creates a new function for theta_hat and loss estimation
# ============================================================================

library(DBI)
library(RSQLite)
library(dplyr)
library(porridge)
# Source your files (adjust paths as needed)
# source("/Users/alexandra/Documents/Werk/Documents/Wessel - project 1/GLMsignpostTest/R/simulationFunctions.R")
# source("/Users/alexandra/Documents/Werk/Documents/Wessel - project 1/GLMsignpostTest/R/signpostTestFunctions.R")
# source("/Users/alexandra/Documents/Werk/Documents/Wessel - project 1/GLMsignpostTest/R/ridgePoisson.R")  # Contains ridge_efficient
# source("/Users/alexandra/Documents/Werk/Documents/Wessel - project 1/GLMsignpostTest/R/prediction_simulation.R")  # Contains ridge_efficient
# source("/Users/alexandra/Documents/Werk/Documents/Wessel - project 1/setup_databases.R")
# 
# 
# setwd("/Users/alexandra/Documents/Werk/Documents/Wessel - project 1/")



# ============================================================================
# POWER SIMULATIONS (using existing functions)
# ============================================================================

#' Run and Store Power Simulations
#' 
#' Uses power_signpost_test_extended_AS() and power_LR_test_extended()
#' from simulationFunctions.R but stores ONLY p-values
#' 
#' @param db_path Path to power database
#' @param n_grid Vector of sample sizes
#' @param p Number of predictors
#' @param lambda Ridge penalty for theta estimation
#' @param test_types c("AS_SW_plugin", "LR")
#' @param model "logistic" or "poisson"
#' @param specifications c("well_specified", "misspecified")
#' @param n_beta_values list(NULL, 300) - NULL=oracle, number=estimated
#' @param lambda_beta Ridge penalty for beta estimation (glmnet scale)
#' @param B_power Number of replicates
#' @param B_X Number of design matrices (for AS test only)
#' @param gammas Gamma grid (default: 0, 0.1, ..., 0.5)
run_and_store_power_simulations <- function(
    db_path = "power_simulations.db",
    n_grid = c(10, 20, 40, 60, 80, 100, 150),
    p = 100,
    lambda = Inf,
    test_types = c("AS_SW_plugin", "LR"),
    model = "logistic",
    specifications = c("well_specified", "misspecified"),
    n_beta_values = list(NULL, 300),
    lambda_beta = 2,
    B_power = 10,
    B_X = 10,
    gammas = seq(0, 0.5, by = 0.1)
) {
  
  con <- dbConnect(SQLite(), db_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  # Ensure schema exists
  if (length(dbListTables(con)) == 0) {
    create_power_database_schema(con)
  }
  
  total_sims <- length(n_grid) * length(test_types) * length(specifications) * length(n_beta_values)
  sim_count <- 0
  
  for (n in n_grid) {
    for (test_type in test_types) {
      for (spec in specifications) {
        for (n_beta in n_beta_values) {
          
          sim_count <- sim_count + 1
          message(sprintf("\n[%d/%d] Running: n=%d, test=%s, spec=%s, n_beta=%s",
                          sim_count, total_sims, n, test_type, spec, 
                          ifelse(is.null(n_beta), "oracle", as.character(n_beta))))
          
          # Determine beta source for database storage
          beta_source <- if (is.null(n_beta)) "oracle" else "estimated"

          # Get or insert parameter ID
          param_id <- get_or_insert_power_param_id(
            con, n, p, lambda, test_type, model, spec, n_beta, lambda_beta,
            beta_source = beta_source
          )

          # Determine parameters for simulation
          estimate_beta <- !is.null(n_beta)
          misspecified <- (spec == "misspecified")
          
          # Run simulation based on test type
          if (test_type == "AS_SW_plugin") {
            result <- power_signpost_test_extended_AS(
              n = n, 
              p = p, 
              gammaH0 = 0,
              model = model, 
              lambda = lambda,  # glmnet scale
              B_X = B_X,
              B_power = B_power,
              n_beta = n_beta,
              lambda_beta = lambda_beta,  # glmnet scale
              estimate_beta = estimate_beta,
              use_theta_hat_in_variance = TRUE,
              fisher = FALSE,
              term1ONLY = TRUE,
              lambda_est = lambda ,  # Not used for power, but required parameter
              misspecified = misspecified,
              scaled = TRUE
            )
          } else if (test_type == "LR") {
            result <- power_LR_test_extended(
              n = n,
              p = p,
              model = model,
              lambda = lambda ,  # glmnet scale
              B_power = B_power,
              n_beta = n_beta,
              lambda_beta = lambda_beta,  # glmnet scale
              estimate_beta = estimate_beta,
              lambda_est = lambda ,  # Not used for power, but required parameter
              misspecified = misspecified,
              scaled = TRUE
            )
          }
          
          # Extract ONLY p-values and insert
          insert_power_results_batch(con, param_id, result$p_values)
          
          message(sprintf("  ✓ Stored %d replicates", nrow(result$p_values)))
        }
      }
    }
  }
  
  message("\n✓ All power simulations completed!")
}


# ============================================================================
# ESTIMATION SIMULATIONS (theta_hat and loss)
# ============================================================================

#' Estimation Simulation Function for Theta_hat and Loss
#' 
#' Similar to power functions but with gamma = 0, 0.1, ..., 1.0
#' and WITHOUT p-value calculation
#' 
#' @param n Sample size
#' @param p Number of predictors
#' @param lambda Ridge penalty for theta estimation
#' @param model "logistic" or "poisson"
#' @param misspecified Logical
#' @param n_beta Sample size for beta estimation (NULL for oracle)
#' @param lambda_beta Ridge penalty for beta estimation (glmnet scale)
#' @param lambda_est Ridge penalty for loss calculation (glmnet scale)
#' @param B_power Number of replicates
#' @param gammas Gamma grid (default: 0, 0.1, ..., 1.0)
#' @return List with theta_hats, target_loss, null_loss matrices
estimation_simulation <- function(n, p, lambda, model,
                                  misspecified = FALSE,
                                  n_beta = NULL, 
                                  lambda_beta = 2,
                                  lambda_est = 0.5,
                                  B_power = 10,
                                  gammas = seq(0, 1, by = 0.1),
                                  scaled = TRUE) {
  
  G <- length(gammas)
  estimate_beta <- !is.null(n_beta)
  
  # Preallocate matrices
  theta_hat_matrix    <- matrix(NA, nrow = B_power, ncol = G)
  target_loss_matrix  <- matrix(NA, nrow = B_power, ncol = G)
  null_loss_matrix    <- matrix(NA, nrow = B_power, ncol = G)
  
  progress <- create_progress_reporter(B_power, "Estimation simulation")
  on.exit(progress$close(), add = TRUE)
  
  for (i in seq_len(B_power)) {
    progress$update(i)
    
    # Generate true betas
    BETA_A <- rnorm(p, sd = 0.06)
    BETA_0 <- rnorm(p, sd = 0.06)
    scl <- sqrt(sum((BETA_A - BETA_0)^2))
    BETA_A <- 2 * BETA_A / scl
    BETA_0 <- 2 * BETA_0 / scl
    
    # Working beta_0 (misspecified or not)
    if (misspecified) {
      BETA_0_info <- 0 * BETA_0
    } else {
      BETA_0_info <- BETA_0
    }
    
    # Estimate beta_a if needed
    if (estimate_beta) {
      if (is.null(n_beta) || is.null(lambda_beta)) {
        stop("n_beta and lambda_beta must be specified when estimate_beta = TRUE.")
      }
      X_a <- generate_X(n_beta, p)
      BETA_A_hat <- generate_beta_est(BETA_A, n_beta, X_a, model, 
                                      lambda_beta, scaled = scaled)
    } else {
      BETA_A_hat <- BETA_A
    }
    
    # Generate design matrix for test
    X <- generate_X(n, p)
    
    # Loop over gamma values
    for (j in seq_len(G)) {
      # Generate data from TRUE betas
      eta_j <- calculate_eta(X, gammas[j], BETA_0, BETA_A)
      Y <- generate_Y(eta_j, model, n)
      
      # Estimate theta_hat using WORKING betas
      if (is.infinite(lambda)) {
        theta_hat <- theta_inf_hat(Y, X, matrix(ncol = 0, nrow = n),
                                   BETA_0_info, BETA_A_hat, model)
      } else {
        theta_hat <- theta_hat_lambda(Y, X, BETA_0_info, BETA_A_hat, 
                                      lambda , model)
      }
      
      # Calculate losses using TRUE betas for target
      true_beta <- (1 - gammas[j]) * BETA_0 + gammas[j] * BETA_A
      
      loss_target <- squared_loss_beta(
        theta_hat, Y, X, BETA_0_info, BETA_A_hat,
        lambda_est, true_beta, model
      )
      
      loss_null <- squared_loss_beta(
        0, Y, X, BETA_0_info, BETA_A_hat,
        lambda_est, true_beta, model
      )
      
      # Store results
      theta_hat_matrix[i, j] <- theta_hat
      target_loss_matrix[i, j] <- loss_target
      null_loss_matrix[i, j] <- loss_null
    }
  }
  
  return(list(
    theta_hats = theta_hat_matrix,
    target_loss = target_loss_matrix,
    null_loss = null_loss_matrix
  ))
}


#' Run and Store Estimation Simulations
#' 
#' @param db_path Path to estimation database
#' @param n_grid Vector of sample sizes
#' @param p Number of predictors
#' @param lambda Ridge penalty for theta estimation
#' @param model "logistic" or "poisson"
#' @param specifications c("well_specified", "misspecified")
#' @param n_beta_values list(NULL, 300) - NULL=oracle, number=estimated
#' @param lambda_beta Ridge penalty for beta estimation (glmnet scale)
#' @param lambda_est Ridge penalty for loss calculation (glmnet scale)
#' @param B_power Number of replicates
#' @param gammas Gamma grid (default: 0, 0.1, ..., 1.0)
run_and_store_estimation_simulations <- function(
    db_path = "estimation_simulations.db",
    n_grid = c(10, 50, 150),
    p = 100,
    lambda = Inf,
    model = "logistic",
    specifications = c("well_specified", "misspecified"),
    n_beta_values = list(NULL, 300),
    lambda_beta = 2,
    lambda_est = 0.5,
    B_power = 10,
    gammas = seq(0, 1, by = 0.1)
) {
  
  con <- dbConnect(SQLite(), db_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  # Ensure schema exists
  if (length(dbListTables(con)) == 0) {
    create_estimation_database_schema(con)
  }
  
  # Get or insert estimation settings ID
  estimation_id <- get_or_insert_estimation_settings_id(con, lambda_est)
  
  total_sims <- length(n_grid) * length(specifications) * length(n_beta_values)
  sim_count <- 0
  
  for (n in n_grid) {
    for (spec in specifications) {
      for (n_beta in n_beta_values) {
        
        sim_count <- sim_count + 1
        message(sprintf("\n[%d/%d] Running: n=%d, spec=%s, n_beta=%s",
                        sim_count, total_sims, n, spec, 
                        ifelse(is.null(n_beta), "oracle", as.character(n_beta))))
        
        # Determine beta source for database storage
        beta_source <- if (is.null(n_beta)) "oracle" else "estimated"

        # Get or insert parameter ID
        param_id <- get_or_insert_estimation_param_id(
          con, n, p, lambda, model, spec, n_beta, lambda_beta,
          beta_source = beta_source
        )
        
        # Run simulation
        misspecified <- (spec == "misspecified")
        
        result <- estimation_simulation(
          n = n, 
          p = p, 
          lambda = lambda,  # Will be multiplied by n inside function
          model = model,
          misspecified = misspecified,
          n_beta = n_beta,
          lambda_beta = lambda_beta,  # Will be multiplied by n_beta inside function
          lambda_est = lambda_est,    # Will be multiplied by n inside function
          B_power = B_power,
          gammas = gammas,
          scaled = TRUE
        )
        
        # Insert results
        insert_estimation_results_batch(
          con, param_id, estimation_id, 
          result$theta_hats, result$target_loss, result$null_loss
        )
        
        message(sprintf("  ✓ Stored %d replicates", nrow(result$theta_hats)))
      }
    }
  }
  
  message("\n✓ All estimation simulations completed!")
}
