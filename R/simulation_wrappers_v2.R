# ============================================================================
# Simulation Wrappers Using Existing Functions
# ============================================================================
# This uses the existing functions from simulationFunctions.R:
# - power_signpost_test_extended_AS()
# - power_LR_test_extended()
# Plus creates a new function for theta_hat and loss estimation
# ============================================================================
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
#' Wraps [`power_signpost_test_extended_AS()`] and
#' [`power_LR_test_extended()`] to execute batches of simulations and persist
#' the resulting p-values to the power simulation database.
#'
#' @param db_path Path to the SQLite database that stores power simulations.
#' @param n_grid Integer vector of sample sizes.
#' @param p Integer number of predictors.
#' @param lambda Numeric ridge penalty for theta estimation.
#' @param test_types Character vector of test identifiers (e.g.
#'   `c("AS_SW_plugin", "LR")`).
#' @param model Character string selecting the GLM family (either
#'   `"logistic"` or `"poisson"`).
#' @param specifications Character vector describing whether each run should be
#'   well specified or misspecified.
#' @param n_beta_values List containing `NULL` (oracle) and/or integers for the
#'   working coefficient sample sizes.
#' @param lambda_beta Numeric ridge penalty used for coefficient estimation on
#'   the glmnet scale.
#' @param B_power Integer number of Monte Carlo replicates.
#' @param B_X Integer number of design matrices used by the AS test.
#' @param gammas Numeric vector of gamma values. Defaults to `seq(0, 0.5, by = 0.1)`.
#'
#' @return Invisibly returns `NULL`. Called for its side effects.
#' @export
#'
#' @importFrom DBI dbConnect dbDisconnect dbListTables
#' @importFrom RSQLite SQLite
#' @examples
#' \dontrun{
#' tmp_db <- tempfile(fileext = ".sqlite")
#' initialize_power_database(tmp_db)
#' run_and_store_power_simulations(db_path = tmp_db, B_power = 2, B_X = 2,
#'                                 n_grid = c(10, 20), gammas = c(0, 0.1))
#' }
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
          
          # Get or insert parameter ID
          param_id <- get_or_insert_power_param_id(
            con, n, p, lambda, test_type, model, spec, n_beta, lambda_beta
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

#' Run Estimation Simulations for Theta-hat and Loss Summaries
#'
#' Generates Monte Carlo samples to approximate the distribution of theta-hat
#' as well as target and null losses over a supplied gamma grid.
#'
#' @param n Integer sample size.
#' @param p Integer number of predictors.
#' @param lambda Numeric ridge penalty for theta estimation.
#' @param model Character string selecting the GLM family (either
#'   `"logistic"` or `"poisson"`).
#' @param misspecified Logical flag indicating whether the working model is
#'   misspecified.
#' @param n_beta Optional integer sample size used to estimate working
#'   coefficients. Use `NULL` for oracle coefficients.
#' @param lambda_beta Numeric ridge penalty (glmnet scale) for coefficient
#'   estimation.
#' @param lambda_est Numeric ridge penalty (glmnet scale) used when computing
#'   losses.
#' @param B_power Integer number of Monte Carlo replicates.
#' @param gammas Numeric vector of gamma values. Defaults to `seq(0, 1, by = 0.1)`.
#' @param scaled Logical, forwarded to the helper routines to signal whether
#'   the ridge penalties are provided on the scaled or unscaled scale.
#'
#' @return A list containing matrices `theta_hats`, `target_loss`, and
#'   `null_loss`.
#' @export
#'
#' @importFrom porridge create_progress_reporter
#' @examples
#' \dontrun{
#' sim <- estimation_simulation(
#'   n = 10, p = 20, lambda = Inf, model = "logistic",
#'   n_beta = NULL, lambda_beta = 2, lambda_est = 0.5,
#'   B_power = 2, gammas = c(0, 0.5)
#' )
#' str(sim)
#' }
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
#' Executes [estimation_simulation()] across grids of sample sizes, model
#' specifications, and coefficient estimation settings, storing the resulting
#' summaries in the estimation simulation database.
#'
#' @param db_path Path to the SQLite database that stores estimation
#'   simulations.
#' @param n_grid Integer vector of sample sizes.
#' @param p Integer number of predictors.
#' @param lambda Numeric ridge penalty for theta estimation.
#' @param model Character string selecting the GLM family (either
#'   `"logistic"` or `"poisson"`).
#' @param specifications Character vector describing whether each run should be
#'   well specified or misspecified.
#' @param n_beta_values List containing `NULL` (oracle) and/or integers for the
#'   working coefficient sample sizes.
#' @param lambda_beta Numeric ridge penalty (glmnet scale) used when estimating
#'   working coefficients.
#' @param lambda_est Numeric ridge penalty (glmnet scale) used in the loss
#'   calculation.
#' @param B_power Integer number of Monte Carlo replicates.
#' @param gammas Numeric vector of gamma values. Defaults to `seq(0, 1, by = 0.1)`.
#'
#' @return Invisibly returns `NULL`. Called for its side effects.
#' @export
#'
#' @importFrom DBI dbConnect dbDisconnect dbListTables
#' @importFrom RSQLite SQLite
#' @examples
#' \dontrun{
#' tmp_db <- tempfile(fileext = ".sqlite")
#' initialize_estimation_database(tmp_db)
#' run_and_store_estimation_simulations(
#'   db_path = tmp_db, B_power = 2, n_grid = c(10, 20), gammas = c(0, 0.5)
#' )
#' }
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
        
        # Get or insert parameter ID
        param_id <- get_or_insert_estimation_param_id(
          con, n, p, lambda, model, spec, n_beta, lambda_beta
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
