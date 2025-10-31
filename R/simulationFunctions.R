#---- generation functions ----


#' Sample/generate normal design matrix
#'
#' @param n positive integer, sample size
#' @param p positive integer, dimension
#' @param sd positive numeric, standard deviation of columns, DEFAULT = 1
#'
#' @return X design matrix of the penalized covariates.
#' @export
generate_X = function(n,p, sd =1){
  X = matrix( rnorm(n*p,mean=0,sd=sd), n, p)
  X[,1] = 1
  return(X)
}

#' Calculate eta
#'
#' this function calculates X * ((1-theta) beta_0 + theta . beta_a)
#' @param X design matrix of the penalized covariates.
#' @param theta numeric, in [0,1] weight parameter of targets beta_0 and beta_a
#' @param beta_0 numeric, The number of elements should match the number of columns of X.
#' @param beta_a numeric, The number of elements should match the number of columns of X.
#'
#' @return eta, numeric with same number of elements as rows in X
#' @export
calculate_eta = function(X,theta, beta_0,beta_a){
  beta_true = (1-theta)*beta_0 + theta*beta_a
  eta = X%*%beta_true #Canonical so theta = eta = X*beta_true
  return(eta)
}

#' Generate Response Vector Y
#'
#' This function generates a response vector Y from a given linear predictor (eta)
#' according to the specified model. For a logistic model, it generates binary responses
#' using a logistic link; for a Poisson model, it generates count data; and for a linear model,
#' it adds normally distributed error to eta.
#'
#' @param eta A numeric vector representing the linear predictor.
#' @param model A character string specifying the model type. Must be one of
#'   \code{"linear"}, \code{"logistic"}, or \code{"poisson"}.
#' @param n A positive integer indicating the sample size. Must be equal to \code{length(eta)}.
#' @param sigma Optional. A positive numeric value specifying the standard deviation of the
#'   noise for the linear model. Default is 1. Only used when \code{model = "linear"}.
#'
#' @return A numeric vector containing the generated responses.
#'
#' @examples
#' # Logistic model example:
#' eta <- rnorm(100)
#' Y_logistic <- generate_Y(eta, model = "logistic", n = 100)
#'
#' # Poisson model example:
#' eta <- rnorm(100)
#' Y_poisson <- generate_Y(eta, model = "poisson", n = 100)
#'
#' # Linear model example:
#' eta <- rnorm(100)
#' Y_linear <- generate_Y(eta, model = "linear", n = 100, sigma = 2)
#'
#' @export
generate_Y <- function(eta, model, n, sigma = 1) {
  # Input validation
  if (!is.numeric(eta))
    stop("eta must be a numeric vector.")
  if (!is.numeric(n) || length(n) != 1 || n <= 0)
    stop("n must be a positive integer.")
  if (length(eta) != n)
    stop("Length of eta must equal n.")

  # Ensure model is one of the supported types (case-insensitive)
  model <- match.arg(tolower(model), choices = c("linear", "logistic", "poisson"))

  if (model == "logistic") {
    # Compute probabilities using the logistic function
    mu <- exp(eta) / (1 + exp(eta))
    # Generate binary responses vectorized over mu
    Y <- rbinom(n, size = 1, prob = mu)
  } else if (model == "poisson") {
    # Compute mean counts
    mu <- exp(eta)
    # Generate count responses vectorized over mu
    Y <- rpois(n, lambda = mu)
  } else if (model == "linear") {
    # Generate continuous responses with additive Gaussian noise
    Y <- eta + rnorm(n, mean = 0, sd = sigma)
  } else {
    stop("Unsupported model specified. Choose 'linear', 'logistic', or 'poisson'.")
  }

  return(Y)
}








#---- simulation to obtain power ----



#' #' Generate Estimated Beta Coefficients
#' #'
#' #' This function generates an estimate of the beta coefficients using ridge regression.
#' #' It first generates a design matrix `X`, computes the linear predictor `eta`,
#' #' and then simulates the response variable `Y` based on the given model.
#' #' Finally, it estimates the beta coefficients using `ridgeGLM`.
#' #'
#' #' @param true_beta Numeric vector, the true beta coefficients.
#' #' @param n_beta Integer, the number of samples used for estimation.
#' #' @param p Integer, the number of predictors.
#' #' @param model Character, the type of model, e.g., "linear", "logistic", or "poisson".
#' #' @param lambda_beta Numeric, the regularization parameter for ridge regression.
#' #' @importFrom porridge ridgeGLM
#' #' @return A numeric vector of estimated beta coefficients.
#' generate_beta_est <- function(true_beta, n_beta, p, model, lambda_beta) {
#'   X <- generate_X(n_beta, p)  # Generate design matrix
#'   eta <- calculate_eta(X, 0, true_beta, true_beta)  # Compute linear predictor
#'   Y <- generate_Y(eta, model, n_beta)  # Generate response variable
#'   return(ridgeGLM(Y = Y, X = X, lambda = lambda_beta, model = model))  # Estimate beta coefficients
#' }


#' #' Generate Estimated Beta Coefficients (Updated)
#' #'
#' #' Updated version that uses ridge_efficient instead of ridgeGLM directly.
#' #'
#' #' @param true_beta Numeric vector, the true beta coefficients (length p)
#' #' @param n_beta Integer, sample size for estimation
#' #' @param X Numeric matrix, design matrix (n_beta x p)
#' #' @param model Character, one of "linear", "logistic", or "poisson"
#' #' @param lambda_beta Numeric, ridge penalty parameter (glmnet scale)
#' #'
#' #' @return Numeric vector, normalized estimated beta coefficients
#' #' @export
#' generate_beta_est <- function(true_beta, n_beta, X, model, lambda_beta) {
#'   eta <- calculate_eta(X, 0, true_beta, true_beta)  # Compute linear predictor
#'   Y <- generate_Y(eta, model, n_beta)  # Generate response variable
#'   
#'   # Use ridge_efficient instead of ridgeGLM
#'   if (model %in% c("logistic", "linear")) {
#'     beta_est <- ridge_efficient(Y = Y, X = X, lambda = lambda_beta, model = model)
#'   } else if (model == "poisson") {
#'     result <- ridge_efficient(Y = Y, X = X, lambda = lambda_beta, model = model)
#'     beta_est <- result$penalized
#'   }
#'   
#'   return(beta_est / sqrt(sum(beta_est^2)))  # Return SCALED beta estimate
#' }




#' Generate Estimated Beta Coefficients for given X
#'
#' This function generates an estimate of the beta coefficients using ridge regression.
#' It first computes the linear predictor `eta`, simulates the response variable `Y` 
#' based on the given model, and estimates the beta coefficients using `glmnet`.
#'
#' @param true_beta Numeric vector, the true beta coefficients.
#' @param n_beta Integer, the number of samples used for estimation.
#' @param X design matrix of the penalized covariates
#' @param model Character, the type of model, e.g., "linear", "logistic", or "poisson".
#' @param lambda_beta Numeric, the regularization parameter for ridge regression.
#' @return A numeric vector of estimated beta coefficients (scaled to unit norm).
#' @export
generate_beta_est <- function(true_beta, n_beta, X, model, lambda_beta, scaled = TRUE) {
  # Generate response data
  eta <- calculate_eta(X, 0, true_beta, true_beta)  # Compute linear predictor
  Y <- generate_Y(eta, model, n_beta)  # Generate response variable
  
  # Use glmnet for ridge regression
  if (model == "logistic") {
    fit <- glmnet::glmnet(X, Y, 
                          family = "binomial", 
                          lambda = lambda_beta, 
                          alpha = 0, 
                          intercept = FALSE,
                          standardize = FALSE)
  } else if (model == "poisson") {
    fit <- glmnet::glmnet(X, Y, 
                          family = "poisson", 
                          lambda = lambda_beta, 
                          alpha = 0, 
                          intercept = FALSE,
                          standardize = FALSE)
  } else if (model == "linear") {
    fit <- glmnet::glmnet(X, Y, 
                          family = "gaussian", 
                          lambda = lambda_beta, 
                          alpha = 0, 
                          intercept = FALSE,
                          standardize = FALSE)
  } else {
    stop("Unsupported model type: ", model)
  }
  
  # Extract coefficients (excluding intercept)
  beta_est <- as.vector(coef(fit))[-1]
  
  # Return scaled beta estimate
  if(scaled){
    norm_beta <- sqrt(sum(beta_est^2))
    if (norm_beta > 0) {
      return(beta_est / norm_beta)
    } else {
      warning("Estimated beta has zero norm")
      return(beta_est)
    }
  }else {
    return(beta_est)
  }
  

}


# ---- Power Analysis Functions ----


#' Calculates power of signpost test - RESAMPLE (Updated)
#'
#' Updated version using ridge_efficient with glmnet scaling.
#'
#' @param n Positive integer, sample size.
#' @param p Positive integer, dimension (length of beta_0 and beta_a).
#' @param gammaH0 Numeric in [0,1], specifying the null hypothesis value.
#' @param model Character, specifying "linear", "logistic", or "poisson".
#' @param lambda Positive numeric, ridge penalty parameter (glmnet scale).
#' @param B_X Numeric, number of X's to be generated.
#' @param B_theta_null Numeric, size of the null distribution (per test).
#' @param B_power Numeric, number of p-values per gamma per generated X.
#' @param n_cores Optional integer, number of cores for parallel programming.
#' @param n_beta Positive integer, sample size to estimate beta_0 and beta_a.
#' @param lambda_beta Positive numeric, ridge penalty parameter when estimating beta_0 and beta_a (glmnet scale).
#' @param estimate_beta Logical, if TRUE beta_a are estimated; if FALSE, true values are used.
#'
#' @return A matrix of p-values with `B_X * B_power` rows and `length(gammas)` columns.
#' @export
power_signpost_test <- function(n, p, gammaH0, model, lambda, B_X, B_theta_null, B_power, n_cores = 8, n_beta = NULL, lambda_beta = NULL, estimate_beta = FALSE) {
  # Define gamma values
  gammas <- seq(0, 0.5, 0.1)
  num_gammas <- length(gammas)
  
  # Preallocate matrix for efficiency
  p_val_matrix <- matrix(NA, nrow = B_X * B_power, ncol = num_gammas)
  row_index <- 1  # Track row position
  
  progress <- create_progress_reporter(B_X, "Power calculation")
  
  for (i in seq_len(B_X)) {
    progress$update(i)
    
    # Define beta vectors
    BETA_A <- rnorm(p, sd = 0.06)
    BETA_0 <- rnorm(p, sd = 0.06)
    scl <- sqrt(sum((BETA_A - BETA_0)^2))
    BETA_A <- 2 * BETA_A / scl
    BETA_0 <- 2 * BETA_0 / scl
    
    if (estimate_beta) {
      if (is.null(n_beta) || is.null(lambda_beta)) {
        stop("n_beta and lambda_beta must be specified when estimate_beta = TRUE.")
      }
      X_a <- generate_X(n_beta, p)
      beta_a_hat <- generate_beta_est(BETA_A, n_beta, X_a, model, lambda_beta)
    } else {
      beta_a_hat <- BETA_A
    }
    
    X <- generate_X(n, p)
    eta <- calculate_eta(X, gammaH0, BETA_0, beta_a_hat) #for null distr the estimated betas used
    thetasH0 <- unlist(mclapply(1:B_theta_null, function(i) {
      if (estimate_beta) {
        beta_a_hat_k <- generate_beta_est(beta_a_hat, n_beta, X_a, model, lambda_beta)
      } else{
        beta_a_hat_k <- beta_a_hat
      }
      theta_Y(X, BETA_0, beta_a_hat_k, eta, model, lambda)
    }, mc.cores = n_cores))
    
    for (j in seq_len(num_gammas)) {
      eta <- calculate_eta(X, gammas[j], BETA_0, BETA_A) #Data is generated with true BETA_0, BETA_A
      theta_obs_gamma_j <- replicate(B_power, theta_Y(X, BETA_0, beta_a_hat, eta, model, lambda)) #for estimation the estimated beta used
      
      # Store p-values in the matrix
      p_val_matrix[row_index:(row_index + B_power - 1), j] <-
        sapply(theta_obs_gamma_j, calc_p_val, gammaH0 = gammaH0, thetasH0 = thetasH0)
    }
    
    row_index <- row_index + B_power  # Move to next set of rows
  }
  
  progress$close()
  
  return(p_val_matrix)
}

#' Calculates power of signpost test with rotation (Updated)
#'
#' Updated version using ridge_efficient with glmnet scaling.
#'
#' @param n Positive integer, sample size.
#' @param p Positive integer, dimension.
#' @param model Character, specifying one of "linear", "logistic", or "poisson".
#' @param B_power Positive integer, the number of p-values per gamma per generated design X.
#' @param B_theta_null Positive integer, the size of the null distribution (per test).
#' @param lambda Positive numeric, ridge penalty parameter (glmnet scale).
#' @param gammas Optional numeric vector of gamma values. Defaults to seq(0, 0.5, by = 0.1).
#' @param verbose Logical, if TRUE displays a progress bar. Default is FALSE.
#' @param seed Optional integer. If provided, sets the seed for reproducibility.
#' @param n_cores Optional integer, number of cores for parallel programming.
#' @param n_beta Positive integer, sample size to estimate beta_0 and beta_a.
#' @param lambda_beta Positive numeric, ridge penalty parameter when estimating beta_0 and beta_a (glmnet scale).
#' @param estimate_beta Logical, if TRUE beta_0 and beta_a are estimated; if FALSE, true values are used.
#'
#' @return A matrix of p-values with `B_power` rows and length(gammas) columns.
#' @export
power_ROTATION_signpost_test <- function(n, p, model, B_power, B_theta_null,
                                         lambda = Inf,
                                         gammas = seq(0, 0.5, by = 0.1),
                                         verbose = FALSE, seed = NULL, n_cores = 8, n_beta = NULL, lambda_beta = NULL,estimate_beta = FALSE) {
  # Input validation
  if (!is.numeric(n) || length(n) != 1 || n <= 0)
    stop("n must be a positive integer.")
  if (!is.numeric(p) || length(p) != 1 || p <= 0)
    stop("p must be a positive integer.")
  if (!is.numeric(B_power) || length(B_power) != 1 || B_power <= 0)
    stop("B_power must be a positive integer.")
  if (!is.numeric(B_theta_null) || length(B_theta_null) != 1 || B_theta_null <= 0)
    stop("B_theta_null must be a positive integer.")
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda <= 0)
    stop("lambda must be a positive number (use Inf for infinite lambda).")
  
  # Validate model: allowed models "linear", "logistic", "poisson"
  model <- match.arg(tolower(model), choices = c("linear", "logistic", "poisson"))
  
  # Set seed for reproducibility if provided
  if (!is.null(seed)) set.seed(seed)
  
  # Preallocate matrix for p-values
  num_gammas <- length(gammas)
  p_val_matrix <- matrix(NA, nrow = B_power, ncol = num_gammas)
  
  # Optionally set up a progress bar
  if (verbose) {
    progress <- create_progress_reporter(B_power, "Rotation power calculation")
  }
  
  # Outer loop over B_power iterations (each with a new design matrix X)
  for (i in seq_len(B_power)) {
    if (verbose) progress$update(i)
    # Define beta vectors
    BETA_A <- rnorm(p, sd = 0.06)
    BETA_0 <- rnorm(p, sd = 0.06)
    scl <- sqrt(sum((BETA_A - BETA_0)^2))
    BETA_A <- 2 * BETA_A / scl
    BETA_0 <- 2 * BETA_0 / scl
    
    X <- generate_X(n, p)
    if (estimate_beta) {
      if (is.null(n_beta) || is.null(lambda_beta)) {
        stop("n_beta and lambda_beta must be specified when estimate_beta = TRUE.")
      }
      beta_a_hat <- generate_beta_est(BETA_A, n_beta, p, model, lambda_beta)
      beta_0_hat <- generate_beta_est(BETA_0, n_beta, p, model, lambda_beta)
    } else {
      beta_a_hat <- BETA_A
      beta_0_hat <- BETA_0
    }
    
    # Loop over gamma values
    for (j in seq_along(gammas)) {
      eta <- calculate_eta(X, gammas[j], BETA_0, BETA_A)
      Y <- generate_Y(eta, model, n)
      
      if (is.infinite(lambda)) {
        # For infinite lambda, use theta_inf_hat
        theta_obs <- theta_inf_hat(Y, X, matrix(ncol = 0, nrow = n),
                                   beta_0_hat, beta_a_hat, model, do_rotation = FALSE)
        # Generate null distribution with rotation
        null_thetas <- unlist(mclapply(1:B_theta_null, function(i) {
          theta_inf_hat(Y, X, matrix(ncol = 0, nrow = n),
                        beta_0_hat, beta_a_hat, model, do_rotation = TRUE)
        }, mc.cores = n_cores))
      } else {
        # For finite lambda, use theta_hat_lambda
        theta_obs <- theta_hat_lambda(Y, X, beta_0_hat, beta_a_hat, lambda, model, do_rotation = FALSE)
        # Generate null distribution with rotation
        null_thetas <- unlist(mclapply(1:B_theta_null, function(i) {
          theta_hat_lambda(Y, X, beta_0_hat, beta_a_hat, lambda, model, do_rotation = TRUE)
        }, mc.cores = n_cores))
      }
      
      # Compute the p-value for gamma = gammas[j]
      p_val_matrix[i, j] <- calc_p_val(theta_obs, 0, null_thetas)
    }
  }
  
  if (verbose) progress$close()
  
  return(p_val_matrix)
}


#' Calculates raw p-values for the Likelihood Ratio test (Updated)
#'
#' Updated version using ridge_efficient with glmnet scaling.
#'
#' @param n Positive integer, sample size.
#' @param p Positive integer, dimension.
#' @param model Character specifying one of "linear", "logistic", or "poisson".
#' @param lambda Positive numeric, ridge penalty parameter (glmnet scale).
#' @param B_power Positive integer, the number of simulation iterations.
#' @param n_beta Positive integer, sample size to estimate beta_0 and beta_a.
#' @param lambda_beta Positive numeric, ridge penalty parameter when estimating beta_0 and beta_a (glmnet scale).
#' @param estimate_beta Logical, if TRUE beta_0 and beta_a are estimated; if FALSE, true values are used.
#'
#' @return A matrix of raw p-values with `B_power` rows and one column per gamma value.
#' @export
power_LR_test <- function(n, p, model, lambda, B_power, n_beta = NULL, lambda_beta = NULL,estimate_beta = FALSE) {
  # Input validation
  if (!is.numeric(n) || length(n) != 1 || n <= 0)
    stop("n must be a positive integer.")
  if (!is.numeric(p) || length(p) != 1 || p <= 0)
    stop("p must be a positive integer.")
  if (!model %in% c("linear", "logistic", "poisson"))
    stop("model must be one of 'linear', 'logistic', or 'poisson'.")
  if (!is.numeric(lambda) || length(lambda) != 1 || (lambda <= 0 && !is.infinite(lambda)))
    stop("lambda must be a positive number (or Inf).")
  if (!is.numeric(B_power) || length(B_power) != 1 || B_power <= 0)
    stop("B_power must be a positive integer.")
  
  # Define gamma values (0, 0.1, ..., 0.5)
  gammas <- seq(0, 0.5, by = 0.1)
  num_gammas <- length(gammas)
  
  # Preallocate a matrix to store raw p-values
  p_val_matrix <- matrix(NA, nrow = B_power, ncol = num_gammas)
  
  # Loop over simulation iterations
  for (i in seq_len(B_power)) {
    message("Iteration ", i, " of ", B_power)
    # Define beta vectors
    BETA_A <- rnorm(p, sd = 0.06)
    BETA_0 <- rnorm(p, sd = 0.06)
    scl <- sqrt(sum((BETA_A - BETA_0)^2))
    BETA_A <- 2 * BETA_A / scl
    BETA_0 <- 2 * BETA_0 / scl
    
    X <- generate_X(n, p)
    if (estimate_beta) {
      if (is.null(n_beta) || is.null(lambda_beta)) {
        stop("n_beta and lambda_beta must be specified when estimate_beta = TRUE.")
      }
      beta_a_hat <- generate_beta_est(BETA_A, n_beta, p, model, lambda_beta)
      beta_0_hat <- generate_beta_est(BETA_0, n_beta, p, model, lambda_beta)
    } else {
      beta_a_hat <- BETA_A
      beta_0_hat <- BETA_0
    }
    
    for (j in seq_along(gammas)) {
      eta <- calculate_eta(X, gammas[j], BETA_0, BETA_A)
      Y <- generate_Y(eta, model, n)
      
      # Compute the observed test statistic
      if (is.infinite(lambda)) {
        theta_obs <- theta_inf_hat(Y, X, matrix(ncol = 0, nrow = n),
                                   beta_0_hat, beta_a_hat, model)
      } else {
        theta_obs <- theta_hat_lambda(Y, X, beta_0_hat, beta_a_hat, lambda, model)
      }
      
      # Calculate the log-likelihoods at theta_obs and at 0
      l_hat <- loglik_theta(theta_obs, Y, X, beta_0_hat, beta_a_hat, lambda, model)
      l_0   <- loglik_theta(0, Y, X, beta_0_hat, beta_a_hat, lambda, model)
      
      # Compute the LR test p-value using the chi-squared approximation
      p_val_matrix[i, j] <- pchisq(-2 * (l_0 - l_hat), df = 1, lower.tail = FALSE)
    }
  }
  
  return(p_val_matrix)
}

#' Calculates power of signpost test with asymptotic standard errors (Updated)
#'
#' Updated version using ridge_efficient with glmnet scaling.
#'
#' @param n Positive integer, sample size.
#' @param p Positive integer, dimension.
#' @param gammaH0 Numeric in [0,1], specifying the null hypothesis value.
#' @param model Character, specifying "linear", "logistic", or "poisson".
#' @param lambda Positive numeric, ridge penalty parameter (glmnet scale).
#' @param B_X Numeric, number of X's to be generated.
#' @param B_power Numeric, number of p-values per gamma per generated X.
#' @param n_cores Optional integer, number of cores for parallel programming.
#' @param n_beta Positive integer, sample size to estimate beta_0 and beta_a.
#' @param lambda_beta Positive numeric, ridge penalty parameter when estimating beta_0 and beta_a (glmnet scale).
#' @param estimate_beta Logical, if TRUE beta_0 and beta_a are estimated; if FALSE, true values are used.
#' @param use_theta_hat_in_variance Logical, whether to use theta_hat in variance calculation.
#' @param fisher Logical, whether to use Fisher information for standard error.
#'
#' @return A matrix of p-values with `B_X * B_power` rows and `length(gammas)` columns.
#' @export
power_signpost_test_AS <- function(n, p, gammaH0, model, lambda, B_X, B_power, n_cores = 8, n_beta = NULL, lambda_beta = NULL, estimate_beta = FALSE, use_theta_hat_in_variance = FALSE, fisher = TRUE) {
  # Define beta vectors
  betas <- seq(-2.5, 2.5, length.out = p)
  betas <- betas / sqrt(sum(betas^2))
  BETA_A <- betas[p:1]
  BETA_0 <- betas
  
  # Define gamma values
  gammas <- seq(0, 0.5, 0.1)
  num_gammas <- length(gammas)
  
  # Preallocate matrix for efficiency
  p_val_matrix_AS <- matrix(NA, nrow = B_X * B_power, ncol = num_gammas)
  
  row_index <- 1  # Track row position
  
  progress <- create_progress_reporter(B_power, "AS power calculation")
  
  for (i in seq_len(B_X)) {
    progress$update(i)
    
    X <- generate_X(n, p)
    
    for (j in seq_len(num_gammas)) {
      eta_j <- calculate_eta(X, gammas[j], BETA_0, BETA_A) #Data is generated with true BETA_0, BETA_A
      for (k in seq_len(B_power)) {
        
        Yk <- generate_Y(eta_j, model, n)
        
        # Estimate theta_hat
        if (is.infinite(lambda)) {
          theta_hat <- theta_inf_hat(Yk, X, matrix(ncol = 0, nrow = n), BETA_0, BETA_A, model)
        } else {
          theta_hat <- theta_hat_lambda(Yk, X, BETA_0, BETA_A, lambda, model)
        }
        theta2_for_variance <- if (use_theta_hat_in_variance) theta_hat else 0
        
        stdef = sqrt(AS_stdef(BETA_0, BETA_A, theta2 = theta2_for_variance, y2 = Yk, X2 = X, fisher = fisher))
        #Asymptotic p-value
        row <- row_index + k - 1
        if (gammaH0 == 0) {
          # one‐sided right tail
          p_val_matrix_AS[row, j] <-
            pnorm(theta_hat, mean = 0, sd = stdef, lower.tail = FALSE)
        } else {
          # two‐sided
          p_val_matrix_AS[row, j] <-
            2 * pnorm(-abs(theta_hat - gammaH0) / stdef)
        }
      }
      
    }
    row_index <- row_index + B_power
  }
  
  progress$close()
  
  return(p_val_matrix_AS)
}

# ---- Extended Power Analysis Functions ----
#' Extended power analysis for signpost test with resampling (Updated)
#'
#' Updated version using ridge_efficient with glmnet scaling.
#'
#' @param n Positive integer, sample size.
#' @param p Positive integer, dimension.
#' @param gammaH0 Numeric in [0,1], specifying the null hypothesis value.
#' @param model Character, specifying "linear", "logistic", or "poisson".
#' @param lambda Positive numeric, ridge penalty parameter (glmnet scale).
#' @param B_X Numeric, number of X's to be generated.
#' @param B_theta_null Numeric, size of the null distribution (per test).
#' @param B_power Numeric, number of p-values per gamma per generated X.
#' @param n_cores Optional integer, number of cores for parallel programming.
#' @param n_beta Positive integer, sample size to estimate beta_0 and beta_a.
#' @param lambda_beta Positive numeric, ridge penalty parameter when estimating beta_0 and beta_a (glmnet scale).
#' @param estimate_beta Logical, if TRUE beta_a are estimated; if FALSE, true values are used.
#' @param lambda_est Numeric, regularization parameter for loss estimation (glmnet scale).
#'
#' @return List containing p_values, theta_hats, theta_hat_loss, theta_null_loss matrices.
#' @export
power_signpost_test_extended_resample <- function(n, p, gammaH0, model, lambda, B_X, B_theta_null, B_power, 
                                                  n_cores = 8, n_beta = NULL, lambda_beta = NULL, estimate_beta = FALSE, lambda_est = lambda) {
  # Define gamma values (grid)
  gammas <- seq(0, 0.5, 0.1)
  num_gammas <- length(gammas)
  
  # Preallocate matrices for each metric:
  p_val_matrix     <- matrix(NA, nrow = B_X * B_power, ncol = num_gammas)
  theta_hat_matrix <- matrix(NA, nrow = B_X * B_power, ncol = num_gammas)
  theta_hat_loss     <- matrix(NA, nrow = B_X * B_power, ncol = num_gammas)
  theta_null_loss     <- matrix(NA, nrow = B_X * B_power, ncol = num_gammas)
  
  row_index <- 1  # to keep track of where to store the replicates
  progress <- create_progress_reporter(B_X, "Extended resample power")
  
  for (i in seq_len(B_X)) {
    progress$update(i)
    
    # Define beta vectors for this iteration
    BETA_A <- rnorm(p, sd = 0.06)
    BETA_0 <- rnorm(p, sd = 0.06)
    scl <- sqrt(sum((BETA_A - BETA_0)^2))
    BETA_A <- 2 * BETA_A / scl
    BETA_0 <- 2 * BETA_0 / scl
    
    # If beta needs to be estimated, do so
    if (estimate_beta) {
      if (is.null(n_beta) || is.null(lambda_beta)) {
        stop("n_beta and lambda_beta must be specified when estimate_beta = TRUE.")
      }
      X_a <- generate_X(n_beta, p)
      beta_a_hat <- generate_beta_est(BETA_A, n_beta, X_a, model, lambda_beta)
    } else {
      beta_a_hat <- BETA_A
    }
    
    # Generate design matrix X and null distribution using gammaH0
    X <- generate_X(n, p)
    eta_null <- calculate_eta(X, gammaH0, BETA_0, beta_a_hat)
    thetasH0 <- unlist(mclapply(1:B_theta_null, function(idx) {
      if (estimate_beta) {
        beta_a_hat_k <- generate_beta_est(beta_a_hat, n_beta, X_a, model, lambda_beta)
      } else {
        beta_a_hat_k <- beta_a_hat
      }
      theta_Y(X, BETA_0, beta_a_hat_k, eta_null, model, lambda)
    }, mc.cores = n_cores))
    
    # Loop over gamma grid values
    for (j in seq_len(num_gammas)) {
      eta_j <- calculate_eta(X, gammas[j], BETA_0, BETA_A)  # using true BETA_0, BETA_A for data generation
      
      # For each gamma, run B_power replicates:
      for (k in 1:B_power) {
        # Generate response vector Y for current eta_j
        Yk <- generate_Y(eta_j, model, n)
        
        # Compute theta_hat (depending on lambda)
        if (is.infinite(lambda)) {
          theta_hat <- theta_inf_hat(Yk, X, matrix(ncol = 0, nrow = n), BETA_0, beta_a_hat, model)
        } else {
          theta_hat <- theta_hat_lambda(Yk, X, BETA_0, beta_a_hat, lambda, model)
        }
        
        # Calculate p-value based on the null distribution
        p_val <- calc_p_val(theta_hat, gammaH0, thetasH0)
        
        # Compute two loss metrics
        true_beta <- (1-gammas[j])*(BETA_0) + gammas[j]*BETA_A
        # Loss 1: squared loss with true beta using theta_hat
        loss1 <- squared_loss_beta(theta_hat, Yk, X, BETA_0, beta_a_hat, lambda_est, true_beta,model)
        # Loss 2: squared loss with true beta using theta = 0
        loss2 <- squared_loss_beta(0, Yk, X, BETA_0, beta_a_hat, lambda_est, true_beta,model)
        
        # Store results in the preallocated matrices
        current_row <- row_index + k - 1
        p_val_matrix[current_row, j]     <- p_val
        theta_hat_matrix[current_row, j] <- theta_hat
        theta_hat_loss[current_row, j]     <- loss1
        theta_null_loss[current_row, j]     <- loss2
      }
    }
    
    # Update row index for next set of replicates
    row_index <- row_index + B_power
  }
  
  progress$close()
  
  #Return a list of result matrices
  return(list(p_values   = p_val_matrix,
              theta_hats = theta_hat_matrix,
              theta_hat_loss      = theta_hat_loss,
              theta_null_loss      = theta_null_loss))
}


#' Extended power analysis for signpost test with asymptotic standard errors (Updated)
#'
#' Updated version using ridge_efficient with glmnet scaling.
#'
#' @param n Positive integer, sample size.
#' @param p Positive integer, dimension.
#' @param gammaH0 Numeric in [0,1], specifying the null hypothesis value.
#' @param model Character, specifying "linear", "logistic", or "poisson".
#' @param lambda Positive numeric, ridge penalty parameter (glmnet scale).
#' @param B_X Numeric, number of X's to be generated.
#' @param B_power Numeric, number of p-values per gamma per generated X.
#' @param n_cores Optional integer, number of cores for parallel programming.
#' @param n_beta Positive integer, sample size to estimate beta_0 and beta_a.
#' @param lambda_beta Positive numeric, ridge penalty parameter when estimating beta_0 and beta_a (glmnet scale).
#' @param estimate_beta Logical, if TRUE beta_a are estimated; if FALSE, true values are used.
#' @param use_theta_hat_in_variance Logical, whether to use theta_hat in variance calculation.
#' @param fisher Logical, whether to use Fisher information for standard error.
#' @param term1ONLY Logical, whether to use only first term in variance calculation.
#' @param lambda_est Numeric, regularization parameter for loss estimation (glmnet scale).
#' @param misspecified Logical, whether to use misspecified beta_0 (all zeros).
#' @param scaled Logical, whether to scale estimated beta_a.
#'
#' @return List containing p_values, theta_hats, theta_hat_loss, theta_null_loss matrices.
#' @export
power_signpost_test_extended_AS <- function(n, p, gammaH0, model, lambda,
                                            B_X, B_power, n_cores = 8,
                                            n_beta = NULL, lambda_beta = NULL,
                                            estimate_beta = FALSE,
                                            use_theta_hat_in_variance = FALSE,
                                            fisher = FALSE, term1ONLY = FALSE, lambda_est = lambda,
                                            misspecified = FALSE, scaled = TRUE) {
  # Define gamma values
  gammas <- seq(0, 0.5, 0.1)
  num_gammas <- length(gammas)
  
  # Preallocate matrices
  p_val_matrix_AS     <- matrix(NA, nrow = B_X * B_power, ncol = num_gammas)
  theta_hat_matrix_AS <- matrix(NA, nrow = B_X * B_power, ncol = num_gammas)
  theta_hat_loss_mat  <- matrix(NA, nrow = B_X *B_power, ncol = num_gammas)
  theta_null_loss_mat <- matrix(NA, nrow = B_X * B_power, ncol = num_gammas)
  
  row_index <- 1
  progress <- create_progress_reporter(B_X, "Extended AS power")
  on.exit(progress$close(), add = TRUE)
  
  for (i in seq_len(B_X)) {
    progress$update(i)
    
    # Generate random beta vectors and normalize
    BETA_A <- rnorm(p, sd = 0.06)
    BETA_0 <- rnorm(p, sd = 0.06)
    scl <- sqrt(sum((BETA_A - BETA_0)^2))
    BETA_A <- 2 * BETA_A / scl
    BETA_0 <- 2 * BETA_0 / scl
    
    if(misspecified){
      BETA_0_info = 0*BETA_0
    }else{
      BETA_0_info = BETA_0
    }
    
    # Optionally estimate beta_a
    if (estimate_beta) {
      if (is.null(n_beta) || is.null(lambda_beta)) {
        stop("n_beta and lambda_beta must be specified when estimate_beta = TRUE.")
      }
      X_a <- generate_X(n_beta, p)
      eta_a <- calculate_eta(X_a, 0, BETA_A, BETA_A)  # Compute linear predictor
      Y_a <- generate_Y(eta_a, model, n_beta)  # Generate response variable
      
      # Use ridge_efficient for estimation
      if (model %in% c("logistic", "linear")) {
        beta_a_hat <- ridge_efficient(Y = Y_a, X = X_a, lambda = lambda_beta, model = model)
      } else if (model == "poisson") {
        result <- ridge_efficient(Y = Y_a, X = X_a, lambda = lambda_beta, model = model)
        beta_a_hat <- result$penalized
      }
      
      if(scaled){
        beta_a_hat = beta_a_hat / sqrt(sum(beta_a_hat^2))
      }
    } else {
      beta_a_hat <- BETA_A
    }
    
    # Generate design matrix for test
    X <- generate_X(n, p)
    
    for (j in seq_len(num_gammas)) {
      # Generate data from true beta with gamma_j
      eta_j <- calculate_eta(X, gammas[j], BETA_0, BETA_A)
      
      for (k in seq_len(B_power)) {
        Yk <- generate_Y(eta_j, model, n)
        
        # Estimate theta_hat
        if (is.infinite(lambda)) {
          theta_hat <- theta_inf_hat(Yk, X, matrix(ncol = 0, nrow = n), BETA_0_info, beta_a_hat, model)
        } else {
          theta_hat <- theta_hat_lambda(Yk, X, BETA_0_info, beta_a_hat, lambda, model)
        }
        
        if(fisher){
          if(use_theta_hat_in_variance){
            diff_beta <- X %*% (beta_a_hat - BETA_0_info)
            lp0       <- ((1-theta_hat)*BETA_0_info + theta_hat *beta_a_hat) %*% t(X)
            ddg       <- diag(c(exp(-lp0) / (1 + exp(-lp0))^2))
            Fisher_info <- as.numeric(t(diff_beta) %*% ddg %*% diff_beta)
            stdef <- sqrt(1 / Fisher_info)
          }else{
            diff_beta <- X %*% (beta_a_hat - BETA_0_info)
            lp0       <- BETA_0_info %*% t(X)
            ddg       <- diag(c(exp(-lp0) / (1 + exp(-lp0))^2))
            Fisher_info <- as.numeric(t(diff_beta) %*% ddg %*% diff_beta)
            stdef <- sqrt(1 / Fisher_info)
          }
        }else{
          theta2_for_variance <- if (use_theta_hat_in_variance) theta_hat else 0
          v22 = calc_V22(beta_0 = BETA_0_info, beta_a = beta_a_hat, theta2 = theta2_for_variance, y1 = Y_a, y2 = Yk, X1= X_a, X2 = X, lambda_beta = lambda_beta, term1ONLY = term1ONLY)
          
          if (!is.finite(v22) || v22 <= 0) {
            warning("Invalid variance estimate; skipping")
            next
          }
          stdef <- sqrt(v22)
        }
        
        # Compute asymptotic p-value
        if (gammaH0 == 0) {
          p_val <- pnorm(theta_hat, mean = 0, sd = stdef, lower.tail = FALSE)
        } else {
          p_val <- 2 * pnorm(-abs(theta_hat - gammaH0) / stdef)
        }
        
        row <- row_index + k - 1
        p_val_matrix_AS[row, j]     <- p_val
        theta_hat_matrix_AS[row, j] <- theta_hat
        
        true_beta <- (1 - gammas[j]) * BETA_0 + gammas[j] * BETA_A
        loss_hat  <- squared_loss_beta(theta_hat, Yk, X, BETA_0_info, beta_a_hat,
                                       lambda_est, true_beta, model)
        loss_null <- squared_loss_beta(0,         Yk, X, BETA_0_info, beta_a_hat,
                                       lambda_est, true_beta, model)
        # store ----------------------------------------------------------------
        theta_hat_loss_mat [row, j] <- loss_hat
        theta_null_loss_mat[row, j] <- loss_null
      }
    }
    
    row_index <- row_index + B_power
  }
  
  return(list(
    p_values   = p_val_matrix_AS,
    theta_hats = theta_hat_matrix_AS,
    theta_hat_loss      = theta_hat_loss_mat,
    theta_null_loss      = theta_null_loss_mat
  ))
}

#' Extended power analysis for likelihood ratio test (Updated)
#'
#' Updated version using ridge_efficient with glmnet scaling.
#'
#' @param n Positive integer, sample size
#' @param p Positive integer, dimension
#' @param model Character specifying "linear", "logistic", or "poisson"
#' @param lambda Positive numeric, ridge penalty parameter (glmnet scale)
#' @param B_power Positive integer, number of simulation iterations
#' @param n_beta Positive integer, sample size for beta estimation
#' @param lambda_beta Positive numeric, penalty for beta estimation (glmnet scale)
#' @param estimate_beta Logical, whether to estimate beta values
#' @param lambda_est Numeric, regularization parameter for loss estimation (glmnet scale)
#' @param misspecified Logical, whether to use misspecified beta_0
#' @param scaled Logical, whether to scale estimated beta_a
#'
#' @return List containing p_values matrix
#' @export
power_LR_test_extended <- function(n, p, model, lambda,
                                   B_power,
                                   n_beta = NULL, lambda_beta = NULL,
                                   estimate_beta = FALSE, lambda_est = lambda, misspecified = FALSE, scaled = TRUE) {
  
  ## ---- input checks --------------------------------------------------------
  stopifnot(length(n)       == 1, n  > 0,
            length(p)       == 1, p  > 0,
            model %in% c("linear", "logistic", "poisson"),
            length(lambda)  == 1, (lambda > 0 | is.infinite(lambda)),
            length(B_power) == 1, B_power > 0)
  
  gammas <- seq(0, 0.5, by = 0.1)
  G      <- length(gammas)
  
  # wide result matrices ------------------------------------------------------
  p_val_matrix        <- matrix(NA, nrow = B_power, ncol = G)
  
  for (i in seq_len(B_power)) {
    message("Iteration ", i, "/", B_power)
    
    ## draw fresh (β₀, βₐ) for this replicate --------------------------------
    BETA_A <- rnorm(p, sd = 0.06)
    BETA_0 <- rnorm(p, sd = 0.06)
    scl     <- sqrt(sum((BETA_A - BETA_0)^2))
    BETA_A  <- 2 * BETA_A / scl
    BETA_0  <- 2 * BETA_0 / scl
    
    if(misspecified){
      BETA_0_info = 0*BETA_0
    }else{
      BETA_0_info = BETA_0
    }
    ## optionally re‑estimate ONLY βₐ ---------------------------------------
    if (estimate_beta) {
      if (is.null(n_beta) || is.null(lambda_beta))
        stop("n_beta and lambda_beta must be supplied when estimate_beta = TRUE.")
      X_a <- generate_X(n_beta, p)
      BETA_A_hat <- generate_beta_est(BETA_A, n_beta, X_a, model, lambda_beta, scaled = scaled)
    } else {
      BETA_A_hat <- BETA_A
    }
    # β₀ is **never** estimated
    
    ## design matrix for this replicate --------------------------------------
    X <- generate_X(n, p)
    
    for (j in seq_len(G)) {
      gamma_j <- gammas[j]
      eta_j   <- calculate_eta(X, gamma_j, BETA_0, BETA_A)
      Y       <- generate_Y(eta_j, model, n)
      
      # θ̂ under λ
      theta_hat <- if (is.infinite(lambda)) {
        theta_inf_hat(Y, X, matrix(ncol = 0, nrow = n),
                      BETA_0_info, BETA_A_hat, model)
      } else {
        theta_hat_lambda(Y, X, BETA_0_info, BETA_A_hat, lambda, model)
      }
      
      # LR p‑value -----------------------------------------------------------
      l_hat <- loglik_theta(theta_hat, Y, X, BETA_0_info, BETA_A_hat, lambda, model)
      l_0   <- loglik_theta(0,         Y, X, BETA_0_info, BETA_A_hat, lambda, model)
      p_val <- pchisq(-2 * (l_0 - l_hat), df = 1, lower.tail = FALSE)
      
      # store ----------------------------------------------------------------
      p_val_matrix       [i, j] <- p_val
    }
  }
  
  return(list(p_values        = p_val_matrix))
}


#' Extended power analysis for rotation test (Updated)
#'
#' Updated version using ridge_efficient with glmnet scaling.
#'
#' @param n Positive integer, sample size.
#' @param p Positive integer, dimension.
#' @param model Character, specifying one of "linear", "logistic", or "poisson".
#' @param B_power Positive integer, the number of p-values per gamma per generated design X.
#' @param B_theta_null Positive integer, the size of the null distribution (per test).
#' @param lambda Positive numeric, ridge penalty parameter (glmnet scale).
#' @param gammas Optional numeric vector of gamma values. Defaults to seq(0, 0.5, by = 0.1).
#' @param verbose Logical, if TRUE displays a progress bar. Default is TRUE.
#' @param seed Optional integer. If provided, sets the seed for reproducibility.
#' @param n_cores Optional integer, number of cores for parallel programming.
#' @param n_beta Positive integer, sample size to estimate beta_0 and beta_a.
#' @param lambda_beta Positive numeric, ridge penalty parameter when estimating beta_0 and beta_a (glmnet scale).
#' @param estimate_beta Logical, if TRUE beta_0 and beta_a are estimated; if FALSE, true values are used.
#' @param lambda_est Numeric, regularization parameter for loss estimation (glmnet scale).
#'
#' @return List containing p_values, theta_hats, theta_hat_loss, theta_null_loss matrices.
#' @export
power_signpost_test_extended_ROTATION <- function(n, p, model, B_power, B_theta_null,
                                                  lambda = Inf,
                                                  gammas = seq(0, 0.5, by = 0.1),
                                                  verbose = TRUE, seed = NULL, n_cores = 8, n_beta = NULL, lambda_beta = NULL, estimate_beta = FALSE, lambda_est = lambda) {
  if (!is.numeric(n) || length(n) != 1 || n <= 0)
    stop("n must be a positive integer.")
  if (!is.numeric(p) || length(p) != 1 || p <= 0)
    stop("p must be a positive integer.")
  if (!is.numeric(B_power) || length(B_power) != 1 || B_power <= 0)
    stop("B_power must be a positive integer.")
  if (!is.numeric(B_theta_null) || length(B_theta_null) != 1 || B_theta_null <= 0)
    stop("B_theta_null must be a positive integer.")
  if (!is.numeric(lambda) || length(lambda) != 1 || (lambda <= 0 && !is.infinite(lambda)))
    stop("lambda must be a positive number (use Inf for infinite lambda).")
  
  model <- match.arg(tolower(model), choices = c("linear", "logistic", "poisson"))
  if (!is.null(seed)) set.seed(seed)
  
  gammas <- seq(0, 0.5, by = 0.1)
  G <- length(gammas)
  p_val_matrix        <- matrix(NA, nrow = B_power, ncol = G)
  theta_hat_matrix    <- matrix(NA, nrow = B_power, ncol = G)
  theta_hat_loss_mat  <- matrix(NA, nrow = B_power, ncol = G)
  theta_null_loss_mat <- matrix(NA, nrow = B_power, ncol = G)
  
  if (verbose) progress <- create_progress_reporter(B_power, "Extended rotation power")
  
  for (i in seq_len(B_power)) {
    if (verbose) progress$update(i)
    
    BETA_A <- rnorm(p, sd = 0.06)
    BETA_0 <- rnorm(p, sd = 0.06)
    scl     <- sqrt(sum((BETA_A - BETA_0)^2))
    BETA_A  <- 2 * BETA_A / scl
    BETA_0  <- 2 * BETA_0 / scl
    
    if (estimate_beta) {
      if (is.null(n_beta) || is.null(lambda_beta)) {
        stop("n_beta and lambda_beta must be specified when estimate_beta = TRUE.")
      }
      X_beta     <- generate_X(n_beta, p)
      BETA_A_hat <- generate_beta_est(BETA_A, n_beta, X_beta, model, lambda_beta)
    } else {
      BETA_A_hat <- BETA_A
    }
    
    X <- generate_X(n, p)
    
    for (j in seq_len(G)) {
      gamma <- gammas[j]
      eta <- calculate_eta(X, gamma, BETA_0, BETA_A)
      Y <- generate_Y(eta, model, n)
      
      if (is.infinite(lambda)) {
        theta_obs <- theta_inf_hat(Y, X, matrix(ncol = 0, nrow = n),
                                   BETA_0, BETA_A_hat, model, do_rotation = FALSE)
        null_thetas <- unlist(mclapply(1:B_theta_null, function(i) {
          if (estimate_beta) {
            beta_a_hat_k <- generate_beta_est(BETA_A_hat, n_beta, X_beta, model, lambda_beta)
          } else {
            beta_a_hat_k <- BETA_A_hat
          }
          eta_null <- calculate_eta(X, 0, BETA_0, BETA_A_hat)
          Yk <- generate_Y(eta_null, model, n)
          theta_inf_hat(Yk, X, matrix(ncol = 0, nrow = n),
                        BETA_0, beta_a_hat_k, model, do_rotation = TRUE)
        }, mc.cores = n_cores))
      } else {
        theta_obs <- theta_hat_lambda(Y, X, BETA_0, BETA_A_hat, lambda, model, do_rotation = FALSE)
        null_thetas <- unlist(mclapply(1:B_theta_null, function(i) {
          if (estimate_beta) {
            beta_a_hat_k <- generate_beta_est(BETA_A_hat, n_beta, X_beta, model, lambda_beta)
          } else {
            beta_a_hat_k <- BETA_A_hat
          }
          eta_null <- calculate_eta(X, 0, BETA_0, BETA_A_hat)
          Yk <- generate_Y(eta_null, model, n)
          
          theta_hat_lambda(Yk, X, BETA_0, beta_a_hat_k, lambda, model, do_rotation = TRUE)
        }, mc.cores = n_cores))
      }
      
      p_val <- calc_p_val(theta_obs, 0, null_thetas)
      true_beta <- (1 - gamma) * BETA_0 + gamma * BETA_A
      loss_hat <- squared_loss_beta(theta_obs, Y, X, BETA_0, BETA_A_hat, lambda_est, true_beta, model)
      loss_null <- squared_loss_beta(0, Y, X, BETA_0, BETA_A_hat, lambda_est, true_beta, model)
      
      p_val_matrix[i, j]        <- p_val
      theta_hat_matrix[i, j]    <- theta_obs
      theta_hat_loss_mat[i, j]  <- loss_hat
      theta_null_loss_mat[i, j] <- loss_null
    }
  }
  
  if (verbose) progress$close()
  
  list(
    p_values        = p_val_matrix,
    theta_hats      = theta_hat_matrix,
    theta_hat_loss  = theta_hat_loss_mat,
    theta_null_loss = theta_null_loss_mat
  )
}




