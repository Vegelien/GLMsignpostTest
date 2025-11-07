
# ============================================================================
#                        SIGNPOST TEST FUNCTIONS
# ============================================================================

# ---- Utility Functions ----

#' Validate common inputs for signpost test functions
#' @keywords internal
validate_signpost_inputs <- function(Y, X, beta_0, beta_a, model) {
  if (length(Y) != nrow(X)) {
    stop("Length of Y must equal number of rows in X")
  }
  if (length(beta_0) != ncol(X) || length(beta_a) != ncol(X)) {
    stop("Length of beta_0 and beta_a must equal number of columns in X")
  }
  if (!model %in% c("linear", "logistic", "poisson")) {
    stop("Model must be one of: 'linear', 'logistic', 'poisson'")
  }
}

#' Create a progress reporter
#' @keywords internal
create_progress_reporter <- function(total, name = "Progress") {
  if (interactive()) {
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    list(
      update = function(i) setTxtProgressBar(pb, i),
      close = function() close(pb)
    )
  } else {
    list(
      update = function(i) if (i %% 10 == 0) message(name, ": ", i, "/", total),
      close = function() invisible(NULL)
    )
  }
}

#' Resolve GLM family utilities for variance calculations
#'
#' Returns link-inverse and variance functions for supported models.
#'
#' @param model Character string identifying the model ("logistic" or "poisson").
#' @keywords internal
resolve_glm_family <- function(model) {
  switch(
    model,
    logistic = list(
      linkinv = function(eta) 1 / (1 + exp(-eta)),
      variance = function(mu) mu * (1 - mu)
    ),
    linear = list(
      linkinv = function(eta) eta,
      variance = function(mu) rep(1, length(mu))
    ),
    poisson = list(
      linkinv = function(eta) exp(eta),
      variance = function(mu) mu
    ),
    stop(sprintf("Unsupported model '%s'", model))
  )
}

#' Compute GLM mean/variance given penalized and unpenalized blocks
#'
#' @param X Numeric matrix for penalized coefficients.
#' @param beta Numeric vector of coefficients for X.
#' @param model Character string ("logistic" or "poisson").
#' @param U Optional numeric matrix for unpenalized covariates.
#' @param delta Optional numeric vector of coefficients for U.
#' @param offset Optional numeric vector to add to the linear predictor.
#' @param mu Optional pre-computed mean vector.
#' @param variance Optional pre-computed variance vector.
#' @keywords internal
compute_eta_mu_variance <- function(X, beta, model, U = NULL, delta = NULL,
                                    offset = NULL, mu = NULL, variance = NULL) {
  family_utils <- resolve_glm_family(model)
  
  if (!is.null(mu)) {
    if (is.null(variance)) {
      variance <- family_utils$variance(mu)
    }
    return(list(eta = NULL, mu = mu, variance = variance))
  }
  
  if (is.null(beta)) {
    stop("'beta' must be provided when 'mu' is NULL")
  }
  
  eta <- drop(X %*% beta)
  if (!is.null(offset)) {
    eta <- eta + offset
  }
  if (!is.null(U)) {
    if (is.null(delta)) {
      stop("'delta' must be supplied when U is provided and mean is not pre-computed")
    }
    eta <- eta + drop(U %*% delta)
  }
  
  mu <- family_utils$linkinv(eta)
  variance <- family_utils$variance(mu)
  
  list(eta = eta, mu = mu, variance = variance)
}

#' Refit unpenalized coefficients for a fixed theta
#'
#' @param Y Response vector.
#' @param U Matrix of unpenalized covariates.
#' @param offset Linear predictor offset (typically X %*% beta(theta)).
#' @param model Model family ("logistic" or "poisson").
#' @keywords internal
estimate_unpenalized_delta <- function(Y, U, offset, model) {
  if (is.null(U) || ncol(U) == 0) {
    return(numeric(0))
  }
  
  family <- switch(
    model,
    logistic = stats::binomial(),
    linear = stats::gaussian(),
    poisson = stats::poisson(),
    stop(sprintf("Unsupported model '%s'", model))
  )
  
  fit <- try(stats::glm.fit(x = U, y = Y, family = family, offset = offset), silent = TRUE)
  if (inherits(fit, "try-error") || is.null(fit$coefficients)) {
    return(rep(0, ncol(U)))
  }
  
  delta_hat <- as.numeric(fit$coefficients)
  delta_hat[!is.finite(delta_hat)] <- 0
  delta_hat
}

# ---- Core Theta Calculation Functions ----


# thetaInf <- function(theta, Y, X, U, lpTarget0, lpTargetD, model) {
#   linPred <- lpTarget0 + lpTargetD * theta
#   
#   if (ncol(U) > 0) {
#     # CLIP just for the ridge step to avoid exp()/weights overflow in IRLS
#     linPred_clip <- pmax(pmin(linPred, 35), -35)
#     
#     adj <- ridge_efficient(
#       Y = Y,
#       X = matrix(linPred_clip, ncol = 1),
#       U = U,
#       lambda = 1e10,          # keep as-is for now (we can reduce after verifying)
#       target = c(1),
#       model = model
#     )[1:ncol(U)]
#     
#     # if adj has any non-finite, skip the adjustment (rare with clip)
#     if (!all(is.finite(adj))) adj[] <- 0
#     
#     linPred <- linPred + as.vector(U %*% adj)
#   }
#   
#   # ... (rest identical)
# }


#' Helper function to theta_inf_hat, calculates RHS of estimating equation
#'
#' @param theta Numeric in [0,1], shrinkage parameter.
#' @param Y Numeric vector, response variable.
#' @param X Numeric matrix, design matrix of penalized covariates. Rows should match length of Y.
#' @param U Numeric matrix, design matrix of unpenalized covariates. Rows should match length of Y.
#' @param lpTarget0 Precomputed numeric vector, linear predictor for beta_0.
#' @param lpTargetD Precomputed numeric vector, linear predictor difference for beta_a - beta_0.
#' @param model String, one of "linear", "logistic", or "poisson".
#' @importFrom porridge ridgeGLM
#' @return Numeric, right-hand side of estimating equation.
thetaInf <- function(theta, Y, X, U, lpTarget0, lpTargetD, model) {
  linPred <- lpTarget0 + lpTargetD * theta
  
  if (ncol(U) > 0) {
    linPred_clip <- pmax(pmin(linPred, 35), -35)
    dHat <- ridge_efficient(Y = Y, X = matrix(linPred_clip, ncol = 1), U = U, lambda = 10^10, target = c(1), model = model)[1:ncol(U)]
    linPred <- linPred + dHat %*% t(U)

    
    # fit_unpen <- ridge_efficient(
    #   Y = Y,
    #   X = matrix(linPred, ncol = 1),
    #   U = U,
    #   lambda = 10^10,
    #   target = c(1),
    #   model = model
    # )
    # 
    # dHat <- extract_unpenalized_coefficients(fit_unpen, ncol(U))
    # linPred <- linPred + tcrossprod(dHat, U)
    # 
  }
  
  if (model == "linear") {
    Yexp <- linPred
  } else if (model == "logistic") {
    Yexp <- 1 / (1 + exp(-linPred))
  } else if (model == "poisson") {
    Yexp <- exp(linPred)
  }
  
  return(sum((Y - Yexp) * lpTargetD))
}

#' Extract unpenalized coefficients from a ridge solver result
#'
#' @param fit Result returned by `ridge_efficient()`.
#' @param n_unpenalized Integer, number of unpenalized coefficients expected.
#' @keywords internal
extract_unpenalized_coefficients <- function(fit, n_unpenalized) {
  if (n_unpenalized <= 0L) {
    return(numeric(0))
  }
  
  if (is.null(fit)) {
    return(rep(0, n_unpenalized))
  }
  
  coeffs <- NULL
  
  if (is.list(fit)) {
    if (!is.null(fit$unpenalized)) {
      coeffs <- fit$unpenalized
    } else {
      alpha <- if (!is.null(fit$alpha)) fit$alpha else numeric(0)
      gamma <- if (!is.null(fit$gamma)) fit$gamma else numeric(0)
      coeffs <- c(alpha, gamma)
      
      if (length(coeffs) == 0 && !is.null(fit$coef)) {
        coeffs <- fit$coef
      }
      
      if (length(coeffs) == 0) {
        coeffs <- unlist(fit, use.names = FALSE)
      }
    }
  } else {
    coeffs <- as.numeric(fit)
  }
  
  coeffs <- as.numeric(coeffs)
  if (!length(coeffs)) {
    coeffs <- numeric(0)
  }
  
  if (length(coeffs) < n_unpenalized) {
    coeffs <- c(coeffs, rep(0, n_unpenalized - length(coeffs)))
  }
  
  coeffs <- coeffs[seq_len(n_unpenalized)]
  coeffs[!is.finite(coeffs)] <- 0
  
  coeffs
}


#' Generalized theta_inf_hat with custom endpoints
#'
#' Computes theta_hat for lambda = infinity with custom interval endpoints.
#' Allows testing hypotheses beyond the standard [0,1] interval.
#'
#' @param Y Numeric vector, response variable.
#' @param X Numeric matrix, design matrix of penalized covariates. Rows should match the length of Y.
#' @param U Numeric matrix, design matrix of unpenalized covariates. Rows should match the length of Y.
#' @param beta_0 Numeric vector, shrinkage target (same length as columns of X).
#' @param beta_a Numeric vector, alternative hypothesis shrinkage target (same length as columns of X).
#' @param model String, one of "linear", "logistic", or "poisson".
#' @param theta_min Numeric, lower bound for theta search (default: 0).
#' @param theta_max Numeric, upper bound for theta search (default: 1).
#' @param do_rotation Logical, whether to apply random rotation.
#'
#' @return Numeric, the estimated theta_hat within [theta_min, theta_max].
#' @importFrom pracma randortho
#' @export
theta_inf_hat <- function(Y, X, U, beta_0, beta_a, model, 
                          theta_min = 0, theta_max = 1, do_rotation = FALSE) {
  
  lpTarget0 <- beta_0 %*% t(X)
  
  if (do_rotation) {
    Rp <- randortho(length(beta_a), type = "orthonormal")
    lpTargetD <- as.numeric(Rp %*% as.numeric(beta_a - beta_0))
    lpTargetD <- tcrossprod(lpTargetD, X)
  } else {
    lpTargetD <- tcrossprod(beta_a - beta_0, X)
  }
  
  # Evaluate at endpoints
  thetaIlow <- thetaInf(theta_min, Y, X, U, lpTarget0, lpTargetD, model)
  thetaIup <- thetaInf(theta_max, Y, X, U, lpTarget0, lpTargetD, model)
  
  # Check if root exists in interval
  if (thetaIup * thetaIlow < 0) {
    thetaI <- uniroot(thetaInf, c(theta_min, theta_max), 
                      Y = Y, X = X, U = U, 
                      lpTarget0 = lpTarget0, lpTargetD = lpTargetD, model = model)$root
  } else {
    # Choose endpoint with smaller absolute value
    print("no root in interval")
    thetaI <- ifelse(abs(thetaIup) > abs(thetaIlow), theta_min, theta_max)
  }
  
  return(thetaI)
}

#' Direct θ̂∞ when U ≠ ∅ via re-parameterization (robust intercept handling)
#'
#' Estimates theta for the λ = ∞ case by fixing Xβ0 at coefficient 1
#' (via a huge ridge penalty with target=1) and reading off θ as the
#' coefficient on X(βa−β0) included in the unpenalized block.
#'
#' Robustness features:
#' - Drops constant (≈zero-variance) columns from U to avoid backend dropping.
#' - If a constant column (intercept proxy) is removed, the solver's intercept
#'   is enabled so we still fit exactly one intercept.
#' - Optionally bumps lambda once if coef(base) deviates noticeably from 1.
#'
#' @param Y Numeric vector (length n)
#' @param X Numeric matrix (n × p) of penalized covariates
#' @param U Numeric matrix (n × q) of unpenalized covariates (q can be 0)
#' @param beta_0 Numeric length-p shrinkage target
#' @param beta_a Numeric length-p alternative target
#' @param model "linear", "logistic", or "poisson"
#' @param theta_min,theta_max numeric bounds used for optional clipping
#' @param lambda_fix large ridge penalty to nail coef(base) ≈ 1
#' @param intercept logical; if TRUE, solver adds an intercept (may be flipped on
#'                  internally when U contains/loses a constant column)
#' @param standardize logical; FALSE means target=1 is on raw scale
#' @param check_boundaries logical; clip to [theta_min, theta_max]
#' @param tol_var tiny variance guard for identifiability of diff
#' @param const_tol tolerance for "constant column" detection in U
#' @param bump_lambda_tol threshold for |coef(base) - 1| to trigger one bump
#' @param bump_factor how much to multiply lambda_fix by when bumping
#' @param ...
#' @return numeric scalar θ̂∞
#' @export
theta_inf_hat_direct <- function(
    Y, X, U, beta_0, beta_a,
    model = c("linear","logistic","poisson"),
    method = c("glmoffset", "ridge"),
    theta_min = 0, theta_max = 1,
    lambda_fix = 1e10,
    intercept = TRUE,
    standardize = FALSE,
    check_boundaries = TRUE,
    tol_var = 1e-10,
    const_tol = 1e-12,
    bump_lambda_tol = 1e-6,
    bump_factor = 100,
    ...
) {
  model <- match.arg(model)
  method <- match.arg(method)
  validate_signpost_inputs(Y, X, beta_0, beta_a, model)
  
  n <- length(Y)
  if (is.null(U)) U <- matrix(nrow = n, ncol = 0)
  U <- as.matrix(U)
  
  # Fallback to root-based when U is empty.
  if (ncol(U) == 0) {
    return(theta_inf_hat(Y, X, U, beta_0, beta_a, model,
                         theta_min = theta_min, theta_max = theta_max,
                         do_rotation = FALSE))
  }
  
  # base = X beta0 ; diff = X (beta_a - beta0)
  base <- drop(X %*% beta_0)
  diff <- drop(X %*% (beta_a - beta_0))
  
  # Identifiability guard: if diff has ~0 variance, choose boundary matching the smaller score
  if (stats::sd(diff) < tol_var) {
    if (!check_boundaries) return(theta_min)
    lpTarget0 <- base
    lpTargetD <- diff
    g_low <- thetaInf(theta_min, Y, X, U, lpTarget0, lpTargetD, model)
    g_up  <- thetaInf(theta_max, Y, X, U, lpTarget0, lpTargetD, model)
    return(if (abs(g_up) > abs(g_low)) theta_min else theta_max)
  }
  
  # === Intercept & constant-column canonicalization ===
  # Drop constant columns from U (prevents backend "Dropped zero-variance" surprises)
  if (ncol(U) > 0) {
    sds <- apply(U, 2L, function(col) stats::sd(col))
    const_idx <- which(!is.finite(sds) | sds < const_tol)
    if (length(const_idx)) {
      # If we removed a constant column (very likely an intercept proxy) and caller asked
      # for no intercept, we still enable solver intercept to keep exactly one intercept.
      if (!intercept) intercept <- TRUE
      U <- U[, setdiff(seq_len(ncol(U)), const_idx), drop = FALSE]
    }
  }
  
  if (method == "glmoffset") {
    if (standardize) {
      stop("standardize = TRUE is not supported for method = 'glmoffset'")
    }

    family <- switch(
      model,
      linear = stats::gaussian(),
      logistic = stats::binomial(),
      poisson = stats::poisson(),
      stop("Unsupported model: ", model)
    )

    df <- data.frame(Y = Y)
    if (ncol(U) > 0) {
      U_df <- as.data.frame(U)
      if (is.null(colnames(U_df))) {
        colnames(U_df) <- paste0("U", seq_len(ncol(U_df)))
      }
      df <- cbind(df, U_df)
      u_names <- colnames(U_df)
    } else {
      u_names <- character(0)
    }
    df$diff <- diff
    df$base <- base

    rhs_terms <- c(u_names, "diff")
    if (length(rhs_terms) == 0) {
      formula_str <- if (intercept) "Y ~ 1" else "Y ~ 0"
    } else {
      rhs <- paste(rhs_terms, collapse = " + ")
      formula_str <- if (intercept) paste("Y ~", rhs) else paste("Y ~ 0 +", rhs)
    }

    fit <- stats::glm(
      stats::as.formula(formula_str),
      family = family,
      data = df,
      offset = df$base,
      ...
    )

    coef_diff <- stats::coef(fit)["diff"]
    if (is.na(coef_diff)) {
      stop("Failed to estimate coefficient for 'diff' in glm fit")
    }
    theta_hat <- unname(coef_diff)
  } else {
    # One-shot fit: penalized block is base (target=1, huge lambda),
    # unpenalized block is [U, diff]; θ is the last gamma.
    do_fit <- function(lam) {
      ridge_complete(
        Y = Y,
        X = matrix(base, ncol = 1),
        lambda = lam,
        target = 1,
        model = if (model == "linear") "gaussian" else model,
        intercept = intercept,
        U_extra = cbind(U, diff),
        standardize = standardize,
        ...
      )
    }

    fit <- do_fit(lambda_fix)

    # If coef(base) deviates meaningfully from 1, bump lambda once.
    if (!is.null(fit$beta) && length(fit$beta) == 1L && is.finite(fit$beta[1L])) {
      if (abs(fit$beta[1L] - 1) > bump_lambda_tol) {
        fit <- do_fit(lambda_fix * bump_factor)
      }
    }

    theta_hat <- utils::tail(fit$gamma, 1L)
  }
  
  if (check_boundaries) {
    if (theta_hat < theta_min) theta_hat <- theta_min
    if (theta_hat > theta_max) theta_hat <- theta_max
  }
  theta_hat
}



#' Calculates theta_hat for given lambda < Inf with one-time random rotation and optimization
#'
#' @param Y Numeric, response vector.
#' @param X Design matrix of the penalized covariates; the number of rows should match the length of Y.
#' @param beta_0 Numeric, shrinkage target; length must match the number of columns of X.
#' @param beta_a Numeric, alternative target; length must match the number of columns of X.
#' @param lambda Numeric, penalty parameter (glmnet scale).
#' @param model Character, one of "linear", "logistic", or "poisson". Defaults to "logistic".
#' @param do_rotation Logical, default FALSE. If TRUE, applies a random orthonormal rotation to (beta_a - beta_0) once.
#' @param grid_length Optional integer; if provided, a grid search with that many points will be performed instead of using optimize().
#'
#' @return theta_hat Numeric in [0,1], the estimated optimal theta.
#' @importFrom pracma randortho
#' @export
theta_hat_lambda <- function(Y, X, beta_0, beta_a, lambda,
                             model = c("logistic", "poisson", "linear"),
                             do_rotation = FALSE, grid_length = NULL) {
  model <- match.arg(model)
  
  # Validate dimensions
  validate_signpost_inputs(Y, X, beta_0, beta_a, model)
  if (!is.finite(lambda) && lambda < 0) {
    stop("lambda must be positive (or Inf).")
  }
  
  # Compute rotated beta_a if requested
  if (do_rotation) {
    Rp <- randortho(length(beta_a), type = "orthonormal")
    beta_a <- beta_0 + as.numeric(Rp %*% as.numeric(beta_a - beta_0))
  } 
  
  # Define objective: negative log-likelihood (for minimization)
  obj_fun <- function(theta) {
    -loglik_theta(theta, Y, X, beta_0, beta_a, lambda, model)
  }
  
  # Use either a grid search or optimize() based on grid_length
  if (is.null(grid_length)) {
    opt_res <- optimize(obj_fun, interval = c(0, 1))
    thetaObs <- opt_res$minimum
  } else {
    theta_seq <- seq(0, 1, length.out = grid_length)
    vals <- sapply(theta_seq, obj_fun)
    thetaObs <- theta_seq[which.min(vals)]
  }
  
  return(thetaObs)
}

#' Generates Y and calculates theta_hat for given lambda (possibly inf)
#'
#' @param X design matrix of the penalized covariates, The number of rows should
#'               match the number of elements of Y.
#' @param beta_0 numeric, towards which the estimate is shrunken. The number of
#'               elements should match the number of columns of X.
#' @param beta_a numeric, towards which the estimate is shrunken. The number of
#'               elements should match the number of columns of X.
#' @param eta numeric linear predictor with same number of elements as rows in X
#' @param model string specifying "linear", "logistic" or "poisson"
#' @param lambda positive numeric, the regular ridge penalty parameter
#'
#' @return theta_hat numeric in [0,1]
#' @export
theta_Y = function(X, beta_0, beta_a, eta, model, lambda, theta_min = 0, theta_max = 1){
  n = length(eta)
  
  Yk = generate_Y(eta,model,n)
  
  if(is.infinite(lambda)){
    return(theta_inf_hat(Yk, X, matrix(ncol=0, nrow=n), beta_0, beta_a, model, theta_min = theta_min, theta_max = theta_max))
  }else{
    return( theta_hat_lambda(Yk,X,beta_0, beta_a, lambda, model))
  }
  
}

# ---- Likelihood and Penalty Functions ----

#' Calculate loglikelihood for a given theta using ridge_efficient
#' 
#' Updated version that uses ridge_efficient with glmnet scaling convention.
#' 
#' @param theta Numeric in [0,1], weight parameter for the targets beta_0 and beta_a.
#' @param Y Numeric, response vector.
#' @param X Design matrix of the penalized covariates; the number of rows should match the length of Y.
#' @param beta_0 Numeric, shrinkage target; length must match the number of columns of X.
#' @param beta_a Numeric, alternative target; if rotated, pass the rotated version.
#' @param lambda Numeric, penalty parameter (glmnet scale).
#' @param model Character, one of "linear", "logistic", or "poisson". Defaults to "logistic".
#'
#' @return Numeric value of the loglikelihood.
#' @export
loglik_theta <- function(theta, Y, X, beta_0, beta_a, lambda,
                         model = c("logistic", "poisson", "linear")) {
  model <- match.arg(model)
  n <- length(Y)
  
  # Compute the target coefficients based on theta
  target <- (1 - theta) * beta_0 + theta * beta_a
  
  # Compute the penalized estimate using ridge_efficient
  if (model %in% c("logistic", "linear")) {
    beta_hat <- ridge_efficient(Y = Y, X = X, lambda = lambda, target = target, model = model)
  } else if (model == "poisson") {
    result <- ridge_efficient(Y = Y, X = X, lambda = lambda, target = target, model = model)
    beta_hat <- result$penalized
  }
  
  # Compute the predicted response and corresponding log-likelihood
  if (model == "logistic") {
    mu_hat <- plogis(X %*% beta_hat)
    loglik_hat <- sum(dbinom(Y, size = 1, prob = mu_hat, log = TRUE))
  } else if (model == "poisson") {
    eta_hat <- X %*% beta_hat
    mu_hat <- exp(eta_hat)
    loglik_hat <- sum(dpois(Y, lambda = mu_hat, log = TRUE))
  } else if (model == "linear") {
    residuals <- Y - X %*% beta_hat
    loglik_hat <- -sum(residuals^2)
  }
  
  return(loglik_hat)
}

#' Calculates penalty
#'
#' @param theta numeric, in [0,1] weight parameter of targets beta_0 and beta_a
#' @param Y numeric, response vector
#' @param X design matrix of the penalized covariates, The number of rows should 
#'               match the number of elements of Y.
#' @param beta_0 numeric, towards which the estimate is shrunken. The number of  
#'               elements should match the number of columns of X.
#' @param beta_a numeric, towards which the estimate is shrunken. The number of  
#'               elements should match the number of columns of X.
#' @param lambda positive numeric, the regular ridge penalty parameter
#' @param model string specifying "linear", "logistic" or "poisson"
#'
#' @return penalty numeric
#' @importFrom porridge ridgeGLM
#' @export
penalty_theta <- function(theta,Y,X,beta_0, beta_a, lambda, model = "logistic"){
  print("penalty_theta uses ridgeGLM")
  target = (1-theta)*beta_0 + theta*beta_a
  
  if(model == "logistic"){
    beta_hat = ridgeGLM(Y=Y, X=X,  lambda=lambda, target=target, model="logistic")
  }else if(model == "poisson"){
    beta_hat = AridgePLM(Y=Y, X=X,  lambda=lambda, target=target, model="poisson")$penalized
  }
  
  penalty = 
  
  return(penalty)
}

#' Calculate Squared Loss Between Estimated and True Beta (Updated)
#'
#' Updated version using ridge_efficient.
#'
#' @param theta Numeric in [0,1], weight parameter for convex combination of targets
#' @param Y Numeric vector, response variable
#' @param X Numeric matrix, design matrix
#' @param beta_0 Numeric vector, null hypothesis target
#' @param beta_a Numeric vector, alternative hypothesis target
#' @param lambda Numeric, ridge penalty parameter (glmnet scale)
#' @param true_beta Numeric vector, true coefficients for loss calculation
#' @param model Character, one of "linear", "logistic", or "poisson"
#'
#' @return Numeric, squared loss between estimated and true beta
#' @export
squared_loss_beta <- function(theta, Y, X, beta_0, beta_a, lambda, true_beta, model = "logistic") {
  target <- (1 - theta) * beta_0 + theta * beta_a
  
  if (model %in% c("logistic", "linear")) {
    beta_hat <- ridge_efficient(Y = Y, X = X, lambda = lambda, target = target, model = model)
  } else if (model == "poisson") {
    result <- ridge_efficient(Y = Y, X = X, lambda = lambda, target = target, model = model)
    beta_hat <- result$penalized
  } else {
    stop("Model must be 'logistic', 'poisson', or 'linear'")
  }
  
  sum((beta_hat - true_beta)^2)
}

# ---- P-value Calculation ----

#' Calculates p value for null distr and observed value
#'
#' @param thetaObs numeric in [0,1] test statistic
#' @param gammaH0 numeric in [0,1] the true value of gamma under H0
#' @param thetasH0 numeric in [0,1] the generated thetas under H0
#'
#' @return p value numeric in [0,1]
#' @export
calc_p_val = function(thetaObs,gammaH0,  thetasH0){
  nRot = length(thetasH0)
  
  if(gammaH0 == 0){
    pvalR  = sum((c(thetaObs, thetasH0) >= thetaObs))/ (nRot + 1)
  }else if(gammaH0 ==1){
    pvalR  = sum((c(thetaObs, thetasH0) <= thetaObs))/ (nRot + 1)
  }else if(thetaObs > gammaH0){
    pvalR  = 2* sum((c(thetaObs, thetasH0) >= thetaObs))/ (nRot + 1)
  }else{
    pvalR  = 2* sum((c(thetaObs, thetasH0) <= thetaObs))/ (nRot + 1)
  }
  return(pvalR)
}

#' Maps array of p_values to the power
#'
#' @param pval Array of p_values (numeric)
#' @param alpha numeric, default = 0.05
#'
#' @return pval
#' @export
pval2power = function(pval, alpha = 0.05){
  ## ---------------------------------------------------------------------
  ## Maps array of p_values to the power
  ## ---------------------------------------------------------------------
  ## Arguments
  ## pval        : Array of p_values (numeric)
  ## alpha       : numeric, default = 0.05
  ## ---------------------------------------------------------------------
  
  return(sum(pval < alpha)/length(pval))
}

# ---- Asymptotic Variance Calculation ----

#' Calculate Asymptotic Standard Error for Signpost Test
#'
#' Computes the asymptotic variance for the signpost test statistic theta_hat.
#' Can compute either the sandwich estimator or Fisher information based variance.
#'
#' @param beta_0 Numeric vector, null hypothesis coefficients (length p)
#' @param beta_a Numeric vector, alternative hypothesis coefficients (length p)
#' @param theta2 Numeric scalar, value of theta for variance calculation
#' @param y1 Numeric vector, response from beta estimation (length n_beta)
#' @param y2 Numeric vector, response from test data (length n)
#' @param X1 Numeric matrix, design matrix for beta estimation (n_beta x p)
#' @param X2 Numeric matrix, design matrix for test data (n x p)
#' @param lambda_beta Numeric scalar, ridge penalty for beta estimation
#' @param variance_type Character, either "sandwich" or "fisher"
#' @param include_estimation_uncertainty Logical, whether to include uncertainty from beta estimation
#' @param model Character string indicating the GLM family ("logistic" or "poisson")
#' @param U1 Optional numeric matrix of unpenalized covariates used in beta estimation
#' @param U2 Optional numeric matrix of unpenalized covariates for the test data
#' @param delta Optional numeric vector of unpenalized coefficients for the test data
#' @param delta_est Optional numeric vector of unpenalized coefficients for the estimation data
#' @param offset1 Optional numeric vector offset for the estimation data linear predictor
#' @param offset2 Optional numeric vector offset for the test data linear predictor
#' @param mu2 Optional numeric vector of pre-computed means for the test data
#' @param variance2 Optional numeric vector of pre-computed variances for the test data
#' @return Numeric scalar, the variance estimate
#' @export
calc_asymptotic_variance <- function(beta_0, beta_a, theta2, y1 = NULL, y2,
                                     X1 = NULL, X2, lambda_beta = NULL,
                                     variance_type = c("sandwich", "fisher"),
                                     include_estimation_uncertainty = TRUE,
                                     model = c("logistic", "poisson"),
                                     U1 = NULL, U2 = NULL,
                                     delta = NULL, delta_est = NULL,
                                     offset1 = NULL, offset2 = NULL,
                                     mu2 = NULL, variance2 = NULL) {

  variance_type <- match.arg(variance_type)
  #sigmoid <- function(x) 1 / (1 + exp(-x))
  model <- match.arg(model)

  if (is.null(lambda_beta)) {
    lambda_beta <- 0
  }
  
  n <- nrow(X2)
  p <- length(beta_0)
  # delta <- beta_a - beta_0
  # 
  # # Calculate components for test data (always needed)
  # eta2 <- X2 %*% (beta_0 + theta2 * delta)
  # p2 <- sigmoid(eta2)
  # R22 <- mean(p2 * (1 - p2) * (X2 %*% delta)^2)
  
  delta_pen <- beta_a - beta_0
  g2 <- as.vector(X2 %*% delta_pen)
  
  beta_theta <- beta_0 + theta2 * delta_pen
  test_moments <- compute_eta_mu_variance(
    X = X2,
    beta = beta_theta,
    model = model,
    U = U2,
    delta = delta,
    offset = offset2,
    mu = mu2,
    variance = variance2
  )
  mu2 <- test_moments$mu
  variance2 <- test_moments$variance
  
  R22 <- mean(variance2 * g2^2)
  
  if (variance_type == "fisher") {
    # Fisher information approach
    return(1 / (R22 * n))
  }
  
  # Sandwich estimator approach
  # psi2 <- (y2 - p2) * as.vector(X2 %*% delta)
  psi2 <- (y2 - mu2) * g2
  
  S22 <- mean(psi2^2)
  term1 <- S22 / (R22^2 * n)
  
  if (!include_estimation_uncertainty || is.null(X1) || is.null(y1)) {
    return(term1)
  }
  #--------- usually until here -------
  
  
  # Include uncertainty from beta estimation
  n_beta <- nrow(X1)
  
  # Components from estimation step
  # p1 <- sigmoid(X1 %*% beta_a)
  # W1 <- diag(as.vector(p1 * (1 - p1)))
  # R11 <- -crossprod(X1, W1 %*% X1) / n_beta - lambda_beta * diag(p)
  # 
  # # Cross terms
  # R12 <- t(X2) %*% (y2 - p2) - theta2 * t(X2) %*% diag(c(p2 * (1 - p2))) %*% (X2 %*% delta)
  # 
  # # S11 calculation
  # psi1 <- t(X1) %*% (y1 - p1) - lambda_beta * beta_a
  
  est_moments <- compute_eta_mu_variance(
    X = X1,
    beta = beta_a,
    model = model,
    U = U1,
    delta = delta_est,
    offset = offset1
  )
  mu1 <- est_moments$mu
  variance1 <- est_moments$variance
  
  X1_weighted <- X1 * as.vector(variance1)
  R11 <- -crossprod(X1, X1_weighted) / n_beta - lambda_beta * diag(p)
  
  R12 <- crossprod(X2, (y2 - mu2) - theta2 * variance2 * g2)
  
  psi1 <- crossprod(X1, y1 - mu1) - lambda_beta * beta_a
  
  S11 <- matrix(0, ncol = p, nrow = p)
  for (i in 1:n_beta) {
    psi1_i <- (y1[i] - mu1[i]) * X1[i, ]
    S11 <- S11 + tcrossprod(psi1_i)
  }
  S11 <- S11 / n_beta
  
  # Combine terms
  R11_inv <- solve(R11)
  term2_middle <- -R11_inv %*% R12 / R22
  term2 <- as.numeric(t(term2_middle) %*% S11 %*% term2_middle) / n_beta
  
  return(term1 + term2)
}

# Wrapper for backward compatibility
calc_V22 <- function(beta_0, beta_a, theta2, y1, y2, X1, X2, lambda_beta,
                     term1ONLY = FALSE,
                     model = c("logistic", "poisson"),
                     U1 = NULL, U2 = NULL,
                     delta = NULL, delta_est = NULL,
                     offset1 = NULL, offset2 = NULL,
                     mu2 = NULL, variance2 = NULL) {
  calc_asymptotic_variance(
    beta_0 = beta_0, beta_a = beta_a, theta2 = theta2,
    y1 = y1, y2 = y2, X1 = X1, X2 = X2,
    lambda_beta = lambda_beta,
    variance_type = "sandwich",
    include_estimation_uncertainty = !term1ONLY, 
    model = model,
    U1 = U1, U2 = U2,
    delta = delta, delta_est = delta_est,
    offset1 = offset1, offset2 = offset2,
    mu2 = mu2, variance2 = variance2
  )
}

AS_stdef <- function(beta_0, beta_a, theta2, y2, X2, fisher = FALSE,
                     model = c("logistic", "poisson"),
                     U2 = NULL, delta = NULL,
                     offset2 = NULL, mu2 = NULL, variance2 = NULL) {
  model <- match.arg(model)
  sqrt(calc_asymptotic_variance(
    beta_0 = beta_0, beta_a = beta_a, theta2 = theta2,
    y2 = y2, X2 = X2,
    variance_type = if (fisher) "fisher" else "sandwich",
    include_estimation_uncertainty = FALSE,
    model = model,
    U2 = U2, delta = delta,
    offset2 = offset2, mu2 = mu2, variance2 = variance2
  ))
}
# ---- Main Signpost Test Functions ----

#' Signpost Test for Real Data with resampling (Updated)
#'
#' Updated version using ridge_efficient with glmnet scaling and generalized theta_inf_hat.
#'
#' @param Y Numeric vector of response variables (length n)
#' @param X Numeric matrix of predictors (n x p)
#' @param gammaH0 Numeric scalar. Null hypothesis value for gamma parameter (default: 0)
#' @param model Character string. Model type, either "logistic" or "poisson" (default: "logistic")
#' @param lambda Numeric scalar. Regularization parameter (glmnet scale)
#' @param beta_0 Numeric vector. Null hypothesis beta coefficients (length p)
#' @param beta_a Numeric vector. Alternative hypothesis beta coefficients (length p)
#' @param B numeric, specifying the size of the null distribution (per test)
#' @param lambda_beta Numeric scalar. Regularization parameter for beta estimation (glmnet scale)
#' @param Y_a Numeric vector. Response variables for beta estimation (required if estimate_beta = TRUE)
#' @param X_a Numeric matrix. Predictors for beta estimation (required if estimate_beta = TRUE)
#' @param estimate_beta Logical. Whether to estimate beta_a from auxiliary data (default: FALSE)
#' @param scaled Logical. Whether to scale estimated beta_a (default: TRUE)
#' @param n_cores Integer. Number of cores for parallel processing (default: 8)
#' @param theta_min Numeric. Lower bound for theta search (default: 0)
#' @param theta_max Numeric. Upper bound for theta search (default: 1)
#'
#' @return A list containing p_value and theta_hat
#' @importFrom glmnet glmnet coef
#' @export
Signpost_test_resample <- function(Y, X, gammaH0 = 0, model = "logistic", lambda,
                                   beta_0, beta_a, B = 1000,
                                   lambda_beta = NULL, Y_a = NULL, X_a = NULL,
                                   estimate_beta = FALSE,
                                   scaled = TRUE, n_cores = 8, 
                                   theta_min = 0, theta_max = 1) {
  # Input validation
  if (estimate_beta && (is.null(lambda_beta) || is.null(Y_a) || is.null(X_a))) {
    stop("When estimate_beta = TRUE, you must provide lambda_beta, Y_a, and X_a")
  }
  
  n <- nrow(X)
  p <- ncol(X)
  n_beta <- if (!is.null(X_a)) nrow(X_a) else NULL
  
  # Estimate beta_a if requested
  if (estimate_beta) {
    # Use ridge_efficient for beta estimation
    if (model %in% c("logistic")) {
      #beta_a_hat <- ridge_efficient(Y = Y_a, X = X_a, lambda = lambda_beta, model = model)
      fit <- glmnet::glmnet(X_a, Y_a, 
                            family = "binomial",
                            alpha = 0,           # Ridge (not lasso)
                            lambda = lambda_beta,
                            intercept = FALSE,   # Match your existing setup
                            standardize = FALSE) # Don't standardize (match ridge_efficient)
      
      # Extract coefficients (excluding intercept)
      beta_a_hat <- as.vector(coef(fit))[-1]
    } else if (model == "poisson") {
      # result <- ridge_efficient(Y = Y_a, X = X_a, lambda = lambda_beta, model = model)
      # beta_a_hat <- result$penalized
      fit <- glmnet::glmnet(X_a, Y_a,
                            family = "poisson",
                            alpha = 0, 
                            lambda = lambda_beta,
                            intercept = FALSE,
                            standardize = FALSE)
      
      beta_a_hat <- as.vector(coef(fit))[-1]
    }
    
    if (scaled) {
      beta_a_hat <- beta_a_hat / sqrt(sum(beta_a_hat^2))
    }
  } else {
    beta_a_hat <- beta_a
  }
  
  # Estimate theta_hat with custom endpoints
  if (is.infinite(lambda)) {
    theta_hat <- theta_inf_hat(Y, X, matrix(ncol = 0, nrow = n), beta_0, beta_a_hat, 
                               model, theta_min = theta_min, theta_max = theta_max)
  } else {
    theta_hat <- theta_hat_lambda(Y, X, beta_0, beta_a_hat, lambda, model)
  }
  
  # Generate null distribution using gammaH0
  eta_null <- calculate_eta(X, gammaH0, beta_0, beta_a_hat)
  thetasH0 <- unlist(mclapply(1:B, function(idx) {
    if (estimate_beta) {
      beta_a_hat_k <- generate_beta_est(beta_a_hat, n_beta, X_a, model, lambda_beta)
    } else {
      beta_a_hat_k <- beta_a_hat
    }
    theta_Y(X, beta_0, beta_a_hat_k, eta_null, model, lambda, theta_min = theta_min, theta_max = theta_max)
  }, mc.cores = n_cores))
  
  # Calculate p-value based on the null distribution
  p_val <- calc_p_val(theta_hat, gammaH0, thetasH0)
  
  # Return results
  return(list(
    p_value = p_val,
    theta_hat = theta_hat,
    thetasH0 = thetasH0
  ))
}




#' Signpost Test for Real Data (Updated with ridge_efficient)
#'
#' Updated version using ridge_efficient with glmnet scaling convention.
#'
#' @param Y Numeric vector of response variables (length n)
#' @param X Numeric matrix of predictors (n x p)
#' @param U Numeric matrix of predictors (n x q) that will not be penalized (default: empty)
#' @param gammaH0 Numeric scalar. Null hypothesis value for gamma parameter (default: 0)
#' @param model Character string. Model type, either "logistic" or "poisson" (default: "logistic")
#' @param lambda Numeric scalar. Regularization parameter (glmnet scale)
#' @param beta_0 Numeric vector. Null hypothesis beta coefficients (length p)
#' @param beta_a Numeric vector. Alternative hypothesis beta coefficients (length p)
#' @param n_beta Integer. Sample size for beta estimation (required if estimate_beta = TRUE)
#' @param lambda_beta Numeric scalar. Regularization parameter for beta estimation (glmnet scale)
#' @param Y_a Numeric vector. Response variables for beta estimation (required if estimate_beta = TRUE)
#' @param X_a Numeric matrix. Predictors for beta estimation (required if estimate_beta = TRUE)
#' @param estimate_beta Logical. Whether to estimate beta_a from auxiliary data (default: FALSE)
#' @param use_theta_hat_in_variance Logical. Whether to use theta_hat in variance calculation (default: FALSE)
#' @param fisher Logical. Whether to use Fisher information for standard error (default: FALSE)
#' @param term1ONLY Logical. Whether to use only the first term in variance calculation (default: FALSE)
#' @param scaled Logical. Whether to scale estimated beta_a (default: TRUE)
#'
#' @return A list containing p_value, theta_hat, std_error, test_statistic
#' @export
signpost_test_AS <- function(Y, X, U = matrix(ncol = 0, nrow = n), gammaH0 = 0, model = "logistic", lambda,
                             beta_0, beta_a, 
                             n_beta = NULL, lambda_beta = NULL, Y_a = NULL, X_a = NULL,
                             estimate_beta = FALSE,
                             use_theta_hat_in_variance = FALSE,
                             fisher = FALSE, term1ONLY = FALSE,
                             scaled = TRUE) {
  
  # Input validation
  if (estimate_beta && (is.null(n_beta) || is.null(lambda_beta) || is.null(Y_a) || is.null(X_a))) {
    stop("When estimate_beta = TRUE, you must provide n_beta, lambda_beta, Y_a, and X_a")
  }
  
  n <- nrow(X)
  p <- ncol(X)
  
  # Estimate beta_a if requested
  if (estimate_beta) {
    # Use ridge_efficient for beta estimation
    if (model %in% c("logistic", "linear")) {
      beta_a_hat <- ridge_efficient(Y = Y_a, X = X_a, lambda = lambda_beta,U = U, model = model)
    } else if (model == "poisson") {
      result <- ridge_efficient(Y = Y_a, X = X_a, lambda = lambda_beta, U = U, model = model)
      beta_a_hat <- result$penalized
    }
    
    if (scaled) {
      beta_a_hat <- beta_a_hat / sqrt(sum(beta_a_hat^2))
    }
  } else {
    beta_a_hat <- beta_a
  }
  
  # Estimate theta_hat
  if (is.infinite(lambda)) {
    theta_hat <- theta_inf_hat(Y, X, U, beta_0, beta_a_hat, model)
  } else {
    theta_hat <- theta_hat_lambda(Y, X, U, beta_0, beta_a_hat, lambda, model)
  }
  
  # Calculate standard error
  if (fisher) {
    # Fisher information approach
    if (use_theta_hat_in_variance) {
      diff_beta <- X %*% (beta_a_hat - beta_0)
      lp0 <- ((1 - theta_hat) * beta_0 + theta_hat * beta_a_hat) %*% t(X)
      ddg <- diag(c(exp(-lp0) / (1 + exp(-lp0))^2))
      Fisher_info <- as.numeric(t(diff_beta) %*% ddg %*% diff_beta)
      stdef <- sqrt(1 / Fisher_info)
    } else {
      diff_beta <- X %*% (beta_a_hat - beta_0)
      lp0 <- beta_0 %*% t(X)
      ddg <- diag(c(exp(-lp0) / (1 + exp(-lp0))^2))
      Fisher_info <- as.numeric(t(diff_beta) %*% ddg %*% diff_beta)
      stdef <- sqrt(1 / Fisher_info)
    }
  } else {
    # Asymptotic variance approach
    theta2_for_variance <- if (use_theta_hat_in_variance) theta_hat else 0
    v22 <- calc_V22(beta_0 = beta_0, beta_a = beta_a_hat, theta2 = theta2_for_variance, 
                    y1 = Y_a, y2 = Y, X1 = X_a, X2 = X, lambda_beta = lambda_beta, 
                    term1ONLY = term1ONLY)
    
    if (!is.finite(v22) || v22 <= 0) {
      warning("Invalid variance estimate")
      return(list(p_value = NA, theta_hat = theta_hat, std_error = NA))
    }
    stdef <- sqrt(v22)
  }
  
  # Compute p-value
  if (gammaH0 == 0) {
    p_value <- pnorm(theta_hat, mean = 0, sd = stdef, lower.tail = FALSE)
  } else {
    p_value <- 2 * pnorm(-abs(theta_hat - gammaH0) / stdef)
  }
  
  # Return results
  return(list(
    p_value = p_value,
    theta_hat = theta_hat,
    std_error = stdef,
    test_statistic = theta_hat / stdef
  ))
}

