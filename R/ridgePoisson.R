# ============================================================================
#                    RIDGE EFFICIENT - UNIFIED SCALING & EFFICIENCY
# ============================================================================
#' Targeted ridge wrapper with local standardization, intercept, and near-zero drop
#'
#' Fits (targeted) ridge for GLMs in a standardized feature space and returns
#' coefficients on the raw/original scale. Supports an unpenalized intercept
#' and optional user-supplied unpenalized covariates (`U_extra`).
#'
#' @param Y numeric vector of responses (length n).
#' @param X numeric matrix of penalized covariates (n x p) on the RAW scale.
#' @param lambda nonnegative scalar regularization parameter (glmnet scale; see `lambda_scale`).
#' @param target optional numeric vector of length p (target for \eqn{\beta}); if NULL, untargeted ridge.
#' @param model character, one of c("logistic","gaussian","poisson").
#' @param intercept logical; if TRUE, include an unpenalized intercept internally (recommended).
#' @param intercept_policy character; strategy for handling zero-variance penalized columns.
#'   One of:
#'   - "auto" (default): drop zero-variance columns and, if `intercept = FALSE`, flip the
#'     intercept on when any are found.
#'   - "keep_one": when `intercept = FALSE`, retain exactly one zero-variance column (preferring
#'     the one with the largest absolute mean) instead of adding an intercept. Drops all when
#'     `intercept = TRUE`.
#'   - "drop_all": drop every zero-variance penalized column without modifying the intercept
#'     setting; emits a warning if this removes all such columns while `intercept = FALSE`.
#' @param U_extra optional (n x q) matrix of extra unpenalized covariates (RAW scale; not centered/scaled).
#' @param standardize logical; if TRUE, center & scale penalized columns locally (see `center`, `scale`).
#' @param center logical; if TRUE, subtract column means of X (penalized part only).
#' @param scale logical; if TRUE, divide columns by `s` (penalized part only).
#' @param scale_convention character, "glmnet" (default; divisor n) or "sd_unbiased" (divisor n-1).
#' @param zero_var_policy list with elements:
#'   - tol_abs: absolute floor on post-centering scale (default 1e-8)
#'   - tol_rel: relative floor vs median scale (default 1e-6)
#'   - rule: "abs_or_rel" (default) — drop if below either threshold
#' @param lambda_scale character, "glmnet" (default) or "n_times".
#'   - "glmnet": internal solver expects lambda as given
#'   - "n_times": internal solver expects lambda*n; wrapper multiplies by n
#' @param solver_fun function(Y, X, lambda, target, model, U, ...) -> either:
#'   - numeric vector of length ncol(U) + ncol(X) (unpenalized first, then penalized), or
#'   - list with elements `alpha` (scalar intercept if present), `gamma` (U_extra coefs),
#'     and `beta` (penalized coefs in the STANDARDIZED space).
#'   Default is `ridgeGLM`.
#' @param ... passed on to `solver_fun` (e.g., maxIter, tolerance).
#'
#' @return list with elements:
#'   - beta: length-p vector (RAW scale), zeros in dropped positions
#'   - alpha: scalar intercept on RAW scale (if intercept=TRUE), else 0
#'   - gamma: coefficients for U_extra (RAW scale), if provided
#'   - kept, dropped: integer indices of columns kept/dropped
#'   - mu, s: centering/scaling stats for penalized X (length p)
#'   - diagnostics: list with KKT residuals (norm2, max_abs), objective pieces, etc.
#'   - meta: list with settings (model, lambda_in, lambda_internal, scale options, policy)
#'
#' @examples
#' # beta_hat <- ridge_complete(Y, X, lambda = 0.5, target = rep(0, ncol(X)),
#' #                             model = "logistic", intercept = TRUE)
ridge_complete <- function(
    Y, X, lambda, target = NULL,
    model = c("logistic","gaussian","poisson"),
    intercept = TRUE,
    U_extra = NULL,
    standardize = TRUE,
    center = standardize,
    scale = standardize,
    scale_convention = c("glmnet","sd_unbiased"),
    zero_var_policy = list(tol_abs = 1e-8, tol_rel = 1e-6, rule = "abs_or_rel"),
    lambda_scale = "glmnet", #c("glmnet","n_times"),
    solver_fun = ridgeGLM,
    intercept_policy = c("auto","keep_one","drop_all"),
    ...
) {
  model <- match.arg(model)
  scale_convention <- match.arg(scale_convention)
  lambda_scale <- match.arg(lambda_scale)
  intercept_policy <- match.arg(intercept_policy)
  intercept   <- isTRUE(intercept)
  intercept_requested <- intercept
  intercept_auto_flipped <- FALSE
  kept_constant_idx <- integer(0)
  standardize <- isTRUE(standardize)
  center      <- isTRUE(center)
  scale       <- isTRUE(scale)
  
  
  # Basic dims
  Y <- as.numeric(Y)
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  if (length(Y) != n) stop("Y and X have incompatible lengths.")
  if (!is.null(target)) {
    if (length(target) != p) stop("target must have length equal to ncol(X).")
    target <- as.numeric(target)
  }
  
  # ---- Build unpenalized design U (intercept + U_extra) on RAW scale ----
  if (!is.null(U_extra)) U_extra <- as.matrix(U_extra)
  if (!is.null(U_extra) && nrow(U_extra) != n) stop("U_extra has incompatible number of rows.")
  
  # ---- Compute centering and scale for penalized X (local; for dropping & std) ----
  # For dropping, we compute post-centering scale with divisor "glmnet" (n) regardless of scale_convention.
  mu_drop <- colMeans(X)
  s_drop <- sqrt(colMeans( sweep(X, 2L, mu_drop, FUN = "-")^2 ))
  
  # Apply near-zero policy
  tol_abs <- zero_var_policy$tol_abs %||% 1e-8
  tol_rel <- zero_var_policy$tol_rel %||% 1e-6
  rule    <- zero_var_policy$rule    %||% "abs_or_rel"
  med_s   <- stats::median(s_drop[is.finite(s_drop)], na.rm = TRUE)
  if (!is.finite(med_s) || med_s <= 0) med_s <- 1.0
  
  # ---- Check for zero-variance columns in U_extra (before building U) ----
  if (!is.null(U_extra) && ncol(U_extra) > 0) {
    # Compute column means and standard deviations for U_extra
    mu_U <- colMeans(U_extra)
    s_U <- sqrt(colMeans(sweep(U_extra, 2L, mu_U, FUN = "-")^2))
    
    # Apply near-zero policy to U_extra (using same policy as X)
    med_s_U <- stats::median(s_U[is.finite(s_U)], na.rm = TRUE)
    if (!is.finite(med_s_U) || med_s_U <= 0) med_s_U <- 1.0
    
    small_abs_U <- s_U < tol_abs
    small_rel_U <- s_U < (tol_rel * med_s_U)
    drop_mask_U <- switch(rule,
                          abs_or_rel = (small_abs_U | small_rel_U),
                          abs_and_rel = (small_abs_U & small_rel_U),
                          abs_only = small_abs_U,
                          rel_only = small_rel_U,
                          stop("Unknown zero_var_policy$rule")
    )
    keep_idx_U <- which(!drop_mask_U)
    
    # Filter U_extra to keep only non-zero-variance columns
    if (length(keep_idx_U) < ncol(U_extra)) {
      warning(sprintf("Dropped %d zero-variance columns from U_extra", 
                      ncol(U_extra) - length(keep_idx_U)))
      U_extra <- U_extra[, keep_idx_U, drop = FALSE]
    }
  }
  
  small_abs <- s_drop < tol_abs
  small_rel <- s_drop < (tol_rel * med_s)
  drop_mask <- switch(rule,
                      abs_or_rel = (small_abs | small_rel),
                      abs_and_rel = (small_abs & small_rel),
                      abs_only = small_abs,
                      rel_only = small_rel,
                      stop("Unknown zero_var_policy$rule")
  )
  zero_var_idx <- which(drop_mask)

  if (length(zero_var_idx) > 0L) {
    if (identical(intercept_policy, "auto")) {
      if (!intercept) {
        intercept <- TRUE
        intercept_auto_flipped <- TRUE
      }
    } else if (identical(intercept_policy, "keep_one")) {
      if (!intercept) {
        ord <- zero_var_idx[order(-abs(mu_drop[zero_var_idx]))]
        keep_const <- ord[1L]
        drop_mask[keep_const] <- FALSE
        kept_constant_idx <- as.integer(keep_const)
      }
    } else if (identical(intercept_policy, "drop_all")) {
      if (!intercept) {
        warning(sprintf(
          "Dropping %d zero-variance penalized columns while intercept = FALSE.",
          length(zero_var_idx)
        ))
      }
    }
  }

  keep_idx <- which(!drop_mask)
  drop_idx <- which(drop_mask)

  # Edge case: if nothing kept, we still want an intercept-only (and U_extra) fit.
  # We'll pass a 0-column matrix to the solver and handle fallbacks if needed.
  X_keep <- if (length(keep_idx)) X[, keep_idx, drop = FALSE] else X[, FALSE, drop = FALSE]
  p_keep <- ncol(X_keep)
  
  # Now build U after U_extra has been filtered
  if (intercept) {
    U <- cbind(`(Intercept)` = rep(1.0, n), U_extra)
  } else {
    U <- U_extra
  }
  pU <- if (is.null(U)) 0L else ncol(U)
  if(is.null(U)){
    U =  matrix(ncol = 0, nrow = n)
  }
  
  
  # Prepare centering/scaling statistics for the penalized block (length p)
  mu <- rep(0, p); s <- rep(1, p)    # RAW-to-STD transform stats
  if (p_keep > 0) {
    if (center) {
      mu[keep_idx] <- colMeans(X_keep)
      if (length(kept_constant_idx)) {
        mu[kept_constant_idx] <- 0
      }
    } else {
      mu[keep_idx] <- 0
    }
    if (scale) {
      # choose divisor
      if (scale_convention == "glmnet") {
        s[keep_idx] <- sqrt(colMeans( sweep(X_keep, 2, mu[keep_idx], "-")^2 ))
      } else {
        # unbiased: divisor n-1
        cm <- sweep(X_keep, 2, mu[keep_idx], "-")
        s[keep_idx] <- apply(cm, 2, function(z) stats::sd(z))  # sd uses n-1
      }
      # Protect against zeros after policy (shouldn't happen if drop policy used s_drop)
      s[keep_idx][s[keep_idx] == 0] <- 1
      if (length(kept_constant_idx)) {
        s[kept_constant_idx] <- 1
      }
    } else {
      s[keep_idx] <- 1
    }
  }
  
  # ---- Standardize penalized X (kept columns only) ----
  Xc <- if (p_keep > 0) sweep(X_keep, 2L, mu[keep_idx], FUN = "-") else X_keep
  Xs <- if (p_keep > 0) sweep(Xc, 2L, s[keep_idx], FUN = "/") else Xc
  
  # ---- Transform target to standardized space (kept columns only) ----
  if (!is.null(target) && p_keep > 0) {
    t_keep <- target[keep_idx]
    t_star <- as.numeric(s[keep_idx] * t_keep)
  } else if (!is.null(target) && p_keep == 0) {
    t_keep <- numeric(0); t_star <- numeric(0)
  } else {
    t_keep <- NULL; t_star <- NULL
  }
  # ---- Lambda scaling for the solver ----
  lambda_internal <- switch(lambda_scale,
                            glmnet = lambda,
                            n_times = lambda * n
  )
  # ---- Call the core solver on standardized space ----
  fit_raw <- tryCatch(
    ridge_efficient(
      Y = Y,
      X = Xs,                         # standardized penalized matrix
      lambda = lambda_internal,
      target = t_star,                # target in standardized space (or NULL)
      model = model,
      U = U,                          # unpenalized design (intercept + U_extra)
      ...
    ),
    error = function(e) e
  )

  if (inherits(fit_raw, "error")) {
    # Fallback if solver cannot handle zero columns: do an intercept-only solution
    if (p_keep == 0L) {
      # Closed-form-ish intercept for canonical GLMs given offset 0
      alpha_star <- switch(model,
                           gaussian = mean(Y),
                           logistic = qlogis( mean(Y) ),
                           poisson  = log( mean(Y) ),
                           stop("Unknown model for fallback.")
      )
      gamma <- if (pU > as.integer(intercept)) rep(0, pU - as.integer(intercept)) else numeric(0)
      beta_star <- numeric(0)
      solver_status <- list(note = "solver fallback: intercept-only")
    } else {
      stop(fit_raw$message)
    }
  } else {
    # Parse solver output
    if (is.list(fit_raw)) {
      # Expect: list(alpha, gamma, beta) in STANDARDIZED space
      alpha_star <- if (!is.null(fit_raw$alpha)) as.numeric(fit_raw$alpha) else 0
      gamma      <- if (!is.null(fit_raw$gamma)) as.numeric(fit_raw$gamma) else {
        if (pU > as.integer(intercept)) rep(0, pU - as.integer(intercept)) else numeric(0)
      }
      beta_star  <- if (!is.null(fit_raw$beta)) as.numeric(fit_raw$beta) else {
        # try a generic 'coef' name
        if (!is.null(fit_raw$coef)) as.numeric(fit_raw$coef) else stop("solver_fun list output missing $beta.")
      }
    } else if (is.numeric(fit_raw)) {
      # Numeric vector: [unpenalized (pU), then penalized (p_keep)]
      if (length(fit_raw) != (pU + p_keep)) {
        stop("solver_fun returned a numeric vector of unexpected length.")
      }
      if (pU > 0) {
        u_coef <- fit_raw[seq_len(pU)]
        alpha_star <- if (intercept) u_coef[1L] else 0
        gamma <- if (intercept && pU > 1L) u_coef[-1L] else if (!intercept) u_coef else numeric(0)
      } else {
        alpha_star <- 0; gamma <- numeric(0)
      }
      beta_star <- if (p_keep > 0) fit_raw[(pU + 1L):(pU + p_keep)] else numeric(0)
    } else {
      stop("solver_fun returned neither list nor numeric vector.")
    }
    solver_status <- list(note = "ok")
  }
  
  # ---- Back-transform to RAW scale ----
  beta <- numeric(p)
  if (p_keep > 0) {
    beta_keep <- beta_star / s[keep_idx]
    beta[keep_idx] <- beta_keep
  }
  # Intercept mapping: alpha = alpha* - mu^T beta_keep
  alpha <- if (intercept) {
    if (p_keep > 0) alpha_star - sum(mu[keep_idx] * beta[keep_idx]) else alpha_star
  } else {
    0
  }
  
  # ---- Diagnostics ---------------------------------------------------------
  
  # 1) KKT residuals (standardized space) for penalized block
  kkt <- list(norm2 = NA_real_, max_abs = NA_real_)
  obj <- list(nll_mean = NA_real_, penalty = NA_real_, obj_mean = NA_real_)
  loglik_raw <- NA_real_
  
  # Only compute if we actually had a penalized block (or even if empty, g=0, penalty=0)
  # Build linear predictor in standardized space:
  eta_star <- {
    offs <- if (intercept) alpha_star else 0
    offs <- if (pU > as.integer(intercept) && !is.null(U_extra)) {
      # add U_extra contribution on standardized space (U not standardized; coefs already on RAW scale)
      # Note: solver returns gamma already applicable to U (unstandardized), so eta includes U %*% gamma
      # We add both intercept and U terms:
      offs + drop(U %*% c(if (intercept) alpha_star else NULL, gamma)) - if (intercept) alpha_star else 0
    } else offs
    if (p_keep > 0) {
      # Prefer consistent computation: eta_star = (intercept+U part) + Xs %*% beta_star
      offs + drop(Xs %*% beta_star)
    } else {
      offs
    }
  }
  
  # Mean negative loglik (STANDARDIZED space; U terms are already included via alpha_star/gamma)
  nll_mean <- switch(model,
                     logistic = mean( log1p(exp(eta_star)) - Y * eta_star ),
                     gaussian = 0.5 * mean( (Y - eta_star)^2 ),
                     poisson  = mean( exp(eta_star) - Y * eta_star )
  )
  
  # Penalty (STANDARDIZED space) on penalized block
  if (p_keep > 0) {
    if (!is.null(t_star)) {
      pen <- 0.5 * lambda * sum((beta_star - t_star)^2)
    } else {
      pen <- 0.5 * lambda * sum(beta_star^2)
    }
  } else {
    pen <- 0
  }
  
  obj$nll_mean <- nll_mean
  obj$penalty  <- pen
  obj$obj_mean <- nll_mean + pen
  
  # KKT: gradient w.r.t. penalized coefficients in standardized space
  if (p_keep > 0) {
    mu_hat <- switch(model,
                     logistic = plogis(eta_star),
                     gaussian = eta_star,
                     poisson  = exp(eta_star)
    )
    # gradient of mean loss: (1/n) Xs^T (mu_hat - Y) for logistic/poisson;
    # for gaussian with nll=0.5*mean((Y-eta)^2), grad = (1/n) Xs^T(eta - Y) == (1/n) Xs^T(mu_hat - Y)
    g_beta <- as.numeric(crossprod(Xs, (mu_hat - Y)) / n)
    # penalty gradient: lambda * (beta_star - t_star)  [glmnet-style scale]
    # This is correct regardless of internal lambda scaling, so we use user-supplied `lambda`.
    pen_grad <- if (!is.null(t_star)) (beta_star - t_star) else beta_star
    kkt_vec <- g_beta + lambda * pen_grad
    kkt$norm2  <- sqrt(sum(kkt_vec^2))
    kkt$max_abs <- max(abs(kkt_vec))
  } else {
    kkt$norm2 <- 0; kkt$max_abs <- 0
  }
  
  # Unpenalized log-likelihood on RAW scale (for interpretability)
  # eta_raw = alpha + U_extra%*%gamma + X%*%beta
  eta_raw <- {
    offs <- if (intercept) alpha else 0
    if (!is.null(U_extra) && length(gamma)) offs <- offs + drop(U_extra %*% gamma)
    offs + drop(X %*% beta)
  }
  loglik_raw <- switch(model,
                       logistic = sum( Y * eta_raw - log1p(exp(eta_raw)) ),
                       gaussian = -0.5 * sum( (Y - eta_raw)^2 ),
                       poisson  = sum( Y * eta_raw - exp(eta_raw) )
  )
  
  # ---- Return --------------------------------------------------------------
  out <- list(
    beta = beta,              # RAW scale; zeros in dropped columns
    alpha = alpha,            # RAW intercept
    gamma = if (!is.null(U_extra)) gamma else numeric(0),
    kept = keep_idx,
    dropped = drop_idx,
    mu = mu,                  # centering stats for penalized block (length p)
    s  = s,                   # scaling stats for penalized block (length p)
    diagnostics = list(
      kkt = kkt,
      objective = obj,
      loglik_raw = loglik_raw,
      solver_status = solver_status
    ),
    meta = list(
      model = model,
      intercept = intercept,
      p = p,
      n = n,
      lambda_in = lambda,
      lambda_internal = lambda_internal,
      lambda_scale = lambda_scale,
      standardize = standardize,
      center = center,
      scale = scale,
      scale_convention = scale_convention,
      intercept_requested = intercept_requested,
      intercept_policy = intercept_policy,
      intercept_auto_flipped = intercept_auto_flipped,
      kept_constant = kept_constant_idx,
      zero_var_policy = list(tol_abs = tol_abs, tol_rel = tol_rel, rule = rule),
      solver_fun = deparse(substitute(solver_fun))
    )
  )
  class(out) <- "ridge_efficient_fit"
  out
}

# Helper for default list elements
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Efficient Ridge Regression with Standardized Scaling
#'
#' Unified function that handles both scaling consistency (glmnet convention) and
#' computational efficiency for high-dimensional problems. Automatically chooses
#' between ridgeGLM for smaller problems and efficient IRWLS for large problems.
#'
#' @param Y Numeric vector, response variable
#' @param X Numeric matrix, design matrix  
#' @param lambda Numeric, penalty parameter (glmnet scale)
#' @param target Numeric vector, shrinkage target (default: zero vector)
#' @param model Character, one of "linear", "logistic", "poisson"
#' @param high_dim_threshold Integer, threshold for switching to efficient method (default: 1000)
#' @param ... Additional arguments
#'
#' @return Numeric vector of estimated coefficients (for logistic/linear) or 
#'         list with 'penalized' component (for poisson)
#' @importFrom porridge ridgeGLM
#' @importFrom pracma randortho
#' @export
ridge_efficient <- function(Y, X, lambda,U = matrix(ncol = 0, nrow = length(Y)), target = rep(0, ncol(X)), model = "logistic", 
                            high_dim_threshold = 10000, ...) {
  
  n <- length(Y)
  p <- ncol(X)
  
  # Default target to zero vector if not provided
  if (is.null(target)) {
    target <- rep(0, p)
  }
  
  # Handle infinite lambda case
  if (is.infinite(lambda)) {
    if (model == "poisson") {
      return(list(penalized = target, unpenalized = numeric(0))) #TODO this does not seem good yet.
    } else {
      return(target)
    }
  }
  
  # Choose method based on dimensions
  use_efficient <- (p > high_dim_threshold) || (p >= n && p > 10000)
 if (use_efficient && model == "logistic") {
    # Use efficient high-dimensional method with glmnet scaling
   print("efficient ridge method")
    return(ridgeGLM_noblk_memsafe(Y = Y, X = X,U = U, lambda = (lambda * n), target = target,  ...))
    
  } else if (model %in% c("logistic", "linear")) {
    # Use standard ridgeGLM with scaling conversion
    return(ridgeGLM(Y = Y, X = X, U= U, lambda = lambda * n,  target = target, model = model, ...))
    
  } else if (model == "poisson") {
    # Use AridgePLM with scaling conversion
    return(AridgePLM(Y = Y, X = X,U = U, lambda = lambda * n, target = target, model = model, ...))
  }
}

# if (p > 50000) {
#   return(ridge_targeted_chunked_exact(Y, X, lambda, target, model))
# }else


#' Fit glmnet-style scaling stats
#' @param X matrix or dgCMatrix
#' @param intercept logical; if TRUE center by column means
#' @param standardize logical; if TRUE scale to unit variance (denom n)
#' @param eps numeric; floor for tiny scales to detect zero variance
#' @return list with mu, s, keep_idx, intercept, standardize, and helpers
fit_glmnet_scaler <- function(X, intercept = TRUE, standardize = TRUE, eps = 1e-12) {
  n <- nrow(X); p <- ncol(X)
  stopifnot(n > 0L, p > 0L)
  
  # column means (dense or sparse)
  mu <- if (intercept) {
    if (inherits(X, "dgCMatrix")) {
      Matrix::colMeans(X)
    } else {
      colMeans(X)
    }
  } else {
    rep(0, p)
  }
  
  # column scales with denominator n
  if (standardize) {
    if (inherits(X, "dgCMatrix")) {
      # Variance with denominator n: E[(X - mu)^2] or E[X^2] if not centered
      Ex2 <- Matrix::colSums(X@x^2 * 0 + 1) # placeholder; we'll compute properly below
      # Proper sparse-safe computation:
      # E[X^2] = (colSums(X^2))/n; works whether intercept is TRUE or FALSE
      Ex2 <- Matrix::colSums(X^2) / n
      s <- if (intercept) {
        # Var = E[X^2] - (E[X])^2
        sqrt(pmax(Ex2 - mu^2, 0))
      } else {
        sqrt(pmax(Ex2, 0))
      }
    } else {
      if (intercept) {
        # center-then-scale with denom n
        xm <- sweep(X, 2L, mu, FUN = "-")
        s <- sqrt(colSums(xm * xm) / n)
      } else {
        s <- sqrt(colSums(X * X) / n)
      }
    }
  } else {
    s <- rep(1, p)
  }
  
  # detect zero-variance
  keep_idx <- which(s > eps)
  dropped  <- setdiff(seq_len(p), keep_idx)
  
  # protect against zeros
  s_eff <- s
  s_eff[s_eff <= eps] <- 1   # avoid divide-by-zero; those cols are dropped anyway
  
  # helpers
  transform_X <- function(Xnew) {
    if (!length(keep_idx)) return(Xnew[, integer(0L), drop = FALSE])
    Xk <- Xnew[, keep_idx, drop = FALSE]
    if (intercept) {
      Xk <- sweep(Xk, 2L, mu[keep_idx], FUN = "-")
    }
    if (standardize) {
      Xk <- sweep(Xk, 2L, s_eff[keep_idx], FUN = "/")
    }
    Xk
  }
  
  inverse_beta <- function(gamma_std) {
    # gamma_std is coef on standardized X (excluding intercept)
    beta <- numeric(length(s))
    beta[keep_idx] <- gamma_std / s_eff[keep_idx]
    beta
  }
  
  adjust_intercept <- function(gamma0, gamma_std) {
    if (!intercept) return(gamma0)
    # gamma0 is intercept on standardized problem
    # map to original: beta0 = gamma0 - sum_j gamma_j * mu_j / s_j  (over kept j)
    adj <- sum((gamma_std / s_eff[keep_idx]) * mu[keep_idx])
    gamma0 - adj
  }
  
  list(
    mu = mu, s = s, s_eff = s_eff,
    keep_idx = keep_idx, dropped_idx = dropped,
    intercept = intercept, standardize = standardize,
    transform_X = transform_X,
    inverse_beta = inverse_beta,
    adjust_intercept = adjust_intercept
  )
}

#same as ridgeGLM but chunks large matrix multiplications
ridgeGLM_noblk_memsafe <- function(
    Y, X, lambda, target = rep(0, ncol(X)), U = matrix(ncol = 0, nrow = length(Y)),  minSuccDiff = 1e-8, maxIter = 100, chunk_size = 5000
) {
  # Assumes: no U, no Dg, no lambdaG. Returns bHat (length p).
  # Works for dense X with n <= ~600 and p >> n (e.g., 60k+).
  # Uses chunked ops to avoid large temporaries.
  print(paste("lambda ", lambda))
  # -- light checks (avoid copying X) ----------------------------------------
  Y <- as.numeric(Y)
  if (!is.double(X)) storage.mode(X) <- "double"  # no as.matrix() to avoid a copy
  n <- nrow(X); p <- ncol(X)
  if (length(target) != p) stop("length(target) must equal ncol(X).")
  target <- as.numeric(target)
  
  # numerically-stable loglik(y, lp) for logistic
  # log1pexp <- function(x) ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
  # loglikBLMlp <- function(y, lp) sum(y * lp - log1pexp(lp))
  
  loglikBLMlp <- function(Y, lp){
    loglik1 <- exp(Y * lp)/(1 + exp(lp))
    loglik2 <- exp(((Y - 1) * lp))/(1 + exp(-lp))
    loglik1[!is.finite(log(loglik1))] <- NA
    loglik2[!is.finite(log(loglik2))] <- NA
    loglik <- sum(log(apply(cbind(loglik1, loglik2), 1, mean,na.rm = TRUE)))
    return(loglik)
  }
  
  
  
  # -- helpers: chunked ops over columns of X --------------------------------
  iter_cols <- function(p, by) {
    starts <- seq.int(1L, p, by = by)
    ends   <- pmin(starts + by - 1L, p)
    Map(seq.int, starts, ends)
  }
  
  # K = X X^T / lambda (n x n), built in chunks over columns
  build_K <- function(X, lambda, chunk) {
    n <- nrow(X)
    K <- matrix(0.0, n, n)
    for (jj in iter_cols(ncol(X), chunk)) {
      Xj <- X[, jj, drop = FALSE]
      # K += Xj %*% t(Xj)
      K <- K + tcrossprod(Xj)  # n x n, cheap
    }
    K / lambda
  }
  
  # Xtarget = X %*% target (length n), chunked
  build_Xtarget <- function(X, target, chunk) {
    n <- nrow(X)
    xt <- numeric(n)
    for (jj in iter_cols(length(target), chunk)) {
      xt <- xt + X[, jj, drop = FALSE] %*% target[jj]
    }
    as.numeric(xt)
  }
  
  # v = crossprod(X, s) (length p), chunked
  chunked_crossprod_X_mat <- function(X, S, chunk) {
    if (is.vector(S)) S <- matrix(S, ncol = 1L)
    if (nrow(S) != nrow(X)) stop("nrow(S) must equal nrow(X).")
    p <- ncol(X); c <- ncol(S)
    V <- matrix(0.0, nrow = p, ncol = c)
    for (jj in iter_cols(p, chunk)) {
      Xb <- X[, jj, drop = FALSE]
      V[jj, ] <- as.matrix(crossprod(Xb, S))   # k×c -> assign to rows jj
    }
    V
  }
  
  
  # -- precomputations --------------------------------------------------------
  bigp <- (p >= n)
  if (!bigp) {
    stop("This memsafe version is intended for p >> n; for p < n use a small-p path.")
  }
  
  K <- build_K(X, lambda, chunk_size)           # n x n, small
  Xtarget <- build_Xtarget(X, target, chunk_size) # length n
  diagK_orig <- diag(K)
  
  # -- IRLS loop (identical algebra to your big-p branch) --------------------
  lp <- numeric(n)
  loglikPrev <- -1e10
  
  
  
  if(ncol(U)==0){
    
    for (iter in 1:maxIter) {
      Ypred <- 1.0 / (1.0 + exp(-lp))
      W0 <- Ypred * (1.0 - Ypred)
      if (min(W0) <= .Machine$double.eps) {
        idx <- which(W0 < .Machine$double.eps)
        W0[idx] <- .Machine$double.eps
      }
      
      Z <- lp + (Y - Ypred) / W0
      
      # M = K + diag(1/W0)   (update diag only; no big allocs)
      diag(K) <- diagK_orig + 1.0 / W0
      
      # slh = solve(M, Z - Xtarget) via Cholesky (stable & fast)
      # Cf <- chol(K)  # n x n
      # rhs <- Z - Xtarget
      # slh <- backsolve(Cf, forwardsolve(t(Cf), rhs, upper.tri = TRUE, transpose = TRUE))
      slh <- solve(K, Z - Xtarget)
      
      # restore diag(K)
      diag(K) <- diagK_orig
      
      # lp_temp = crossprod(K, slh) using ORIGINAL K
      lp <- as.numeric(crossprod(K, slh))
      
      # penalty = sum(lp_temp * slh)/2
      penalty <- sum(lp * slh) / 2.0
      
      # lp = Xtarget + lp_temp
      lp <- Xtarget + lp
      
      loglik <- loglikBLMlp(Y, lp) - penalty
      if (abs(loglik - loglikPrev) < minSuccDiff) {
        break
      } else {
        loglikPrev <- loglik
      }
    }
    
    # -- Recover beta: bHat = target + (1/lambda) X' slh (chunked)
    Xt_s <- chunked_crossprod_X_mat(X, slh, chunk_size)
    bHat <- target + Xt_s / lambda
    bHat
  }else{
    tUTX <- chunked_crossprod_X_mat(X, U, chunk_size)
    
    for (iter in 1:maxIter) {
      Ypred <- 1.0 / (1.0 + exp(-lp))
      W0 <- as.numeric(Ypred * (1.0 - Ypred))
      if (min(W0) <= .Machine$double.eps) {
        idx <- which(W0 < .Machine$double.eps)
        W0[idx] <- .Machine$double.eps
      }
      
      Z <- lp + (Y - Ypred) / W0
      
      # M = K + diag(1/W0)   (update diag only; no big allocs)
      diag(K) <- diagK_orig + 1.0 / W0
      
      # slh = solve(M, Z - Xtarget) via Cholesky (stable & fast)
      # Cf <- chol(K)  # n x n
      # rhs <- Z - Xtarget
      # slh <- backsolve(Cf, forwardsolve(t(Cf), rhs, upper.tri = TRUE, transpose = TRUE))
      slh <- solve(K, Z - Xtarget)
      
      gHat <- solve(crossprod(U, solve(K, U)), crossprod(U,  slh))
      slh <- solve(K, Z - Xtarget - U %*% gHat)
      
      # restore diag(K)
      diag(K) <- diagK_orig
      
      # lp_temp = crossprod(K, slh) using ORIGINAL K
      lp <- as.numeric(crossprod(K, slh))
      
      # penalty = sum(lp_temp * slh)/2
      penalty <- sum(lp * slh) / 2.0
      
      # lp = Xtarget + lp_temp
      lp <- Xtarget + lp+ U %*% gHat
      
      loglik <- loglikBLMlp(Y, lp) - penalty
      if (abs(loglik - loglikPrev) < minSuccDiff) {
        break
      } else {
        loglikPrev <- loglik
      }
    }
    
    # -- Recover beta: bHat = target + (1/lambda) X' slh (chunked)
    Xt_s <- chunked_crossprod_X_mat(X, slh, chunk_size)
    bHat <- target + Xt_s / lambda
    return(c(gHat, bHat))
  }
  
}






#' Targeted Ridge Estimation for Poisson Regression
#'
#' Performs targeted ridge estimation for a Poisson regression model using an iteratively
#' reweighted least squares (IRWLS) algorithm. The estimator shrinks the penalized
#' regression coefficients towards a specified target. Note that although the parameter
#' \code{model} defaults to "logistic", the current implementation only supports Poisson regression.
#'
#' @param Y Numeric vector. The response variable (observed counts), assumed to follow a Poisson distribution.
#' @param X Numeric matrix. The design matrix for penalized covariates; the number of rows must equal the length of Y.
#' @param U Numeric matrix. The design matrix for unpenalized covariates. Default is an empty matrix with \code{nrow = length(Y)}.
#' @param lambda Positive numeric. The ridge penalty parameter.
#' @param lambdaG Numeric. The generalized ridge penalty parameter. Default is 0.
#' @param Dg Numeric matrix. A matrix of dimension \code{ncol(X)} x \code{ncol(X)} used for generalized ridge penalization. Default is a zero matrix.
#' @param target Numeric vector. The target towards which the penalized estimates are shrunk.
#'   Its length must equal the number of columns of \code{X}. Default is a zero vector.
#' @param model Character. The model type. Only Poisson regression is currently supported.
#' @param minSuccDiff Numeric. The minimum absolute change in log-likelihood between iterations to declare convergence. Default is \code{1e-10}.
#' @param maxIter Integer. The maximum number of IRWLS iterations allowed. Default is 100.
#' @param verbose Logical. If TRUE, prints the log-likelihood at each iteration. Default is FALSE.
#'
#' @return A named list with components:
#'   \item{unpenalized}{Numeric vector of estimated unpenalized coefficients (empty if none).}
#'   \item{penalized}{Numeric vector of estimated penalized coefficients.}
#'
#' @examples
#' set.seed(123)
#' n <- 100; p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' Y <- rpois(n, lambda = exp(X %*% rep(0.5, p)))
#' result <- AridgePLM(Y, X, lambda = 1, target = rep(0, p), verbose = TRUE)
#' result$penalized
#'
#' @export
AridgePLM <- function(Y, X, U = matrix(ncol = 0, nrow = length(Y)),
                      lambda, lambdaG = 0,
                      Dg = matrix(0, ncol = ncol(X), nrow = ncol(X)),
                      target = rep(0, ncol(X)), model = "poisson",
                      minSuccDiff = 1e-10, maxIter = 100, verbose = FALSE) {

  # Input validation
  if (!is.numeric(Y) || !is.matrix(X))
    stop("Y must be a numeric vector and X must be a numeric matrix.")
  if (length(Y) != nrow(X))
    stop("The number of rows in X must equal the length of Y.")
  if (length(target) != ncol(X))
    stop("The length of target must equal the number of columns of X.")

  # Currently, generalized ridge (Dg, lambdaG) and unpenalized covariates (U) are not supported.
  if ((max(abs(Dg)) != 0 || lambdaG != 0) || (ncol(U) != 0))
    stop("Generalized ridge penalization and unpenalized covariates are not supported in the current implementation.")

  # Initialize variables
  loglikPrev <- -1e10
  eta <- rep(0, length(Y))  # initial linear predictor
  high_dim <- (ncol(X) >= nrow(X))

  # Pre-calculate for high-dimensional case
  if (high_dim) {
    XXT <- tcrossprod(X) / lambda
    Xtarget <- t(tcrossprod(target, X))
  }

  # IRWLS iteration
  for (iter in 1:maxIter) {
    # Calculate the predicted mean (for Poisson, mu = exp(eta))
    mu <- exp(eta)
    # Calculate weights; for Poisson, variance equals the mean
    W <- mu
    # Avoid extremely small weights
    W[W < .Machine$double.eps] <- .Machine$double.eps

    if (high_dim) {
      # Compute the adjusted response (working variable)
      Z <- eta + (Y - mu) / W

      # Adjust the XXT matrix diagonal by adding 1/W
      diag(XXT) <- diag(XXT) + 1 / W
      # Solve for the update (efficient IRWLS update)
      slh <- solve(XXT, Z - Xtarget)
      diag(XXT) <- diag(XXT) - 1 / W

      eta_update <- crossprod(XXT, slh)
      penalty <- sum(eta_update * slh) / 2
      eta <- Xtarget + eta_update
      loglik <- AloglikPLMlp(Y, eta) - penalty
    } else {
      # Low-dimensional case
      XWZpT <- lambda * target + as.numeric(crossprod(X, W * eta + Y - mu))
      XTX <- crossprod(sweep(X, 1, sqrt(W), FUN = "*"))
      diag(XTX) <- diag(XTX) + lambda

      bHat <- solve(XTX, XWZpT)
      eta <- tcrossprod(X, t(bHat))
      loglik <- AloglikPLMlp(Y, eta) - 0.5 * lambda * sum((bHat - target)^2)
    }

    if (verbose) {
      cat("Iteration", iter, ": log-likelihood =", loglik, "\n")
    }

    # Check convergence
    if (is.infinite(loglik)) {
      loglikPrev <- -1e10
    } else if (abs(loglik - loglikPrev) < minSuccDiff) {
      break
    } else {
      loglikPrev <- loglik
    }
  }

  # In the high-dimensional case, recover the penalized coefficients
  if (high_dim) {
    bHat <- target + crossprod(X, slh) / lambda
  }

  # gHat is reserved for unpenalized covariates, which are not supported here.
  gHat <- numeric()

  # Return a named list for clarity
  return(list(unpenalized = gHat, penalized = bHat))
}
