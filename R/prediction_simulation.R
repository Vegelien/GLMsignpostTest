simulate_target_gain_logistic <- function(
    X_A, Y_A, X_B, Y_B,
    n_sims        = 100L,
    test_frac     = 0.5,
    lambda_A        = 2,
    lambda_B        = 2,
    theta_max     = 1,
    stratify      = TRUE,
    parallel      = FALSE,
    seed          = 1,
    return_linear_preds = FALSE,
    verbose       = TRUE,
    intercept = TRUE,
    standardize = TRUE,
    unpenalized_cols = NULL, 
    unpen_as_signpost= FALSE,
    ...
) {
  stopifnot(is.matrix(X_A), is.matrix(X_B))
  stopifnot(length(Y_A) == nrow(X_A), length(Y_B) == nrow(X_B))
  stopifnot(all(Y_A %in% c(0,1)), all(Y_B %in% c(0,1)))
  stopifnot(ncol(X_A) == ncol(X_B))
  
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for plain ridge.")
  }
  # --- map unpenalized columns ------------------------------------------------
  p <- ncol(X_B)
  if (is.null(unpenalized_cols) || length(unpenalized_cols) == 0L) {
    idx_unp <- integer(0)
  } else if (is.character(unpenalized_cols)) {
    if (is.null(colnames(X_B))) stop("unpenalized_cols are names but X has no colnames.")
    idx_unp <- match(unpenalized_cols, colnames(X_B))
    if (anyNA(idx_unp)) stop("Some unpenalized column names not found in X.")
  } else if (is.numeric(unpenalized_cols)) {
    idx_unp <- as.integer(unpenalized_cols)
    if (any(idx_unp < 1 | idx_unp > p)) stop("unpenalized_cols indices out of range.")
  } else {
    stop("unpenalized_cols must be NULL, character (names), or integer (indices).")
  }
  idx_unp <- sort(unique(idx_unp))
  print(idx_unp)
  idx_pen <- setdiff(seq_len(p), idx_unp)
  
  make_penalty_factor <- function(p, idx_unp) {
    pf <- rep(1, p)
    if (length(idx_unp)) pf[idx_unp] <- 0
    pf
  }
  # --- helpers (package-light, stable) ---------------------------------------
  eps <- 1e-15
  log_loss <- function(y, p) {
    p <- pmin(pmax(p, eps), 1 - eps)
    mean(-(y * log(p) + (1 - y) * log(1 - p)))
  }
  brier <- function(y, p) mean((p - y)^2)
  
  # AUC via Mann–Whitney U (fast, no deps)
  auc_fast <- function(y, p) {
    o <- order(p) # ascending
    y <- y[o]; p <- p[o]
    n1 <- sum(y == 1L); n0 <- sum(y == 0L)
    if (n1 == 0 || n0 == 0) return(NA_real_)
    # Rank sum for positives
    r <- rank(p, ties.method = "average")
    RS <- sum(r[y == 1L])
    (RS - n1 * (n1 + 1) / 2) / (n0 * n1)
  }
  
  # Average Precision (area under PR curve, step-wise)
  average_precision <- function(y, p) {
    o <- order(p, decreasing = TRUE)
    y <- y[o]; p <- p[o]
    tp <- cumsum(y == 1L)
    fp <- cumsum(y == 0L)
    n_pos <- sum(y == 1L)
    if (n_pos == 0) return(NA_real_)
    precision <- tp / (tp + fp)
    recall    <- tp / n_pos
    # AP = sum over threshold increments of precision * delta_recall
    # Handle duplicated scores implicitly (already monotone by cumsum)
    idx <- which(diff(c(0, recall)) > 0)
    sum(precision[idx] * diff(c(0, recall))[idx])
  }
  
  acc_at <- function(y, p, thr = 0.5) {
    mean((p >= thr) == (y == 1L))
  }
  
  # Calibration intercept & slope: glm(y ~ logit(p), binomial())
  cal_params <- function(y, p) {
    p <- pmin(pmax(p, eps), 1 - eps)
    lp <- qlogis(p)
    df <- data.frame(y = y, lp = lp)
    fit <- try(stats::glm(y ~ lp, data = df, family = stats::binomial()), silent = TRUE)
    if (inherits(fit, "try-error")) return(c(intercept = NA_real_, slope = NA_real_))
    co <- coef(fit)
    c(intercept = unname(co[1L]), slope = unname(co[2L]))
  }
  
  # Train/test split, optionally stratified
  split_once <- function(y, test_frac, stratify) {
    n <- length(y); n_test <- floor(test_frac * n)
    if (!stratify) {
      idx <- sample.int(n, n_test)
      list(test = idx, train = setdiff(seq_len(n), idx))
    } else {
      pos <- which(y == 1L); neg <- which(y == 0L)
      n_test_pos <- floor(length(pos) * test_frac)
      n_test_neg <- floor(length(neg) * test_frac)
      idx <- c(sample(pos, n_test_pos), sample(neg, n_test_neg))
      list(test = idx, train = setdiff(seq_len(n), idx))
    }
  }
  
  
  # --- precompute beta_a_hat on A with glmnet (unpenalized via penalty.factor=0)
  pf_A <- make_penalty_factor(ncol(X_A), idx_unp)
  fit_A <- glmnet::glmnet(
    x = X_A, y = Y_A,
    family = "binomial", alpha = 0,
    lambda = lambda_A,
    intercept = intercept, standardize = standardize,
    penalty.factor = pf_A
  )
  co_A <- as.numeric(stats::coef(fit_A, s = lambda_A))[-1L]  # drop intercept
  # target direction for PENALIZED block only, normalized on that block
  beta_a_hat_pen <- if(unpen_as_signpost){
    print("here")
    print(length(co_A))
    co_A
  }else if(length(idx_pen)) co_A[idx_pen] else numeric(0)
  if (length(beta_a_hat_pen)) {
    norm_pen <- sqrt(sum(beta_a_hat_pen^2))
    print(norm_pen)
    if (norm_pen > 0) beta_a_hat_pen <- beta_a_hat_pen / norm_pen
  }
  print(head(beta_a_hat_pen))

  
  
  
  # --- simulation loop --------------------------------------------------------
  #set.seed(seed)
  
  runner <- function(rep_id) {
    # Fresh split
    sp <- split_once(Y_B, test_frac = test_frac, stratify = stratify)
    tr <- sp$train; te <- sp$test
    
    Xtr <- X_B[tr, , drop = FALSE]
    Ytr <- Y_B[tr]
    Xte <- X_B[te, , drop = FALSE]
    Yte <- Y_B[te]
    
    ntr <- nrow(Xtr)
    
    
    # 1) Plain ridge on train (glmnet) with unpenalized via penalty.factor
    pf_B <- make_penalty_factor(ncol(Xtr), idx_unp)
    fit_plain <- glmnet::glmnet(
      x = Xtr, y = Ytr,
      family = "binomial", alpha = 0,
      lambda = lambda_B,
      intercept = intercept, standardize = standardize,
      penalty.factor = pf_B
    )
    beta0_hat <- as.numeric(stats::coef(fit_plain, s = lambda_B))  # (no intercept)

    # 2) theta_hat on train 
    if(unpen_as_signpost){
      Xtr_pen = Xtr
      Utr_unp = matrix(nrow = nrow(Xtr), ncol = 0)
      U_theta =  as.matrix(data.frame("intercept" = rep(1, nrow(Xtr))) )
    }else{
      Xtr_pen <- if (length(idx_pen)) Xtr[, idx_pen, drop = FALSE] else Xtr[, FALSE, drop = FALSE]
      Utr_unp <- if (length(idx_unp)) Xtr[, idx_unp, drop = FALSE] else matrix(nrow = nrow(Xtr), ncol = 0)
      
      U_theta = if (isTRUE(intercept)) cbind(`(Intercept)` = rep(1, nrow(Xtr)), Utr_unp) else Utr_unp
      U_theta = as.matrix(U_theta)
    }

    print(dim(U_theta))
    beta_00_pen <- if (length(beta_a_hat_pen)) rep(0, length(beta_a_hat_pen)) else numeric(0)
    print(dim(Xtr_pen))
    print(length(beta_00_pen))
    print(length(beta_a_hat_pen))
    
    theta_hat <- theta_inf_hat(
      Y = Ytr,
      X = Xtr_pen,
      U = U_theta,
      beta_0 = beta_00_pen,
      beta_a = beta_a_hat_pen,
      model  = "logistic",
      theta_max = theta_max
    )
    print(theta_hat)
    
    if(theta_hat < 1e-5) {
      beta_tgt_hat = beta0_hat
    }else{
      # 3) Target vector for penalized block
      target_vec_pen <- if (length(idx_pen)) {
         theta_hat * beta_a_hat_pen
      } else numeric(0)
      
      # 4) Targeted ridge on train (ridge_complete) with U_extra = unpenalized block
      fit_tgt <- ridge_complete(
        Y = Ytr, X = Xtr_pen,
        lambda = lambda_B,
        target = if (length(idx_pen)) target_vec_pen else NULL,
        model = "logistic",
        intercept = intercept, standardize = standardize,
        U_extra = Utr_unp,
        ...
      )
      # Reconstruct full-length coefficient vector: put gamma at unpenalized spots, beta at penalized spots
      beta_tgt_hat <- numeric(1 + ncol(Xtr))   # 1 for intercept
      beta_tgt_hat[1] <- fit_tgt$alpha
      if(unpen_as_signpost){
        beta_tgt_hat[-1] <- fit_tgt$beta
      }else{
        if (length(idx_pen))  beta_tgt_hat[1 + idx_pen] <- fit_tgt$beta
        if (length(idx_unp))  beta_tgt_hat[1 + idx_unp] <- fit_tgt$gamma
        
      }
    }
    

    # 5) Evaluate on test (paired metrics)
    lp_plain <- as.vector(cbind(1, Xte) %*% beta0_hat)
    lp_tgt   <- as.vector(cbind(1, Xte) %*% beta_tgt_hat)
    p_plain  <- plogis(lp_plain)
    p_tgt    <- plogis(lp_tgt)
    
    # metrics
    ll_plain <- log_loss(Yte, p_plain)
    ll_tgt   <- log_loss(Yte, p_tgt)
    
    br_plain <- brier(Yte, p_plain)
    br_tgt   <- brier(Yte, p_tgt)
    
    auc_p    <- auc_fast(Yte, p_plain)
    auc_t    <- auc_fast(Yte, p_tgt)
    
    ap_p     <- average_precision(Yte, p_plain)
    ap_t     <- average_precision(Yte, p_tgt)
    
    acc_p    <- acc_at(Yte, p_plain, 0.5)
    acc_t    <- acc_at(Yte, p_tgt, 0.5)
    
    cal_p    <- cal_params(Yte, p_plain)
    cal_t    <- cal_params(Yte, p_tgt)
    
    out <- data.frame(
      rep = rep_id,
      n_train = ntr, n_test = length(Yte),
      lambda = lambda_B,
      theta_hat = as.numeric(theta_hat),
      
      logloss_null = ll_plain,
      logloss_tgt  = ll_tgt,
      delta_logloss = ll_plain - ll_tgt,
      
      brier_null = br_plain,
      brier_tgt  = br_tgt,
      delta_brier = br_plain - br_tgt,
      
      auc_null = auc_p,
      auc_tgt  = auc_t,
      delta_auc = auc_t - auc_p,          # higher is better
      
      prauc_null = ap_p,
      prauc_tgt  = ap_t,
      delta_prauc = ap_t - ap_p,          # higher is better
      
      acc05_null = acc_p,
      acc05_tgt  = acc_t,
      delta_acc05 = acc_t - acc_p,        # higher is better
      
      cal_int_null = cal_p["intercept"],
      cal_slope_null = cal_p["slope"],
      cal_int_tgt = cal_t["intercept"],
      cal_slope_tgt = cal_t["slope"],
      
      stringsAsFactors = FALSE
    )
    
    if (return_linear_preds) {
      out$lp_null_mean <- mean(lp_plain)
      out$lp_tgt_mean  <- mean(lp_tgt)
      # You can also stash full vectors externally if desired
    }
    
    if (verbose && (rep_id %% max(1L, floor(n_sims / 10))) == 0) {
      message(sprintf("[rep %d/%d] Δlogloss=%.6f, ΔAUC=%.4f",
                      rep_id, n_sims, out$delta_logloss, out$delta_auc))
    }
    
    out
  }
  
  # choose apply engine
  use_future <- isTRUE(parallel) && requireNamespace("future.apply", quietly = TRUE)
  reps <- seq_len(n_sims)
  res_list <- if (use_future) {
    future.apply::future_lapply(reps, runner, future.seed = TRUE)
  } else {
    lapply(reps, runner)
  }
  
  do.call(rbind, res_list)
}


simulate_target_gain_logistic_LEAKAGE <- function(
    X_A, Y_A, X_B, Y_B,
    id_A = NULL, id_B = NULL,              # NEW: patient IDs
    leakage_guard = TRUE,                  # NEW: exclude A rows that appear in B-test for each sim
    n_sims        = 100L,
    test_frac     = 0.5,
    lambda_A      = 2,
    lambda_B      = 2,
    theta_max     = 1,
    stratify      = TRUE,
    parallel      = FALSE,
    seed          = 1,
    return_linear_preds = FALSE,
    verbose       = TRUE,
    intercept     = TRUE,
    standardize   = TRUE,
    unpenalized_cols = NULL,
    unpen_as_signpost = FALSE,
    ...
) {
  stopifnot(is.matrix(X_A), is.matrix(X_B))
  stopifnot(length(Y_A) == nrow(X_A), length(Y_B) == nrow(X_B))
  stopifnot(all(Y_A %in% c(0,1)), all(Y_B %in% c(0,1)))
  
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for plain ridge.")
  }
  
  # ---------- IDs (explicit first; else rownames) ----------
  if (is.null(id_A)) {
    id_A <- rownames(X_A)
    if (is.null(id_A)) {
      if (isTRUE(leakage_guard)) stop("Provide id_A or rownames(X_A) when leakage_guard=TRUE.")
      id_A <- seq_len(nrow(X_A))  # dummy IDs (only used if leakage_guard=FALSE)
      message("No id_A provided and no rownames(X_A); using sequence IDs (leakage guard off or ineffective).")
    } else {
      message("Using rownames(X_A) as id_A.")
    }
  }
  if (is.null(id_B)) {
    id_B <- rownames(X_B)
    if (is.null(id_B)) {
      if (isTRUE(leakage_guard)) stop("Provide id_B or rownames(X_B) when leakage_guard=TRUE.")
      id_B <- seq_len(nrow(X_B))
      message("No id_B provided and no rownames(X_B); using sequence IDs (leakage guard off or ineffective).")
    } else {
      message("Using rownames(X_B) as id_B.")
    }
  }
  stopifnot(length(id_A) == nrow(X_A), length(id_B) == nrow(X_B))
  
  # Optional sanity check: uniqueness within each set (warn only)
  if (anyDuplicated(id_A)) warning("id_A contains duplicate patient IDs.")
  if (anyDuplicated(id_B)) warning("id_B contains duplicate patient IDs.")
  
  # ---------- Align columns across A and B by name if available ----------
  if (!is.null(colnames(X_A)) && !is.null(colnames(X_B))) {
    if (!setequal(colnames(X_A), colnames(X_B))) {
      diffA <- setdiff(colnames(X_A), colnames(X_B))
      diffB <- setdiff(colnames(X_B), colnames(X_A))
      stop(sprintf("Column name mismatch between X_A and X_B.\nOnly in A: %s\nOnly in B: %s",
                   paste(diffA, collapse=", "), paste(diffB, collapse=", ")))
    }
    # Reorder B to match A (or vice-versa; choose one consistently)
    X_B <- X_B[, colnames(X_A), drop = FALSE]
  } else {
    # No names: we can only assume the columns already correspond
    if (ncol(X_A) != ncol(X_B)) stop("X_A and X_B must have same number of columns when no colnames present.")
  }
  
  # ---------- Map unpenalized columns after alignment ----------
  p <- ncol(X_B)
  if (is.null(unpenalized_cols) || length(unpenalized_cols) == 0L) {
    idx_unp <- integer(0)
  } else if (is.character(unpenalized_cols)) {
    if (is.null(colnames(X_B))) stop("unpenalized_cols are names but X has no colnames.")
    idx_unp <- match(unpenalized_cols, colnames(X_B))
    if (anyNA(idx_unp)) stop("Some unpenalized column names not found in X.")
  } else if (is.numeric(unpenalized_cols)) {
    idx_unp <- as.integer(unpenalized_cols)
    if (any(idx_unp < 1 | idx_unp > p)) stop("unpenalized_cols indices out of range.")
  } else {
    stop("unpenalized_cols must be NULL, character (names), or integer (indices).")
  }
  idx_unp <- sort(unique(idx_unp))
  idx_pen <- setdiff(seq_len(p), idx_unp)
  
  make_penalty_factor <- function(p, idx_unp) {
    pf <- rep(1, p); if (length(idx_unp)) pf[idx_unp] <- 0; pf
  }
  
  # ---------- Metrics helpers ----------
  eps <- 1e-15
  log_loss <- function(y, pr) { pr <- pmin(pmax(pr, eps), 1 - eps); mean(-(y*log(pr) + (1-y)*log(1-pr))) }
  brier   <- function(y, pr) mean((pr - y)^2)
  auc_fast <- function(y, s) {
    o <- order(s); y <- y[o]; n1 <- sum(y==1L); n0 <- sum(y==0L)
    if (n1==0 || n0==0) return(NA_real_)
    r <- rank(s[o], ties.method="average"); RS <- sum(r[y==1L])
    (RS - n1*(n1+1)/2) / (n0*n1)
  }
  average_precision <- function(y, s) {
    o <- order(s, decreasing=TRUE); y <- y[o]
    tp <- cumsum(y==1L); fp <- cumsum(y==0L); n_pos <- sum(y==1L)
    if (n_pos==0) return(NA_real_)
    precision <- tp / (tp + fp); recall <- tp / n_pos
    idx <- which(diff(c(0, recall)) > 0)
    sum(precision[idx] * diff(c(0, recall))[idx])
  }
  acc_at <- function(y, pr, thr=0.5) mean((pr >= thr) == (y==1L))
  cal_params <- function(y, pr) {
    pr <- pmin(pmax(pr, eps), 1 - eps); lp <- qlogis(pr)
    fit <- try(stats::glm(y ~ lp, family = stats::binomial()), silent = TRUE)
    if (inherits(fit, "try-error")) return(c(intercept=NA_real_, slope=NA_real_))
    co <- coef(fit); c(intercept = unname(co[1L]), slope = unname(co[2L]))
  }
  
  # ---------- B-driven split ----------
  split_once <- function(y, test_frac, stratify) {
    n <- length(y); n_test <- floor(test_frac * n)
    if (!stratify) {
      idx <- sample.int(n, n_test)
      list(test=idx, train=setdiff(seq_len(n), idx))
    } else {
      pos <- which(y==1L); neg <- which(y==0L)
      n_test_pos <- floor(length(pos)*test_frac); n_test_neg <- floor(length(neg)*test_frac)
      idx <- c(sample(pos, n_test_pos), sample(neg, n_test_neg))
      list(test=idx, train=setdiff(seq_len(n), idx))
    }
  }
  
  # ---------- Simulation loop ----------
  #set.seed(seed)
  use_future <- isTRUE(parallel) && requireNamespace("future.apply", quietly = TRUE)
  reps <- seq_len(n_sims)
  
  runner <- function(rep_id) {
    # Split on B
    sp <- split_once(Y_B, test_frac = test_frac, stratify = stratify)
    tr_B <- sp$train; te_B <- sp$test
    ids_teB <- id_B[te_B]
    
    # Compose A-training indices for beta_a (exclude B-test overlaps if requested)
    if (isTRUE(leakage_guard)) {
      tr_A_for_beta <- which(!(id_A %in% ids_teB))
    } else {
      tr_A_for_beta <- seq_len(nrow(X_A))
    }
    print("length of A:")
    print(length(tr_A_for_beta))
    
    # If A becomes empty after exclusion, skip this repetition gracefully
    if (length(tr_A_for_beta) == 0L) {
      if (verbose) message(sprintf("[rep %d/%d] Skipped: no A rows left after excluding B-test overlaps.", rep_id, n_sims))
      return(data.frame(
        rep = rep_id, n_train = length(tr_B), n_test = length(te_B),
        lambda = lambda_B, theta_hat = NA_real_,
        logloss_null = NA_real_, logloss_tgt = NA_real_, delta_logloss = NA_real_,
        brier_null = NA_real_, brier_tgt = NA_real_, delta_brier = NA_real_,
        auc_null = NA_real_, auc_tgt = NA_real_, delta_auc = NA_real_,
        prauc_null = NA_real_, prauc_tgt = NA_real_, delta_prauc = NA_real_,
        acc05_null = NA_real_, acc05_tgt = NA_real_, delta_acc05 = NA_real_,
        cal_int_null = NA_real_, cal_slope_null = NA_real_,
        cal_int_tgt = NA_real_,  cal_slope_tgt = NA_real_,
        skipped = TRUE,
        stringsAsFactors = FALSE
      ))
    }
    
    # ---- Fit beta_a on A (minus B-test overlaps), extract penalized-direction ----
    pf_A <- make_penalty_factor(ncol(X_A), idx_unp)
    fit_A <- glmnet::glmnet(
      x = X_A[tr_A_for_beta, , drop=FALSE], y = Y_A[tr_A_for_beta],
      family = "binomial", alpha = 0, lambda = lambda_A,
      intercept = intercept, standardize = standardize,
      penalty.factor = pf_A
    )
    co_A <- as.numeric(stats::coef(fit_A, s = lambda_A))[-1L]  # drop intercept
    beta_a_hat_pen <- if (unpen_as_signpost) {
      co_A
    } else if (length(idx_pen)) {
      co_A[idx_pen]
    } else numeric(0)
    
    if (length(beta_a_hat_pen)) {
      norm_pen <- sqrt(sum(beta_a_hat_pen^2))
      if (is.finite(norm_pen) && norm_pen > 0) beta_a_hat_pen <- beta_a_hat_pen / norm_pen
    }
    
    # ---- Prepare B train/test matrices ----
    Xtr <- X_B[tr_B, , drop = FALSE]; Ytr <- Y_B[tr_B]
    Xte <- X_B[te_B, , drop = FALSE]; Yte <- Y_B[te_B]
    ntr <- nrow(Xtr)
    
    # ---- Null model on B-train ----
    pf_B <- make_penalty_factor(ncol(Xtr), idx_unp)
    fit_plain <- glmnet::glmnet(
      x = Xtr, y = Ytr,
      family = "binomial", alpha = 0, lambda = lambda_B,
      intercept = intercept, standardize = standardize,
      penalty.factor = pf_B
    )
    beta0_hat <- as.numeric(stats::coef(fit_plain, s = lambda_B))  # includes intercept as first element
    
    # ---- Theta and targeted fit on B-train ----
    if (unpen_as_signpost) {
      Xtr_pen <- Xtr
      Utr_unp <- matrix(nrow=nrow(Xtr), ncol=0)
      U_theta <- data.frame(`(Intercept)` = rep(1, nrow(Xtr)))
      U_theta <- as.matrix(U_theta)
    } else {
      Xtr_pen <- if (length(idx_pen)) Xtr[, idx_pen, drop = FALSE] else Xtr[, FALSE, drop = FALSE]
      Utr_unp <- if (length(idx_unp)) Xtr[, idx_unp, drop = FALSE] else matrix(nrow = nrow(Xtr), ncol = 0)
      U_theta <- if (isTRUE(intercept)) cbind(`(Intercept)` = rep(1, nrow(Xtr)), Utr_unp) else Utr_unp
      U_theta <- as.matrix(U_theta)
    }
    
    beta_00_pen <- if (length(beta_a_hat_pen)) rep(0, length(beta_a_hat_pen)) else numeric(0)
    
    theta_hat <- theta_inf_hat(
      Y = Ytr,
      X = Xtr_pen,
      U = U_theta,
      beta_0 = beta_00_pen,
      beta_a = beta_a_hat_pen,
      model  = "logistic",
      theta_max = theta_max
    )
    print(theta_hat)
    # Targeted vs null selection
    if (!is.finite(theta_hat) || theta_hat < 1e-5) {
      beta_tgt_hat <- beta0_hat
    } else {
      target_vec_pen <- if (length(idx_pen)) theta_hat * beta_a_hat_pen else numeric(0)
      fit_tgt <- ridge_complete(
        Y = Ytr, X = Xtr_pen,
        lambda = lambda_B,
        target = if (length(idx_pen)) target_vec_pen else NULL,
        model = "logistic",
        intercept = intercept, standardize = standardize,
        U_extra = Utr_unp,
        ...
      )
      # Reconstruct full-length coef vector: (Intercept) + all columns in original order
      beta_tgt_hat <- numeric(1 + ncol(Xtr))
      beta_tgt_hat[1] <- fit_tgt$alpha
      if (unpen_as_signpost) {
        beta_tgt_hat[-1] <- fit_tgt$beta
      } else {
        if (length(idx_pen)) beta_tgt_hat[1 + idx_pen] <- fit_tgt$beta
        if (length(idx_unp)) beta_tgt_hat[1 + idx_unp] <- fit_tgt$gamma
      }
    }
    
    # ---- Evaluate on B-test ----
    lp_plain <- as.vector(cbind(1, Xte) %*% beta0_hat)
    lp_tgt   <- as.vector(cbind(1, Xte) %*% beta_tgt_hat)
    p_plain  <- plogis(lp_plain)
    p_tgt    <- plogis(lp_tgt)
    
    ll_plain <- log_loss(Yte, p_plain); ll_tgt <- log_loss(Yte, p_tgt)
    br_plain <- brier(Yte, p_plain);    br_tgt <- brier(Yte, p_tgt)
    auc_p <- auc_fast(Yte, p_plain);    auc_t <- auc_fast(Yte, p_tgt)
    ap_p  <- average_precision(Yte, p_plain); ap_t <- average_precision(Yte, p_tgt)
    acc_p <- acc_at(Yte, p_plain, 0.5);       acc_t <- acc_at(Yte, p_tgt, 0.5)
    cal_p <- cal_params(Yte, p_plain);        cal_t <- cal_params(Yte, p_tgt)
    
    out <- data.frame(
      rep = rep_id,
      n_train = ntr, n_test = length(Yte),
      lambda = lambda_B,
      theta_hat = as.numeric(theta_hat),
      logloss_null = ll_plain, logloss_tgt = ll_tgt, delta_logloss = ll_plain - ll_tgt,
      brier_null = br_plain, brier_tgt = br_tgt, delta_brier = br_plain - br_tgt,
      auc_null = auc_p, auc_tgt = auc_t, delta_auc = auc_t - auc_p,
      prauc_null = ap_p, prauc_tgt = ap_t, delta_prauc = ap_t - ap_p,
      acc05_null = acc_p, acc05_tgt = acc_t, delta_acc05 = acc_t - acc_p,
      cal_int_null = cal_p["intercept"], cal_slope_null = cal_p["slope"],
      cal_int_tgt  = cal_t["intercept"], cal_slope_tgt  = cal_t["slope"],
      skipped = FALSE,
      stringsAsFactors = FALSE
    )
    
    if (return_linear_preds) {
      out$lp_null_mean <- mean(lp_plain)
      out$lp_tgt_mean  <- mean(lp_tgt)
    }
    
    if (verbose && (rep_id %% max(1L, floor(n_sims / 10))) == 0) {
      message(sprintf("[rep %d/%d] Δlogloss=%.6f, ΔAUC=%.4f",
                      rep_id, n_sims, out$delta_logloss, out$delta_auc))
    }
    out
  }
  
  res_list <- if (use_future) {
    future.apply::future_lapply(reps, runner, future.seed = TRUE)
  } else {
    lapply(reps, runner)
  }
  do.call(rbind, res_list)
}


summarize_target_gain <- function(
    res,
    metrics = c("logloss","brier","auc","prauc","acc05"),
    conf_level = 0.95,
    use_ggplot = FALSE,
    file = NULL,
    width = 8, height = 6,
    show_calibration = TRUE
) {
  stopifnot(is.data.frame(res))
  # Check presence of expected columns
  needed <- unlist(lapply(metrics, function(m) {
    paste0(c("", "_null","_tgt","delta_"), if (m=="logloss") "logloss" else m)
  }))
  # Only deltas & paired columns matter; be lenient
  deltas <- paste0("delta_", ifelse(metrics=="logloss","logloss",metrics))
  missing <- setdiff(deltas, names(res))
  if (length(missing)) {
    stop("Result is missing expected columns: ", paste(missing, collapse = ", "))
  }
  
  # ----- helpers -----
  ci_mean <- function(x, level = 0.95) {
    x <- x[is.finite(x)]
    n <- length(x)
    m <- mean(x)
    s <- stats::sd(x)
    if (n <= 1 || is.na(s)) return(c(lower = NA_real_, upper = NA_real_))
    alpha <- 1 - level
    tcrit <- stats::qt(1 - alpha/2, df = n - 1)
    half <- tcrit * s / sqrt(n)
    c(lower = m - half, upper = m + half)
  }
  
  # ----- paired summaries for deltas -----
  make_row <- function(m) {
    core <- if (m == "logloss") "logloss" else m
    dname <- paste0("delta_", core)
    x <- res[[dname]]
    x <- x[is.finite(x)]
    n <- length(x)
    m_hat <- mean(x)
    med  <- stats::median(x)
    sdv  <- stats::sd(x)
    ap   <- mean(x > 0)    # share of reps where targeting "wins"
    ci   <- ci_mean(x, conf_level)
    data.frame(
      metric = core,
      n = n,
      mean_delta = m_hat,
      median_delta = med,
      sd_delta = sdv,
      frac_positive = ap,
      ci_lower = ci[1], ci_upper = ci[2],
      stringsAsFactors = FALSE
    )
  }
  summary_df <- do.call(rbind, lapply(metrics, make_row))
  
  # ----- (optional) calibration summary -----
  calib_summary <- NULL
  if (isTRUE(show_calibration)) {
    # slope ideally ~1, intercept ~0. Summarize both arms and deltas.
    for (arm in c("null","tgt")) {
      ci_slope <- ci_mean(res[[paste0("cal_slope_", arm)]], conf_level)
      ci_int   <- ci_mean(res[[paste0("cal_int_", arm)]], conf_level)
      tmp <- data.frame(
        what = c(paste0("cal_slope_", arm), paste0("cal_int_", arm)),
        n = sum(is.finite(res[[paste0("cal_slope_", arm)]])),
        mean = c(mean(res[[paste0("cal_slope_", arm)]], na.rm=TRUE),
                 mean(res[[paste0("cal_int_", arm)]],   na.rm=TRUE)),
        sd   = c(sd(res[[paste0("cal_slope_", arm)]], na.rm=TRUE),
                 sd(res[[paste0("cal_int_", arm)]],   na.rm=TRUE)),
        ci_lower = c(ci_slope[1], ci_int[1]),
        ci_upper = c(ci_slope[2], ci_int[2]),
        stringsAsFactors = FALSE
      )
      calib_summary <- rbind(calib_summary, tmp)
    }
    # deltas (tgt - null for slopes/intercepts)
    res$delta_cal_slope <- res$cal_slope_tgt - res$cal_slope_null
    res$delta_cal_int   <- res$cal_int_tgt   - res$cal_int_null
    ci_ds <- ci_mean(res$delta_cal_slope, conf_level)
    ci_di <- ci_mean(res$delta_cal_int,   conf_level)
    calib_summary <- rbind(
      calib_summary,
      data.frame(
        what = c("delta_cal_slope (tgt-null)","delta_cal_int (tgt-null)"),
        n = c(sum(is.finite(res$delta_cal_slope)), sum(is.finite(res$delta_cal_int))),
        mean = c(mean(res$delta_cal_slope, na.rm=TRUE), mean(res$delta_cal_int, na.rm=TRUE)),
        sd   = c(sd(res$delta_cal_slope, na.rm=TRUE),   sd(res$delta_cal_int,   na.rm=TRUE)),
        ci_lower = c(ci_ds[1], ci_di[1]),
        ci_upper = c(ci_ds[2], ci_di[2]),
        stringsAsFactors = FALSE
      )
    )
  }
  
  # ----- plots -----
  plot_handles <- list()
  
  # device handling (optional save)
  if (!is.null(file)) {
    ext <- tolower(tools::file_ext(file))
    if (ext %in% c("pdf","png")) {
      if (ext == "pdf") grDevices::pdf(file, width = width, height = height)
      if (ext == "png") grDevices::png(file, width = width, height = height, units = "in", res = 150)
      on.exit(grDevices::dev.off(), add = TRUE)
    } else {
      warning("Unknown file extension; plotting to current device instead.")
    }
  }
  
  # Plot set A: distributions of deltas (one panel per metric)
  if (!use_ggplot) {
    old_par <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(old_par), add = TRUE)
    nM <- length(metrics); nrow <- ceiling(sqrt(nM)); ncol <- ceiling(nM / nrow)
    graphics::par(mfrow = c(nrow, ncol), mar = c(4,4,2,1))
    for (m in metrics) {
      core <- if (m=="logloss") "logloss" else m
      d <- res[[paste0("delta_", core)]]
      d <- d[is.finite(d)]
      if (!length(d)) { plot.new(); title(main = paste("delta", core)); next }
      rng <- range(d, 0, na.rm = TRUE)
      h <- graphics::hist(d, breaks = "FD", col = "grey85", border = "white",
                          main = paste0("Δ ", core, " (null - tgt)"),
                          xlab = "paired delta", yaxt = "n")
      graphics::axis(2, las = 1)
      graphics::abline(v = 0, lty = 2)
      graphics::abline(v = mean(d), lwd = 2)
      # 95% CI
      ci <- ci_mean(d, conf_level)
      graphics::segments(ci[1], 0, ci[2], 0, lwd = 4, lend = 1)
      graphics::legend("topright",
                       legend = c(sprintf("mean = %.4g", mean(d)),
                                  sprintf("%d%% CI", as.integer(conf_level*100))),
                       bty = "n")
    }
  } else if (use_ggplot && requireNamespace("ggplot2", quietly = TRUE)) {
    gg <- lapply(metrics, function(m) {
      core <- if (m=="logloss") "logloss" else m
      dname <- paste0("delta_", core)
      df <- data.frame(delta = res[[dname]])
      df <- df[is.finite(df$delta), , drop = FALSE]
      ci <- ci_mean(df$delta, conf_level)
      ggplot2::ggplot(df, ggplot2::aes(x = delta)) +
        ggplot2::geom_histogram(bins = max(10, floor(sqrt(nrow(df)))), fill = "grey80", color = "white") +
        ggplot2::geom_vline(xintercept = 0, linetype = 2) +
        ggplot2::geom_vline(xintercept = mean(df$delta), linewidth = 0.6) +
        ggplot2::geom_segment(x = ci[1], xend = ci[2], y = 0, yend = 0, linewidth = 1.2) +
        ggplot2::labs(title = paste0("Δ ", core, " (null - tgt)"),
                      x = "paired delta", y = "count")
    })
    names(gg) <- paste0("hist_delta_", ifelse(metrics=="logloss","logloss",metrics))
    plot_handles <- c(plot_handles, gg)
    # Print to device if not saving; otherwise user can print them later
    if (is.null(file)) for (g in plot_handles) print(g)
  } else {
    warning("ggplot2 not available; falling back to base plots.")
    return(summarize_target_gain(res, metrics, conf_level, use_ggplot=FALSE, file=file,
                                 width=width, height=height, show_calibration=show_calibration))
  }
  
  # Plot set B: paired scatter (null vs. targeted) with 45° line (first 3 metrics only for space)
  m_show <- head(metrics, 3)
  if (!use_ggplot) {
    graphics::par(mfrow = c(1, length(m_show)), mar = c(4,4,2,1))
    for (m in m_show) {
      core <- if (m=="logloss") "logloss" else m
      x <- res[[paste0(core, "_null")]]
      y <- res[[paste0(core, "_tgt")]]
      ok <- is.finite(x) & is.finite(y)
      x <- x[ok]; y <- y[ok]
      if (!length(x)) { plot.new(); title(main = paste(core, "(paired)")); next }
      lim <- range(c(x, y))
      graphics::plot(x, y, xlab = paste(core, "null"), ylab = paste(core, "targeted"),
                     main = paste(core, "(paired)"), xlim = lim, ylim = lim, pch = 16, col = "#00000040")
      graphics::abline(0, 1, lty = 2)
    }
  } else if (use_ggplot) {
    gg2 <- lapply(m_show, function(m) {
      core <- if (m=="logloss") "logloss" else m
      df <- data.frame(null = res[[paste0(core, "_null")]],
                       tgt  = res[[paste0(core, "_tgt")]])
      df <- df[is.finite(df$null) & is.finite(df$tgt), , drop = FALSE]
      lim <- range(c(df$null, df$tgt))
      ggplot2::ggplot(df, ggplot2::aes(null, tgt)) +
        ggplot2::geom_point(alpha = 0.4) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
        ggplot2::coord_equal(xlim = lim, ylim = lim, expand = TRUE) +
        ggplot2::labs(title = paste0(core, " (paired)"),
                      x = paste(core, "null"), y = paste(core, "targeted"))
    })
    names(gg2) <- paste0("scatter_", ifelse(m_show=="logloss","logloss",m_show))
    plot_handles <- c(plot_handles, gg2)
    if (is.null(file)) for (g in gg2) print(g)
  }
  
  # Optional: quick calibration comparison scatter (slopes)
  if (isTRUE(show_calibration)) {
    if (!use_ggplot) {
      graphics::par(mfrow = c(1,2), mar = c(4,4,2,1))
      # Slope
      s_ok <- is.finite(res$cal_slope_null) & is.finite(res$cal_slope_tgt)
      sl <- range(c(res$cal_slope_null[s_ok], res$cal_slope_tgt[s_ok], 1))
      graphics::plot(res$cal_slope_null[s_ok], res$cal_slope_tgt[s_ok],
                     xlab = "cal slope null", ylab = "cal slope targeted",
                     main = "Calibration slope", pch = 16, col = "#00000040",
                     xlim = sl, ylim = sl)
      graphics::abline(0, 1, lty = 2); graphics::abline(v = 1, h = 1, col = "#00000040")
      # Intercept
      i_ok <- is.finite(res$cal_int_null) & is.finite(res$cal_int_tgt)
      il <- range(c(res$cal_int_null[i_ok], res$cal_int_tgt[i_ok], 0))
      graphics::plot(res$cal_int_null[i_ok], res$cal_int_tgt[i_ok],
                     xlab = "cal intercept null", ylab = "cal intercept targeted",
                     main = "Calibration intercept", pch = 16, col = "#00000040",
                     xlim = il, ylim = il)
      graphics::abline(0, 1, lty = 2); graphics::abline(v = 0, h = 0, col = "#00000040")
    } else {
      if (requireNamespace("ggplot2", quietly = TRUE)) {
        df_s <- data.frame(null = res$cal_slope_null, tgt = res$cal_slope_tgt)
        df_s <- df_s[is.finite(df_s$null) & is.finite(df_s$tgt), , drop = FALSE]
        sl <- range(c(df_s$null, df_s$tgt, 1))
        p_s <- ggplot2::ggplot(df_s, ggplot2::aes(null, tgt)) +
          ggplot2::geom_point(alpha = 0.4) +
          ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
          ggplot2::geom_vline(xintercept = 1, linetype = 3) +
          ggplot2::geom_hline(yintercept = 1, linetype = 3) +
          ggplot2::coord_equal(xlim = sl, ylim = sl, expand = TRUE) +
          ggplot2::labs(title = "Calibration slope", x = "null", y = "targeted")
        df_i <- data.frame(null = res$cal_int_null, tgt = res$cal_int_tgt)
        df_i <- df_i[is.finite(df_i$null) & is.finite(df_i$tgt), , drop = FALSE]
        il <- range(c(df_i$null, df_i$tgt, 0))
        p_i <- ggplot2::ggplot(df_i, ggplot2::aes(null, tgt)) +
          ggplot2::geom_point(alpha = 0.4) +
          ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
          ggplot2::geom_vline(xintercept = 0, linetype = 3) +
          ggplot2::geom_hline(yintercept = 0, linetype = 3) +
          ggplot2::coord_equal(xlim = il, ylim = il, expand = TRUE) +
          ggplot2::labs(title = "Calibration intercept", x = "null", y = "targeted")
        plot_handles$cal_slope <- p_s
        plot_handles$cal_int   <- p_i
        if (is.null(file)) { print(p_s); print(p_i) }
      }
    }
  }
  
  out <- list(
    summary = summary_df,
    calibration_summary = calib_summary,
    by_rep = res,
    plots = if (length(plot_handles)) plot_handles else NULL
  )
  return(out)
}



#' Split-Based Empirical Validity with Zero Overlap (intercept + unpenalized support)
#'
#' Splits A and B in disjoint halves per bootstrap. First halves estimate
#' the targets (β_0 from A, β_a from B). Second halves are used to form
#' mixed test sets at proportions in `replacement_props`, and θ̂ is obtained
#' by solving the θ estimating equation with U = [Intercept, X_unpen].
#'
#' @param Y_A,Y_B Numeric {0,1} outcomes for datasets A (reference) and B (alternative)
#' @param X_A,X_B Numeric matrices with identical columns (same p)
#' @param lambda Numeric test penalty (glmnet scale). Use Inf for θ̂ at λ = ∞.
#' @param lambda_beta Numeric penalty for target estimation (glmnet scale)
#' @param model "logistic" | "poisson" | "linear"
#' @param n_bootstrap Integer, number of random disjoint 50/50 splits
#' @param n_reps Integer, repetitions per replacement proportion per split
#' @param replacement_props numeric vector in [0,1]
#' @param scaled Logical; if TRUE, normalize the target direction on the **penalized block only**
#' @param verbose Logical
#' @param theta_min,theta_max numeric endpoints for θ search
#' @param n_cores Integer, parallel workers
#' @param output_file NULL or path; if not NULL, writes a single CSV at the end
#' @param enhancement FALSE or numeric factor (bootstrap within disjoint halves)
#' @param intercept Logical; if TRUE include an unpenalized intercept in U (θ step)
#' @param unpenalized_cols NULL | integer | character, columns treated as unpenalized in θ step
#' @param unpen_as_signpost Logical; if TRUE, treat all columns as “signpost” (no U except intercept); mirrors simulator.
#'
#' @return data.frame with columns: split_id, prop_B, rep_id, theta_hat
#' @export
signpost_empirical_validity_split <- function(
    Y_A, X_A, Y_B, X_B,
    lambda, lambda_beta,
    model = "logistic",
    n_bootstrap = 50L,
    n_reps = 10L,
    replacement_props = seq(0, 1, by = 0.1),
    scaled = TRUE,
    verbose = TRUE,
    theta_min = 0, theta_max = 1,
    n_cores = 8L,
    output_file = NULL,
    enhancement = FALSE,
    intercept = TRUE,
    unpenalized_cols = NULL,
    unpen_as_signpost = FALSE
) {
  cat("LETSGOOO")
  # ---- basic checks ---------------------------------------------------------
  stopifnot(length(Y_A) == nrow(X_A), length(Y_B) == nrow(X_B))
  stopifnot(ncol(X_A) == ncol(X_B))
  n_A <- nrow(X_A); n_B <- nrow(X_B); p <- ncol(X_A)
  
  if (is.numeric(enhancement) && enhancement <= 0) {
    stop("enhancement must be FALSE or a positive number")
  }
  
  if (verbose) {
    cat("=== SPLIT-BASED EMPIRICAL VALIDITY (ZERO OVERLAP) ===\n")
    cat("Dataset A:", n_A, "samples | Dataset B:", n_B, "samples | p =", p, "\n")
  }
  
  # ---- map unpenalized columns (mirrors simulator) --------------------------
  # unpenalized_cols: NULL | indices | names
  if (is.null(unpenalized_cols) || length(unpenalized_cols) == 0L) {
    idx_unp <- integer(0)
  } else if (is.character(unpenalized_cols)) {
    if (is.null(colnames(X_A))) stop("unpenalized_cols are names but X has no colnames.")
    idx_unp <- match(unpenalized_cols, colnames(X_A))
    if (anyNA(idx_unp)) stop("Some unpenalized column names not found in X.")
  } else if (is.numeric(unpenalized_cols)) {
    idx_unp <- as.integer(unpenalized_cols)
    if (any(idx_unp < 1 | idx_unp > p)) stop("unpenalized_cols indices out of range.")
  } else {
    stop("unpenalized_cols must be NULL, character (names), or integer (indices).")
  }
  idx_unp <- sort(unique(idx_unp))
  idx_pen <- setdiff(seq_len(p), idx_unp)
  
  make_penalty_factor <- function(p, idx_unp) {
    pf <- rep(1, p); if (length(idx_unp)) pf[idx_unp] <- 0; pf
  }
  
  # glmnet family mapping
  fam <- switch(model,
                logistic = "binomial",
                linear   = "gaussian",
                poisson  = "poisson",
                stop("model must be 'logistic','poisson','linear'"))
  
  # sizes for 50/50 disjoint halves (use the smaller to match test sizes)
  min_size    <- min(n_A, n_B)
  test_size   <- floor(min_size / 2)
  target_size_A <- n_A - test_size
  target_size_B <- n_B - test_size
  
  if (verbose) {
    cat("Disjoint split sizes:\n")
    cat("  A: target", target_size_A, "| test", test_size, "\n")
    cat("  B: target", target_size_B, "| test", test_size, "\n")
    if (is.numeric(enhancement)) {
      cat("Bootstrap enhancement factor:", enhancement, "\n")
    } else {
      cat("No bootstrap enhancement (pure 50/50 split)\n")
    }
    cat("Workers:", n_cores, "| Splits:", n_bootstrap, "\n")
  }
  
  # ---- worker ---------------------------------------------------------------
  worker_function <- function(split_id) {
    cat("LETSGOOO")
    set.seed(split_id * 100000L)
    
    # 1) Disjoint random halves
    A_idx <- sample.int(n_A)
    B_idx <- sample.int(n_B)
    A_target_idx <- A_idx[ 1:target_size_A ]
    A_test_idx   <- A_idx[ (target_size_A + 1):(target_size_A + test_size) ]
    B_target_idx <- B_idx[ 1:target_size_B ]
    B_test_idx   <- B_idx[ (target_size_B + 1):(target_size_B + test_size) ]
    
    # Optional bootstrap enhancement within halves
    if (is.numeric(enhancement)) {
      tsA <- round(enhancement * target_size_A)
      tsB <- round(enhancement * target_size_B)
      tsT <- round(enhancement * test_size)
      A_target_idx <- sample(A_target_idx, tsA, replace = TRUE)
      B_target_idx <- sample(B_target_idx, tsB, replace = TRUE)
      A_test_idx   <- sample(A_test_idx,   tsT, replace = TRUE)
      B_test_idx   <- sample(B_test_idx,   tsT, replace = TRUE)
    }
    
    # Split design/response
    Y_A_tg <- Y_A[A_target_idx]; X_A_tg <- X_A[A_target_idx, , drop = FALSE]
    Y_B_tg <- Y_B[B_target_idx]; X_B_tg <- X_B[B_target_idx, , drop = FALSE]
    Y_A_te <- Y_A[A_test_idx];   X_A_te <- X_A[A_test_idx,   , drop = FALSE]
    Y_B_te <- Y_B[B_test_idx];   X_B_te <- X_B[B_test_idx,   , drop = FALSE]
    
    # 2) Estimate target directions on target halves (glmnet ridge, no standardization)
    #    Penalty factors drop the unpenalized block if present.
    pf_A <- make_penalty_factor(p, idx_unp)
    pf_B <- make_penalty_factor(p, idx_unp)
    
    fit_A <- glmnet::glmnet(x = X_A_tg, y = Y_A_tg,
                            family = fam, alpha = 0,
                            lambda = lambda_beta,
                            intercept = intercept, standardize = FALSE,
                            penalty.factor = pf_A)
    co_A <- as.numeric(stats::coef(fit_A, s = lambda_beta))[-1L]  # drop intercept
    
    fit_B <- glmnet::glmnet(x = X_B_tg, y = Y_B_tg,
                            family = fam, alpha = 0,
                            lambda = lambda_beta,
                            intercept = intercept, standardize = FALSE,
                            penalty.factor = pf_B)
    co_B <- as.numeric(stats::coef(fit_B, s = lambda_beta))[-1L]  # drop intercept
    
    # 3) Build penalized targets for θ (mirror simulator behavior)
    beta_0_pen <- if (unpen_as_signpost) co_A else if (length(idx_pen)) co_A[idx_pen] else numeric(0)
    beta_a_pen <- if (unpen_as_signpost) co_B else if (length(idx_pen)) co_B[idx_pen] else numeric(0)
    
    # Scale only the penalized block when requested
    if (isTRUE(scaled)) {
      n0 <- sqrt(sum(beta_0_pen^2)); if (is.finite(n0) && n0 > 0) beta_0_pen <- beta_0_pen / n0
      nA <- sqrt(sum(beta_a_pen^2)); if (is.finite(nA) && nA > 0) beta_a_pen <- beta_a_pen / nA
    }
    
    # 4) Loop over replacement proportions and reps; compute θ̂
    out_rows <- vector("list", length(replacement_props) * n_reps)
    ridx <- 1L
    
    for (prop_B in replacement_props) {
      base_n <- length(Y_A_te)          # == length(Y_B_te)
      n_from_A <- round((1 - prop_B) * base_n)
      n_from_B <- round(prop_B * base_n)
      
      for (rep in seq_len(n_reps)) {
        set.seed(split_id * 100000L + rep * 1000L + which(replacement_props == prop_B))
        
        if (n_from_A == 0L) {
          selB <- sample.int(base_n, base_n, replace = FALSE)
          Y_mix <- Y_B_te[selB]; X_mix <- X_B_te[selB, , drop = FALSE]
        } else if (n_from_B == 0L) {
          Y_mix <- Y_A_te;      X_mix <- X_A_te
        } else {
          selA <- sample.int(base_n, n_from_A, replace = FALSE)
          selB <- sample.int(base_n, n_from_B, replace = FALSE)
          Y_mix <- c(Y_A_te[selA], Y_B_te[selB])
          X_mix <- rbind(X_A_te[selA, , drop = FALSE], X_B_te[selB, , drop = FALSE])
        }
        
        # Build (X_pen, U) for θ (mirrors simulate_target_gain_logistic)
        if (unpen_as_signpost) {
          X_pen  <- X_mix
          U_unp  <- matrix(nrow = nrow(X_mix), ncol = 0)
          U      <- if (isTRUE(intercept)) cbind(`(Intercept)` = 1.0, U_unp) else U_unp
          b0_pen <- beta_0_pen
          ba_pen <- beta_a_pen
        } else {
          X_pen  <- if (length(idx_pen)) X_mix[, idx_pen, drop = FALSE] else X_mix[, FALSE, drop = FALSE]
          U_unp  <- if (length(idx_unp)) X_mix[, idx_unp, drop = FALSE] else matrix(nrow = nrow(X_mix), ncol = 0)
          U      <- if (isTRUE(intercept)) cbind(`(Intercept)` = 1.0, U_unp) else U_unp
          b0_pen <- beta_0_pen
          ba_pen <- beta_a_pen
        }
        
        theta_hat <- tryCatch({
          
          if (is.infinite(lambda)) {
            theta_inf_hat(Y_mix, X_pen, U, b0_pen, ba_pen, model, theta_min = theta_min, theta_max = theta_max)
          } else {
            theta_hat_lambda(Y_mix, X_pen, b0_pen, ba_pen, lambda, model)
          }
        }, error = function(e) {
          if (isTRUE(verbose)) message(sprintf(
            "[split %d, prop=%.2f, rep %d] θ failed: %s",
            split_id, prop_B, rep, conditionMessage(e)
          ))
          NA_real_
        })
        
        # theta_hat <- tryCatch({
        #   if (is.infinite(lambda)) {
        #     theta_inf_hat(Y_mix, X_pen, U, b0_pen, ba_pen, model, theta_min = theta_min, theta_max = theta_max)
        #   } else {
        #     # λ < ∞ path (kept for completeness)
        #     theta_hat_lambda(Y_mix, X_pen, b0_pen, ba_pen, lambda, model)
        #   }
        # }, error = function(e) NA_real_)
        # 
        out_rows[[ridx]] <- data.frame(
          split_id = split_id,
          prop_B   = prop_B,
          rep_id   = rep,
          theta_hat = if (is.finite(theta_hat)) theta_hat else NA_real_,
          stringsAsFactors = FALSE
        )
        ridx <- ridx + 1L
      }
    }
    
    do.call(rbind, out_rows)
  }
  
  # ---- parallel work --------------------------------------------------------
  if (verbose) cat("\n=== PARALLEL SPLIT PROCESSING ===\n")
  cl <- parallel::makeCluster(n_cores)
  
  
  
  # After: cl <- parallel::makeCluster(n_cores)
  parallel::clusterEvalQ(cl, {
    library(glmnet)
    library(porridge)  # <- needed by ridge_efficient() / ridgeGLM
    # library(pracma)  # optional; for theta rotation if ever enabled
  })
  
  # Export all functions used by the worker, directly or indirectly
  parallel::clusterExport(cl, varlist = c(
    # core theta functions
    "theta_inf_hat", "thetaInf", "theta_hat_lambda",
    # ridge wrappers + helpers
    "ridge_complete", "ridge_efficient", "ridgeGLM_noblk_memsafe",
    "ridgeGLM", "AridgePLM", "fit_glmnet_scaler", "%||%",
    # any small utilities the worker calls
    "make_penalty_factor"
  ), envir = environment())
  
  
  
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  
  parallel::clusterExport(cl, varlist = c(
    # data and settings
    "Y_A","X_A","Y_B","X_B","n_A","n_B","p",
    "lambda","lambda_beta","model","fam",
    "n_bootstrap","n_reps","replacement_props",
    "idx_unp","idx_pen","intercept","unpen_as_signpost",
    "target_size_A","target_size_B","test_size","enhancement",
    "scaled","theta_min","theta_max","make_penalty_factor"
  ), envir = environment())
  
  parallel::clusterEvalQ(cl, { library(glmnet) })
  
  res_list <- parallel::parLapply(cl, X = seq_len(n_bootstrap), fun = worker_function)
  results  <- do.call(rbind, res_list)
  
  # ---- write once (no early header) -----------------------------------------
  if (!is.null(output_file)) {
    utils::write.csv(results, file = output_file, row.names = FALSE)
    if (verbose) cat("Results written to:", output_file, "\n")
  }
  
  results
}








#' Split-Based Empirical Validity with Zero Overlap
#'
#' Implements empirical validity by splitting datasets in half:
#' - First half used for target estimation (zero overlap)
#' - Second half used for testing 
#' - Optional bootstrap enhancement for larger effective samples
#' - Full dataset warm start for faster convergence
#'
#' @param Y_A Numeric vector, responses from dataset A (reference group)
#' @param X_A Numeric matrix, predictors from dataset A (reference group)
#' @param Y_B Numeric vector, responses from dataset B (alternative group)
#' @param X_B Numeric matrix, predictors from dataset B (alternative group)
#' @param lambda Numeric, regularization parameter for test
#' @param lambda_beta Numeric, regularization parameter for target estimation
#' @param model Character, one of "logistic", "poisson", "linear"
#' @param n_bootstrap Integer, number of random splits to evaluate
#' @param n_reps Integer, number of repetitions per replacement proportion
#' @param replacement_props Numeric vector, proportions of B samples
#' @param scaled Logical, whether to scale estimated targets
#' @param verbose Logical, whether to show progress
#' @param theta_min Numeric, lower bound for theta search
#' @param theta_max Numeric, upper bound for theta search
#' @param n_cores Integer, number of cores for parallel processing
#' @param output_file Character, path to output CSV file
#' @param enhancement Logical or numeric. If FALSE, use pure 50/50 split. 
#'                   If numeric (e.g., 0.75), bootstrap that fraction from each half.
#'
#' @return Data frame with empirical validity results
#' @export
# signpost_empirical_validity_split <- function(Y_A, X_A, Y_B, X_B,
#                                               lambda, lambda_beta, 
#                                               model = "logistic",
#                                               n_bootstrap = 50, 
#                                               n_reps = 10,
#                                               replacement_props = seq(0, 1, by = 0.1),
#                                               scaled = TRUE,
#                                               verbose = TRUE,
#                                               theta_min = 0, 
#                                               theta_max = 1,
#                                               n_cores = 8,
#                                               output_file = "empirical_validity_split.csv",
#                                               enhancement = FALSE) {
#   
#   # Input validation
#   if (length(Y_A) != nrow(X_A) || length(Y_B) != nrow(X_B)) {
#     stop("Y and X dimensions must match within each dataset")
#   }
#   if (ncol(X_A) != ncol(X_B)) {
#     stop("X_A and X_B must have same number of columns")
#   }
#   if (is.numeric(enhancement) && enhancement <= 0) {
#     stop("enhancement must be FALSE or a positive number")
#   }
#   
#   n_A <- nrow(X_A)
#   n_B <- nrow(X_B)
#   p <- ncol(X_A)
#   
#   # Calculate split sizes (use smaller dataset to ensure equal test set sizes)
#   min_size <- min(n_A, n_B)
#   test_size <- floor(min_size / 2)
#   target_size_A <- n_A - test_size
#   target_size_B <- n_B - test_size
#   
#   if (verbose) {
#     cat("=== SPLIT-BASED EMPIRICAL VALIDITY (ZERO OVERLAP) ===\n")
#     cat("Dataset A:", n_A, "samples, Dataset B:", n_B, "samples\n")
#     cat("Split strategy: Target A =", target_size_A, ", Test A =", test_size, "\n")
#     cat("                Target B =", target_size_B, ", Test B =", test_size, "\n")
#     
#     if (is.numeric(enhancement)) {
#       target_enhanced_size_A <- round(enhancement * target_size_A)
#       target_enhanced_size_B <- round(enhancement * target_size_B)
#       test_enhanced_size <- round(enhancement * test_size)
#       cat("Bootstrap enhancement factor:", enhancement, "\n")
#       cat("Enhanced sizes: Target A =", target_enhanced_size_A, 
#           ", Target B =", target_enhanced_size_B, 
#           ", Test =", test_enhanced_size, "\n")
#     } else {
#       cat("No bootstrap enhancement (pure 50/50 split)\n")
#     }
#     
#     cat("Workers:", n_cores, ", Splits:", n_bootstrap, "\n")
#     cat("Zero overlap guaranteed between target and test sets\n")
#   }
#   
#   # ========================================================================
#   # PHASE 1: COMPUTE WARM START FROM FULL DATASETS
#   # ========================================================================
#   
#   if (verbose) cat("\n=== COMPUTING WARM START FROM FULL DATASETS ===\n")
#   
#   # Robust target estimation function
#   estimate_target_robust <- function(Y_data, X_data, lambda_beta, model, scaled = TRUE, 
#                                      warm_start = NULL, max_attempts = 3) {
#     
#     for (attempt in 1:max_attempts) {
#       target_est <- tryCatch({
#         if (model == "logistic") {
#           fit <- glmnet::glmnet(X_data, Y_data, 
#                                 family = "binomial", alpha = 0,
#                                 lambda = lambda_beta, intercept = FALSE,
#                                 standardize = FALSE)
#           beta_est <- as.vector(coef(fit))[-1]
#         } else if (model == "poisson") {
#           fit <- glmnet::glmnet(X_data, Y_data,
#                                 family = "poisson", alpha = 0,
#                                 lambda = lambda_beta, intercept = FALSE,
#                                 standardize = FALSE)
#           beta_est <- as.vector(coef(fit))[-1]
#         } else if (model == "linear") {
#           fit <- glmnet::glmnet(X_data, Y_data,
#                                 family = "gaussian", alpha = 0,
#                                 lambda = lambda_beta, intercept = FALSE,
#                                 standardize = FALSE)
#           beta_est <- as.vector(coef(fit))[-1]
#         }
#         
#         if (any(!is.finite(beta_est)) || all(abs(beta_est) < 1e-15)) {
#           if (attempt < max_attempts) {
#             lambda_beta <- lambda_beta * (attempt + 1)
#             next
#           } else {
#             stop("All coefficients are zero or non-finite")
#           }
#         }
#         
#         if (scaled && sum(beta_est^2) > 1e-12) {
#           beta_est / sqrt(sum(beta_est^2))
#         } else {
#           beta_est
#         }
#         
#       }, error = function(e) {
#         if (attempt == max_attempts) {
#           stop("Target estimation failed after ", max_attempts, " attempts: ", e$message)
#         }
#         NULL
#       })
#       
#       if (!is.null(target_est)) {
#         return(target_est)
#       }
#     }
#   }
#   
#   # Compute warm start targets from full datasets
#   if (verbose) cat("Computing beta_0 warm start from full dataset A...\n")
#   beta_0_warm <- estimate_target_robust(Y_A, X_A, lambda_beta, model, scaled)
#   
#   if (verbose) cat("Computing beta_a warm start from full dataset B...\n")
#   beta_a_warm <- estimate_target_robust(Y_B, X_B, lambda_beta, model, scaled)
#   
#   if (verbose) cat("Warm start computation completed\n")
#   
#   # ========================================================================
#   # PHASE 2: PARALLEL WORKER PROCESSING
#   # ========================================================================
#   
#   if (verbose) cat("\n=== PARALLEL SPLIT PROCESSING ===\n")
#   
#   # Initialize output file
#   header_df <- data.frame(
#     split_id = integer(0), 
#     prop_B = numeric(0),
#     rep_id = integer(0),
#     theta_hat = numeric(0),
#     stringsAsFactors = FALSE
#   )
#   write.csv(header_df, output_file, row.names = FALSE)
#   
#   # Set up cluster
#   cl <- makeCluster(n_cores)
#   
#   # Export necessary objects
#   clusterExport(cl, c(
#     # Data
#     "Y_A", "X_A", "Y_B", "X_B",
#     "lambda", "lambda_beta", "model", 
#     "n_A", "n_B", "p", "test_size", "target_size_A", "target_size_B",
#     "replacement_props", "n_reps", "scaled",
#     "theta_min", "theta_max", "enhancement",
#     # Warm start
#     "beta_0_warm", "beta_a_warm",
#     # Functions
#     "theta_inf_hat", "theta_hat_lambda", "generate_Y", "calculate_eta",
#     "thetaInf", "estimate_target_robust"
#   ), envir = environment())
#   
#   clusterEvalQ(cl, {
#     library(glmnet)
#   })
#   
#   # Worker function - each worker processes one split independently
#   worker_function <- function(split_id) {
#     tryCatch({
#       # Each worker gets unique seed
#       #set.seed(split_id * 100000)
#       
#       # STEP 1: Create 50/50 split (disjoint halves)
#       # Split A: first target_size_A for target, last test_size for test
#       A_indices <- sample(n_A)  # Random permutation
#       A_target_indices <- A_indices[1:target_size_A]
#       A_test_indices <- A_indices[(target_size_A + 1):(target_size_A + test_size)]
#       
#       # Split B: first target_size_B for target, last test_size for test  
#       B_indices <- sample(n_B)  # Random permutation
#       B_target_indices <- B_indices[1:target_size_B]
#       B_test_indices <- B_indices[(target_size_B + 1):(target_size_B + test_size)]
#       
#       # STEP 2: Apply bootstrap enhancement if requested
#       if (is.numeric(enhancement)) {
#         # Calculate enhanced sizes
#         target_enhanced_size_A <- round(enhancement * target_size_A)
#         target_enhanced_size_B <- round(enhancement * target_size_B)
#         test_enhanced_size <- round(enhancement * test_size)
#         
#         # Bootstrap from each split (oversample from disjoint halves)
#         A_target_enhanced <- sample(A_target_indices, target_enhanced_size_A, replace = TRUE)
#         B_target_enhanced <- sample(B_target_indices, target_enhanced_size_B, replace = TRUE)
#         A_test_enhanced <- sample(A_test_indices, test_enhanced_size, replace = TRUE)
#         B_test_enhanced <- sample(B_test_indices, test_enhanced_size, replace = TRUE)
#         
#         # Use enhanced samples
#         Y_A_target <- Y_A[A_target_enhanced]
#         X_A_target <- X_A[A_target_enhanced, , drop = FALSE]
#         Y_B_target <- Y_B[B_target_enhanced]
#         X_B_target <- X_B[B_target_enhanced, , drop = FALSE]
#         
#         Y_A_test <- Y_A[A_test_enhanced]
#         X_A_test <- X_A[A_test_enhanced, , drop = FALSE]
#         Y_B_test <- Y_B[B_test_enhanced]
#         X_B_test <- X_B[B_test_enhanced, , drop = FALSE]
#         
#       } else {
#         # Use pure splits (no enhancement)
#         Y_A_target <- Y_A[A_target_indices]
#         X_A_target <- X_A[A_target_indices, , drop = FALSE]
#         Y_B_target <- Y_B[B_target_indices]
#         X_B_target <- X_B[B_target_indices, , drop = FALSE]
#         
#         Y_A_test <- Y_A[A_test_indices]
#         X_A_test <- X_A[A_test_indices, , drop = FALSE]
#         Y_B_test <- Y_B[B_test_indices]
#         X_B_test <- X_B[B_test_indices, , drop = FALSE]
#       }
#       
#       # Estimate targets for this split using warm start
#       beta_0_hat <- estimate_target_robust(Y_A_target, X_A_target, lambda_beta, 
#                                            model, scaled, warm_start = beta_0_warm)
#       beta_a_hat <- estimate_target_robust(Y_B_target, X_B_target, lambda_beta, 
#                                            model, scaled, warm_start = beta_a_warm)
#       
#       # Test evaluation on this split's test sets
#       results_this_split <- list()
#       result_idx <- 1
#       
#       # Test each replacement proportion
#       for (prop_B in replacement_props) {
#         base_size <- length(Y_A_test)  # Both test sets have same size
#         n_from_A <- round((1 - prop_B) * base_size)
#         n_from_B <- round(prop_B * base_size)
#         
#         # Multiple repetitions for each prop_B
#         for (rep in 1:n_reps) {
#           theta_hat_value <- tryCatch({
#             set.seed(split_id * 100000 + rep * 1000 + which(replacement_props == prop_B))
#             
#             # Create mixed dataset
#             if (n_from_A == 0) {
#               # Pure B
#               sample_B_idx <- sample(length(Y_B_test), base_size, replace = FALSE)
#               Y_mixed <- Y_B_test[sample_B_idx]
#               X_mixed <- X_B_test[sample_B_idx, , drop = FALSE]
#             } else if (n_from_B == 0) {
#               # Pure A
#               Y_mixed <- Y_A_test
#               X_mixed <- X_A_test
#             } else {
#               # Mixed
#               sample_A_idx <- sample(length(Y_A_test), n_from_A, replace = FALSE)
#               sample_B_idx <- sample(length(Y_B_test), n_from_B, replace = FALSE)
#               
#               Y_mixed <- c(Y_A_test[sample_A_idx], Y_B_test[sample_B_idx])
#               X_mixed <- rbind(X_A_test[sample_A_idx, , drop = FALSE], 
#                                X_B_test[sample_B_idx, , drop = FALSE])
#             }
#             
#             # Validate mixed dataset
#             if (any(!is.finite(Y_mixed)) || any(!is.finite(X_mixed))) {
#               return(NA_real_)
#             }
#             
#             # Compute theta_hat
#             if (is.infinite(lambda)) {
#               result <- theta_inf_hat(Y_mixed, X_mixed, 
#                                       matrix(ncol = 0, nrow = length(Y_mixed)), 
#                                       beta_0_hat, beta_a_hat, model, 
#                                       theta_min = theta_min, theta_max = theta_max)
#             } else {
#               result <- theta_hat_lambda(Y_mixed, X_mixed, 
#                                          beta_0_hat, beta_a_hat, lambda, model)
#             }
#             
#             if (is.finite(result)) result else NA_real_
#             
#           }, error = function(e) {
#             NA_real_
#           })
#           
#           # Store result
#           results_this_split[[result_idx]] <- data.frame(
#             split_id = split_id,
#             prop_B = prop_B,
#             rep_id = rep,
#             theta_hat = theta_hat_value,
#             stringsAsFactors = FALSE
#           )
#           result_idx <- result_idx + 1
#         }
#       }
#       
#       # Return results for this split
#       do.call(rbind, results_this_split)
#       
#     }, error = function(e) {
#       # Return empty data frame on complete failure
#       data.frame(
#         split_id = split_id,
#         prop_B = replacement_props[1],
#         rep_id = 1,
#         theta_hat = NA_real_,
#         stringsAsFactors = FALSE
#       )[0, ]
#     })
#   }
#   
#   # Process all splits in parallel
#   if (verbose) {
#     pb <- txtProgressBar(min = 0, max = n_bootstrap, style = 3)
#     cat("Processing", n_bootstrap, "splits in parallel...\n")
#   }
#   
#   all_results_list <- parLapply(cl, 1:n_bootstrap, worker_function)
#   
#   stopCluster(cl)
#   
#   if (verbose) {
#     close(pb)
#     cat("\nCombining and saving results...\n")
#   }
#   
#   # Combine all results
#   final_results <- do.call(rbind, all_results_list)
#   
#   # Write to file
#   write.csv(final_results, output_file, row.names = FALSE)
#   
#   # Analysis and reporting
#   if (verbose) {
#     na_count <- sum(is.na(final_results$theta_hat))
#     total_count <- nrow(final_results)
#     cat("Final results - Total:", total_count, ", NA:", na_count, 
#         "(", round(100 * na_count / total_count, 1), "%)\n")
#     
#     # Analyze monotonicity
#     if (nrow(final_results) > 0) {
#       avg_by_prop <- aggregate(theta_hat ~ prop_B, final_results, mean, na.rm = TRUE)
#       avg_by_prop <- avg_by_prop[order(avg_by_prop$prop_B), ]
#       cat("\nAverage theta by prop_B:\n")
#       print(avg_by_prop)
#       
#       # Check for monotonicity
#       if (nrow(avg_by_prop) > 1) {
#         diffs <- diff(avg_by_prop$theta_hat)
#         is_monotonic <- all(diffs >= -1e-10)
#         cat("Monotonic increase:", is_monotonic, "\n")
#         
#         if (!is_monotonic) {
#           decreases <- which(diffs < -1e-10)
#           cat("Non-monotonic transitions at prop_B:", 
#               avg_by_prop$prop_B[decreases], "to", avg_by_prop$prop_B[decreases + 1], "\n")
#         }
#         
#         # Check for jumps
#         abs_diffs <- abs(diffs)
#         max_jump_idx <- which.max(abs_diffs)
#         max_jump <- diffs[max_jump_idx]
#         cat("Largest change:", round(max_jump, 4), 
#             "from prop_B =", avg_by_prop$prop_B[max_jump_idx], 
#             "to", avg_by_prop$prop_B[max_jump_idx + 1], "\n")
#       }
#     }
#     
#     cat("Zero overlap guaranteed - results reflect true generalization\n")
#     cat("Results saved to:", output_file, "\n")
#   }
#   
#   return(final_results)
# }
# 
# 
