augmented_data_poisson_lognormal_noloop <- function(Y_obs,
                                                    eta_draw_current,
                                                    eta_tilde_mean,
                                                    Ht_variance,
                                                    pswitch_Y_current,
                                                    ctuning_val) {
  dim_Y <- dim(Y_obs)
  N_spatial <- dim_Y[1] # Number of spatial units
  T_time <- dim_Y[2]    # Number of time points
  
  eta_draw_new <- eta_draw_current # Initialize new eta with current values
  pswitch_Y_new <- pswitch_Y_current
  
  
  current_variance <- Ht_variance # Will be scalar if Ht_variance is scalar
  
  for (i in 1:N_spatial) {
    # Extract row-specific data
    eta_draw_i_current <- eta_draw_current[i, ]
    eta_tilde_mean_i <- eta_tilde_mean[i, ]
    Y_obs_i <- Y_obs[i, ]
    
  
    variance_i <- current_variance # This will be a scalar repeated for all t if Ht_variance is scalar
    
    sigma_proposal_i <- sqrt(variance_i) # Standard deviation for proposal
    
    
    if (length(ctuning_val) == T_time) {
      
      proposal_scaled_random_part <- ctuning_val * rnorm(T_time) * sigma_proposal_i
    } else { # Scalar ctuning_val
      proposal_scaled_random_part <- ctuning_val[1] * rnorm(T_time) * sigma_proposal_i
    }
    eta_candidate_i <- eta_draw_i_current + proposal_scaled_random_part
    
    # Calculate log-posterior for current and candidate eta values for row i
    log_post_eta_PLN <- function(eta_t_values, theta_t_values, Y_t_values, variance_eta_t) {
      log_lik_t <- Y_t_values * eta_t_values - exp(eta_t_values) 
      log_prior_t <- -0.5 * ((eta_t_values - theta_t_values)^2) / variance_eta_t
      return(log_lik_t + log_prior_t) # Returns a vector of length T
    }
    
    log_post_vec_current_i <- log_post_eta_PLN(eta_draw_i_current, eta_tilde_mean_i, Y_obs_i, variance_i)
    log_post_vec_candidate_i <- log_post_eta_PLN(eta_candidate_i, eta_tilde_mean_i, Y_obs_i, variance_i)
    
    log_acc_prob_i_vector <- log_post_vec_candidate_i - log_post_vec_current_i
    
    accepted_indices <- log(runif(T_time)) < min(0, log_acc_prob_i_vector)
    
    
    # Update eta_draw and pswitch_Y where accepted
    if (any(accepted_indices, na.rm = TRUE)) {
      eta_draw_new[i, accepted_indices] <- eta_candidate_i[accepted_indices]
      pswitch_Y_new[i, accepted_indices] <- pswitch_Y_new[i, accepted_indices] + 1
    }
  } # End of loop over spatial units i
  
  return(list(eta_draw = eta_draw_new, pswitch_Y = pswitch_Y_new))
}


posterior_conditional_variance <- function(data_slice, invCorrZ, a1, b1, nstar, tstar) {
  
  # Get dimensions (N=p, T=t_slice)
  dims <- dim(data_slice)
  N <- dims[1]
  T_slice <- dims[2]
  
  # Posterior shape parameter for sigma^2 ~ IG(shape, rate)
  # Equivalent to shape parameter for 1/sigma^2 ~ Gamma(shape, rate)
  g1_pos <- (nstar * tstar / 2) + a1
  
  quad_form_sum <- 0
  for (t in 1:T_slice) {
    z_t <- matrix(data_slice[, t], ncol = 1) 
    invCorrZ_zt <- invCorrZ %*% z_t
    quad_form_sum <- quad_form_sum + crossprod(z_t, invCorrZ_zt)[1, 1]
  }
  
  # Posterior rate parameter for sigma^2 ~ IG(shape, rate)
  g2_pos <- b1 + 0.5 * quad_form_sum
  
  # Sample precision (h_eta_star = 1/rho1_star) 
  h_eta_star <- rgamma(1, shape = g1_pos, rate = g2_pos)
  
  rho1_star <- 1 / h_eta_star
  
  
  return(rho1_star)
}

MH_spatial_correlation_CAR <- function(W, data_slice, rho1, min_tau, max_tau, n_grid = 15) {
  
  # Ensure W is sparse
  W <- as(W, "sparseMatrix")
  
  # Get dimensions
  # N is spatial dimension (p), T is time dimension (t or 1)
  dims <- dim(data_slice)
  N <- dims[1]
  T_slice <- dims[2] # Number of 'time' points in the slice
  
  # Create the grid for the spatial parameter rho2 
  allowed_range <- seq(min_tau, max_tau, length.out = n_grid)
  nr <- length(allowed_range)
  lpos <- numeric(nr) # Vector to store log posterior densities
  
  D <- Matrix::Diagonal(N, Matrix::rowSums(W))
  
  # Calculate log posterior for each candidate rho2 value
  for (i in 1:nr) {
    range_val <- allowed_range[i]
    
    # Calculate precision matrix for this rho2 value
    # invS = (D - range_val * W) * (1 / rho1)
    invS <- Matrix::forceSymmetric((D - range_val * W) * (1 / rho1))
    
    quad_form_sum <- 0
    for (t in 1:T_slice) {
      z_t <- matrix(data_slice[, t], ncol = 1) 
      invS_z_t <- invS %*% z_t
      quad_form_sum <- quad_form_sum + Matrix::crossprod(z_t, invS_z_t)[1,1]
    }
    ee = eigen(invS)$values
    # Calculate log determinant term
    log_det_val = sum(log(1/ee[ee>0]))
    
    # Log posterior density (up to constant)
    lpos[i] <- -(0.5 * T_slice * log_det_val) - (0.5 * quad_form_sum)
  }
  
  # Normalize probabilities using log-sum-exp trick
  max_lpos <- max(lpos, na.rm = TRUE) # Find max log density
  # Check if all lpos are -Inf
  if(is.infinite(max_lpos) && max_lpos < 0) {
    stop("All grid points resulted in -Inf log-posterior.")
  } else {
    # Calculate probabilities: exp(lpos - max_lpos) / sum(exp(lpos - max_lpos))
    log_sum_exp <- max_lpos + log(sum(exp(lpos - max_lpos)))
    probx <- exp(lpos - log_sum_exp)
    # Handle potential NaN if all lpos were identical before exp
    probx[is.na(probx)] <- 0 
    # Renormalize if any NAs were removed or if probabilities don't sum to 1
    sum_probx <- sum(probx)
    if (abs(sum_probx - 1) > 1e-9 || sum_probx == 0) {
      stop("Probabilities do not sum to 1 after log-sum-exp, possibly due to numerical issues or all -Inf.")
      
    }
  }
  
  
  # Sample one value from the allowed range based on calculated probabilities
  range_draw <- sample(allowed_range, size = 1, prob = probx)
  
  return(range_draw)
}




fit_SDGLMC <- function(Y, X_cov, Xmean_val, Zmeasconf_list, 
                       W_adj, offset_term,
                       interaction_type,
                       nrep, nburn, thin,
                       ctuning_poisson, priors,
                       ave_restart = NULL, out_restart = NULL
) {
  
  # --- Initial Setup ---
  dim_Y <- dim(Y)
  p <- dim_Y[1]
  t <- dim_Y[2]
  
  ST_interaction_X <- 0 # For the X covariate's space-time term
  
  a0_prior = priors$a0
  a1_prior = priors$a1
  b0_prior = priors$b0
  b1_prior = priors$b1
  s2_a = priors$s2_a
  s2_b = priors$s2_b
  
  if (is.null(out_restart)) {
    restart <- FALSE
    out <- list()
  } else {
    restart <- TRUE
    out <- out_restart
    message("Restarting MCMC from provided 'out_restart' state.")
  }
  if (is.null(ave_restart)){
    ave <- list()
  } else {
    ave <- ave_restart
    message("Using provided 'ave_restart' structure.")
  }
  
  
  # Initial variance/correlation parameters (if not restarting)
  if (!restart) {
    rho1_space_draw <- 1e-6
    rho1_space_time_draw <- 1e-4       # For Bspacetimedraw (X related)
    rho1_0_space_time_draw <- 1e-7     # For B0spacetimedraw (baseline ST)
    if (interaction_type == 5) {
      rho1_0_space_time_draw <- matrix(1e-7, nrow = p, ncol = 1) # Vector for interaction 5
    }
    rho1_0_space_draw <- 1e-6        # For B0spacedraw (baseline S)
    
    rho2_space_draw <- 0.5
    rho2_space_time_draw <- 0.5      # For Bspacetimedraw
    rho2_0_space_draw <- 0.5         # For B0spacedraw
    rho2_0_space_time_draw <- 0.5    # For B0spacetimedraw
    
    Q1invdraw_time <- 100000       # For Btimedraw (X related)
    Q0invdraw_time <- 100000       # For B0timedraw (baseline T)
  }
  

  # Pre-calculate some stats
  sum_zeros_obs_ <- sum(Y == 0, na.rm = TRUE)
  p95_obs_ <- apply(Y, 1, quantile, probs = 0.95, na.rm = TRUE)
  
  W_adj <- as(W_adj, "sparseMatrix")
  D_adj <- Matrix::Diagonal(p, Matrix::rowSums(W_adj))
  
  # --- Process Measured Confounders (Zmeasconf_list) ---
  # Zmeasconf is a list of matrices (nzmconf x t for each p)
  # Add intercept to Zmeasconf if not present for each site i
  for (i in 1:p) {
    if (!is.null(Zmeasconf_list[[i]]) && ncol(Zmeasconf_list[[i]]) > 0) {
      sumZmeas_i <- rowSums(Zmeasconf_list[[i]], na.rm = TRUE)
      # Check if a row sums to t (indicating a constant intercept)
      if (!any(abs(sumZmeas_i - t) < 1e-9)) { # If no intercept-like row
        Zmeasconf_list[[i]] <- rbind(Zmeasconf_list[[i]], matrix(1, nrow = 1, ncol = t))
      }
    } else if (is.null(Zmeasconf_list[[i]]) || ncol(Zmeasconf_list[[i]]) == 0) {
      # If empty for this site, create an intercept-only term
      Zmeasconf_list[[i]] <- matrix(1, nrow = 1, ncol = t)
    }
  }
  
  if (is.null(Zmeasconf_list[[1]]) || nrow(Zmeasconf_list[[1]]) == 0) {
    measured_confounders <- FALSE
    nzmconf <- 0
    MConf_mat <- NULL # Matrix (pt x nzmconf)
  } else {
    measured_confounders <- TRUE
    nzmconf <- nrow(Zmeasconf_list[[1]])
    MConf_mat <- matrix(0, nrow = p * t, ncol = nzmconf)
    for (i in 1:p) {
      MConf_mat[seq(i, p * t, by = p), ] <- t(Zmeasconf_list[[i]])
    }
  }
  
  m0 <- p # Number of elements for spatial baseline (B0spacedraw, B0spacetimedraw rows)
  m1 <- 1 # Number of elements for X-related temporal (Btimedraw) - scalar X effect
  
  offset_mat <- offset_term
  # --- Initialize Latent Data (eta_tilde) based on distribution ---
  if (!restart) {
    eta_draw <- Y
    
    eta_draw <- log(Y + 0.5) # Initial eta for Poisson
    eta_tilde <- eta_draw + matrix(rnorm(p * t, 0, 0.00005), nrow = p, ncol = t)
    Ytemp <- Y
    Otemp <- offset_mat
  }
  
  
  # --- Preliminaries for X covariate ---
  if (is.null(Xmean_val) || length(Xmean_val)==0) { # Check if Xmean_val is truly empty or NULL
    X_centered <- X_cov - matrix(rowMeans(X_cov), p, t) # Center by site-specific mean
  } else {
    X_centered <- X_cov - Xmean_val
  }
  
  pswitch_Y <- matrix(0, nrow = p, ncol = t)
  # --- Initial GLM fit for starting values (if not restarting) ---
  if (!restart) {
    meanZ_mat <- matrix(rnorm(p*t), nrow = p, ncol = t)
    
    # Basis for initial smooth fit of baseline
    nbasis <- max(4 * round(t / 365), 3)
    basis_time <- ns(1:t, nbasis, intercept = T)
    
    # Construct Design Matrix for GLM
    DesignX_list <- list()
    if (measured_confounders && nzmconf > 1) { # If MConf has more than just an intercept
      DesignX_list$MConf_no_intercept <- MConf_mat[, 1:(nzmconf - 1), drop = FALSE]
    } else if (measured_confounders && nzmconf == 1) { # MConf is just an intercept
      # Don't add MConf if it's only the intercept, which is covered by basis_time for each p
    }
    DesignX_list$X_vec <- as.vector(X_centered)
    DesignX_list$SplineBasis_kron <- kronecker(diag(p), basis_time)
    
    DesignX_mat <- do.call(cbind, DesignX_list)
    
    
    target_poisson_init <- as.vector(eta_tilde) - as.vector(Otemp)
    fit <- lm(target_poisson_init ~ DesignX_mat - 1)
    bml <- coefficients(fit)
    
    idx_start_spline <- ncol(DesignX_mat) - ncol(DesignX_list$SplineBasis_kron) + 1
    BB_vec <- DesignX_list$SplineBasis_kron %*% bml[idx_start_spline:ncol(DesignX_mat)]
    BB_mat <- matrix(BB_vec, nrow = p, ncol = t)
    
    idx_start_MConf <- 1
    if (!is.null(DesignX_list$MConf_no_intercept)) {
      len_MConf_no_int <- ncol(DesignX_list$MConf_no_intercept)
      meanZ_mat <- matrix(DesignX_list$MConf_no_intercept %*% bml[idx_start_MConf:(idx_start_MConf + len_MConf_no_int - 1)],
                          nrow = p, ncol = t) + mean(BB_mat)
    } else {
      meanZ_mat <- matrix(mean(BB_mat), nrow = p, ncol = t)
    }
    
    idx_X_coeff <- if (!is.null(DesignX_list$MConf_no_intercept)) ncol(DesignX_list$MConf_no_intercept) + 1 else 1
    Bdraw <- bml[idx_X_coeff]
    
    mean_BB_overall <- mean(BB_mat)
    B0timedraw_init <- colMeans(BB_mat) - mean_BB_overall   # 1 x t
    B0spacedraw_init <- rowMeans(BB_mat) - mean_BB_overall  # p x 1
    B0spacetimedraw_init <- BB_mat -
      matrix(colMeans(BB_mat), nrow = p, ncol = t, byrow = TRUE) -
      matrix(rowMeans(BB_mat), nrow = p, ncol = t, byrow = FALSE) +
      mean_BB_overall                 # p x t
    
    s2_err_mis_draw_init <- 0.01
    thetay_draw_init <- offset_mat + meanZ_mat + BB_mat + (X_centered * Bdraw) - mean_BB_overall
    aug_data_result <- augmented_data_poisson_lognormal_noloop(Y, eta_tilde, thetay_draw_init, s2_err_mis_draw_init, pswitch_Y, ctuning_poisson)
    eta_tilde <- aug_data_result$eta_draw
      
    
    # Assign initial values to sampler variables
    B0timedraw <- matrix(B0timedraw_init, nrow=1) # Ensure 1xt
    B0spacedraw <- matrix(B0spacedraw_init, ncol=1) # Ensure px1
    B0spacetimedraw <- B0spacetimedraw_init
    s2_err_mis_draw <- s2_err_mis_draw_init
    # Bdraw already set
    
    Btimedraw <- matrix(0, nrow = m1, ncol = t)
    Bspacedraw <- matrix(0, nrow = p, ncol = 1) # For X related effect
    Bspacetimedraw <- matrix(0, nrow = p, ncol = t) # For X related effect
    gamma_draw <- if(measured_confounders) matrix(0, nrow=nzmconf, ncol=1) else NULL
    if(measured_confounders && !is.null(DesignX_list$MConf_no_intercept)){
      gamma_draw[1:len_MConf_no_int,1] <- bml[idx_start_MConf:(idx_start_MConf + len_MConf_no_int - 1)]
    }
    
    
  }
  
  # --- Setup for MCMC (Precision matrices, H matrices, BigG matrices) ---
  VB_OLS_prior <- diag(m1) * 1
  B_0_prvar <- diag(m1) * 4
  
  q_val <- m0 + m1
  Tq <- t * q_val
  
  # H matrix for random walk (common structure)
  H_mat <- Diagonal(Tq) - sparseMatrix(i = (q_val + 1):Tq, j = 1:((t - 1) * q_val), x = 1, dims = c(Tq, Tq))
  
  # H2 matrix for interaction_type == 3
  H2_mat <- NULL
  if (interaction_type == 3) {
    H2_mat <- Diagonal(Tq) - sparseMatrix(i = seq(q_val + 1 + m0, Tq, by=m0+1), 
                                          j = seq(m0+1, (t - 1) * q_val, by=m0+1),
                                          x = 1, dims = c(Tq, Tq))
  }
  
  
  # Construct BigG matrices
  # bigG: for (B0spacetimedraw, Btimedraw)
  # bigG2: for (B0timedraw, Bspacetimedraw)
  # bigG3: for (B0spacedraw, Bspacedraw)
  
  Xtemp_list_G <- vector("list", t) # For bigG
  Xtemp_list_G2 <- vector("list", t) # For bigG2
  
  for (i in 1:t) {
    Xtemp_list_G[[i]] <- cbind(Diagonal(p), X_centered[, i])
    
    Xtemp_list_G2[[i]] <- cbind(1, Diagonal(p, X_centered[, i]))
  }
  bigG_mat <- Matrix::bdiag(Xtemp_list_G)  # pt x (p+1)t = pt x Tq
  bigG2_mat <- Matrix::bdiag(Xtemp_list_G2) # pt x (1+p)t = pt x Tq
  
  # bigG3: for (B0spacedraw, Bspacedraw) - state is [B0spacedraw (p x 1); Bspacedraw (p x 1)] = 2p x 1
  
  row_idx_G3 <- rep(0, p * t * 2)
  col_idx_G3 <- rep(0, p * t * 2)
  val_G3 <- numeric(p * t * 2)
  current_g3_ptr <- 1
  for (i in 1:t) {
    time_block_rows <- ((i - 1) * p + 1):(i * p)
    row_idx_G3[current_g3_ptr:(current_g3_ptr + p - 1)] <- time_block_rows
    col_idx_G3[current_g3_ptr:(current_g3_ptr + p - 1)] <- 1:p
    val_G3[current_g3_ptr:(current_g3_ptr + p - 1)] <- 1
    current_g3_ptr <- current_g3_ptr + p
    
    row_idx_G3[current_g3_ptr:(current_g3_ptr + p - 1)] <- time_block_rows
    col_idx_G3[current_g3_ptr:(current_g3_ptr + p - 1)] <- (p + 1):(2 * p)
    val_G3[current_g3_ptr:(current_g3_ptr + p - 1)] <- X_centered[, i]
    current_g3_ptr <- current_g3_ptr + p
  }
  bigG3_mat <- sparseMatrix(i = row_idx_G3, j = col_idx_G3, x = val_G3, dims = c(p * t, 2 * p))
  
  
  # --- Load from restart if applicable ---
  if (restart) {
    message("Loading parameters from restart object 'out'...")
    last_idx <- NCOL(out$S2_err_mis_)
    
    Bspacedraw <- matrix(out$Bspacedraw[, last_idx], ncol=1)
    Btimedraw <- matrix(out$Btimedraw, nrow=1)
    Bspacetimedraw <- out$Bspacetimedraw # Should be p x t from last state save
    B0timedraw <- matrix(out$B0timedraw, nrow=1)     # Should be 1 x t
    B0spacedraw <- matrix(out$B0spacedraw, ncol=1)  # Should be p x 1
    if(!is.null(out$B0spacetimedraw)) B0spacetimedraw <- out$B0spacetimedraw else B0spacetimedraw <- matrix(0,p,t)
    
    
    rho1_space_draw <- out$RHO1_space_[last_idx]
    rho1_space_time_draw <- out$RHO1_space_time_[last_idx]
    
    rho1_0_space_time_draw <- out$RHO1_0_space_time_[last_idx]
    if (interaction_type == 5) {
      rho1_0_space_time_draw <- matrix(out$RHO1_0_space_time_[,last_idx], ncol=1)
    }
    
    rho1_0_space_draw <- out$RHO1_0_space_[last_idx]
    rho2_space_draw <- out$RHO2_space_[last_idx]
    rho2_space_time_draw <- out$RHO2_space_time_[last_idx]
    rho2_0_space_time_draw <- out$RHO2_0_space_time_[last_idx]
    rho2_0_space_draw <- out$RHO2_0_space_[last_idx]
    
    Q1invdraw_time <- out$Q1inv_time[last_idx]
    Q0invdraw_time <- out$Q0inv_time[last_idx]
    s2_err_mis_draw <- out$S2_err_mis_[last_idx]
    
    eta_tilde <- out$eta_tilde # p x t matrix
    if (measured_confounders && !is.null(MConf_mat)) {
      gamma_draw <- matrix(out$G_[, last_idx], ncol=1)
      meanZ_mat <- matrix(MConf_mat %*% gamma_draw, nrow = p, ncol = t)
    } else {
      gamma_draw <- NULL
      meanZ_mat <- matrix(0, nrow = p, ncol = t)
    }
    Bdraw <- out$Bdraw[last_idx]
  }
  s2_err_mis_draw = 1e-6
  
  # Initialize CAR precision matrices (Qinv...)
  Q1invdraw_space <- (1 / rho1_space_draw) * forceSymmetric(D_adj - rho2_space_draw * W_adj)
  Q1invdraw_spacetime <- (1 / rho1_space_time_draw) * forceSymmetric(D_adj - rho2_space_time_draw * W_adj)
  Q0invdraw_space <- (1 / rho1_0_space_draw) * forceSymmetric(D_adj - rho2_0_space_draw * W_adj)
  
  # Q0invdraw_spacetime depends on interaction_type
  if (interaction_type == 1) {
    rho1_0_space_time_draw <- 1e-9 # Fixed small variance, no CAR
    Q0invdraw_spacetime <- (1 / rho1_0_space_time_draw) * Diagonal(p)
  } else if (interaction_type == 2) { # Diagonal variance, no spatial corr
    Q0invdraw_spacetime <- (1 / rho1_0_space_time_draw) * Diagonal(p) # rho1_0_space_time_draw is scalar
  } else if (interaction_type == 3 || interaction_type == 4) { # Full CAR
    Q0invdraw_spacetime <- (1 / rho1_0_space_time_draw) * forceSymmetric(D_adj - rho2_0_space_time_draw * W_adj)
  } else if (interaction_type == 5) { # Heterogeneous diagonal variance
    Q0invdraw_spacetime <- Diagonal(p, x = 1 / as.vector(rho1_0_space_time_draw)) # rho1_0_space_time_draw is p x 1
  } else {
    stop("Invalid interaction_type.")
  }
  
  # Initial composite state vectors (used in subtractions in samplers)
  # Bdraw1_vec: B0spacetimedraw (pxt); Btimedraw (1xt) -> stacked (p+1)t x 1
  Bdraw1_vec <- as.vector(rbind(B0spacetimedraw, Btimedraw))
  # Bdraw2_vec: B0timedraw (1xt); Bspacetimedraw (pxt) -> stacked (1+p)t x 1
  Bdraw2_vec <- as.vector(rbind(B0timedraw, Bspacetimedraw))
  # Bdraw3_vec: B0spacedraw (px1); Bspacedraw (px1) -> stacked 2p x 1
  Bdraw3_vec <- c(as.vector(B0spacedraw), as.vector(Bspacedraw))
  
  
  # --- Storage Initialization ---
  collections <- floor(nrep / thin)
  it_print <- 1
  
  ave$Eta_tilde_mean <- matrix(0, p, t); ave$B_postmean <- 0; ave$B2_postmean <- 0
  ave$Btime_postmean <- matrix(0, m1, t); ave$Btime2_postmean <- matrix(0, m1, t)
  ave$Bspace_postmean <- matrix(0, p, 1); ave$Bspace2_postmean <- matrix(0, p, 1)
  ave$Bspacetime_postmean <- matrix(0, p, t); ave$Bspacetime2_postmean <- matrix(0, p, t)
  ave$B0spacetime_postmean <- matrix(0, m0, t); ave$B0spacetime2_postmean <- matrix(0, m0, t)
  ave$B0time_postmean <- matrix(0, 1, t); ave$B0time2_postmean <- matrix(0, 1, t)
  ave$B0space_postmean <- matrix(0, p, 1); ave$B0space2_postmean <- matrix(0, p, 1)
  ave$meanZmean <- matrix(0, p, t); ave$meanY1mean <- matrix(0, p, t); ave$meanY0mean <- matrix(0, p, t)
  
  ave$B2_c_t <- matrix(0,m1,t); ave$B2_c_s <- matrix(0,p,1)
  ave$B2_c_t_s_st <- matrix(0,p,t); ave$B2_c_t_st <- matrix(0,p,t)
  ave$B2_t_st <- matrix(0,p,t); ave$B2_t_s <- matrix(0,p,t)
  ave$B2_t_s_st <- matrix(0,p,t); ave$B2_c_t_s <- matrix(0,p,t)
  ave$B2_s_st <- matrix(0,p,t); ave$B2_c_s_st <- matrix(0,p,t)
  ave$E_t_s <- matrix(0,p,t); ave$E_t_st <- matrix(0,p,t); ave$E_s_st <- matrix(0,p,t)
  
  
  out$Bdraw <- numeric(collections)
  out$Bspacedraw <- matrix(0, p, collections)
  out$S2_err_mis_ <- numeric(collections)
  out$RHO1_space_ <- numeric(collections)
  out$RHO1_space_time_ <- numeric(collections)
  if (interaction_type == 5) {
    out$RHO1_0_space_time_ <- matrix(0, p, collections)
  } else {
    out$RHO1_0_space_time_ <- numeric(collections)
  }
  out$RHO1_0_space_ <- numeric(collections)
  out$RHO2_space_ <- numeric(collections)
  out$RHO2_space_time_ <- numeric(collections)
  out$RHO2_0_space_time_ <- numeric(collections)
  out$RHO2_0_space_ <- numeric(collections)
  out$Q1inv_time <- matrix(0, m1, collections)
  out$Q0inv_time <- numeric(collections)
  if (measured_confounders) {
    out$G_ <- matrix(0, nzmconf, collections)
  }
  out$store_llike <- numeric(collections)
  YPRED_storage <- array(0, dim = c(p, t, collections))
  out$RMSE_ <- numeric(collections)
  out$MAE_ <- numeric(collections)
  
  # Diagnostic storage
  store_LS_sum_per_iter <- matrix(0, p, t) # For log score
  store_CRPS_1_sum_per_iter <- numeric(collections)
  store_CRPS_2_sum_per_iter <- numeric(collections)
  
  sum_zeros_pred_ <- numeric(collections)
  p95_pred_ <- matrix(0, p, collections)
  
  chi_sq_pred_ <- numeric(collections)
  chi_sq_etapred_ <- numeric(collections)
  chi_sq_obs_ <- numeric(collections)
  chi_sq_obs_etapred_ <- numeric(collections)
  Freeman_Tukey_obs_ <- numeric(collections)
  Freeman_Tukey_obs_etapred_ <- numeric(collections)
  Freeman_Tukey_pred_ <- numeric(collections)
  Freeman_Tukey_etapred_ <- numeric(collections)
  pvalue_ResgrRespred_sum <- matrix(0, p, t)
  
  
  YPRED_storage <- array(0, dim = c(p, t, collections))
  pvalue_ResgrRespred_sum <- matrix(0,p,t)
  
  
  tick <- 1 # Storage index
  
  # --- MCMC Loop ---
  message("Starting MCMC...")
  total_iter <- nrep + nburn
  
  # Pre-calculate for efficiency in loop
  eta_tilde_vec <- as.vector(eta_tilde)
  offset_vec <- as.vector(offset_mat)
  X_centered_vec <- as.vector(X_centered)
  
  t_bigG_mat <- Matrix::t(bigG_mat)
  t_bigG2_mat <- Matrix::t(bigG2_mat)
  t_bigG3_mat <- Matrix::t(bigG3_mat)
  
  MConf_X_bind <- if(measured_confounders) cbind(MConf_mat, X_centered_vec) else matrix(X_centered_vec, ncol=1)
  t_MConf_X_bind <- Matrix::t(MConf_X_bind)
  
  
  for (irep in 1:total_iter) {
    if (irep %% it_print == 0) {
      message(paste("Iteration:", irep, "/", total_iter))
    }
    
    eta_tilde_vec <- as.vector(eta_tilde)
    meanZ_vec <- as.vector(meanZ_mat)
    
    # --- Step I.a: Sample (B0spacetimedraw, Btimedraw) ---
    residual_a <- eta_tilde_vec - offset_vec -
      as.vector(bigG2_mat %*% Bdraw2_vec) -
      as.vector(bigG3_mat %*% Bdraw3_vec) -
      meanZ_vec - (X_centered_vec * Bdraw)
    
    # Prior precision for state evolution
    invOmega22_a <- bdiag(Q0invdraw_spacetime, Diagonal(m1, x = Q1invdraw_time))
    invS_a_list <- list(Diagonal(p) * 4, solve(B_0_prvar)) # Priors for t=1
    if (t > 1) {
      invS_a_list_kron <- rep(list(invOmega22_a), t-1)
      invS_a <- Matrix::bdiag(c(invS_a_list, invS_a_list_kron))
    } else {
      invS_a <- Matrix::bdiag(invS_a_list)
    }
    
    
    H_current <- if(interaction_type == 3 && !is.null(H2_mat)) H2_mat else H_mat
    K_a <- Matrix::t(H_current) %*% invS_a %*% H_current
    
    GinvOmega11_a <- t_bigG_mat / s2_err_mis_draw
    GinvOmega11G_a <- GinvOmega11_a %*% bigG_mat
    invP_a <- K_a + GinvOmega11G_a
    invP_a <- forceSymmetric(invP_a)
    
    chol_P_a <- as(Cholesky(invP_a, perm = F), 'Matrix')
    mean_param_a <- solve(t(chol_P_a), 
                          solve(chol_P_a, GinvOmega11_a %*% residual_a))
    Bdraw1_sampled_vec <- mean_param_a + solve(t(chol_P_a), rnorm(Tq))
    
    Bdraw1_vec <- Bdraw1_sampled_vec
    
    # Reshape and apply constraints
    Bdrawc_a <- matrix(Bdraw1_vec, nrow = q_val, ncol = t) # (p+1) x t
    B0spacetimedraw <- Bdrawc_a[1:m0, , drop = FALSE]      # p x t
    Btimedraw <- Bdrawc_a[(m0 + 1):q_val, , drop = FALSE] # 1 x t
    
    if (interaction_type == 3 || interaction_type == 4) { # Sum-to-zero over space for each time
      B0spacetimedraw <- B0spacetimedraw - matrix(colMeans(B0spacetimedraw), nrow = p, ncol = t, byrow = TRUE)
    }
    if (interaction_type == 2 || interaction_type == 4 || interaction_type == 5) { # Sum-to-zero over time for each space
      B0spacetimedraw <- B0spacetimedraw - matrix(rowMeans(B0spacetimedraw), nrow = p, ncol = t, byrow = FALSE)
    }
    if (interaction_type == 1) B0spacetimedraw <- B0spacetimedraw * 0 # Set to zero if interaction 1
    
    Btimedraw <- Btimedraw - mean(Btimedraw) # Sum-to-zero for Btimedraw (X-effect)
    
    Bdraw1_vec <- as.vector(rbind(B0spacetimedraw, Btimedraw))
    
    
    # --- Step I.b: Sample (B0timedraw, Bspacetimedraw) ---
    residual_b <- eta_tilde_vec - offset_vec -
      as.vector(bigG_mat %*% Bdraw1_vec) -
      as.vector(bigG3_mat %*% Bdraw3_vec) -
      meanZ_vec - (X_centered_vec * Bdraw)
    
    invOmega22_b <- bdiag(Diagonal(1, x=Q0invdraw_time), Q1invdraw_spacetime)
    invS_b_list <- list(matrix(1/4,1,1), Diagonal(p) * (1/10))
    if (t > 1) {
      invS_b_list_kron <- rep(list(invOmega22_b), t-1)
      invS_b <- Matrix::bdiag(c(invS_b_list, invS_b_list_kron))
    } else {
      invS_b <- Matrix::bdiag(invS_b_list)
    }
    
    K_b <- Matrix::t(H_mat) %*% invS_b %*% H_mat
    GinvOmega11_b <- t_bigG2_mat / s2_err_mis_draw
    GinvOmega11G_b <- GinvOmega11_b %*% bigG2_mat
    invP_b <- K_b + GinvOmega11G_b
    invP_b <- forceSymmetric(invP_b)
    
    chol_P_b <- as(Cholesky(invP_b, perm = F), 'Matrix')
    mean_param_b <- solve(t(chol_P_b), 
                          solve(chol_P_b, GinvOmega11_b %*% residual_b))
    Bdraw2_sampled_vec <- mean_param_b + solve(t(chol_P_b), rnorm(Tq))
    
    Bdraw2_vec <- Bdraw2_sampled_vec
    
    # Reshape and apply constraints
    Bdrawc_b <- matrix(Bdraw2_vec, nrow = q_val, ncol = t) # (1+p) x t
    B0timedraw <- Bdrawc_b[1, , drop = FALSE]                 # 1 x t
    Bspacetimedraw <- Bdrawc_b[2:q_val, , drop = FALSE]    # p x t
    
    B0timedraw <- B0timedraw - mean(B0timedraw) # Sum-to-zero
    
    Bspacetimedraw <- Bspacetimedraw - matrix(colMeans(Bspacetimedraw), nrow=p, ncol=t, byrow=TRUE)
    Bspacetimedraw <- Bspacetimedraw - matrix(rowMeans(Bspacetimedraw), nrow=p, ncol=t, byrow=FALSE)
    if (ST_interaction_X == 0) Bspacetimedraw <- Bspacetimedraw * 0
    
    Bdraw2_vec <- as.vector(rbind(B0timedraw, Bspacetimedraw))
    
    
    # --- Step I.c: Sample (B0spacedraw, Bspacedraw) ---
    residual_c <- eta_tilde_vec - offset_vec -
      as.vector(bigG_mat %*% Bdraw1_vec) -
      as.vector(bigG2_mat %*% Bdraw2_vec) -
      meanZ_vec - (X_centered_vec * Bdraw)
    
    K_c <- bdiag(Q0invdraw_space, Q1invdraw_space) # 2p x 2p
    GinvOmega11_c <- t_bigG3_mat / s2_err_mis_draw
    GinvOmega11G_c <- GinvOmega11_c %*% bigG3_mat
    invP_c <- K_c + GinvOmega11G_c
    invP_c <- forceSymmetric(invP_c)
    
    chol_P_c <- as(Cholesky(invP_c, perm = F), 'Matrix')
    mean_param_c <- solve(t(chol_P_c), 
                          solve(chol_P_c, GinvOmega11_c %*% residual_c))
    Bdraw3_sampled_vec <- mean_param_c + solve(t(chol_P_c), rnorm(2 * p))
    
    Bdraw3_vec <- Bdraw3_sampled_vec
    
    # Reshape and apply constraints
    Bdrawc_c <- matrix(Bdraw3_vec, nrow = p, ncol = 2) # p x 2
    B0spacedraw <- Bdrawc_c[, 1, drop = FALSE]          # p x 1
    Bspacedraw <- Bdrawc_c[, 2, drop = FALSE]           # p x 1
    
    B0spacedraw <- B0spacedraw - mean(B0spacedraw)
    Bspacedraw <- Bspacedraw - mean(Bspacedraw)
    
    Bdraw3_vec <- c(as.vector(B0spacedraw), as.vector(Bspacedraw))
    
    
    # --- Step I.d: Sample constant coefficients (gamma_draw for MConf, Bdraw for X) ---
    residual_d <- eta_tilde_vec - offset_vec -
      as.vector(bigG_mat %*% Bdraw1_vec) -
      as.vector(bigG2_mat %*% Bdraw2_vec) -
      as.vector(bigG3_mat %*% Bdraw3_vec)
    
    num_const_coeffs <- ncol(MConf_X_bind)
    K_d_list <- list()
    if (measured_confounders) K_d_list$MConf <- Diagonal(nzmconf) * 1e-12 # Prior for gamma_draw
    K_d_list$Xcoeff <- matrix(1/1, 1, 1)
    K_d <- Matrix::bdiag(K_d_list)
    
    GinvOmega11_d <- t_MConf_X_bind / s2_err_mis_draw
    GinvOmega11G_d <- GinvOmega11_d %*% MConf_X_bind
    invP_d <- K_d + GinvOmega11G_d
    invP_d <- forceSymmetric(invP_d)
    
    chol_P_d <- t(chol(invP_d))
    mean_param_d <- solve(t(chol_P_d), 
                          solve(chol_P_d, GinvOmega11_d %*% residual_d))
    gamma_B_draw_vec <- mean_param_d + solve(t(chol_P_d), rnorm(num_const_coeffs))
    
    Bdraw <- gamma_B_draw_vec[num_const_coeffs] # Last one is Bdraw for X
    if (measured_confounders) {
      gamma_draw <- matrix(gamma_B_draw_vec[1:nzmconf], ncol = 1)
      meanZ_mat <- matrix(MConf_mat %*% gamma_draw, nrow = p, ncol = t)
    } else {
      gamma_draw <- NULL # Ensure NULL if no measured confounders
      meanZ_mat <- matrix(0,p,t) # Reset if no confounders
    }
    
    
    # --- Calculate Mean Components for thetay_draw ---
    meanY1 <- X_centered * (Bdraw +
                              matrix(Btimedraw, nrow = p, ncol = t, byrow = TRUE) +
                              matrix(Bspacedraw, nrow = p, ncol = t, byrow = FALSE) +
                              Bspacetimedraw)
    
    meanY0 <- B0spacetimedraw +
      matrix(B0timedraw, nrow = p, ncol = t, byrow = TRUE) +
      matrix(B0spacedraw, nrow = p, ncol = t, byrow = FALSE)
    
    
    # --- Step III: Sample Variances of State Vectors & CAR parameters ---
    
    # Variance for B0spacetimedraw (Q0invdraw_spacetime related rho1_0_space_time_draw, rho2_0_space_time_draw)
    if (t > 1) {
      B0st_diff <- t(B0spacetimedraw[, 2:t, drop=FALSE]) - t(B0spacetimedraw[, 1:(t - 1), drop=FALSE]) # (t-1) x p
      B0st_temp_for_var <- rbind(matrix(0, nrow=1, ncol=p), B0st_diff) # t x p
    } else {
      B0st_temp_for_var <- matrix(0, nrow=1, ncol=p)
    }
    
    if (interaction_type == 2) { # Diagonal variance, no spatial correlation
      rho1_0_space_time_draw <- posterior_conditional_variance(t(B0st_temp_for_var), Diagonal(p), 0.01, 0.01, p, max(1,t-1))
      Q0invdraw_spacetime <- (1 / rho1_0_space_time_draw) * Diagonal(p)
    } else if (interaction_type == 3) { # Full CAR on B0spacetimedraw (not differenced version)
      # rho2 for B0spacetimedraw itself
      rho2_0_space_time_draw <- MH_spatial_correlation_CAR(W_adj, B0spacetimedraw, rho1_0_space_time_draw, 0.8, 1.0)
      pstar_0st <- ifelse(abs(rho2_0_space_time_draw - 1) < 1e-9, p - 1, p)
      Q0st_struct <- forceSymmetric(D_adj - rho2_0_space_time_draw * W_adj)
      rho1_0_space_time_draw <- posterior_conditional_variance(B0spacetimedraw, Q0st_struct, 0.01, 0.01, pstar_0st, t)
      Q0invdraw_spacetime <- (1 / rho1_0_space_time_draw) * Q0st_struct
    } else if (interaction_type == 4) { # Full CAR on differenced B0st_temp_for_var
      rho2_0_space_time_draw <- MH_spatial_correlation_CAR(W_adj, t(B0st_temp_for_var), rho1_0_space_time_draw, 0.8, 1.0)
      pstar_0st <- ifelse(abs(rho2_0_space_time_draw - 1) < 1e-9, p - 1, p)
      Q0st_struct <- forceSymmetric(D_adj - rho2_0_space_time_draw * W_adj)
      rho1_0_space_time_draw <- posterior_conditional_variance(t(B0st_temp_for_var), Q0st_struct, 0.01, 0.01, pstar_0st, max(1, t-1))
      Q0invdraw_spacetime <- (1 / rho1_0_space_time_draw) * Q0st_struct
    } else if (interaction_type == 5) { # Heterogeneous diagonal variances (p x 1)
      nu02_0st <- matrix(0.01, nrow = p, ncol = 1)
      S02_0st <- matrix(0.01, nrow = p, ncol = 1)
      newnu2_0st <- nu02_0st + (max(1,t-1) / 2)
      newS2_0st <- S02_0st + rowSums(t(B0st_temp_for_var)^2) / 2
      rho1_0_space_time_draw <- 1 / rgamma(p, shape = newnu2_0st, rate = newS2_0st)
      Q0invdraw_spacetime <- Diagonal(p, x = 1 / rho1_0_space_time_draw)
    }
    # If interaction_type == 1, rho1_0_space_time_draw is fixed small
    
    # Variance for B0timedraw (Q0invdraw_time) - scalar
    if (t > 1) {
      e2_0t <- B0timedraw[1, 2:t, drop=FALSE] - B0timedraw[1, 1:(t - 1), drop=FALSE]
      sum_sq_e2_0t <- sum(e2_0t^2)
    } else { sum_sq_e2_0t <- 0 }
    Q0invdraw_time <- rgamma(1, shape = (0.01 + max(1,t-1)/2), rate = (0.01 + sum_sq_e2_0t/2))
    
    # Variance for B0spacedraw (Q0invdraw_space related rho1_0_space_draw, rho2_0_space_draw)
    rho2_0_space_draw <- MH_spatial_correlation_CAR(W_adj, B0spacedraw, rho1_0_space_draw, 0.8, 1.0)
    pstar_0s <- ifelse(abs(rho2_0_space_draw-1)<1e-9, p-1, p)
    Q0s_struct <- forceSymmetric(D_adj - rho2_0_space_draw * W_adj)
    rho1_0_space_draw <- posterior_conditional_variance(B0spacedraw, Q0s_struct, 0.01, 0.01, pstar_0s, 1)
    Q0invdraw_space <- (1/rho1_0_space_draw) * Q0s_struct
    
    # Variance for Btimedraw (Q1invdraw_time) - scalar
    if (t > 1) {
      e2_1t <- Btimedraw[1, 2:t, drop=FALSE] - Btimedraw[1, 1:(t - 1), drop=FALSE]
      sum_sq_e2_1t <- sum(e2_1t^2)
    } else { sum_sq_e2_1t <- 0 }
    Q1invdraw_time <- rgamma(1, shape = (a1_prior + max(1,t-1)/2), rate = (b1_prior + sum_sq_e2_1t/2))
    
    # Variance for Bspacetimedraw (Q1invdraw_spacetime related rho1_space_time_draw, rho2_space_time_draw)
    if (ST_interaction_X == 1) {
      if (t > 1) {
        Bst_diff <- t(Bspacetimedraw[, 2:t, drop=FALSE]) - t(Bspacetimedraw[, 1:(t - 1), drop=FALSE])
        Bst_temp_for_var <- rbind(matrix(0,1,p), Bst_diff) # t x p
      } else { Bst_temp_for_var <- matrix(0,1,p) }
      rho2_space_time_draw <- MH_spatial_correlation_CAR(W_adj, t(Bst_temp_for_var), rho1_space_time_draw, 0.8, 1.0)
      pstar_st <- ifelse(abs(rho2_space_time_draw-1)<1e-9, p-1, p)
      Qst_struct <- forceSymmetric(D_adj - rho2_space_time_draw * W_adj)
      rho1_space_time_draw <- posterior_conditional_variance(t(Bst_temp_for_var), Qst_struct, 0.01, 0.01, pstar_st, max(1,t-1))
    } else {
      rho1_space_time_draw <- 1e-9
    }
    Q1invdraw_spacetime <- (1/rho1_space_time_draw) * forceSymmetric(D_adj - rho2_space_time_draw * W_adj)
    
    # Variance for Bspacedraw (Q1invdraw_space related rho1_space_draw, rho2_space_draw)
    rho2_space_draw <- MH_spatial_correlation_CAR(W_adj, Bspacedraw, rho1_space_draw, 0.8, 1.0)
    pstar_s <- ifelse(abs(rho2_space_draw-1)<1e-9, p-1, p)
    Qs_struct <- forceSymmetric(D_adj - rho2_space_draw * W_adj)
    rho1_space_draw <- posterior_conditional_variance(Bspacedraw, Qs_struct, a0_prior, b0_prior, pstar_s, 1)
    Q1invdraw_space <- (1/rho1_space_draw) * Qs_struct
    
    
    # --- Step V: Sample latent process eta_tilde (for Poisson) and measurement error variance ---
    thetay_draw <- meanY1 + offset_mat + meanY0 + meanZ_mat
    
    aug_data_result <- augmented_data_poisson_lognormal_noloop(Y, eta_tilde, thetay_draw, s2_err_mis_draw, pswitch_Y, ctuning_poisson)
    eta_tilde <- aug_data_result$eta_draw
    pswitch_Y <- aug_data_result$pswitch_Y # Accumulates acceptance rates
   
    yhat_residuals <- eta_tilde - thetay_draw
    sse_2 <- sum(yhat_residuals^2)
    g1_pos_s2err <- (p * t) * 0.5 + s2_a
    g2_pos_s2err <- s2_b + 0.5 * sse_2
    s2_err_mis_draw <- 1 / rgamma(1, shape = g1_pos_s2err, rate = g2_pos_s2err)
    
    
    # --- Store Results After Burn-in and Thinning ---
    if (irep > nburn && (irep - nburn) %% thin == 0) {
      # Accumulate posterior means in 'ave' list
      ave$Eta_tilde_mean <- ave$Eta_tilde_mean + eta_tilde
      ave$B_postmean <- ave$B_postmean + Bdraw
      ave$B2_postmean <- ave$B2_postmean + Bdraw^2
      ave$Btime_postmean <- ave$Btime_postmean + Btimedraw
      ave$Btime2_postmean <- ave$Btime2_postmean + Btimedraw^2
      ave$Bspace_postmean <- ave$Bspace_postmean + Bspacedraw
      ave$Bspace2_postmean <- ave$Bspace2_postmean + Bspacedraw^2
      ave$Bspacetime_postmean <- ave$Bspacetime_postmean + Bspacetimedraw
      ave$Bspacetime2_postmean <- ave$Bspacetime2_postmean + Bspacetimedraw^2
      
      ave$B2_c_t = ave$B2_c_t + (Bdraw + Btimedraw)^2;
      ave$B2_c_s = ave$B2_c_s + (Bdraw + Bspacedraw)^2;
      ave$B2_c_t_s_st = ave$B2_c_t_s_st + (Bdraw + 
                                             matrix(Btimedraw, nrow = p, ncol = t, byrow = TRUE) +
                                             matrix(Bspacedraw, nrow = p, ncol = t, byrow = FALSE) +
                                             Bspacetimedraw)^2;
      ave$B2_c_t_st = ave$B2_c_t_st + (Bdraw + 
                                         matrix(Btimedraw, nrow = p, ncol = t, byrow = TRUE) +
                                         Bspacetimedraw)^2;
      ave$B2_t_st = ave$B2_t_st + (matrix(Btimedraw, nrow = p, ncol = t, byrow = TRUE) +
                                     Bspacetimedraw)^2;
      ave$B2_t_s = ave$B2_t_s + (matrix(Btimedraw, nrow = p, ncol = t, byrow = TRUE) +
                                   matrix(Bspacedraw, nrow = p, ncol = t, byrow = FALSE))^2;
      ave$B2_t_s_st = ave$B2_t_s_st + (matrix(Btimedraw, nrow = p, ncol = t, byrow = TRUE) +
                                         matrix(Bspacedraw, nrow = p, ncol = t, byrow = FALSE) +
                                         Bspacetimedraw)^2;
      ave$B2_c_t_s = ave$B2_c_t_s + (Bdraw + 
                                       matrix(Btimedraw, nrow = p, ncol = t, byrow = TRUE) +
                                       matrix(Bspacedraw, nrow = p, ncol = t, byrow = FALSE))^2;
      ave$B2_s_st = ave$B2_s_st + (matrix(Bspacedraw, nrow = p, ncol = t, byrow = FALSE) + 
                                     Bspacetimedraw)^2;
      ave$B2_c_s_st = ave$B2_c_s_st + (Bdraw + 
                                         matrix(Bspacedraw, nrow = p, ncol = t, byrow = FALSE) + 
                                         Bspacetimedraw)^2;
      
      ave$B0spacetime_postmean <- ave$B0spacetime_postmean + B0spacetimedraw
      ave$B0spacetime2_postmean <- ave$B0spacetime2_postmean + B0spacetimedraw^2
      ave$B0time_postmean <- ave$B0time_postmean + B0timedraw
      ave$B0time2_postmean <- ave$B0time2_postmean + B0timedraw^2
      ave$B0space_postmean <- ave$B0space_postmean + B0spacedraw
      ave$B0space2_postmean <- ave$B0space2_postmean + B0spacedraw^2
      
      ave$meanZmean <- ave$meanZmean + meanZ_mat
      ave$meanY1mean <- ave$meanY1mean +  meanY1
      ave$meanY0mean <- ave$meanY0mean +  meanY0
      
      # Calculate likelihood and predictions for this iteration
      llike_iter <- 0
      Ypred_iter <- matrix(0, p, t)
      mu_hat_iter <- matrix(0,p,t) # E[Y | params] or E[eta | params]
      V_hat_iter <- matrix(0,p,t)  # Var[Y | params] or Var[eta | params]
      
      curr_llike = matrix(dpois(Y, lambda = exp(eta_tilde), log = TRUE), p, t)
      llike_iter <- sum( curr_llike)
      
      # For prediction Ypred and diagnostics mu_hat, V_hat
      mu_hat_iter <- exp(eta_tilde) # E[Y | eta_tilde] for Poisson
      V_hat_iter <- mu_hat_iter      # Var[Y | eta_tilde] for Poisson
      Ypred_iter <- matrix(rpois(p * t, lambda = as.vector(mu_hat_iter)), nrow = p, ncol = t)
      Ypred2_iter <- matrix(rpois(p * t, lambda = as.vector(mu_hat_iter)), nrow = p, ncol = t) # For CRPS
      
      store_LS_sum_per_iter = store_LS_sum_per_iter + exp(curr_llike)
      
      mu_hat_etapred_iter = exp(thetay_draw + 0.5*s2_err_mis_draw)
      V_hat_etapred_iter = mu_hat_etapred_iter + (exp(s2_err_mis_draw)-1) * exp(2*thetay_draw + s2_err_mis_draw)
      
      
      YPRED_storage[, , tick] <- Ypred_iter
      
      # Diagnostics
      # Pearson residuals: (Y - E[Y|params]) / sqrt(Var[Y|params])
      valid_idx <- !is.na(Y)
      pearson_res_iter <- (Y[valid_idx] - mu_hat_iter[valid_idx]) / sqrt(V_hat_iter[valid_idx])
      pearson_res_pred_iter <- (Ypred_iter[valid_idx] - mu_hat_iter[valid_idx]) / sqrt(V_hat_iter[valid_idx])
      pearson_res_etapred_iter <- (Y[valid_idx] - mu_hat_etapred_iter[valid_idx]) / sqrt(V_hat_etapred_iter[valid_idx])
      pearson_res_etapred_pred_iter <- (Ypred_iter[valid_idx] - mu_hat_etapred_iter[valid_idx]) / sqrt(V_hat_etapred_iter[valid_idx])
      
      
      chi_sq_obs_[tick] <- sum(pearson_res_iter^2, na.rm=TRUE)
      chi_sq_pred_[tick] <- sum(pearson_res_pred_iter^2, na.rm=TRUE)
      
      chi_sq_obs_etapred_[tick] <- sum(pearson_res_etapred_iter^2, na.rm=TRUE)
      chi_sq_etapred_[tick] <- sum(pearson_res_etapred_pred_iter^2, na.rm=TRUE)
      
      # Accumulate pvalue_ResgrRespred_sum element-wise
      temp_pval_matrix <- matrix(NA, p,t)
      temp_pval_matrix[valid_idx] <- (pearson_res_iter^2 >= pearson_res_pred_iter^2)
      pvalue_ResgrRespred_sum <- pvalue_ResgrRespred_sum + temp_pval_matrix # Will average later
      
      
      # Freeman-Tukey for Poisson
      # FT = (sqrt(Y_obs) - sqrt(E[Y_obs]))^2
      Freeman_Tukey_obs_[tick] <- sum( (sqrt(Y[valid_idx]) - sqrt(mu_hat_iter[valid_idx]))^2 , na.rm=TRUE)
      Freeman_Tukey_pred_[tick]<- sum( (sqrt(Ypred_iter[valid_idx]) - sqrt(mu_hat_iter[valid_idx]))^2 , na.rm=TRUE)
      Freeman_Tukey_obs_etapred_[tick] <- sum( (sqrt(Y[valid_idx]) - sqrt(mu_hat_etapred_iter[valid_idx]))^2 , na.rm=TRUE)
      Freeman_Tukey_etapred_[tick] <- sum( (sqrt(Ypred_iter[valid_idx]) - sqrt(mu_hat_etapred_iter[valid_idx]))^2 , na.rm=TRUE)
      
      sum_zeros_pred_[tick] <- sum(Ypred_iter == 0, na.rm = TRUE)
      p95_pred_[, tick] <- apply(Ypred_iter, 1, quantile, probs = 0.95, na.rm = TRUE)
      
      # CRPS components (sum over all non-NA Y)
      store_CRPS_1_sum_per_iter[tick] <- sum(abs(Ypred_iter[valid_idx] - Ypred2_iter[valid_idx]), na.rm=TRUE)
      store_CRPS_2_sum_per_iter[tick] <- sum(abs(Ypred_iter[valid_idx] - Y[valid_idx]), na.rm=TRUE)
      
      # Store samples in 'out' list
      out$store_llike[tick] <- llike_iter
      out$Bdraw[tick] <- Bdraw
      out$Bspacedraw[, tick] <- as.vector(Bspacedraw)
      if (measured_confounders && !is.null(gamma_draw)) {
        out$G_[, tick] <- as.vector(gamma_draw)
      }
      out$S2_err_mis_[tick] <- s2_err_mis_draw
      out$RHO1_space_[tick] <- rho1_space_draw
      out$RHO1_space_time_[tick] <- rho1_space_time_draw
      if (interaction_type == 5) {
        out$RHO1_0_space_time_[, tick] <- as.vector(rho1_0_space_time_draw)
      } else {
        out$RHO1_0_space_time_[tick] <- rho1_0_space_time_draw
      }
      out$RHO1_0_space_[tick] <- rho1_0_space_draw
      out$RHO2_space_[tick] <- rho2_space_draw
      out$RHO2_space_time_[tick] <- rho2_space_time_draw
      out$RHO2_0_space_time_[tick] <- rho2_0_space_time_draw
      out$RHO2_0_space_[tick] <- rho2_0_space_draw
      out$Q1inv_time[tick] <- Q1invdraw_time
      out$Q0inv_time[tick] <- Q0invdraw_time
      
      out$RMSE_[tick] <- sqrt(mean((Y[valid_idx] - Ypred_iter[valid_idx])^2, na.rm=TRUE))
      out$MAE_[tick] <- mean(abs(Y[valid_idx] - Ypred_iter[valid_idx]), na.rm=TRUE)
      
      tick <- tick + 1
    } # End of storing results
    
    # Store last state for potential restart
    if (irep == total_iter) {
      out$eta_tilde <- eta_tilde
      out$B0timedraw <- B0timedraw; out$B0spacedraw <- B0spacedraw
      out$Bspacetimedraw <- Bspacetimedraw # pxt
      out$Btimedraw <- Btimedraw
    }
    
  } # End of MCMC loop
  
  # --- Post-processing ---
  num_samples_collected = tick-1
  # Calculate averages from 'ave' list
  fields_to_average <- names(ave)
  for(field in fields_to_average){
    if(is.numeric(ave[[field]]) || is.matrix(ave[[field]]) || is.array(ave[[field]])){
      ave[[field]] <- ave[[field]] / num_samples_collected
    }
  }
  
  # Final calculations for DIC, PMCC, etc.
  mean_Dbar <- -2 * mean(out$store_llike[1:num_samples_collected], na.rm=TRUE)
  
  # D_hat: -2 * logLik(theta_hat)
  # Using posterior means for theta_hat
  S2_ERR_mean <- mean(out$S2_err_mis_[1:num_samples_collected], na.rm=TRUE)
  
  # Likelihood at posterior mean parameters
  llike_at_mean_params <- 0
  mu_at_mean_eta <- exp(ave$Eta_tilde_mean)
  llike_at_mean_params <- sum(dpois(Y, lambda = mu_at_mean_eta, log = TRUE), na.rm=TRUE)
  if(!is.null(ave$rate_acc_Y)) ave$rate_acc_Y <- pswitch_Y / total_iter # Mean acceptance for Poisson eta
  
  D_hat <- -2 * llike_at_mean_params
  
  ave$Dbar <- mean_Dbar
  ave$pd <- mean_Dbar - D_hat
  ave$DIC <- ave$Dbar + ave$pd
  
  # PMCC
  YPRED_mean_overall <- apply(YPRED_storage[,,1:num_samples_collected, drop=FALSE], c(1,2), mean, na.rm=TRUE)
  YPRED_var_overall <- apply(YPRED_storage[,,1:num_samples_collected, drop=FALSE], c(1,2), var, na.rm=TRUE)
  YPRED_var_overall[is.na(YPRED_var_overall)] <- 0 # Handle cases with too few samples for var
  
  valid_idx_pmcc <- !is.na(Y) & !is.na(YPRED_mean_overall) & !is.na(YPRED_var_overall)
  ave$PMCC <- sum( (Y[valid_idx_pmcc] - YPRED_mean_overall[valid_idx_pmcc])^2 + YPRED_var_overall[valid_idx_pmcc] , na.rm=TRUE)
  ave$meanYpredmean <- YPRED_mean_overall
  
  # Other diagnostics
  ave$chisq_pvalue <- mean(chi_sq_obs_[1:num_samples_collected] >= chi_sq_pred_[1:num_samples_collected], na.rm=TRUE)
  # For average p-value over all valid entries:
  ave$pvalue_ResgrRespred_overall <- mean(pvalue_ResgrRespred_sum[valid_idx] / num_samples_collected, na.rm=TRUE)
  
  
  ave$Freeman_Tukey_pvalue <- mean(Freeman_Tukey_obs_[1:num_samples_collected] > Freeman_Tukey_pred_[1:num_samples_collected], na.rm=TRUE)
  ave$Freeman_Tukey_etapred_pvalue <- mean(Freeman_Tukey_obs_etapred_[1:num_samples_collected] > Freeman_Tukey_etapred_[1:num_samples_collected], na.rm=TRUE)
  # Bayesian p-value for Y > Y_pred
  Y_rep <- array(Y, dim=c(p,t,num_samples_collected)) # Replicate Y
  ave$pvalue_YgrYhat <- apply( (Y_rep > YPRED_storage[,,1:num_samples_collected]) + 0.5 * (Y_rep == YPRED_storage[,,1:num_samples_collected]), 1:2, mean)
  
  ave$prop_zeros_pvalue <- mean(sum_zeros_obs_ >= sum_zeros_pred_[1:num_samples_collected], na.rm=TRUE)
  ave$percentile95_pvalue <- rowMeans(matrix(p95_obs_,p,num_samples_collected) >= p95_pred_[,1:num_samples_collected, drop=FALSE], na.rm=TRUE) # p-vector
  
  # CRPS
  num_valid_obs_for_crps <- sum(valid_idx)
  ave$CRPS <- mean( (0.5 * store_CRPS_1_sum_per_iter[1:num_samples_collected] - store_CRPS_2_sum_per_iter[1:num_samples_collected]) / num_valid_obs_for_crps, na.rm=TRUE)
  
  # Log Score (LS
  ave$LS <- mean(log(store_LS_sum_per_iter / num_samples_collected)) 
  
  # Quantile Score (QS)
  percentiles_qs <- c(2.5, 25, 50, 75, 97.5)
  alpha_qs <- percentiles_qs / 100
  
  ave$QS <- numeric(length(alpha_qs))
  for(k_qs in 1:length(alpha_qs)){
    Q_alpha <- apply(YPRED_storage[,,1:num_samples_collected, drop=FALSE], c(1,2), quantile, probs=alpha_qs[k_qs], na.rm=TRUE)
    ave$QS[k_qs] <- mean( (Y - Q_alpha) * ( (Y <= Q_alpha) - alpha_qs[k_qs] ), na.rm=TRUE)
  }
  
  return(list(ave = ave, out = out))
}