library(dplyr)

get_high_res_splines <- function(K = 4, n_points = 100) {
    knots <- quantile(1:18, seq(0, 1, length = (K + 2))[-c(1, (K + 2))])
    
    OMEGA_all <- (abs(outer(knots, knots, "-")))^3
    svd_O <- svd(OMEGA_all)
    sqrt_O <- t(svd_O$v %*% (t(svd_O$u) * sqrt(svd_O$d)))
    
    Z_K_orig <- (abs(outer(1:18, knots, "-")))^3
    Z_18 <- t(solve(sqrt_O, t(Z_K_orig)))
    max_Z <- max(Z_18)
    Z_18 <- Z_18 / max_Z
    
    agr_grid <- seq(1, 18, length.out = n_points)
    Z_K_high <- (abs(outer(agr_grid, knots, "-")))^3
    Z_high <- t(solve(sqrt_O, t(Z_K_high))) / max_Z
    
    list(grid = agr_grid, Z_high = Z_high, Z_18 = Z_18, knots = knots, K = K)
}

# Update this function in R/utils_splines.R

predict_highres_curve <- function(mcmc_mat, spline_obj, b0_col, b1_col, bs_prefix, target_sex, gamma_col = "gamma_sex", baseline_only = FALSE) {
    
    if (!b0_col %in% colnames(mcmc_mat)) return(NULL)
    
    b0 <- mcmc_mat[, b0_col]
    b1 <- mcmc_mat[, b1_col]
    
    # Sex term matches your exact working logic: (target_sex - 1)
    # If baseline_only is TRUE, we force the multiplier to 0, completely ignoring gamma_sex
    gamma_chain <- if(gamma_col %in% colnames(mcmc_mat)) mcmc_mat[, gamma_col] else rep(0, nrow(mcmc_mat))
    sex_multiplier <- if(baseline_only) 0 else (target_sex - 1)
    
    K <- spline_obj$K
    bs <- matrix(0, nrow = nrow(mcmc_mat), ncol = K)
    for (k in 1:K) {
        bs[, k] <- mcmc_mat[, paste0(bs_prefix, k, "]")]
    }
    
    # Core Math 
    base_term   <- b0 + (sex_multiplier * gamma_chain)
    age_term    <- outer(b1, spline_obj$grid, "*")
    spline_term <- bs %*% t(spline_obj$Z_high)
    
    log_rate_mat <- base_term + age_term + spline_term
    
    # Quantiles -> Exp -> 100k
    quants <- apply(log_rate_mat, 2, quantile, probs = c(0.025, 0.5, 0.975))
    
    data.frame(
        agr = spline_obj$grid,
        pred_raw_lo = quants[1, ],
        pred_raw    = quants[2, ],
        pred_raw_hi = quants[3, ],
        pred_rate_100k_lo = exp(quants[1, ]) * 100000,
        pred_rate_100k    = exp(quants[2, ]) * 100000,
        pred_rate_100k_hi = exp(quants[3, ]) * 100000
    )
}

# Add `baseline_only` to the wrapper arguments
generate_hierarchy_curves <- function(mcmc_mat, map_df, level_name, id_col, b0_prefix, b1_prefix, bs_prefix, spline_obj, target_sex, baseline_only = FALSE) {
    results <- lapply(1:nrow(map_df), function(i) {
        id <- map_df[[id_col]][i]
        name <- map_df[[level_name]][i]
        
        curve <- predict_highres_curve(
            mcmc_mat, spline_obj, 
            b0_col = paste0(b0_prefix, "[", id, "]"), 
            b1_col = paste0(b1_prefix, "[", id, "]"), 
            bs_prefix = paste0(bs_prefix, "[", id, ", "),
            target_sex = target_sex,
            baseline_only = baseline_only
        )
        if(is.null(curve)) return(NULL)
        
        curve[[level_name]] <- name
        if("region" %in% names(map_df)) curve$region <- map_df$region[i]
        if("continent" %in% names(map_df)) curve$continent <- map_df$continent[i]
        
        return(curve)
    })
    bind_rows(results)
}

# Add `baseline_only` to the global curve generator
generate_global_curve <- function(mcmc_mat, spline_obj, target_sex, gamma_col = "gamma_sex", baseline_only = FALSE) {
    curve <- predict_highres_curve(
        mcmc_mat, spline_obj, 
        b0_col = "beta0", b1_col = "beta1", bs_prefix = "b0[", 
        target_sex = target_sex, gamma_col = gamma_col, baseline_only = baseline_only
    )
    if(!is.null(curve)) curve$level <- "Global"
    return(curve)
}

# Calculates Standardised ASIR trajectory across Urbanisation
predict_asir_trajectory <- function(b0, b1, b_vec, delta, gamma_sex, 
                                    x_grid_raw, urb_mean, urb_sd, 
                                    Z_18, age_weights, target_sex) {
    
    x_grid_std <- (x_grid_raw - urb_mean) / urb_sd
    b_vec <- matrix(b_vec, ncol = 1)
    
    calc_base_rate <- function(sex_mult) {
        log_rate <- b0 + b1 * (1:18) + (Z_18 %*% b_vec) + (sex_mult * gamma_sex)
        exp(log_rate) * 100000
    }
    
    # Handle Sex logic natively (0=Both, 1=M, 2=F)
    if (target_sex == 0) {
        rate_m <- calc_base_rate(0) # Male multiplier = 0
        rate_f <- calc_base_rate(1) # Female multiplier = 1
        asir_base <- (sum(rate_m * age_weights) + sum(rate_f * age_weights)) / 2
    } else {
        rate <- calc_base_rate(target_sex - 1)
        asir_base <- sum(rate * age_weights)
    }
    
    asir_base * exp(delta * x_grid_std)
}

