
grid_lookup <- function(incidence) {
    grid.lookup <- incidence %>%
        select(hier3.id, hier2.id, hier1.id) %>%
        distinct() %>%
        arrange(hier3.id)
    
    grid.lookup.s2 <- incidence %>%
        select(hier2.id, hier1.id) %>%
        distinct() %>%
        arrange(hier2.id)
    
    grid.lookup <- as.matrix(grid.lookup)
    grid.lookup.s2 <- as.matrix(grid.lookup.s2)
    
    return(list("grid.lookup" = grid.lookup, "grid.lookup.s2" = grid.lookup.s2))
}


init_model <- function(incidence) {
    
    # --------- Run Model ---------
    mod <- glmer(cases ~ offset(log(py)) + cov + (1 | agr) + (1 + cov | hier3.id), family = "poisson", data = incidence)
    
    # --------- Extract Parameters ---------
    # --- Fixed Effects ---
    fixed_table <- coef(summary(mod))
    if ("Estimate" %in% colnames(fixed_table)) {
        vals <- fixed_table[, "Estimate"]
    } else {
        vals <- fixed_table[, 1]
    }
    
    intercept <- vals[1]
    slope     <- vals[2]
    
    # --- Random Effects ---
    n_space <- max(incidence$hier3.id)
    n_ages  <- max(incidence$agr)
    
    space_int_inits   <- rep(0, n_space)
    space_slope_inits <- rep(0, n_space)
    age_inits         <- rep(0, n_ages)
    
    # Space
    re_space <- ranef(mod)$hier3.id
    space_idxs <- as.numeric(rownames(re_space))
    
    space_int_inits[space_idxs] <- re_space$"(Intercept)"
    
    space_slope_inits[space_idxs] <- re_space$"cov"
    
    # Age
    re_age <- ranef(mod)$agr
    age_idxs <- as.numeric(rownames(re_age))
    age_inits[age_idxs] <- re_age$"(Intercept)"
    
    # --------- Package ---------
    inits <- list(intercept, slope, space_int_inits, space_slope_inits, age_inits)
    names(inits) <- c(
        "global.intercept", "global.slope", "space.intercepts",
        "space.slopes", "age.intercepts"
    )
    
    return(inits)
}

thin_plate <- function(x, k) {
    
    knots<-quantile(x, seq(0, 1, length=(k+2))[-c(1,(k+2))])
    
    Z_K <- (abs(outer(x, knots, "-")))^3
    OMEGA_all <- (abs(outer(knots, knots, "-")))^3
    svd.OMEGA_all <- svd(OMEGA_all)
    sqrt.OMEGA_all <- t(svd.OMEGA_all$v %*% (t(svd.OMEGA_all$u)*sqrt(svd.OMEGA_all$d)))
    Z <- t(solve(sqrt.OMEGA_all,t(Z_K)))
    Z <- Z / max(Z)
    
    return(Z)
}


detect_fixed_nodes <- function(data, group_col, target_col, return_col) {
    
    fixed_ids <- data %>%
        group_by(.data[[group_col]]) %>%
        filter(n_distinct(.data[[target_col]]) == 1) %>%
        ungroup() %>%
        pull(all_of(return_col)) %>%
        unique()
    
    max_id <- max(data[[return_col]], na.rm = TRUE)
    fixed_vector <- rep(FALSE, max_id)
    
    fixed_vector[fixed_ids] <- TRUE
    
    return(fixed_vector)
}


generate_mcmc_report <- function(mcmc_raw, out_file, waic_file) {
    
    cat("\n [Report] Extracting parameter samples from nested NIMBLE output...\n")
    params_only <- lapply(mcmc_raw, function(chain) {
        as.mcmc(chain$samples2) 
    })
    
    mcmc_list <- mcmc.list(params_only)
    combined_chains <- as.matrix(mcmc_list)
    
    cat(" [Report] Calculating medians, and 95% Credible Intervals...\n")
    summary_stats <- tibble(
        Parameter = colnames(combined_chains),
        Median    = apply(combined_chains, 2, median),
        CrI_2.5   = apply(combined_chains, 2, quantile, probs = 0.025),
        CrI_97.5  = apply(combined_chains, 2, quantile, probs = 0.975)
    )
    
    cat(" [Report] Calculating convergence diagnostics (R-hat and ESS)...\n")
    rhat_vals <- tryCatch({
        gelman.diag(mcmc_list, multivariate = FALSE)$psrf[, "Point est."]
    }, error = function(e) {
        warning("Gelman diagnostic failed on some parameters. Filling with NA.")
        rep(NA, ncol(combined_chains))
    })
    
    ess_vals <- effectiveSize(mcmc_list)
    
    summary_stats <- summary_stats %>%
        mutate(
            Rhat = rhat_vals[Parameter],
            ESS  = ess_vals[Parameter]
        )
    
    top_level_params <- summary_stats %>%
        filter(!grepl("\\[", Parameter)) %>% 
        mutate(across(where(is.numeric), ~ round(., 4)))
    
    # --- EXTRACT AND SAVE WAIC ---
    waic_vals <- sapply(mcmc_raw, function(chain) {
        if (!is.null(chain$WAIC) && !is.null(chain$WAIC$WAIC)) chain$WAIC$WAIC else NA
    })
    mean_waic <- mean(waic_vals, na.rm = TRUE)
    
    # Create a dataframe for WAIC
    waic_df <- data.frame(
        Source = c(paste("Chain", seq_along(waic_vals)), "Mean"),
        WAIC = c(waic_vals, mean_waic)
    )
    
    # Export WAIC summary
    write.csv(waic_df, waic_file, row.names = FALSE)
    cat(paste0(" [Report] WAIC summary successfully saved to '", waic_file, "'\n"))
    
    # Export the FULL parameter summary
    write.csv(summary_stats, out_file, row.names = FALSE)
    cat(paste0(" [Report] Full parameter summary successfully saved to '", out_file, "'\n"))
    
    # Print Summary to Console
    cat("\n--- Key Metrics ---\n")
    max_rhat <- suppressWarnings(max(summary_stats$Rhat, na.rm = TRUE))
    cat(paste0(" max Rhat = ", ifelse(is.finite(max_rhat), round(max_rhat, 4), NA), "\n"))
    cat(paste0(" waic     = ", round(mean_waic, 4), "\n"))
    
    cat("\n--- Top-Level Parameters Summary ---\n")
    print(as.data.frame(top_level_params), row.names = FALSE)
    
    invisible(summary_stats)
}

