
init_model <- function(inc) {
    
    # --------- Run Model ---------
    mod <- glmer(cases ~ offset(log(py)) + cov + (1 | agr) + (1 + cov | hier3.id), family = "poisson", data = inc)

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
    n_space <- max(inc$hier3.id)
    n_ages  <- max(inc$agr)
    
    space_int_inits   <- rep(0, n_space)
    space_slope_inits <- rep(0, n_space)
    age_inits         <- rep(0, n_ages)
    
    # Space
    re_space <- ranef(mod)$hier3.id
    space_idxs <- as.numeric(rownames(re_space))
    
    space_int_inits[space_idxs] <- re_space$"(Intercept)"
    
    if (method == "complex") {
        space_slope_inits[space_idxs] <- re_space$"cov"
    }
    
    # Age
    re_age <- ranef(mod)$agr
    age_idxs <- as.numeric(rownames(re_age))
    age_inits[age_idxs] <- re_age$"(Intercept)"
    
    # --------- Package ---------
    print(summary(mod))
    inits <- list(intercept, slope, space_int_inits, space_slope_inits, age_inits)
    names(inits) <- c(
        "global.intercept", "global.slope", "space.intercepts",
        "space.slopes", "age.intercepts"
    )

    return(inits)
}
