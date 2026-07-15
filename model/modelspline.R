run_MCMC_spline <- function(seed,
                            model_name,
                            incidence,
                            init_vals,
                            n_iter,
                            n_burn,
                            inputs) {
    
    # ----- BUILD MODEL -----
    code <- nimbleCode({
        
        # ----- HYPERPRIORS -----
        sigma_beta0s1 ~ T(dnorm(0, sd=1.0), 0, Inf)
        sigma_beta0s2 ~ T(dnorm(0, sd=1.0), 0, Inf)
        sigma_beta0s3 ~ T(dnorm(0, sd=1.0), 0, Inf)
        
        sigma_beta1s1 ~ T(dnorm(0, sd=1.0), 0, Inf)
        sigma_beta1s2 ~ T(dnorm(0, sd=1.0), 0, Inf)
        sigma_beta1s3 ~ T(dnorm(0, sd=1.0), 0, Inf)
        
        sigma_bs1 ~ T(dnorm(0, sd=1.0), 0, Inf)
        sigma_bs2 ~ T(dnorm(0, sd=1.0), 0, Inf)
        sigma_bs3 ~ T(dnorm(0, sd=1.0), 0, Inf)
        
        sigma_deltas1 ~ T(dnorm(0, sd=1.0), 0, Inf)
        sigma_deltas2 ~ T(dnorm(0, sd=1.0), 0, Inf)
        sigma_deltas3 ~ T(dnorm(0, sd=sigma), 0, Inf)
        
        # ----- PRIORS -----
        r ~ T(dnorm(0, sd=5.0), 0, Inf)
        
        delta0 ~ dnorm(0, sd=1)
        
        beta0  ~ dnorm(0, sd=10)
        beta1  ~ dnorm(0, sd=10)
        
        for (k in 1:K) {
            b0[k] ~ dnorm(0, sd = 1.0)
        }
        
        gamma_sex ~ dnorm(0, sd = 1.0)
        
        # ----- HIERARCHICAL SPATIAL TERMS -----
        
        # HIERARCHY TIER 1
        for (s1 in 1:N_s1) {
            deltas1[s1] ~ dnorm(delta0, sd = sigma_deltas1)
            beta0s1[s1] ~ dnorm(beta0,  sd = sigma_beta0s1)
            beta1s1[s1] ~ dnorm(beta1,  sd = sigma_beta1s1)
            for (k in 1:K) {
                bs1[s1, k] ~ dnorm(b0[k], sd = sigma_bs1)
            }
        }
        
        # HIERARCHY TIER 2
        # for (i in 1:N_s2) {
        #     deltas2[i] ~ dnorm(deltas1[grid.lookup.s2[i, 2]], sd = sigma_deltas2)
        # 
        #     beta0s2[i] ~ dnorm(beta0s1[grid.lookup.s2[i, 2]], sd = sigma_beta0s2)
        #     beta1s2[i] ~ dnorm(beta1s1[grid.lookup.s2[i, 2]], sd = sigma_beta1s2)
        #     for (k in 1:K) {
        #         bs2[i, k] ~ dnorm(bs1[grid.lookup.s2[i, 2], k], sd = sigma_bs2)
        #     }
        # }
        
        # Fixed indices
        for (i in 1:N_hier2fix) {
            deltas2[hier2fix[i]] <- deltas1[grid.lookup.s2[hier2fix[i], 2]]
            
            beta0s2[hier2fix[i]] <- beta0s1[grid.lookup.s2[hier2fix[i], 2]]
            beta1s2[hier2fix[i]] <- beta1s1[grid.lookup.s2[hier2fix[i], 2]]
            for (k in 1:K) {
                bs2[hier2fix[i], k] <- bs1[grid.lookup.s2[hier2fix[i], 2], k]
            }
        }
        
        # Variable indices
        for (i in 1:N_hier2var) {
            deltas2[hier2var[i]] ~ dnorm(deltas1[grid.lookup.s2[hier2var[i], 2]], sd = sigma_deltas2)
            
            beta0s2[hier2var[i]] ~ dnorm(beta0s1[grid.lookup.s2[hier2var[i], 2]], sd = sigma_beta0s2)
            beta1s2[hier2var[i]] ~ dnorm(beta1s1[grid.lookup.s2[hier2var[i], 2]], sd = sigma_beta1s2)
            for (k in 1:K) {
                bs2[hier2var[i], k] ~ dnorm(bs1[grid.lookup.s2[hier2var[i], 2], k], sd = sigma_bs2)
            }
        }
        
        # HIERARCHY TIER 3
        for (s in 1:N_space) {
            deltas3[s] ~ dnorm(deltas2[grid.lookup[s, 2]], sd = sigma_deltas3)
            
            beta0s3[s] ~ dnorm(beta0s2[grid.lookup[s, 2]], sd = sigma_beta0s3)
            beta1s3[s] ~ dnorm(beta1s2[grid.lookup[s, 2]], sd = sigma_beta1s3)
            for (k in 1:K) {
                bs3[s, k] ~ dnorm(bs2[grid.lookup[s, 2], k], sd = sigma_bs3)
            }
        }
        
        # ----- LIKELIHOOD -----
        for (i in 1:N) {
            
            lograte[i] <- beta0s3[space[i]] + beta1s3[space[i]] * age[i] + 
                inprod(bs3[space[i], 1:K], Z[age[i], 1:K]) +
                gamma_sex * sex[i] + 
                deltas3[space[i]] * predictor[i]
            
            y[i] ~ dnegbin(r / (r + (n[i] * exp(lograte[i]))), r) 
        }
    })
    
    constants <- list(
        N = nrow(incidence),
        N_age_groups = max(incidence$agr),
        N_space = max(incidence$hier3.id),
        N_s1    = max(incidence$hier1.id),
        N_s2    = max(incidence$hier2.id),
        age   = incidence$agr,
        space = incidence$hier3.id,
        sex   = incidence$sex - 1, # 0: Men, 1: Women
        predictor = incidence$cov,
        grid.lookup    = inputs$grid.lookup,
        grid.lookup.s2 = inputs$grid.lookup.s2,
        hier2var       = inputs$var.hier2,
        hier2fix       = inputs$fix.hier2,
        N_hier2var     = length(inputs$var.hier2),
        N_hier2fix     = length(inputs$fix.hier2),
        K = inputs$spline.k,
        Z = inputs$spline.z,
        sigma = inputs$sigma
    )
    
    inits <- list(
        
        gamma_sex = 0,
        
        beta0  = init_vals$global.intercept,
        beta1  = 0,
        delta0 = init_vals$global.slope, 
        
        beta0s1 = rep(0, constants$N_s1),
        beta1s1 = rep(0, constants$N_s1),
        deltas1 = rep(0, constants$N_s1),
        
        beta0s2 = rep(0, constants$N_s2),
        beta1s2 = rep(0, constants$N_s2),
        deltas2 = rep(0, constants$N_s2),
        
        beta0s3 = init_vals$space.intercepts,
        beta1s3 = rep(0, constants$N_space),
        deltas3  = init_vals$space.slopes,
        
        b0 = rep(0, constants$K),
        bs1 = matrix(0, nrow = constants$N_s1, ncol = constants$K),
        bs2 = matrix(0, nrow = constants$N_s2, ncol = constants$K),
        bs3 = matrix(0, nrow = constants$N_space, ncol = constants$K),
        
        r = 10.0,
        
        sigma_beta0s1 = 0.1,
        sigma_beta0s2 = 0.1,
        sigma_beta0s3 = 0.1,
        
        sigma_beta1s1 = 0.1,
        sigma_beta1s2 = 0.1,
        sigma_beta1s3 = 0.1,
        
        sigma_bs1 = 0.1,
        sigma_bs2 = 0.1,
        sigma_bs3 = 0.1,
        
        sigma_deltas1 = 0.1,
        sigma_deltas2 = 0.1,
        sigma_deltas3 = 0.1
        
        )
    
    data  <- list(y = incidence$cases, n = incidence$py)
    
    model <- nimbleModel(code = code, constants = constants, inits = inits, data = data, calculate = FALSE)
    print("----- MODEL BUILT -----")
    
    Cmodel <- compileNimble(model)
    print("----- MODEL COMPILED -----")
    
    sigmas <- c("sigma_beta0s1", "sigma_beta0s2", "sigma_beta0s3", 
                "sigma_beta1s1", "sigma_beta1s2", "sigma_beta1s3",
                "sigma_bs1", "sigma_bs2", "sigma_bs3",
                "sigma_deltas1", "sigma_deltas2", "sigma_deltas3")
    
    monitors <- model$getParents(model$getNodeNames(dataOnly = TRUE), stochOnly = TRUE)
    
    monitors <- c(monitors, sigmas, "beta0", "beta1", "beta0s1", "beta0s2", "beta1s1", "beta1s2", "b0", "bs1", "bs2",
                      "deltas1", "deltas2", "delta0")
    
    mcmcConf <- configureMCMC(
        model     = Cmodel,
        monitors  = c("lograte"),
        monitors2 = monitors,
        thin      = 100,
        thin2     = 100,
        print     = FALSE,
        enableWAIC = TRUE
    )
    
    mcmcConf$removeSamplers(sigmas)
    for (s in sigmas) {
        mcmcConf$addSampler(target = s, type = "RW", control = list(log = TRUE))
    }
    
    Rmcmc <- buildMCMC(mcmcConf)
    Cmcmc <- compileNimble(Rmcmc)
    
    mcmc.out <- runMCMC(
        Cmcmc,
        niter       = n_iter,
        nburnin     = n_burn,
        progressBar = TRUE,
        setSeed     = seed,
        WAIC        = TRUE
    )
    
    return(mcmc.out)
}