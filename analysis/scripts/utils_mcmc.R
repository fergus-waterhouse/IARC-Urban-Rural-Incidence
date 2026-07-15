library(coda)
library(dplyr)

# Formats nested MCMC lists into a flat matrix
load_mcmc_matrix <- function(modfile, sample_target = "samples2") {
    mcmc_raw <- readRDS(modfile)
    as.matrix(mcmc.list(lapply(mcmc_raw, function(x) as.mcmc(x[[sample_target]]))))
}

# Extracts specific parameters as an unflattened mcmc.list
load_mcmc_list <- function(modfile, sample_target = "samples2", param_exact = NULL, param_prefix = NULL) {
    mcmc_raw <- readRDS(modfile)
    
    chain_list <- lapply(mcmc_raw, function(x) {
        mat <- x[[sample_target]]
        if (!is.null(param_exact)) {
            cols <- intersect(param_exact, colnames(mat))
        } else if (!is.null(param_prefix)) {
            cols <- colnames(mat)[startsWith(colnames(mat), param_prefix)]
        } else {
            cols <- colnames(mat)
        }
        
        if (length(cols) == 0) return(NULL)
        coda::as.mcmc(mat[, cols, drop = FALSE])
    })
    
    if (is.null(chain_list[[1]])) return(NULL)
    return(coda::mcmc.list(chain_list))
}

# Calculates Gelman-Rubin R-hat
safe_rhat <- function(mcmc_list_obj) {
    n_params <- ncol(mcmc_list_obj[[1]])
    param_names <- colnames(mcmc_list_obj[[1]])
    
    rhats <- sapply(1:n_params, function(i) {
        tryCatch({
            univ_mcmc <- coda::mcmc.list(lapply(mcmc_list_obj, function(chain) chain[, i, drop = FALSE]))
            coda::gelman.diag(univ_mcmc, autoburnin = FALSE)$psrf[, 1]
        }, error = function(e) NA)
    })
    
    setNames(rhats, param_names)
}

# Extraction function for indexed parameters (e.g., deltas3[1], deltas3[2])
extract_param_summary <- function(mcmc_mat, param_name) {
    search_str <- paste0(param_name, "[")
    cols <- colnames(mcmc_mat)[startsWith(colnames(mcmc_mat), search_str)]
    
    if (length(cols) == 0) return(tibble::tibble())
    
    clean_str <- sub(paste0("^.*", param_name), "", cols)
    indices <- as.integer(gsub("\\D", "", clean_str))
    
    sub_mat <- mcmc_mat[, cols, drop = FALSE]
    
    tibble::tibble(
        id     = indices,
        median = apply(sub_mat, 2, median),
        lo     = apply(sub_mat, 2, quantile, probs = 0.025),
        hi     = apply(sub_mat, 2, quantile, probs = 0.975)
    )
}

# Extracts unindexed, scalar parameters (like global delta0, beta0)
extract_scalar_summary <- function(mcmc_mat, param_name) {
    if (!param_name %in% colnames(mcmc_mat)) return(tibble::tibble())
    
    vec <- mcmc_mat[, param_name]
    tibble::tibble(
        median = median(vec),
        lo     = quantile(vec, 0.025),
        hi     = quantile(vec, 0.975)
    )
}

# Safely extracts spline matrices (e.g., bs2) for all hierarchical IDs
extract_spline_matrix <- function(mcmc_mat, prefix, expected_rows, K = 4) {
    b_mat <- matrix(0, nrow = expected_rows, ncol = K)
    for (i in 1:expected_rows) {
        for (k in 1:K) {
            col_name <- paste0(prefix, "[", i, ", ", k, "]")
            if (col_name %in% colnames(mcmc_mat)) {
                b_mat[i, k] <- median(mcmc_mat[, col_name])
            }
        }
    }
    return(b_mat)
}

# Extracts the literal `lograte` array and cb বিকালinds it safely to the dataset
extract_and_bind_logrates <- function(inc_data, modfile) {
    cat("Extracting literal lograte nodes from $samples...\n")
    
    mcmc_raw <- readRDS(modfile)
    combined_samples <- as.matrix(mcmc.list(lapply(mcmc_raw, function(x) as.mcmc(x$samples))))
    
    n_rows <- nrow(inc_data)
    lograte_cols <- paste0("lograte[", 1:n_rows, "]")
    
    if (!all(lograte_cols %in% colnames(combined_samples))) {
        stop("Mismatch! The rows in inc_data do not match the number of lograte nodes in the model. Check your data filtering.")
    }
    
    # Exponentiate to rate per 100k
    rate_mat <- exp(combined_samples[, lograte_cols, drop = FALSE]) * 100000
    quants <- apply(rate_mat, 2, quantile, probs = c(0.025, 0.5, 0.975))
    
    inc_data %>% mutate(
        pred_rate_100k_lo = quants[1, ],
        pred_rate_100k    = quants[2, ],
        pred_rate_100k_hi = quants[3, ]
    )
}