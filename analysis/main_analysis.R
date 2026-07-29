source("analysis/scripts/config.R")
source("analysis/scripts/utils_mcmc.R")
source("analysis/scripts/utils_splines.R")
source("analysis/scripts/utils_plots.R")

# --- CONFIG ---
DATA_FILE  <- "data/processed/ci5_complete.csv"
MOD_FILE   <- "output/colorectal0.05/mcmc_out.rds"
OUT_DIR    <- "figures/colorectal0.05/"

# --- DATA ---
TRAIN_CANCER <- "Colorectal"
TRAIN_SEX    <- 0      # 0=Both, 1=M, 2=F


dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(OUT_DIR, "preds/country/"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(OUT_DIR, "post/"), recursive = TRUE, showWarnings = FALSE)

predict_sexes <- if (TRAIN_SEX == 0) c(1, 2) else TRAIN_SEX


# ---------------------- LOAD DATA -----------------------------------------
inc <- read.csv(DATA_FILE)

urb_vals <- inc %>% select(registry, urban) %>% unique() %>% pull(urban)
urb_mean <- mean(urb_vals, na.rm = TRUE)
urb_sd   <- sd(urb_vals, na.rm = TRUE)

inc <- inc %>% filter(cancer_lab == TRAIN_CANCER)
if (TRAIN_SEX != 0) inc <- inc %>% filter(sex == TRAIN_SEX)

inc <- inc %>% mutate(
    hier1.id = as.integer(as.factor(continent)),
    hier2.id = as.integer(as.factor(region)),
    hier3.id = as.integer(as.factor(country))
)

map_cont <- inc %>% select(hier1.id, continent) %>% distinct()
map_reg  <- inc %>% select(hier2.id, region, continent) %>% distinct()
map_ctry <- inc %>% select(hier3.id, country, region, continent) %>% distinct()


# ===================== PART I - MODEL VALIDATION =========================

cat("\nLoading MCMC parameters...\n")
mcmc_mat <- load_mcmc_matrix(MOD_FILE, sample_target = "samples2")
spline_env <- get_high_res_splines(K = 4)

inc_with_preds <- extract_and_bind_logrates(inc, MOD_FILE)

# -------------------------------------------------------------------------
# IA: GENERATE POSTERIOR AGE TRENDS
# -------------------------------------------------------------------------
cat("\n--- Generating Posterior Age Trends ---\n")

glob_post <- generate_global_curve(mcmc_mat, spline_env, target_sex = 1, baseline_only = TRUE)
cont_post <- generate_hierarchy_curves(mcmc_mat, map_cont, "continent", "hier1.id", "beta0s1", "beta1s1", "bs1", spline_env, target_sex = 1, baseline_only = TRUE)
reg_post  <- generate_hierarchy_curves(mcmc_mat, map_reg, "region", "hier2.id", "beta0s2", "beta1s2", "bs2", spline_env, target_sex = 1, baseline_only = TRUE)
ctry_post <- generate_hierarchy_curves(mcmc_mat, map_ctry, "country", "hier3.id", "beta0s3", "beta1s3", "bs3", spline_env, target_sex = 1, baseline_only = TRUE)

p_post_glob <- plot_posterior_age_trend(glob_post, "level", "Global", TRAIN_CANCER)
p_post_cont <- plot_posterior_age_trend(cont_post, "continent", "Continent", TRAIN_CANCER)
p_post_reg  <- plot_posterior_age_trend(reg_post, "region", "Region", TRAIN_CANCER)
p_post_ctry <- plot_posterior_age_trend(ctry_post, "country", "Country", TRAIN_CANCER)

ggsave(paste0(OUT_DIR, "post/Global_posterior.png"), p_post_glob, width = 6, height = 5, dpi = 300)
ggsave(paste0(OUT_DIR, "post/Continent_posterior.png"), p_post_cont, width = 10, height = 6, dpi = 300)
ggsave(paste0(OUT_DIR, "post/Region_posterior.png"), p_post_reg, width = 12, height = 8, dpi = 300)
ggsave(paste0(OUT_DIR, "post/Country_posterior.png"), p_post_ctry, width = 16, height = 32, dpi = 300)


# -------------------------------------------------------------------------
# IB:  SEX-SPECIFIC PREDICTED VS OBSERVED ASIR PLOTS
# -------------------------------------------------------------------------
for (sx in predict_sexes) {
    
    sex_label <- ifelse(sx == 1, "Male", "Female")
    cat(sprintf("\nPlotting Country Registry Validation for %s...\n", sex_label))
    
    ctry_curves <- generate_hierarchy_curves(mcmc_mat, map_ctry, "country", "hier3.id", "beta0s3", "beta1s3", "bs3", spline_env, target_sex = sx, baseline_only = FALSE)
    
    obs_df_filtered <- inc %>%
        filter(sex == sx) %>%
        group_by(country, registry, urban, agr) %>%
        summarise(cases = sum(cases, na.rm=TRUE), py = sum(py, na.rm=TRUE), .groups = "drop") %>%
        mutate(obs_rate_100k = (cases / py) * 100000) %>%
        filter(!is.na(obs_rate_100k))
    
    obs_list  <- split(obs_df_filtered, obs_df_filtered$country)
    pred_list <- split(ctry_curves, ctry_curves$country)
    
    for (ctry in unique(ctry_curves$country)) {
        obs_sub <- obs_list[[ctry]]
        pred_df <- pred_list[[ctry]]
        
        if (is.null(obs_sub) || nrow(obs_sub) == 0) next 
        
        p_obs_pred <- plot_obs_vs_pred(obs_sub, pred_df, ctry, TRAIN_CANCER)
        save_path <- paste0(OUT_DIR, "preds/country/", ctry, "_", sex_label, "_obs_vs_pred.png")
        ggsave(save_path, p_obs_pred, width = 9, height = 6, dpi = 300)   
    }
}


# -------------------------------------------------------------------------
# IC: MCMC TRACE PLOTS WITH R-HAT
# -------------------------------------------------------------------------
cat("\n--- Generating MCMC Trace Plots ---\n")

# 1. Trace for gamma_sex 
gamma_list <- load_mcmc_list(MOD_FILE, "samples2", param_exact = "gamma_sex")
if (!is.null(gamma_list)) {
    p_gamma <- plot_mcmc_traces(gamma_list, param_labels = c("gamma_sex" = "Gamma Sex Effect"), ncol = 1)
    ggsave(paste0(OUT_DIR, "post/Trace_gamma_sex.png"), p_gamma, width = 6, height = 4.5, dpi = 300)
} else {
    cat("Note: 'gamma_sex' not found in model samples. Skipping trace plot.\n")
}

# 2. Trace for deltas3 (Country Level IRR)
deltas3_list <- load_mcmc_list(MOD_FILE, "samples2", param_prefix = "deltas3[")
if (!is.null(deltas3_list)) {
    d3_cols <- colnames(deltas3_list[[1]])
    
    d3_idx <- as.integer(sub(".*\\[(\\d+)\\].*", "\\1", d3_cols))

    d3_labels <- map_ctry$country[match(d3_idx, map_ctry$hier3.id)]
    names(d3_labels) <- d3_cols
    
    p_deltas3 <- plot_mcmc_traces(deltas3_list, param_labels = d3_labels, ncol = 5)
    
    d3_height <- max(6, ceiling(length(d3_cols) / 5) * 2.5) 
    ggsave(paste0(OUT_DIR, "post/Trace_deltas3.png"), p_deltas3, width = 16, height = d3_height, dpi = 300)
}

# 3. Trace for 32 Random Logrates
set.seed(123) 
sample_indices <- sample(1:nrow(inc), min(32, nrow(inc)))
lr_exact <- paste0("lograte[", sample_indices, "]")

lr_list <- load_mcmc_list(MOD_FILE, "samples", param_exact = lr_exact)
if (!is.null(lr_list)) {
    lr_sub_data <- inc[sample_indices, ]
    lr_labels <- paste0(
        lr_sub_data$country, " | ", lr_sub_data$registry, 
        "\nAgeGrp: ", lr_sub_data$agr
    )
    names(lr_labels) <- lr_exact
    
    p_lograte <- plot_mcmc_traces(lr_list, param_labels = lr_labels, ncol = 4)
    ggsave(paste0(OUT_DIR, "post/Trace_logrates_sample.png"), p_lograte, width = 14, height = 16, dpi = 300)
}



# ===================== PART II - RESULTS =========================

# -------------------------------------------------------------------------
# IIA: URBAN-RURAL IRR TRACK PLOTS
# -------------------------------------------------------------------------
cat("\n--- Generating Independent Result Plots ---\n")

reg_counts <- inc %>% group_by(country) %>% summarise(n_registries = n_distinct(registry)) %>% ungroup()
cont_counts <- inc %>% group_by(continent) %>% summarise(n_registries = n_distinct(registry)) %>% ungroup()

df_global <- extract_scalar_summary(mcmc_mat, "delta0")
df_cont   <- extract_param_summary(mcmc_mat, "deltas1") %>% inner_join(map_cont, by = c("id" = "hier1.id")) %>% left_join(cont_counts, by = "continent")
df_reg    <- extract_param_summary(mcmc_mat, "deltas2") %>% inner_join(map_reg,  by = c("id" = "hier2.id"))
df_ctry   <- extract_param_summary(mcmc_mat, "deltas3") %>% inner_join(map_ctry, by = c("id" = "hier3.id")) %>% left_join(reg_counts, by = "country") %>% filter(country != "Iran (Islamic Republic of)")

# Replace the continental values with the regional values for those regions which are also continents.
replacements <- df_reg %>% 
    filter(region %in% df_cont$continent) %>% 
    select(continent = region, median, lo, hi)

df_cont <- df_cont %>% 
    rows_update(replacements, by = "continent")

p_unified <- plot_unified_tracks(
    df_global    = df_global, df_cont = df_cont, df_reg = df_reg, df_ctry = df_ctry,
    region_order = REGION_ORDER, sdg_palette = SDG_PALETTE, cont_palette = CONT_PALETTE,
    top_breaks = c(0.8, 0.9, 1, 1.2, 1.4, 1.6, 1.8), top_factor = 4
)

ggsave(paste0(OUT_DIR, TRAIN_CANCER, "_final_unified_tracks.svg"), plot = p_unified, width = 3960, height = 1620, units = "px")


# -------------------------------------------------------------------------
# IIB: URBAN-RURAL IRR WITH CRI95%
# -------------------------------------------------------------------------

cat("\nGenerating 4-Panel Linear Validation Plot...\n")

reg_only_counts <- inc %>% group_by(region) %>% summarise(n_registries = n_distinct(registry)) %>% ungroup()
df_reg_4panel   <- df_reg %>% left_join(reg_only_counts, by = "region")

total_registries <- n_distinct(inc$registry)

p_4panel <- plot_4panel_linear(
    df_global        = df_global,
    df_cont          = df_cont,
    df_reg           = df_reg_4panel,
    df_ctry          = df_ctry,
    sdg_palette      = SDG_PALETTE,
    cont_palette     = CONT_PALETTE,
    total_registries = total_registries
)

sex_tag <- if (TRAIN_SEX == 0) "Both" else ifelse(TRAIN_SEX == 1, "M", "F")
ggsave(paste0(OUT_DIR, TRAIN_CANCER, "_", sex_tag, "_4panel_linear.png"), plot = p_4panel, width = 5000, height = 1500, units = "px")


# -------------------------------------------------------------------------
# IIC: IRR VS HDI Plot
# -------------------------------------------------------------------------
country_irr_df <- df_ctry 

p_hdi <- plot_hdi_regression(country_irr_df = country_irr_df, hdi_vec = HDI_2017, cont_palette = CONT_PALETTE)

ggsave(paste0(OUT_DIR, TRAIN_CANCER, "_hdi_regression.png"), plot = p_hdi, width = 8, height = 6, dpi = 300)

