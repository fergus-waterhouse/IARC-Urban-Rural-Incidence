# 1. Initialize environment
source("analysis/scripts/config.R")
source("analysis/scripts/utils_mcmc.R")
source("analysis/scripts/utils_splines.R")
source("analysis/scripts/utils_plots.R")

# --- PATHS & CONFIG ---
DATA_FILE  <- "data/processed/ci5_complete.csv"
MOD_FILE   <- "output/colorectal0.05/mcmc_out.rds"
OUT_DIR    <- "figures/colorectal0.05/"

# --- MODEL TRAINING STATE ---
TRAIN_CANCER <- "Colorectal"
TRAIN_SEX    <- 0      # 0=Both, 1=M, 2=F



# Create unified directory structure
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(OUT_DIR, "preds/country/"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(OUT_DIR, "post/"), recursive = TRUE, showWarnings = FALSE) # New Posteriors Folder

# Determine which sexes to predict for the validation curves
predict_sexes <- if (TRAIN_SEX == 0) c(1, 2) else TRAIN_SEX


# =========================================================================
# 2. LOAD DATA & RECREATE ID MAPPINGS EXACTLY AS TRAINED
# =========================================================================
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

# 1. Bind exact logrates to dataset (Applies to all data in 'inc')
inc_with_preds <- extract_and_bind_logrates(inc, MOD_FILE)

# -------------------------------------------------------------------------
# GOAL IA: GENERATE POSTERIOR AGE TRENDS
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
# GOAL IB: GENERATE SEX-SPECIFIC VALIDATION PLOTS (AGAINST OBSERVED DATA)
# -------------------------------------------------------------------------
for (sx in predict_sexes) {
    
    sex_label <- ifelse(sx == 1, "Male", "Female")
    cat(sprintf("\nPlotting Country Registry Validation for %s...\n", sex_label))
    
    # Generate curves with the sex effect active
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


# =========================================================================
# GOAL IC: MCMC TRACE PLOTS WITH R-HAT
# =========================================================================
cat("\n--- Generating MCMC Trace Plots ---\n")

# 1. Trace for gamma_sex (if present)
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
    
    # FIX: Use regex to extract only the integer inside the brackets! (Ignores the '3' in 'deltas3')
    d3_idx <- as.integer(sub(".*\\[(\\d+)\\].*", "\\1", d3_cols))
    
    # Map raw indices to clean Country Names
    d3_labels <- map_ctry$country[match(d3_idx, map_ctry$hier3.id)]
    names(d3_labels) <- d3_cols
    
    p_deltas3 <- plot_mcmc_traces(deltas3_list, param_labels = d3_labels, ncol = 5)
    
    # Dynamically scale height based on the number of countries
    d3_height <- max(6, ceiling(length(d3_cols) / 5) * 2.5) 
    ggsave(paste0(OUT_DIR, "post/Trace_deltas3.png"), p_deltas3, width = 16, height = d3_height, dpi = 300)
}

# 3. Trace for 32 Random Logrates
set.seed(123) # For reproducibility
sample_indices <- sample(1:nrow(inc), min(32, nrow(inc)))
lr_exact <- paste0("lograte[", sample_indices, "]")

lr_list <- load_mcmc_list(MOD_FILE, "samples", param_exact = lr_exact)
if (!is.null(lr_list)) {
    # Extract metadata for the sampled rows to create an informative label
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

# =========================================================================
# GOAL IIA: URBAN-RURAL IRR
# =========================================================================
cat("\n--- Generating Independent Result Plots ---\n")

reg_counts <- inc %>% group_by(country) %>% summarise(n_registries = n_distinct(registry)) %>% ungroup()
cont_counts <- inc %>% group_by(continent) %>% summarise(n_registries = n_distinct(registry)) %>% ungroup()

df_global <- extract_scalar_summary(mcmc_mat, "delta0")
df_cont   <- extract_param_summary(mcmc_mat, "deltas1") %>% inner_join(map_cont, by = c("id" = "hier1.id")) %>% left_join(cont_counts, by = "continent")
df_reg    <- extract_param_summary(mcmc_mat, "deltas2") %>% inner_join(map_reg,  by = c("id" = "hier2.id"))
df_ctry   <- extract_param_summary(mcmc_mat, "deltas3") %>% inner_join(map_ctry, by = c("id" = "hier3.id")) %>% left_join(reg_counts, by = "country") %>% filter(country != "Iran (Islamic Republic of)")

# Replace the continental values with the regional values for those regions which are also continents. e.g. N America
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

# =========================================================================
# GOAL IIB: URBAN-RURAL IRR WITH CRI95% (4-Panel Linear Plot)
# =========================================================================
cat("\nGenerating 4-Panel Linear Validation Plot...\n")

# Calculate region-specific n_registries since it wasn't tracked in Goal IIA
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
    total_registries = total_registries,
    cont_order       = CONT_ORDER,
    region_order     = REGION_ORDER
)

# Use dynamic sex labeling if saving by individual sexes
sex_tag <- if (TRAIN_SEX == 0) "Both" else ifelse(TRAIN_SEX == 1, "M", "F")
ggsave(paste0(OUT_DIR, TRAIN_CANCER, "_", sex_tag, "_4panel_linear.png"), plot = p_4panel, width = 5000, height = 1500, units = "px")

# =========================================================================
# 3. GOAL IIB: URBAN-RURAL IRR WITH CRI95%
# =========================================================================

cat("\nGenerating 4-Panel Linear Validation Plot...\n")

# Calculate region-specific n_registries since it wasn't tracked in Goal IIA
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

# Use dynamic sex labeling if saving by individual sexes
sex_tag <- if (TRAIN_SEX == 0) "Both" else ifelse(TRAIN_SEX == 1, "M", "F")
ggsave(paste0(OUT_DIR, TRAIN_CANCER, "_", sex_tag, "_4panel_linear.png"), plot = p_4panel, width = 5000, height = 1500, units = "px")


# =========================================================================
# 4. GOAL IIC: Predicted ASIR with %Urban by Region
# =========================================================================
cat("Calculating Regional ASIR Trajectories...\n")

gamma_sex_med <- ifelse("gamma_sex" %in% colnames(mcmc_mat), median(mcmc_mat[, "gamma_sex"]), 0)
beta0s2 <- extract_param_summary(mcmc_mat, "beta0s2")
beta1s2 <- extract_param_summary(mcmc_mat, "beta1s2")
deltas2 <- extract_param_summary(mcmc_mat, "deltas2")
bs2_mat <- extract_spline_matrix(mcmc_mat, "bs2", expected_rows = max(inc$hier2.id, na.rm=T))

x_grid_raw <- seq(0, 1, length.out = 100)

plot_df_list <- lapply(1:nrow(map_reg), function(i) {
    id   <- map_reg$hier2.id[i]
    if (id %in% beta0s2$id) {
        asir <- predict_asir_trajectory(
            b0 = beta0s2$median[beta0s2$id == id], b1 = beta1s2$median[beta1s2$id == id], 
            b_vec = bs2_mat[id, ], delta = deltas2$median[deltas2$id == id], 
            gamma_sex = gamma_sex_med, x_grid_raw = x_grid_raw, urb_mean = urb_mean, urb_sd = urb_sd, 
            Z_18 = spline_env$Z_18, age_weights = AGE_WEIGHTS, 
            target_sex = TRAIN_SEX # Automatically averages M/F if TRAIN_SEX == 0
        )
        return(data.frame(level = "Region", name = map_reg$region[i], continent = map_reg$continent[i], x = x_grid_raw, y = asir))
    }
    return(NULL)
})

reg_asir_df <- bind_rows(plot_df_list)
region_ranges <- inc %>% group_by(region) %>% summarise(min_urb = min(urban, na.rm=T), max_urb = max(urban, na.rm=T), .groups="drop") %>% rename(name = region)

p_reg_asir <- plot_regional_asir(
    plot_df = reg_asir_df, region_ranges = region_ranges, 
    reg_urban_2017 = REG_URBAN_2017, reg_urban_2050 = REG_URBAN_2050, 
    cont_order = CONT_ORDER, full_palette = c(SDG_PALETTE, CONT_PALETTE)
)

ggsave(paste0(OUT_DIR, TRAIN_CANCER, "_region_asir.png"), plot = p_reg_asir, width = 3500, height = 2100, units = "px")


# =========================================================================
# 5. GOAL IID: HDI vs Incidence Rate Ratio (IRR) Plot
# =========================================================================
country_irr_df <- df_ctry 

p_hdi <- plot_hdi_regression(country_irr_df = country_irr_df, hdi_vec = HDI_2017, cont_palette = CONT_PALETTE)

ggsave(paste0(OUT_DIR, TRAIN_CANCER, "_hdi_regression.png"), plot = p_hdi, width = 8, height = 6, dpi = 300)



# ===================== PART I - HIERACHY VISUALISATION ===================

cat("\n--- Generating Hierarchical Plots ---\n")

# Create the new hierarchy output folder
hier_dir <- paste0(OUT_DIR, "hierarchy/")
dir.create(hier_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Generate curves specifically for Males (target_sex = 1) ONCE to save time
drill_cont <- generate_hierarchy_curves(mcmc_mat, map_cont, "continent", "hier1.id", "beta0s1", "beta1s1", "bs1", spline_env, target_sex = 1, baseline_only = FALSE)
drill_reg  <- generate_hierarchy_curves(mcmc_mat, map_reg,  "region",    "hier2.id", "beta0s2", "beta1s2", "bs2", spline_env, target_sex = 1, baseline_only = FALSE)
drill_ctry <- generate_hierarchy_curves(mcmc_mat, map_ctry, "country",   "hier3.id", "beta0s3", "beta1s3", "bs3", spline_env, target_sex = 1, baseline_only = FALSE)

# Replace the continental values with the regional values for those regions which are also continents. e.g. N America
replacements <- drill_reg %>% 
    filter(region %in% drill_cont$continent) %>% 
    select(continent = region, agr, pred_raw_lo, pred_raw, pred_raw_hi, pred_rate_100k_lo, pred_rate_100k, pred_rate_100k_hi)
drill_cont <- drill_cont %>% 
    rows_update(replacements, by = c("continent", "agr"))

# Get list of all unique countries
all_countries <- unique(map_ctry$country)

# 2. Loop through every country
for (tgt_country in all_countries) {
    
    # Identify the parent region and continent
    tgt_info   <- map_ctry %>% filter(country == tgt_country) %>% slice(1)
    tgt_region <- tgt_info$region
    tgt_cont   <- tgt_info$continent
    
    # --- A. FILTER DATA SETS DYNAMICALLY ---
    # Continents
    cont_bg <- drill_cont %>% filter(continent != tgt_cont)
    cont_hl <- drill_cont %>% filter(continent == tgt_cont)
    
    # Regions
    curr_regs <- drill_reg %>% filter(continent == tgt_cont)
    reg_bg    <- curr_regs %>% filter(region != tgt_region)
    reg_hl    <- curr_regs %>% filter(region == tgt_region)
    
    # Countries
    curr_ctrys <- drill_ctry %>% filter(region == tgt_region)
    ctry_bg    <- curr_ctrys %>% filter(country != tgt_country)
    ctry_hl    <- curr_ctrys %>% filter(country == tgt_country)
    
    # Observed Registry Data
    sex_obs <- 1
    if (TRAIN_SEX == 2) {
        sex_obs <- 2
    }
    obs_data <- inc %>%
        filter(country == tgt_country, sex == sex_obs) %>%
        group_by(country, registry, urban, agr) %>%
        summarise(cases = sum(cases, na.rm=TRUE), py = sum(py, na.rm=TRUE), .groups = "drop") %>%
        mutate(obs_rate_100k = (cases / py) * 100000) %>%
        filter(!is.na(obs_rate_100k))
    
    # --- B. CALCULATE SHARED Y-AXIS MAXIMUM ---
    # Safely extract y-values (handles cases where a region might have no other background countries)
    all_y_vals <- c(
        if (nrow(cont_bg) > 0) cont_bg$pred_rate_100k else 0,
        if (nrow(cont_hl) > 0) cont_hl$pred_rate_100k else 0,
        if (nrow(reg_bg) > 0)  reg_bg$pred_rate_100k else 0,
        if (nrow(reg_hl) > 0)  reg_hl$pred_rate_100k else 0,
        if (nrow(ctry_bg) > 0) ctry_bg$pred_rate_100k else 0,
        if (nrow(ctry_hl) > 0) ctry_hl$pred_rate_100k else 0,
        if (nrow(obs_data) > 0) obs_data$obs_rate_100k else 0
    )
    max_y <- max(all_y_vals, na.rm = TRUE) * 1.05 
    
    # Format labels (adds a line break if the region name is very long)
    lbl_region <- tgt_region
    
    # --- C. GENERATE PLOTS ---
    p_drill_1 <- plot_highlight_spline(
        bg_curves = cont_bg, hl_curve = cont_hl, ref_curve = NULL,
        hl_color = CONT_PALETTE[tgt_cont], hl_label = tgt_cont,
        title = "Continents", group_col = "continent"
    ) + coord_cartesian(ylim = c(0, max_y))
    
    p_drill_2 <- plot_highlight_spline(
        bg_curves = reg_bg, hl_curve = reg_hl, ref_curve = cont_hl,
        hl_color = SDG_PALETTE[tgt_region], hl_label = lbl_region,
        title = paste("Regions in", tgt_cont), group_col = "region"
    ) + coord_cartesian(ylim = c(0, max_y))
    
    # Use dark region color for the country highlight
    ctry_color <- darken_hex(SDG_PALETTE[tgt_region], factor = 1.2)
    
    p_drill_3 <- plot_highlight_spline(
        bg_curves = ctry_bg, hl_curve = ctry_hl, ref_curve = reg_hl,
        hl_color = ctry_color, hl_label = tgt_country,
        title = paste("Countries in\n", tgt_region), group_col = "country"
    ) + coord_cartesian(ylim = c(0, max_y))
    
    p_drill_4 <- ggplot() + 
        geom_line(
            data = obs_data,
            aes(x = agr, y = obs_rate_100k, group = registry), 
            color = "grey80",
            linewidth = 0.5, alpha = 0.5
        ) +
        geom_point(
            data = obs_data,
            aes(x = agr, y = obs_rate_100k, group = registry), 
            color = "grey60",
            size = 1, alpha = 0.5
        ) +
        geom_line(
            data = ctry_hl,
            aes(x = agr, y = pred_rate_100k),
            color = ctry_color, linewidth = 1, linetype = "solid"
        ) +
        labs(
            title = tgt_country,
            x = "Age Group (1-18)", 
            y = "Rate per 100,000"
        ) +
        theme_minimal() +
        theme(
            legend.position = "none",
            panel.grid.minor = element_blank(), 
            plot.title = element_text(face = "bold", size = 12),
            axis.line = element_line(color = "grey30")
        ) +
        scale_x_continuous(limits = c(1, 22), breaks = seq(2, 18, 4)) +
        coord_cartesian(ylim = c(0, max_y)) 
    
    # --- D. COMBINE & SAVE ---
    library(patchwork)
    p_drill_final <- p_drill_1 + p_drill_2 + p_drill_3 + p_drill_4 + plot_layout(ncol = 4)
    
    # Clean up filename (replaces spaces/special chars with underscores)
    safe_name <- gsub("[^A-Za-z0-9]", "_", tgt_country)
    save_path <- paste0(hier_dir, TRAIN_CANCER, "_", safe_name, ".png")
    
    ggsave(save_path, plot = p_drill_final, width = 22, height = 5, dpi = 300)
}

cat("\nAnalysis Complete!\n")
