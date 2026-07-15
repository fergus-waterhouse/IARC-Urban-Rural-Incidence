library(ggplot2)
library(ggrepel)
library(patchwork)
library(tidyr)

# ------ HELPER FUNCTIONS ------
# Helper function to dynamically darken hex colors for borders
darken_hex <- function(color, factor = 1.4) {
    if(is.na(color)) return(NA)
    rgb_vals <- col2rgb(color) / factor
    rgb(t(rgb_vals), maxColorValue = 255)
}


# ------ PLOTS ------
# Main figure
plot_unified_tracks <- function(df_global, df_cont, df_reg, df_ctry, 
                                region_order, sdg_palette, cont_palette, 
                                label_thresh = c(0.95, 1.075), top_breaks = c(0.7, 0.8, 0.9, 1, 1.2, 1.4, 1.6, 1.8, 2.0), top_factor = 4) {
    
    combined_palette <- c(sdg_palette, cont_palette)
    combined_palette <- combined_palette[unique(names(combined_palette))]
    
    # Pre-calculate significance
    df_cont <- df_cont %>% mutate(is_significant = (exp(lo) > 1 | exp(hi) < 1))
    df_ctry <- df_ctry %>% mutate(is_significant = (exp(lo) > 1 | exp(hi) < 1))
    
    # 1. DEFINE NUMERIC Y-AXIS MAPPING
    n_reg <- length(region_order)
    y_gap <- 2.8 
    
    y_map <- data.frame(
        track_id = c(rev(region_order), "Global"),
        y_pos    = c(1:n_reg, n_reg + y_gap) 
    )
    
    # 2. PREPARE TRACK DATA 
    reg_track <- df_reg %>%
        filter(!is.na(region)) %>%
        mutate(track_id = as.character(region)) %>%
        left_join(y_map, by = "track_id") %>%
        rowwise() %>%
        mutate(
            reg_color = sdg_palette[track_id],
            reg_dark  = darken_hex(reg_color, 1.4),
            med_size  = 9
        ) %>% ungroup()
    
    glob_track <- df_global %>%
        mutate(track_id = "Global") %>%
        left_join(y_map, by = "track_id") %>%
        mutate(
            median = log((exp(median) - 1) * top_factor + 1), # Custom scaling from your code
            reg_color = "black",
            reg_dark  = "black",
            med_size  = 11
        )
    
    combined_tracks <- bind_rows(reg_track, glob_track)
    
    # 3. PREPARE POINT DATA 
    cntry_pts <- df_ctry %>%
        filter(!is.na(region), n_registries > 1) %>%
        mutate(track_id = as.character(region)) %>%
        left_join(y_map, by = "track_id") %>%
        rowwise() %>%
        mutate(
            point_fill  = ifelse(is_significant, sdg_palette[track_id], "white"),
            point_size  = ifelse(is_significant, 3.5, 2.5),
            point_color = darken_hex(sdg_palette[track_id], 1.4),
            label_name  = as.character(country)
        ) %>% ungroup()
    
    cont_pts <- df_cont %>%
        filter(continent != "Global") %>%
        mutate(track_id = "Global") %>%
        left_join(y_map, by = "track_id") %>%
        rowwise() %>%
        mutate(
            median      = log((exp(median) - 1) * top_factor + 1), # Custom scaling
            base_color  = combined_palette[as.character(continent)],
            point_fill  = ifelse(is_significant, base_color, "white"),
            point_size  = ifelse(is_significant, 3.5, 2.5),
            point_color = darken_hex(base_color, 1.4),
            label_name  = as.character(continent)
        ) %>% ungroup()
    
    combined_pts <- bind_rows(cntry_pts, cont_pts)
    
    # 4. CALCULATE SHARED X-LIMITS
    all_vals <- c(exp(combined_pts$median), 1.0)
    x_min_all <- min(all_vals, na.rm = TRUE) * 0.8
    x_max_all <- max(all_vals, na.rm = TRUE) * 1.2
    
    # 5. FUNNEL BOUNDARIES AND LABELS
    top_labels <- (top_breaks - 1) / top_factor + 1
    
    y_top_reg  <- n_reg                                      
    y_glob     <- y_map$y_pos[y_map$track_id == "Global"]    
    y_mask_min <- y_top_reg + 0.1                            
    y_mask_max <- y_glob - 0.2                               
    y_text           <- y_glob - 0.5                         
    y_funnel_top     <- y_glob - 0.8                         
    y_funnel_top_mid <- y_glob - 1.1                         
    y_funnel_bot_mid <- y_top_reg + 0.4                      
    y_funnel_bot     <- y_top_reg + 0.1                      
    
    f_left <- data.frame(
        x = c(min(top_breaks), min(top_breaks), min(top_labels), min(top_labels)),
        y = c(y_funnel_top, y_funnel_top_mid, y_funnel_bot_mid, y_funnel_bot)
    )
    f_right <- data.frame(
        x = c(max(top_breaks), max(top_breaks), max(top_labels), max(top_labels)),
        y = c(y_funnel_top, y_funnel_top_mid, y_funnel_bot_mid, y_funnel_bot)
    )
    
    # 6. BUILD THE PLOT
    ggplot() +
        geom_segment(data = combined_tracks, aes(y = y_pos, yend = y_pos, x = x_min_all, xend = x_max_all), color = "grey15", linewidth = 0.5) +
        geom_vline(xintercept = 1, linetype = "solid", color = "grey25", linewidth = 0.8) +
        annotate("rect", xmin = 0.01, xmax = 100, ymin = y_mask_min, ymax = y_mask_max, fill = "white", color = NA) +
        geom_path(data = f_left, aes(x = x, y = y), color = "grey50", linetype = "solid", linewidth = 0.6) +
        geom_path(data = f_right, aes(x = x, y = y), color = "grey50", linetype = "solid", linewidth = 0.6) +
        geom_text(data = data.frame(x = top_breaks, y = rep(y_text, length(top_breaks)), label = top_labels), aes(x = x, y = y, label = label), size = 10 / .pt, color = "grey30", fontface = "plain") +
        geom_point(data = combined_tracks, aes(y = y_pos, x = exp(median), color = reg_dark, size = med_size), shape = 124) +
        ggrepel::geom_text_repel(
            data = cntry_pts %>% filter(exp(median) > label_thresh[2] | exp(median) < label_thresh[1]),
            aes(y = y_pos, x = exp(median), label = label_name),
            size = 3.1, nudge_y = 0.40, segment.size = 0.5, segment.color = "grey10",
            min.segment.length = 0, force = 0.5, force_pull = 1, seed = 2
        ) +
        geom_point(data = combined_pts, aes(y = y_pos, x = exp(median), fill = point_fill, color = point_color, size = point_size), shape = 23, stroke = 1.2) +
        scale_fill_identity() + 
        scale_color_identity() + 
        scale_size_identity() +
        scale_x_log10(expand = c(0, 0), breaks = top_breaks) +
        scale_y_continuous(breaks = y_map$y_pos, labels = y_map$track_id, expand = expansion(add = c(0.5, 0.8))) + 
        coord_cartesian(xlim = c(x_min_all, x_max_all), clip = "off") + 
        labs(x = "IRR", y = "") +
        theme_minimal() +
        theme(
            panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(),
            axis.text.y = element_text(face = "bold", size = 11, color = "black"), axis.text.x = element_text(face = "plain", size = 10, color = "grey30"),
            axis.ticks.x = element_line(color = "grey90"), plot.title = element_text(face = "bold", margin = margin(b = 10))
        )
}


# Generates the 4-Panel Linear IRR plot across Global, Continent, Region, and Country levels
plot_4panel_linear <- function(df_global, df_cont, df_reg, df_ctry, 
                               sdg_palette, cont_palette, total_registries,
                               cont_order = NULL, region_order = NULL) {
    
    combined_palette <- c(sdg_palette, cont_palette)
    
    # Fallback to alphabetical if no explicit geographic orders are passed
    if (is.null(cont_order)) cont_order <- sort(unique(df_cont$continent))
    if (is.null(region_order)) region_order <- sort(unique(df_reg$region))
    
    # --- 1. PREP DATA ---
    df_global <- df_global %>% 
        mutate(continent = "Global",
               n_registries = total_registries, 
               is_significant = (exp(lo) > 1 | exp(hi) < 1), 
               point_alpha = if_else(is_significant, 1, 0.6), 
               base_color = "black", dark_color = "black")
    
    df_cont <- df_cont %>% 
        mutate(n_registries = replace_na(n_registries, 1), 
               is_significant = (exp(lo) > 1 | exp(hi) < 1), 
               point_alpha = if_else(is_significant, 1, 0.6)) %>% 
        rowwise() %>% 
        mutate(base_color = ifelse(continent == "Global", "black", combined_palette[as.character(continent)]), 
               dark_color = darken_hex(base_color, factor = 1.4)) %>% 
        ungroup()
    
    df_reg <- df_reg %>% 
        mutate(n_registries = replace_na(n_registries, 1), 
               is_significant = (exp(lo) > 1 | exp(hi) < 1), 
               point_alpha = if_else(is_significant, 1, 0.6)) %>% 
        rowwise() %>% 
        mutate(base_color = sdg_palette[as.character(region)], 
               dark_color = darken_hex(base_color, factor = 1.4)) %>% 
        ungroup()
    
    df_ctry <- df_ctry %>% 
        mutate(n_registries = replace_na(n_registries, 1), 
               is_significant = (exp(lo) > 1 | exp(hi) < 1), 
               point_alpha = if_else(is_significant, 1, 0.6)) %>% 
        rowwise() %>% 
        mutate(base_color = sdg_palette[as.character(region)], 
               dark_color = darken_hex(base_color, factor = 1.4)) %>% 
        ungroup()
    
    # --- 2. SHARED AESTHETICS ---
    shared_y_scale <- scale_y_log10(breaks = c(0.5, 1.0, 1, 2))
    shared_coord   <- coord_cartesian(ylim = c(0.5, 2))
    shared_size_scale <- scale_size_continuous(range = c(1.5, 4.5), limits = c(1, total_registries), name = "N Registries", transform = "log10")
    shared_theme <- theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9), 
        axis.ticks.x = element_line(colour = "black", linewidth = 1), 
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5), 
        legend.position = "none", 
        plot.title = element_text(hjust = 0, size = 12, face = "bold", margin = margin(b = 1))
    )
    
    # --- 3. BUILD PANELS (WITH STRICT HIERARCHICAL SORTING) ---
    # Panel 1: Global
    p1 <- ggplot(df_global, aes(x = continent, y = exp(median))) + 
        geom_hline(yintercept = 1, linetype = "solid", color = "black", linewidth = 0.8) + 
        geom_errorbar(aes(ymin = exp(lo), ymax = exp(hi), color = dark_color), width = 0.25, linewidth = 0.8) + 
        geom_point(aes(fill = base_color), size = 2.5, shape = 21, stroke = 1, color = "black") + 
        scale_color_identity() + scale_fill_identity() + scale_alpha_identity() + 
        shared_size_scale + shared_y_scale + shared_coord + 
        ylab("IRR") + xlab("") + theme_minimal() + shared_theme + labs(title = "")
    
    # Panel 2: Continents
    conts_filtered <- df_cont %>% 
        filter(continent != "Global") %>% 
        mutate(continent = factor(continent, levels = intersect(cont_order, unique(continent)))) %>% 
        arrange(continent) %>%
        mutate(continent = factor(continent, levels = unique(continent)), cont_idx = as.numeric(continent))
    
    p2 <- ggplot(conts_filtered, aes(x = cont_idx, y = exp(median))) + 
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = exp(df_global$lo), ymax = exp(df_global$hi), fill = "grey60", alpha = 0.15) + 
        geom_hline(yintercept = 1, linetype = "solid", color = "black", linewidth = 0.8) + 
        geom_hline(yintercept = exp(df_global$median), linetype = "dashed", color = "grey30", linewidth = 0.8) + 
        geom_errorbar(aes(ymin = exp(lo), ymax = exp(hi), color = dark_color, alpha = point_alpha), width = 0.25, linewidth = 0.8) + 
        geom_point(aes(fill = base_color, size = n_registries, alpha = point_alpha), shape = 21, stroke = 1, color = "black") + 
        scale_x_continuous(breaks = conts_filtered$cont_idx, labels = levels(conts_filtered$continent), expand = expansion(add = c(0.5, 0.5))) + 
        scale_color_identity() + scale_fill_identity() + scale_alpha_identity() + 
        shared_size_scale + shared_y_scale + shared_coord + 
        ylab("") + xlab("") + theme_minimal() + shared_theme + 
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + labs(title = "Continent")
    
    # Panel 3: Regions (Grouped securely by Continent)
    regs_ordered <- df_reg %>% 
        mutate(
            continent = factor(continent, levels = intersect(cont_order, unique(continent))),
            region = factor(region, levels = intersect(region_order, unique(region)))
        ) %>%
        arrange(continent, region) %>%
        mutate(region = factor(region, levels = unique(region)), region_idx = as.numeric(region))
    
    cont_overlay_for_regions <- regs_ordered %>% 
        group_by(continent) %>% 
        summarise(x_start = min(region_idx) - 0.45, x_end = max(region_idx) + 0.45, .groups="drop") %>% 
        inner_join(df_cont, by = "continent") %>% 
        rowwise() %>% mutate(line_color = combined_palette[as.character(continent)]) %>% ungroup()
    
    p3 <- ggplot(regs_ordered, aes(x = region_idx, y = exp(median))) + 
        geom_rect(data = cont_overlay_for_regions, aes(xmin = x_start, xmax = x_end, ymin = exp(lo), ymax = exp(hi), fill = line_color), alpha = 0.15, inherit.aes = FALSE) + 
        geom_segment(data = cont_overlay_for_regions, aes(x = x_start, xend = x_end, y = exp(median), yend = exp(median), color = line_color), linewidth = 0.8, linetype = "solid") + 
        geom_hline(yintercept = 1, linetype = "solid", color = "black", linewidth = 0.8) + 
        geom_errorbar(aes(ymin = exp(lo), ymax = exp(hi), color = dark_color, alpha = point_alpha), width = 0.25, linewidth = 0.8) + 
        geom_point(aes(fill = base_color, size = n_registries, alpha = point_alpha), shape = 21, stroke = 1, color = "black") + 
        scale_x_continuous(breaks = regs_ordered$region_idx, labels = regs_ordered$region, expand = expansion(add = c(0.5, 0.5))) + 
        scale_color_identity() + scale_fill_identity() + scale_alpha_identity() + 
        shared_size_scale + shared_y_scale + shared_coord + 
        ylab("") + xlab("") + theme_minimal() + shared_theme + 
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + labs(title = "Region")
    
    # Panel 4: Countries (Grouped securely by Continent, then Region)
    countries_filtered <- df_ctry %>% 
        mutate(
            continent = factor(continent, levels = intersect(cont_order, unique(continent))),
            region = factor(region, levels = intersect(region_order, unique(region)))
        ) %>%
        arrange(continent, region, country) %>%
        mutate(country = factor(country, levels = unique(country)), country_idx = as.numeric(country))
    
    region_overlay_for_countries <- countries_filtered %>% 
        group_by(region) %>% 
        summarise(x_start = min(country_idx) - 0.45, x_end = max(country_idx) + 0.45, .groups="drop") %>% 
        inner_join(df_reg, by = "region") %>% 
        rowwise() %>% mutate(line_color = sdg_palette[as.character(region)]) %>% ungroup()
    
    p4 <- ggplot(countries_filtered, aes(x = country_idx, y = exp(median))) + 
        geom_rect(data = region_overlay_for_countries, aes(xmin = x_start, xmax = x_end, ymin = exp(lo), ymax = exp(hi), fill = line_color), alpha = 0.15, inherit.aes = FALSE) + 
        geom_segment(data = region_overlay_for_countries, aes(x = x_start, xend = x_end, y = exp(median), yend = exp(median), color = line_color), linewidth = 0.8, linetype = "solid") + 
        geom_hline(yintercept = 1, linetype = "solid", color = "black", linewidth = 0.8) + 
        geom_errorbar(aes(ymin = exp(lo), ymax = exp(hi), color = dark_color, alpha = point_alpha), width = 0.25, linewidth = 0.8) + 
        geom_point(aes(fill = base_color, size = n_registries, alpha = point_alpha), shape = 21, stroke = 1, color = "black") + 
        scale_x_continuous(breaks = countries_filtered$country_idx, labels = levels(countries_filtered$country), expand = expansion(add = c(0.5, 0.5))) + 
        scale_color_identity() + scale_fill_identity() + scale_alpha_identity() + 
        shared_size_scale + shared_y_scale + shared_coord + 
        ylab("") + xlab("") + theme_minimal() + shared_theme + 
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + labs(title = "Country")
    
    # Combine using Patchwork dynamically distributing space
    final_linear <- p1 + p2 + p3 + p4 + plot_layout(widths = c(1, nrow(conts_filtered), nrow(regs_ordered), nrow(countries_filtered)))
    
    return(final_linear)
}


# Generates the linear regression plot of HDI vs Country IRR
plot_hdi_regression <- function(country_irr_df, hdi_vec, cont_palette) {
    
    # Generate a darkened border palette automatically
    cont_dark_palette <- sapply(cont_palette, darken_hex, factor = 1.4)
    
    # Prep data: map HDI, check significance, and map sizes/colors
    plot_data <- country_irr_df %>%
        mutate(
            hdi = hdi_vec[as.character(country)],
            is_significant = (exp(lo) > 1 | exp(hi) < 1),
            fill_var = if_else(is_significant, as.character(continent), "Non-Significant"),
            point_size = if_else(is_significant, 3.5, 2.5)
        ) %>%
        filter(!is.na(hdi)) # Drop countries without HDI data
    
    # Build and return the ggplot object
    ggplot(plot_data, aes(x = hdi, y = exp(median))) +
        geom_hline(yintercept = 1, linetype = "solid", color = "grey60") +
        geom_point(aes(fill = fill_var, color = continent, size = point_size), 
                   shape = 23, stroke = 1.2) +
        geom_smooth(method = "lm", color = "grey30", linetype = "dashed",
                    linewidth = 1, se = FALSE) +
        scale_fill_manual(
            values = c(cont_palette, "Non-Significant" = "white"),
            name = "Point Fill (Significant)"
        ) +
        scale_color_manual(
            values = cont_dark_palette,
            name = "Continent"
        ) +
        scale_size_identity() +
        labs(
            x = "Human Development Index",
            y = "Incidence Rate Ratio (IRR)"
        ) +
        theme_minimal() +
        theme(
            panel.grid.minor = element_blank(),
            legend.position = "none",
            axis.ticks = element_line(color = "grey30", linewidth = 0.8),
            axis.line = element_line(colour = "grey30", linewidth = 0.8)
        )
}


# Plots predicted rate curve vs observed.
plot_obs_vs_pred <- function(obs_df, pred_df, target_country, cancer) {
    ggplot() + 
        geom_ribbon(
            data = pred_df,
            aes(x = agr, ymin = pred_rate_100k_lo, ymax = pred_rate_100k_hi),
            fill = "black", alpha = 0.4
        ) +
        geom_line(
            data = obs_df,
            aes(x = agr, y = obs_rate_100k, group = registry, color = urban),
            linewidth = 0.5, alpha = 0.4
        ) +
        geom_point(
            data = obs_df,
            aes(x = agr, y = obs_rate_100k, group = registry, color = urban),
            size = 1, alpha = 0.4
        ) +
        geom_line(
            data = pred_df,
            aes(x = agr, y = pred_rate_100k),
            color = "black", linewidth = 1.2
        ) +
        scale_color_continuous(palette = "PiYG") +
        labs(
            title = target_country,
            x = "Age Group (1-18)", 
            y = "Rate per 100,000"
        ) +
        theme_minimal() +
        theme(
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", size = 14),
            legend.position = "right"
        )
}



# Predicted regional ASIR by Region
plot_regional_asir <- function(plot_df, region_ranges, reg_urban_2017, reg_urban_2050, 
                               cont_order, full_palette) {
    
    # Ensure factor ordering aligns with geographic logic
    plot_df$continent <- factor(plot_df$continent, levels = cont_order)
    plot_df <- plot_df %>% arrange(continent, name)
    plot_df$name <- factor(plot_df$name, levels = unique(plot_df$name))
    region_ranges$name <- factor(region_ranges$name, levels = levels(plot_df$name))
    
    # Calculate specific projection points and hex colors natively
    region_points <- plot_df %>%
        group_by(continent, name) %>%
        summarise(
            urban_val = reg_urban_2017[as.character(name[1])], 
            asir_val = approx(x, y, xout = urban_val)$y, 
            .groups = "drop"
        ) %>%
        rowwise() %>% mutate(dark_hex = darken_hex(full_palette[as.character(name)], 1.4)) %>% ungroup()
    
    region_points_2050 <- plot_df %>%
        group_by(continent, name) %>%
        summarise(
            urban_val = reg_urban_2050[as.character(name[1])], 
            asir_val = approx(x, y, xout = urban_val)$y, 
            .groups = "drop"
        ) %>%
        rowwise() %>% mutate(dark_hex = darken_hex(full_palette[as.character(name)], 1.4)) %>% ungroup()
    
    # Build the GGPlot
    ggplot(plot_df) + 
        # Area and trend
        geom_area(aes(x = x, y = y, color = name, fill = name), linewidth = 1, alpha = 0.3) + 
        
        # Extrapolation Masks (Fades areas outside observed data)
        geom_rect(data = region_ranges, aes(xmin = -Inf, xmax = min_urb, ymin = -Inf, ymax = Inf), fill = "white", alpha = 0.4, inherit.aes = FALSE) +
        geom_rect(data = region_ranges, aes(xmin = max_urb, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "white", alpha = 0.4, inherit.aes = FALSE) +
        
        # 2017 Dropdown Line
        geom_segment(data = region_points, aes(x = urban_val, xend = urban_val, y = 0, yend = asir_val), color = "grey40", linetype = "dashed", linewidth = 0.8) +
        
        # 2050 & 2017 Points
        geom_point(data = region_points_2050, aes(x = urban_val, y = asir_val, color = I(dark_hex), fill = I(dark_hex)), shape = 21, stroke = 1.4, size = 1.5) +
        geom_point(data = region_points, aes(x = urban_val, y = asir_val, color = I(dark_hex)), shape = 21, stroke = 1.4, size = 2, fill = "white") +
        
        # Labels
        geom_text(data = region_points_2050 %>% filter(name %in% c("E Asia", "SE Asia", "S-S Africa")), aes(x = urban_val, y = asir_val, label = "2050"), color = "grey30", vjust = -1.2, size = 2.75) +
        geom_text(data = region_points, aes(x = urban_val, y = asir_val, label = "2017"), color = "grey30", vjust = -1.2, size = 2.75) +
        
        scale_color_manual(values = full_palette) + 
        scale_fill_manual(values = full_palette) + 
        facet_wrap(~ name, ncol = 6) + 
        scale_y_continuous(limits = c(0, 1.25 * max(plot_df$y))) +
        scale_x_continuous(limits = c(0, 1), breaks = c(0, 1)) +
        
        labs(y = "ASIR per 100,000", x = "Proportion Urban") +
        theme_minimal() + 
        theme(
            legend.position = "none",
            panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
            strip.text = element_text(face = "plain", size = 9),
            axis.line = element_line(color = "grey30", linewidth = 0.8),
            axis.ticks = element_line(color = "grey30", linewidth = 0.8)
        )
}


# Plots the posterior age trend (Spline + 95% CrI) faceted by the hierarchy level
plot_posterior_age_trend <- function(curve_df, level_col, title_prefix, cancer) {
    p <- ggplot(curve_df, aes(x = agr)) +
        geom_ribbon(aes(ymin = pred_rate_100k_lo / 100000, ymax = pred_rate_100k_hi / 100000), color = "black", alpha = 0.3) +
        geom_line(aes(y = pred_rate_100k / 100000), color = "black", linewidth = 1) +
        labs(
            x = "Age Group (1-18)", 
            y = ""
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(face = "bold"),
            panel.grid.minor = element_blank(),
            strip.text = element_text(face = "bold", size = 9),
            panel.border = element_rect(color = "grey80", fill = NA),
            legend.position = "none"
        )
    
    if (level_col == "country") {
        p <- p + facet_wrap(as.formula(paste0("~", level_col)), ncol = 5, scales = "free_y")
    } else if (level_col == "region") {
        p <- p + facet_wrap(as.formula(paste0("~", level_col)), ncol = 4, scales = "free_y")
    } else if (level_col == "continent") {
        p <- p + facet_wrap(as.formula(paste0("~", level_col)), ncol = 3, scales = "free_y")
    } 
    
    
    return(p)
}



# Generates trace plots for an mcmc.list, showing multiple chains and R-hat
plot_mcmc_traces <- function(mcmc_list, param_labels = NULL, ncol = 4) {
    
    # 1. Calculate R-hats
    rhats <- safe_rhat(mcmc_list)
    param_names <- colnames(mcmc_list[[1]])
    
    # 2. Flatten into long format for ggplot
    df_list <- lapply(1:length(mcmc_list), function(chain_id) {
        df <- as.data.frame(as.matrix(mcmc_list[[chain_id]]))
        df$Iteration <- 1:nrow(df)
        df$Chain <- as.factor(chain_id)
        pivot_longer(df, cols = -c(Iteration, Chain), names_to = "Parameter", values_to = "Value")
    })
    plot_df <- bind_rows(df_list)
    
    # 3. Safely map labels (Clean subtitles)
    if (!is.null(param_labels)) {
        disp_names <- param_labels[param_names]
        # Replace any NAs with the raw parameter name just in case mapping failed
        disp_names[is.na(disp_names)] <- param_names[is.na(disp_names)]
        
        # Smart fallback: If duplicates somehow exist, append the parameter name to disambiguate
        dup_idx <- duplicated(disp_names) | duplicated(disp_names, fromLast = TRUE)
        disp_names[dup_idx] <- paste0(disp_names[dup_idx], " (", param_names[dup_idx], ")")
    } else {
        disp_names <- param_names
    }
    names(disp_names) <- param_names
    
    # Map R-hat into the facet string
    facet_levels <- paste0(disp_names, "\n(R-hat: ", sprintf("%.3f", rhats[param_names]), ")")
    names(facet_levels) <- param_names
    
    # Apply factor levels (Guaranteed strictly unique now)
    plot_df$FacetLabel <- facet_levels[plot_df$Parameter]
    plot_df$FacetLabel <- factor(plot_df$FacetLabel, levels = unique(facet_levels))
    
    # 4. Build Plot
    ggplot(plot_df, aes(x = Iteration, y = Value, color = Chain)) +
        geom_line(alpha = 0.8, linewidth = 0.3) +
        facet_wrap(~ FacetLabel, ncol = ncol, scales = "free_y") +
        scale_color_brewer(palette = "Set1") +
        theme_minimal() +
        theme(
            panel.grid.minor = element_blank(),
            strip.text = element_text(size = 9, face = "bold"),
            legend.position = "bottom",
            panel.border = element_rect(color = "grey80", fill = NA)
        ) +
        labs(x = "Iteration", y = "Parameter Value", color = "Chain")
}


# Generates a drill-down spline plot with background curves, a reference curve, and a highlighted target
plot_highlight_spline <- function(bg_curves, hl_curve, ref_curve = NULL, 
                                  hl_color, hl_label, title, group_col) {
    
    p <- ggplot()
    
    # 1. Background curves (Grey)
    if (!is.null(bg_curves) && nrow(bg_curves) > 0) {
        p <- p + geom_line(
            data = bg_curves, 
            aes(x = agr, y = pred_rate_100k, group = .data[[group_col]]), 
            color = "grey80", linewidth = 0.8, alpha = 0.6
        )
    }
    
    # 2. Reference curve (Black Dashed - Parent Level)
    if (!is.null(ref_curve) && nrow(ref_curve) > 0) {
        p <- p + geom_line(
            data = ref_curve, 
            aes(x = agr, y = pred_rate_100k), 
            color = "black", linetype = "dashed", linewidth = 0.9
        )
    }
    
    # 3. Highlighted target curve (Colored + Label)
    if (!is.null(hl_curve) && nrow(hl_curve) > 0) {
        p <- p + geom_line(
            data = hl_curve, 
            aes(x = agr, y = pred_rate_100k), 
            color = hl_color, linewidth = 1.2
        )
        
        # Add label at the very end of the curve (Age Group 18)
        label_data <- hl_curve %>% filter(agr == max(agr))
        p <- p + ggrepel::geom_text_repel(
            data = label_data, 
            aes(x = agr, y = pred_rate_100k, label = hl_label),
            color = hl_color, fontface = "bold", size = 4,
            nudge_x = 1.5, direction = "y", hjust = 0, segment.color = NA
        )
    }
    
    # Ensure X-axis has room for the ggrepel label at the end
    p + scale_x_continuous(limits = c(1, 22), breaks = seq(2, 18, 4)) + 
        labs(title = title, x = "Age Group (1-18)", y = "Rate per 100,000") +
        theme_minimal() +
        theme(
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", size = 12),
            axis.line = element_line(color = "grey30")
        )
}