#!/usr/bin/env Rscript
# run_pipeline.R

suppressPackageStartupMessages({
  library(tidyverse)
})

# --- SETUP DIRECTORIES ---
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed/description", showWarnings = FALSE, recursive = TRUE)

cat("====================================================\n")
cat("      STARTING FULL DATA PIPELINE \n")
cat("====================================================\n")

# --- 1. RUN CI5 DATA PROCESSING ---
cat("\n>>> RUNNING CI5 PROCESSING (data_process.R)...\n")
system2("Rscript", args = c("scripts/data_step1.R", 
                            "data/raw/CI5XII.csv", 
                            "data/raw/urban.csv", 
                            "data/processed/ci5_step1.csv"))

# --- 2. RUN NORDCAN DATA PROCESSING ---
cat("\n>>> RUNNING NORDCAN PROCESSING (data_nordcan.R)...\n")
system2("Rscript", args = c("scripts/data_nordcan_step1.R", 
                            "data/raw/nordcan", 
                            "data/processed/nordcan_step1.csv"))

# --- 3. COMBINE AND HARMONIZE DATA ---
cat("\n>>> COMBINING DATASETS...\n")
ci5 <- read.csv("data/processed/ci5_step1.csv")
nordcan <- read.csv("data/processed/nordcan_step1.csv")

cat("\n--- BEFORE COMBINING ---\n")
cat(sprintf("CI5 Rows: %s | Registries: %s\n", format(nrow(ci5), big.mark=","), length(unique(ci5$registry))))
cat(sprintf("NORDCAN Rows: %s | Registries: %s\n", format(nrow(nordcan), big.mark=","), length(unique(nordcan$registry))))

# Harmonize NORDCAN columns to match CI5 structure
# CI5 has: period, continent, region, country, registry, n_agr, cancer_lab, sex, age, cases, py, urban
# NORDCAN has: continent, region, country, registry, years, cancer, sex, n_agr, agr, urban, cases, py
nordcan_aligned <- nordcan %>% select(all_of(names(ci5)))

# Combine
# First remove the existing country level data for nordic countries from CI5
nordic_countries <- nordcan_aligned %>% pull(country) %>% unique()
ci5_wonord <- filter(ci5, !(country %in% nordic_countries))
# The combine:
combined_data <- bind_rows(ci5_wonord, nordcan_aligned) %>%
  arrange(continent, region, country, registry, cancer_lab, sex, age)

cat("\n--- AFTER COMBINING ---\n")
cat(sprintf("Total Rows: %s\n", format(nrow(combined_data), big.mark=",")))
cat(sprintf("Total Registries: %s\n", length(unique(combined_data$registry))))
cat(sprintf("Total Cases: %s\n", format(sum(combined_data$cases, na.rm=TRUE), big.mark=",")))

# Save combined data
combined_path <- "data/processed/ci5_step2.csv"
write.csv(combined_data, combined_path, row.names = FALSE)
cat(sprintf("Combined data saved to: %s\n", combined_path))

cat("\n--- FINAL FORMATTING FOR MODEL... ---\n")
system2("Rscript", args = c("scripts/data_step2.R", 
                            "data/processed/ci5_step2.csv", 
                            "data/processed/ci5_complete.csv"))


# --- 4. RUN DESCRIPTIVE FIGURES ---
cat("\n>>> RUNNING DESCRIPTIVE FIGURES (data_description.R)...\n")
system2("Rscript", args = c("scripts/data_description.R", 
                            combined_path,
                            "data/processed/description"))

cat("\n====================================================\n")
cat("      PIPELINE COMPLETE \n")
cat("====================================================\n")