#!/usr/bin/env Rscript
"Data Coordination Pipeline: Prepare CI5 and NORDCAN data.

Usage:
    prepare_data.R <ci5_file> <ci5_urban> <nordcan_dir> <output>
    prepare_data.R (-h | --help)

Options:
    -h --help       Show this help message and exit.

Arguments:
    <ci5_file>      Path to the CI5 cases data CSV (e.g., CI5XII.csv).
    <ci5_urban>     Path to the CI5 urban data CSV (e.g., urban.csv).
    <nordcan_dir>   Path to the directory containing NORDCAN CSV files.
    <output>        Path for the final output CSV.
" -> doc

# ==============================================================================
# SETUP & DEPENDENCIES
# ==============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(docopt)
})

if (!interactive()) {
  args <- docopt(doc)
} else {
  # Mock arguments for interactive testing
  args <- list(ci5_file = "data/CI5XII.csv", ci5_urban = "data/urban.csv", 
               nordcan_dir = "data/raw/nordcan", output = "data/model_ready.csv") 
}

# ==============================================================================
# DICTIONARIES
# ==============================================================================
country_regions <- c(
  "Denmark" = "N Europe", "Finland" = "N Europe", "Iceland" = "N Europe",
  "Norway" = "N Europe", "Sweden" = "N Europe", "UK" = "N Europe",
  "Latvia" = "N Europe", "Lithuania" = "N Europe", "Estonia" = "N Europe",
  "Austria" = "W Europe", "France" = "W Europe", "Germany" = "W Europe",
  "Liechtenstein" = "W Europe", "Switzerland" = "W Europe", "The Netherlands" = "W Europe",
  "Ireland" = "W Europe", "Belgium" = "W Europe", "Belarus" = "E Europe",
  "Bulgaria" = "E Europe", "Czech Republic" = "E Europe", "Poland" = "E Europe",
  "Russian Federation" = "E Europe", "Slovakia" = "E Europe", "Ukraine" = "E Europe",
  "Croatia" = "S Europe", "Italy" = "S Europe", "Malta" = "S Europe",
  "Portugal" = "S Europe", "Slovenia" = "S Europe", "Spain" = "S Europe", "Cyprus" = "S Europe",
  "China" = "E Asia", "Japan" = "E Asia", "Republic of Korea" = "E Asia",
  "Thailand" = "SE Asia", "Singapore" = "SE Asia", "Philippines" = "SE Asia",
  "Brunei Darussalam" = "SE Asia", "Iran (Islamic Republic of)" = "C&S Asia",
  "India" = "C&S Asia", "Canada" = "N America", "USA" = "N America",
  "Argentina" = "S America", "Brazil" = "S America", "Chile" = "S America",
  "Colombia" = "S America", "Ecuador" = "S America", "Peru" = "S America",
  "Uruguay" = "S America", "Costa Rica" = "S America", "Trinidad and Tobago" = "S America",
  "Turkey" = "N Africa & W Asia", "Israel" = "N Africa & W Asia", "Kuwait" = "N Africa & W Asia",
  "Algeria" = "N Africa & W Asia", "Morocco" = "N Africa & W Asia", "Benin" = "N Africa & W Asia",
  "Kenya" = "S-S Africa", "Mauritius" = "S-S Africa", "Seychelles" = "S-S Africa",
  "South Africa" = "S-S Africa", "Uganda" = "S-S Africa", "Zimbabwe" = "S-S Africa",
  "Australia" = "Oceania", "New Zealand" = "Oceania"
)

region_continents <- c(
  "N America" = "N America", "S America" = "S America", "N Europe" = "Europe",
  "W Europe" = "Europe", "E Europe" = "Europe", "S Europe" = "Europe",
  "C&S Asia" = "Asia", "E Asia" = "Asia", "SE Asia" = "Asia",
  "N Africa & W Asia" = "Africa & W Asia", "S-S Africa" = "Africa & W Asia",
  "Oceania" = "Oceania"
)

nordic_countries <- c("Denmark", "Finland", "Iceland", "Norway", "Sweden")

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================
clean_nordcan_names <- function(data) {
  country_codes <- c("DK"="Denmark", "FI"="Finland", "IS"="Iceland", "NO"="Norway", "SE"="Sweden")
  res <- data %>% filter(!(Population %in% c("Faroe Islands", "Greenland")))
  codes <- sapply(res$Population, function(s) substr(s, 1, 2))
  res$country <- country_codes[codes]
  res$registry <- gsub("*", "", substr(res$Population, 5, 100), fixed = TRUE)
  res$Population <- NULL; res$Number <- NULL
  return(res)
}

# ==============================================================================
# [STEP 1] LOAD AND FILTER CI5 DATA
# ==============================================================================
cat("\n[1] PROCESSING CI5 DATA...\n")
ci5 <- read.csv(args$ci5_file) %>%
  rename(reglab = registry_lab, period = registry_years) %>%
  mutate(
    reglab = ifelse(reglab == "Turkey, Eski?ehir", "Turkey, Eskişehir", reglab),
    country = sub(",.*", "", reglab),
    registry = sub("^[^,]+, ", "", reglab),
    region = country_regions[country],
    continent = region_continents[region]
  ) %>%
  # Remove ethnic specific datasets and keep missing mappings intact
  filter(ethnic_group == 99) %>%
  select(-ethnic_group) %>%
  # Remove spurious/duplicate/overseas registries and Nordic countries (to be replaced by NORDCAN)
  filter(
    !(registry %in% c("Vologda Region", "Japan", "Martinique", "Guadeloupe", "Republic of Korea", "NPCR")),
    !(country %in% nordic_countries)
  ) %>%
  select(period, continent, region, country, registry, cancer_lab, sex, age, cases, py)

cat(sprintf(" -> CI5 data subset to %s rows.\n", format(nrow(ci5), big.mark=",")))

# ==============================================================================
# [STEP 2] LOAD AND PREPARE NORDCAN DATA
# ==============================================================================
cat("\n[2] PROCESSING NORDCAN DATA...\n")
agr_names <- c(paste0("X", seq(0, 85, by=5)))

# Calculate Population / Person-Years
all_pop <- list()
for (s in 1:2) {
  pop_file <- file.path(args$nordcan_dir, paste0("pop_", s, ".csv"))
  if(file.exists(pop_file)) {
    all_pop[[s]] <- clean_nordcan_names(read.csv(pop_file)) %>%
      pivot_longer(cols = all_of(agr_names), names_to = "agecat", values_to = "pop") %>%
      mutate(py = 5 * pop, age = match(agecat, agr_names)) %>%
      select(country, registry, age, py)
  }
}

# Process Cases
sites <- c("Kidney", "Rectum", "Stomach", "Breast", "Oesophagus", "Thyroid", "Colon", 
           "Lung (incl. trachea and bronchus)", "Pancreas", "Prostate", "Bladder", "Liver", 
           "Cervix uteri", "Non-Hodgkin lymphoma", "All sites")
nordcan_list <- list()

for (site in sites) {
  for (s in 1:2) {
    file_name <- file.path(args$nordcan_dir, paste0("case_", tolower(site), "_", s, ".csv"))
    if (file.exists(file_name)) {
      raw <- clean_nordcan_names(read.csv(file_name)) %>%
        pivot_longer(cols = all_of(agr_names), names_to = "agecat", values_to = "cases") %>%
        mutate(age = match(agecat, agr_names), cancer_lab = site, sex = s, period = "2013-2017") %>%
        inner_join(all_pop[[s]], by = c("country", "registry", "age")) %>%
        mutate(region = country_regions[country], continent = region_continents[region]) %>%
        select(period, continent, region, country, registry, cancer_lab, sex, age, cases, py)
      
      nordcan_list[[paste0(site, "_", s)]] <- raw
    }
  }
}
nordcan <- bind_rows(nordcan_list)
cat(sprintf(" -> NORDCAN data compiled to %s rows.\n", format(nrow(nordcan), big.mark=",")))

# ==============================================================================
# [STEP 3] COMBINE DATASETS AND HARMONIZE CANCERS
# ==============================================================================
cat("\n[3] MERGING CASES & COMBINING COLORECTAL...\n")
combined_cases <- bind_rows(ci5, nordcan) %>%
  mutate(cancer_lab = if_else(cancer_lab %in% c("Colon", "Rectum"), "Colorectal", cancer_lab)) %>%
  group_by(period, continent, region, country, registry, cancer_lab, sex, age) %>%
  summarise(
    cases = sum(cases, na.rm = TRUE),
    py = first(py), # Py is identical across sites for the same pop/sex/age
    .groups = "drop"
  )
cat(sprintf(" -> Total harmonized cases dataset: %s rows.\n", format(nrow(combined_cases), big.mark=",")))

# ==============================================================================
# [STEP 4] LOAD URBAN DATA AND MERGE
# ==============================================================================
cat("\n[4] PROCESSING URBAN DATA...\n")
# CI5 Urban
ci5_urban <- read.csv(args$ci5_urban) %>%
  rename(urban = Urban) %>%
  mutate(
    reglab = trimws(reglab, "right"),
    country = sub(",.*", "", reglab),
    registry = sub("^[^,]+, ", "", reglab),
    urban = urban / 100
  ) %>% select(country, registry, urban)

# NORDCAN Urban
nordcan_urban_raw <- read.csv(file.path(args$nordcan_dir, "nordcan_urban.csv"))
if("Population" %in% names(nordcan_urban_raw)) {
  nordcan_urban <- clean_nordcan_names(nordcan_urban_raw) %>% rename(urban = Urban)
} else {
  nordcan_urban <- nordcan_urban_raw %>% rename(registry = reglab)
}
# Standardize NORDCAN urban ratio if it exists as percentage
if(max(nordcan_urban$urban, na.rm=TRUE) > 1) nordcan_urban$urban <- nordcan_urban$urban / 100

# Merge urban maps together, then join to cases
all_urban <- bind_rows(ci5_urban, nordcan_urban %>% select(country, registry, urban)) %>% distinct()

final_df <- combined_cases %>%
  inner_join(all_urban, by = c("country", "registry")) %>%
  filter(!is.na(urban))

# ==============================================================================
# [STEP 5] FINAL CLEANING & VARIABLE CREATION
# ==============================================================================
cat("\n[5] FINAL CLEANING & VARIABLE CREATION...\n")
final_df <- final_df %>%
  # Handle impossibilities and zeros
  mutate(cases = ifelse(py == 0, 0, cases)) %>%
  filter(py > 0) %>%
  # Standardize covariates
  mutate(
    agr = as.integer(as.factor(age)),
    urbstd = as.numeric(scale(urban))
  ) %>%
  select(
    cancer_lab, sex,
    continent, region, country, registry,
    agr, urban, urbstd,
    cases, py
  ) %>%
  arrange(sex, cancer_lab, continent, region, country, registry, agr)

# Save Final
write.csv(final_df, file = args$output, row.names = FALSE)
cat(sprintf(" -> Dataset successfully saved to: %s\n", args$output))

# ==============================================================================
# [STEP 6] DESCRIPTIVE SUMMARY
# ==============================================================================
cat("\n==================================================\n")
cat("                FINAL DATA SUMMARY                \n")
cat("==================================================\n")

# Compute base metrics
n_continents <- n_distinct(final_df$continent)
n_regions    <- n_distinct(final_df$region)
n_countries  <- n_distinct(final_df$country)
n_registries <- n_distinct(final_df$registry)

# For person years, we filter by one cancer ('Colorectal') and sum Py 
# to avoid duplicating the population base across multiple cancer types.
total_py <- final_df %>%
  filter(cancer_lab == "Colorectal") %>%
  summarise(tot = sum(py, na.rm = TRUE)) %>%
  pull(tot)

cat(sprintf("Total Observations (Rows): %s\n", format(nrow(final_df), big.mark=",")))
cat(sprintf("Total Person-Years: %s\n", format(total_py, big.mark=",")))
cat(sprintf("\nSpatial Hierarchies:\n"))
cat(sprintf("  - Continents: %d\n", n_continents))
cat(sprintf("  - Regions:    %d\n", n_regions))
cat(sprintf("  - Countries:  %d\n", n_countries))
cat(sprintf("  - Registries: %d\n", n_registries))

cat("\n[ CONTINENTS ]\n")
cat(paste(unique(final_df$continent), collapse = " | "), "\n")

cat("\n[ REGIONS ]\n")
cat(paste(unique(final_df$region), collapse = " | "), "\n")

cat("\n[ REGISTRIES ]\n")
cat(paste(unique(final_df$registry), collapse = ", "), "\n")
cat("==================================================\n\n")
