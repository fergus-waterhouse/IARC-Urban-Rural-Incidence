#!/usr/bin/env Rscript
"Process NORDCAN data.

Usage:
  process_nordcan.R <input_dir> <output_file>
  process_nordcan.R (-h | --help)

Options:
  -h --help       Show this help message and exit.

Arguments:
  <input_dir>     Path to the folder containing NORDCAN csv files (pop, cases, urban).
  <output_file>   Path to save the finalized NORDCAN dataset (e.g., nordcan_clean.csv).
" -> doc

# ------ SETUP ------
suppressPackageStartupMessages({
  library(tidyverse)
  library(docopt)
})

if (!interactive()) {
  args <- docopt(doc)
} else {
  args <- list(input_dir = "data/raw/nordcan", output_file = "data/processed/nordcan_final.csv")
}

# ------ REFERENCE MAPPINGS ------
country_regions <- c(
  "Denmark" = "N Europe",
  "Finland" = "N Europe",
  "Iceland" = "N Europe",
  "Norway"  = "N Europe",
  "Sweden"  = "N Europe"
)

region_continents <- c(
  "N Europe" = "Europe"
)

# ------ PROCESSING FUNCTIONS ------
clean_names <- function(data) {
  country_codes <- c(
    "DK" = "Denmark",
    "FI" = "Finland",
    "IS" = "Iceland",
    "NO" = "Norway",
    "SE" = "Sweden"
  )
  
  result <- data %>% filter(!(Population %in% c("Faroe Islands", "Greenland")))
  codes <- sapply(result$Population, function(s) { substr(s, start = 1, stop = 2) })
  
  result$country <- country_codes[codes]
  result$reglab <- substr(result$Population, start = 5, stop = 100)
  result$reglab <- gsub("*", "", result$reglab, fixed = TRUE)
  result$Population <- NULL
  result$Number <- NULL
  
  return(result)
}

# ------ MAIN EXECUTION ------
cat("\n[1] LOADING NORDCAN DATA...\n")
agr_names <- c("X0", "X5", "X10", "X15", "X20", "X25", "X30", "X35", "X40", 
               "X45", "X50", "X55", "X60", "X65", "X70", "X75", "X80", "X85")

# 1. Load Population
all_pop <- list()
for (sex in 1:2) {
  pop_file <- file.path(args$input_dir, paste0("pop_", sex, ".csv"))
  if(!file.exists(pop_file)) stop(sprintf("Missing population file: %s", pop_file))
  
  pop <- clean_names(read.csv(pop_file)) %>%
    pivot_longer(cols = all_of(agr_names), names_to = "agecat", values_to = "pop") %>%
    mutate(py = 5 * pop, agecat = match(agecat, agr_names)) %>%
    select(country, reglab, agecat, py)
  
  all_pop[[sex]] <- pop
}

# 2. Load Urban
urban_file <- file.path(args$input_dir, "nordcan_urban.csv")
if(!file.exists(urban_file)) stop("Missing urban file: nordcan_urban.csv")
urban <- read.csv(urban_file)

# 3. Load Cases
cat("[2] PROCESSING FILES...\n")
all_results <- list()
sites <- c("Kidney", "Rectum", "Stomach", "Breast", "Oesophagus", "Thyroid", 
           "Colon", "Lung (incl. trachea and bronchus)", "Pancreas", "Prostate", "Bladder", "Liver", 
           "Cervix uteri", "Non-Hodgkin lymphoma", "All sites")

for (site in sites) {
  for (sex in 1:2) {
    file_name <- paste0("case_", tolower(site), "_", sex, ".csv")
    filepath <- file.path(args$input_dir, file_name)
    
    if (file.exists(filepath)) {
      raw <- clean_names(read.csv(filepath)) %>%
        pivot_longer(cols = all_of(agr_names), names_to = "agecat", values_to = "cases") %>%
        mutate(
          n_agr = 18,
          agecat = match(agecat, agr_names),
          cancer_lab = site,
          sex = sex,
          registry_years = "2013-2017"
        )
      
      raw <- merge(raw, all_pop[[sex]], by = c("country", "reglab", "agecat"))
      raw <- merge(raw, urban, by = c("country", "reglab"))
      
      all_results[[paste0(site, "_", sex)]] <- raw
    }
  }
}

final_data <- bind_rows(all_results)

# ------ FINAL FORMATTING TO MATCH MAIN PIPELINE ------
cat("[3] APPLYING FINAL FORMATTING...\n")

final_data <- final_data %>%
  mutate(
    region = country_regions[country],
    continent = region_continents[region]
  ) %>%
  select(period = registry_years, continent, region, country, registry = reglab, 
        cancer_lab, sex, n_agr, 
         age = agecat, cases, py, urban) %>%
  arrange(registry, cancer_lab, sex, age)

# ------ SAVE ------
write.csv(final_data, args$output_file, row.names = FALSE)
cat(sprintf("[4] SUCCESS! Saved %s rows to %s\n", format(nrow(final_data), big.mark=","), args$output_file))


# ------ SUMMARY------
cat("\n[5] SUMMARY\n")
cat(sprintf("# Continents: %s\n", pull(final_data, continent) %>% unique() %>% length()))
cat(sprintf("# Regions: %s\n", pull(final_data, region) %>% unique() %>% length()))
cat(sprintf("# Countries: %s\n", pull(final_data, country) %>% unique() %>% length()))
cat(sprintf("# Populations: %s\n", pull(final_data, registry) %>% unique() %>% length()))
cat(sprintf("# Observations %s\n", final_data %>% nrow()))