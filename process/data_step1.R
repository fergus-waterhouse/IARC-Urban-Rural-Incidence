"Run data processing.

Usage:
		run.R <file> <urban> <output>
		run.R (-h | --help)

Options:
		-h --help                       Show this help message and exit.
		
Arguments:
    <file>          Cases data.
    <urban>         Urban data (as %).
    <output>        Path for output.
" -> doc


# ------ SETUP ------
suppressPackageStartupMessages(library(tidyverse))
if (!interactive()) {
  args <- docopt::docopt(doc)
} else {
  args <- list(file = "data/CI5XII.csv", urban = "data/urban.csv", output = "data/test.csv") 
}

# --- HELPER FUNCTION: TRACK DATA ATTRITION ---
track_data <- function(df, step_name) {
  cat(sprintf("\n--- STATE: %s ---\n", step_name))
  cat(sprintf("Total Rows: %s\n", format(nrow(df), big.mark=",")))
  cat(sprintf("Unique Registries: %s\n", length(unique(df$registry))))
  if("cases" %in% names(df)) {
    cat(sprintf("Total Cases: %s\n", format(sum(df$cases, na.rm = TRUE), big.mark=",")))
  }
}

# ------ DATA/GEOGRAPHICAL STRUCTURE DEFINITION ------
country_regions <- c(
  
  # --- EUROPE ---
  
  # Northern Europe
  "Denmark" = "N Europe",
  "Finland" = "N Europe",
  "Iceland" = "N Europe",
  "Norway" = "N Europe",
  "Sweden" = "N Europe",
  "UK" = "N Europe",
  "Latvia" = "N Europe",
  "Lithuania" = "N Europe",
  "Estonia" = "N Europe",
  
  # Western Europe
  "Austria" = "W Europe",
  "France" = "W Europe",
  "Germany" = "W Europe",
  "Liechtenstein" = "W Europe",
  "Switzerland" = "W Europe",
  "The Netherlands" = "W Europe",
  "Ireland" = "W Europe",
  "Belgium" = "W Europe",
  
  # Eastern Europe
  "Belarus" = "E Europe",
  "Bulgaria" = "E Europe",
  "Czech Republic" = "E Europe",
  "Poland" = "E Europe",
  "Russian Federation" = "E Europe",
  "Slovakia" = "E Europe",
  "Ukraine" = "E Europe",
  
  # Southern Europe
  "Croatia" = "S Europe",
  "Italy" = "S Europe",
  "Malta" = "S Europe",
  "Portugal" = "S Europe",
  "Slovenia" = "S Europe",
  "Spain" = "S Europe",
  "Cyprus" = "S Europe",
  
  # --- ASIA ---
  
  # E.Asia
  "China" = "E Asia",
  "Japan" = "E Asia",
  "Republic of Korea" = "E Asia",
  
  # S.E.Asia
  "Thailand" = "SE Asia",
  "Singapore" = "SE Asia",
  "Philippines" = "SE Asia",
  "Brunei Darussalam" = "SE Asia",
  
  # Central and Southern Asia
  "Iran (Islamic Republic of)" = "C&S Asia",
  "India" = "C&S Asia",
  
  # --- Northern America ---
  "Canada" = "N America",
  "USA" = "N America",
  
  # --- Southern America ---
  "Argentina" = "S America",
  "Brazil" = "S America",
  "Chile" = "S America",
  "Colombia" = "S America",
  "Ecuador" = "S America",
  "Peru" = "S America",
  "Uruguay" = "S America",
  "Costa Rica" = "S America",
  "Trinidad and Tobago" = "S America",
  
  # --- AFRICA ---
  "Turkey" = "N Africa & W Asia",
  "Israel" = "N Africa & W Asia",
  "Kuwait" = "N Africa & W Asia",
  "Algeria" = "N Africa & W Asia",
  "Morocco" = "N Africa & W Asia",
  "Benin" = "N Africa & W Asia",
  
  "Kenya" = "S-S Africa",
  "Mauritius" = "S-S Africa",
  "Seychelles" = "S-S Africa",
  "South Africa" = "S-S Africa",
  "Uganda" = "S-S Africa",
  "Zimbabwe" = "S-S Africa",
  
  # --- OCEANIA ---
  "Australia" = "Oceania",
  "New Zealand" = "Oceania"
)

region_continents <- c(
  # Americas
  "N America" = "N America",
  "S America" = "S America",
  
  # Europe
  "N Europe" = "Europe",
  "W Europe" = "Europe",
  "E Europe" = "Europe",
  "S Europe" = "Europe",
  
  # Asia
  "C&S Asia" = "Asia",
  "E Asia" = "Asia",
  "SE Asia" = "Asia",
  
  # Africa
  "N Africa & W Asia" = "Africa & W Asia",
  "S-S Africa" = "Africa & W Asia",
  
  # Oceania
  "Oceania" = "Oceania"
)

# ------ LOAD FILES AND RENAME VARIABLE ------
cat("\n[1] LOADING DATA...\n")
cases <- read.csv(args$file) %>% 
  rename(reglab = registry_lab, period = registry_years)
urban <- read.csv(args$urban) %>%
  rename(urban = Urban)

# ------ STANDARDISE NAMES, ASSIGN GEOGRAPHY & SELECT VARIABLES ------
cat("\n[2] STANDARDIZING AND ASSIGNING GEOGRAPHY...\n")
cases <- cases %>% 
  mutate(reglab = ifelse(reglab == "Turkey, Eski?ehir", "Turkey, Eskişehir", reglab)) %>%
  mutate(
    country = sub(",.*", "", reglab),
    registry = sub("^[^,]+, ", "", reglab),
    region = country_regions[country],
    continent = region_continents[region]
  ) %>%
  select(period, continent, region, country, registry, n_agr, cancer_lab, sex, age, cases, py, ethnic_group)

urban <- urban %>% 
  mutate(
    reglab = trimws(reglab, "right"),
    country = sub(",.*", "", reglab),
    registry = sub("^[^,]+, ", "", reglab),
    urban = urban / 100 # As a proportion (0, 1)
  ) %>%
  select(country, registry, urban)

# Track initial state
track_data(cases, "Initial Cases Data Loaded")

# ------ REMOVE ETHNIC SPECIFIC POPULATIONS ------  
cat("\n[3] FILTERING FOR 'ALL' ETHNICITIES (ethnic_group = 99)...\n")
cases <- cases %>% filter(ethnic_group == 99) %>% select(-ethnic_group)

track_data(cases, "After filter of ethnic group")

# Check for unmatched geographies:
unmapped <- cases %>% filter(is.na(region)) %>% pull(country) %>% unique()
if(length(unmapped) > 0) {
  warning_text <- sprintf("WARNING: The following countries were not found in the geography dictionary and are assigned NA: %s\n", 
                          paste(unmapped, collapse = ", "))
  cat(cli::col_red(warning_text))
}

# ------ MERGE ------  
cat("\n[3] MERGING CASES WITH URBAN DATA...\n")
inc <- left_join(cases, urban, by = c("country", "registry"))

# Look at registries lost/missing urban
missing_urban_regs <- inc %>% filter(is.na(urban)) %>% pull(registry) %>% unique()
cat(sprintf("Number of distinct registries missing urban data: %d\n", length(missing_urban_regs)))

missing_urban <- sum(is.na(inc$urban))
cat(sprintf("Rows missing 'urban' data: %s (%.1f%%)\n", 
            format(missing_urban, big.mark=","), 
            (missing_urban/nrow(inc))*100))

cat(sprintf("Registries missing urban variable:\n%s\n",
  paste(inc %>% filter(is.na(urban)) %>% pull(registry) %>% unique(), collapse = ", ")
  ))

# Track
inc <- inc %>% filter(!is.na(urban))
track_data(inc, "After Filter of Urban")

# ------ FILTER OUT IMPOSSIBILITIES (cases > 0 & py == 0) ------
cat("\n[4] REMOVING IMPOSSIBILITIES (cases > 0 & py == 0)...\n")

# Count structural errors before fixing
py_zero_cases_not_zero <- inc %>% filter(cases > 0 & py == 0, na.rm = TRUE)
cat(sprintf("Rows where cases > 0 but py == 0 (to be overwritten to 0): %d\n", py_zero_cases_not_zero %>% nrow()))
cat(sprintf("From the Populations: %s\n", 
            paste(py_zero_cases_not_zero %>% pull(registry) %>% unique(), collapse = ", ")))

inc <- mutate(inc, cases = ifelse(py == 0, 0, cases))

# ------ FILTER OUT PY = 0 ------
cat("\n[5] FILTERING PY = 0 POPULATIONS...\n")

inc <- filter(inc, py > 0)

track_data(inc, "py > 0")

# ------ FILTER OUT SPURIOUS REGISTRIES (VOLOGDA REGION) ------
cat("\n[6] FILTERING SPURIOUS REGISTRIES (VOLOGDA REGION)...\n")

inc <- filter(inc, !(registry %in% c("Vologda Region", "Japan", "Martinique", "Guadeloupe")))
# Vologda Region for unreliability in the data
# Japan as it is a national registry to avoid repeated data
# Martinique & Guadeloupe as they are French overseas territories 

track_data(inc, "Final Dataset")

# ------ SAVE ------
cat("\n[7] SAVING OUTPUT...\n")
write.csv(inc, file = args$output, row.names = FALSE)
cat(sprintf("File saved successfully to: %s\n", args$output))

# ------ SUMMARY------
cat("\n[7] SUMMARY OF OUTPUT...\n")
cat(sprintf("# Continents: %s\n", pull(inc, continent) %>% unique() %>% length()))
cat(sprintf("# Regions: %s\n", pull(inc, region) %>% unique() %>% length()))
cat(sprintf("# Countries: %s\n", pull(inc, country) %>% unique() %>% length()))
cat(sprintf("# Populations: %s\n", pull(inc, registry) %>% unique() %>% length()))
cat(sprintf("# Observations %s\n", inc %>% nrow()))
