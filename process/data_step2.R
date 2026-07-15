"Prepare hierarchical IDs and final dataset for the mortality/incidence model.

Usage:
    prepare_model_data.R <input> <output>
    prepare_model_data.R (-h | --help)

Options:
    -h --help                       Show this help message and exit.

Arguments:
    <input>         Processed CSV data (output from your data processing script).
    <output>        Path to save the final model-ready CSV.
" -> doc

# ------ SETUP ------
suppressPackageStartupMessages({
    library(tidyverse)
    library(docopt)
})

if (!interactive()) {
    args <- docopt(doc)
} else {
    args <- list(input = "data/test.csv", output = "data/model_ready.csv") 
}

cat("\n[1] LOADING PROCESSED DATA...\n")
df <- read.csv(args$input)

cat("\n[1.5] COMBINING COLON AND RECTUM INTO COLORECTAL...\n")
df <- df %>%
    mutate(
        cancer_lab = if_else(cancer_lab %in% c("Colon", "Rectum"), "Colorectal", cancer_lab)
    ) %>%
    group_by(cancer_lab, sex, age, continent, region, country, registry, urban) %>%
    summarise(
        cases = sum(cases, na.rm = TRUE),
        py = first(py), 
        .groups = "drop"
    )

cat("\n[2] PREPARING AGE AND URBAN VARIABLES...\n")
df <- df %>%
    mutate(
        agr = as.integer(as.factor(age)),
        urbstd = as.numeric(scale(urban))
    )

cat("\n[3] SELECTING FINAL VARIABLES...\n")
final_df <- df %>%
    select(
        cancer_lab, sex,
        # Spatial IDs
        continent, region, country, registry,
        # Covariates
        agr, urban, urbstd,
        # Outcomes & Exposures
        cases, py
    ) %>%
    arrange(sex, cancer_lab, continent, region, country, registry, agr)

cat("\n[4] SAVING FINAL DATA...\n")
write.csv(final_df, file = args$output, row.names = FALSE)

# ------ SUMMARY ------
cat("\n[5] DATA PREPARATION COMPLETE.\n")
cat(sprintf(" - Total Rows:         %s\n", format(nrow(final_df), big.mark=",")))
cat(sprintf(" - Max Age Groups:     %d\n", max(final_df$agr)))
cat(sprintf("\nFile successfully saved to: %s\n", args$output))