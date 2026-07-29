# Bayesian Hierarchical Spatial Model for Cancer Incidence Rate Ratios

**Author:** Fergus Waterhouse (Early Career Scientist)  
**Institution:** Cancer Surveillance Branch (CSU), International Agency for Research on Cancer (IARC) / World Health Organization (WHO)  

---

## Overview

This repository contains a Bayesian hierarchical modelling framework designed to estimate the Incidence Rate Ratio (IRR) for cancer incidence associated with a specific continuous covariate (e.g., the percentage of the population living in urban areas within a registry catchment). 

The model relies on incidence data derived from the **Cancer Incidence in 5 Continents (CI5)** dataset. By structuring the spatial hierarchy at the global, continental, regional, and country levels, the model leverages information across geographical tiers to produce robust, smoothed posterior estimates of covariate effects on cancer risk. 

## Methodology

The script (`run.R`) implements a spatial hierarchical model using Markov chain Monte Carlo (MCMC) sampling via the `nimble` R package. 

### Model Specification
* **Likelihood:** The observed cancer counts ($y$) are modeled using a Negative Binomial distribution to account for overdispersion relative to a Poisson baseline, utilizing person-years ($n$) as an offset.
* **Linear Predictor:** The log-rate of incidence is defined as:
  $$ \ln(\text{rate}_i) = \beta_{0,s3} + \beta_{1,s3} \times \text{age}_i + \text{splines}(\text{age}_i) + \gamma_{\text{sex}} \times \text{sex}_i + \delta_{s3} \times \text{covariate}_i $$
* **Hierarchical Structure:** Age spline coefficients and the covariate effect ($\delta$) are modeled hierarchically across three spatial tiers:
  * Tier 1: Continent
  * Tier 2: Region
  * Tier 3: Country
* **Age Effects:** Non-linear age effects are captured using thin-plate splines.

Initial values for the MCMC chains are computationally derived using a frequentist mixed-effects Poisson model (`glmer` from `lme4`) to ensure efficient convergence.

## Prerequisites and Dependencies

The model is written in R and designed to be executed via the command line interface using `docopt`. The following R packages are required:

* `nimble` (requires a working C++ compiler on your system)
* `docopt`
* `dplyr`, `tidyverse`, `tibble`
* `lme4`
* `coda`
* `parallel`
* `glue`

## Usage

The model is executed from the terminal using `run.R`. 

```bash
Rscript run.R <inc> <out_dir> <predictor> <cancer> <sex> <num_iter> <num_burn> [options]
```

### Arguments

| Argument | Description |
| :--- | :--- |
| `<inc>` | Path to the input incidence dataset (CSV format). |
| `<out_dir>` | Path to the directory where model outputs will be saved. |
| `<predictor>` | Column name of the covariate/predictor variable in the dataset. |
| `<cancer>` | Cancer site label to subset the data (e.g., `"Lung"`). |
| `<sex>` | Integer indicating the sex to model (`0` = Both, `1` = Male, `2` = Female). |
| `<num_iter>` | Total number of MCMC iterations per chain. |
| `<num_burn>` | Number of iterations to discard as burn-in per chain. |

### Options

| Option | Default | Description |
| :--- | :--- | :--- |
| `--num_chains` | `4` | Number of parallel MCMC chains (must be $> 1$). |
| `--knots` | `2` | Number of internal knots used for the thin-plate splines. |
| `--sigma` | `1.0` | Standard deviation for the country-level variation normal prior. |
| `-h, --help` | | Show the help message and exit. |

### Example Command

```bash
Rscript run.R data/ci5_incidence.csv results/lung_urban urban_percent "Lung" 1 10000 5000 --num_chains=4 --knots=3
```

## Input Data Structure

The input dataset (`<inc>`) must be a CSV file containing at least the following columns:

* `cancer_lab`: Cancer site identifier.
* `sex`: Sex (1 for male, 2 for female).
* `agr`: Age group index.
* `cases`: Observed number of cancer cases.
* `py`: Person-years at risk.
* `continent`: Tier 1 geographic label.
* `region`: Tier 2 geographic label.
* `country`: Tier 3 geographic label.
* `registry`: Registry-level identifier.
* `[predictor]`: The continuous covariate of interest (e.g., `urbstd`).

The **data/process.R** was used in the Urban-Rural study and synthesises data from CI5XII and NORDCAN. See **data/README.md** for more information.

## Output Files

Upon successful completion, the script generates a comprehensive suite of outputs within the specified `<out_dir>`:

1. **`meta_info.csv`**: A summary of the execution arguments and metadata.
2. **`mcmc_out.rds`**: The raw, nested NIMBLE MCMC output object containing posterior samples for all tracked parameters.
3. **`mcmc_summary_report.csv`**: A statistical summary of the MCMC chains, including posterior medians, 95% Credible Intervals (CrI), Gelman-Rubin convergence diagnostics ($\hat{R}$), and Effective Sample Sizes (ESS) for all parameters.
4. **`mcmc_waic_report.csv`**: The Watanabe-Akaike Information Criterion (WAIC) for each chain and the overall mean, utilized for model evaluation and comparison.

## Acknowledgments

This tool was developed at the **Cancer Surveillance Branch (CSU)** of the **International Agency for Research on Cancer (IARC)**. It forms part of an ongoing effort to understand spatial inequalities and the impact of socio-demographic indicators on global cancer incidence patterns.
