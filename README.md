
# Bayesian Hierarchical Spatial Model for Cancer Incidence Rate Ratios
## Used in the study: "Urban-Rural Disparities in Cancer Incidence: A Global Ecological Study 2013-2017"

**Author:** Fergus Waterhouse (Early Career Scientist, CSU Branch)  
**Institution:** Cancer Surveillance Branch (CSU), International Agency for Research on Cancer (IARC) / World Health Organization (WHO)  

---

## Overview

A Bayesian hierarchical model used to estimate the Incidence Rate Ratio (IRR) for cancer incidence associated with a specific continuous covariate (e.g., the standardized percentage of the population living in urban areas within a registry catchment). 

The model was designed to use incidence data derived from the **Cancer Incidence in 5 Continents (CI5)** dataset. By using priors that are spatially structured using a nested hierarchy at the global, continental, regional, and country levels, the model leverages information across geographical tiers to produce robust, smoothed posterior estimates of covariate effects on cancer risk. 

## Methodology

The script (`run.R`) implements a spatial hierarchical model using Markov chain Monte Carlo (MCMC) sampling via the `nimble` R package. 

### Model Formulation

**1. Likelihood**  
The observed number of cancer cases $i$ ($y\_i$) is modeled using a Negative Binomial distribution to account for overdispersion. The expected number of cases $\mu\_i$ is the product of the person-years at risk ($n\_i$) and the incidence rate ($\lambda\_i$):

$$
y_i \sim \text{Negative Binomial}\left(\mu_i, \mu_i + \frac{\mu_i^2}{r}\right) \quad \text{where} \quad \mu_i = n_i \cdot \lambda_i
$$

Here, $r$ represents the global overdispersion parameter.

**2. Linear Predictor & Thin Plate Splines**  
The log-transformed incidence rate is defined by a country-specific continuous baseline function of age, an invariant sex effect ($\gamma$), and the country-specific covariate effect ($\delta\_{s3}$):

$$
\ln(\lambda_i) = \beta_{0,s3} + \beta_{1,s3} \times \text{age}_i + \sum_{k=1}^K b_{s3,k} \cdot z_{i,k} + \gamma \times \text{sex}_i + \delta_{s3} \times \text{cov}_i
$$

Non-linear age effects are captured non-parametrically using thin-plate splines (parameterized by $K$ knots). The radial basis function matrix $z\_{i,k} = \left| \text{age}\_i - \kappa\_k \right|^3$ calculates the absolute distance between the observed age and the knots $\kappa\_k$. The spline weights $b\_{s3,k}$ are assigned normally distributed priors centered at zero, allowing the variance parameter to act as a penalty term that balances curve smoothness with data fit.

**3. Incidence Rate Ratio (IRR)**  
The covariate (e.g., urbanization) is standardized prior to modeling. Consequently, the exponentiated parameter $\exp(\delta\_{s3})$ represents the country-specific **Incidence Rate Ratio (IRR)**, signifying the relative change in the incidence rate for a one-standard-deviation increase in the covariate.

**4. Spatial Hierarchical Structure & Priors for the Covariate Effect ($\delta$)**  
Under the assumption that geographic proximity implies similar incidence patterns, the model employs a nested spatial hierarchy. This is particularly important for the covariate effect ($\delta$), which is modeled recursively from the global level down to the country level:

* **Global Baseline:** $\delta\_0 \sim N(0, 1)$
* **Tier 1 (Continent):** $\delta\_c \sim N(\delta\_0, \sigma\_{\delta\_c})$
* **Tier 2 (Region):** $\delta\_m \sim N(\delta\_c, \sigma\_{\delta\_m})$
* **Tier 3 (Country):** $\delta\_s \sim N(\delta\_m, \sigma\_{\delta\_s})$

**Hyperpriors:** The standard deviation parameters controlling the variance between spatial tiers ($\sigma\_{\delta\_c}$, $\sigma\_{\delta\_m}$, $\sigma\_{\delta\_s}$) are modeled using Half-Normal distributions to regularize the estimates toward the spatial mean (shrinkage). 
* At the continental and regional levels, these are specified as weakly informative: $\sigma \sim N^+(0, 1.0)$. 
* At the country level, the hyperprior is dynamically specified based on cancer-specific informativity: $\sigma\_{\delta\_s} \sim N^+(0, \Sigma)$. For highly informative, well-sampled sites (e.g., Lung, Breast, Colorectal), tighter priors are preferred (e.g., $\Sigma = 0.05$), while sparser sites (e.g., Liver) utilize broader priors (e.g., $\Sigma = 0.5$). This $\Sigma$ value is controlled via the `--sigma` command-line argument.

*Note: Initial values for the MCMC chains are computationally derived using a frequentist mixed-effects Poisson model (`glmer` from `lme4`) to ensure efficient convergence.*

## Prerequisites and Dependencies

The following R packages are required:

* `nimble`
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
| `<predictor>` | Column name of the covariate/predictor variable in the dataset (e.g., `urbstd`). |
| `<cancer>` | Cancer site label. |
| `<sex>` | Integer indicating sex (`0` = Both, `1` = Male, `2` = Female). |
| `<num_iter>` | Total number of MCMC iterations per chain. |
| `<num_burn>` | Number of iterations to discard. |

### Options

| Option | Default | Description |
| :--- | :--- | :--- |
| `--num_chains` | `4` | Number of parallel MCMC chains (must be $> 1$). |
| `--knots` | `2` | Number of internal knots used for the thin-plate splines. |
| `--sigma` | `1.0` | Standard deviation ($\Sigma$) for the country-level variation Half-Normal hyperprior. |
| `-h, --help` | | Show the help message and exit. |

### Example Command

```bash
Rscript run.R data/model_ready.csv results/colorectal urbstd "Colorectal" 0 200000 100000 --num_chains=4 --knots=4 --sigma=0.5
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
* `registry`: Registry-ID.
* `[predictor]`: The continuous covariate of interest (e.g., `urbstd`).

The **`data/process.R`** script was used in the Urban-Rural study and synthesises data from CI5XII and NORDCAN. See **`data/README.md`** for more information on the processing pipeline.

## Output Files

The script generates multiple outputs within the specified `<out_dir>`:

1. **`meta_info.csv`**: A summary of the execution arguments and metadata.
2. **`mcmc_out.rds`**: The raw, nested NIMBLE MCMC output object containing posterior samples for all tracked parameters.
3. **`mcmc_summary_report.csv`**: A statistical summary of the MCMC chains, including posterior medians, 95% Credible Intervals (CrI), Gelman-Rubin convergence diagnostics ($\hat{R}$), and Effective Sample Sizes (ESS) for all parameters.
4. **`mcmc_waic_report.csv`**: The Watanabe-Akaike Information Criterion (WAIC) for each chain and the overall mean, utilized for model evaluation and comparison.

## Acknowledgments

This tool was developed at the **Cancer Surveillance Branch (CSU)** of the **International Agency for Research on Cancer (IARC)**. It forms part of an ongoing effort to understand spatial inequalities and the impact of socio-demographic indicators on global cancer incidence patterns.
