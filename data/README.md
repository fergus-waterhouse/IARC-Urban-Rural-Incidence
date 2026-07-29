**raw/NORCAN** contains all the incidence data derived from NORDCAN (https://nordcan.iarc.fr/en)<br>
**raw/urban.csv** is all the urban data for each registry.

The **CI5XII.csv** file is too large to be stored on the repository. It can be sourced from https://ci5.iarc.fr/ci5-xii/download (CI5-XII summary database).

The final dataset is synthesised from this raw data using the **process.R** script to give a .csv file.<br>
Usage: Rscript process.R <ci5_file> <ci5_urban> <nordcan_dir> <outcsv>

The resulting data takes the form:

| Variable     | Type      | Description |
|:-------------|:----------|:------------|
| `cancer_lab` | String    | Cancer site label. |
| `sex`        | Integer   | Sex indicator (1 = Male, 2 = Female). |
| `continent`  | String    | Spatial Tier 1 (Continent). |
| `region`     | String    | Spatial Tier 2 (Region). |
| `country`    | String    | Spatial Tier 3 (Country). |
| `registry`   | String    | Spatial Tier 4 (Registry). |
| `agr`        | Integer   | Age-group factor standardized as an integer sequence. 5-year age bands (1 = 0-4 years ... 18 = 85+ years). |
| `urban`      | Numeric   | Proportion (0.0 to 1.0) of the registry's population living in an urban area. |
| `urbstd`     | Numeric   | Centered and scaled (Z-score) version of the `urban` variable. |
| `cases`      | Integer   | Observed incidence counts. |
| `py`         | Numeric   | Total person-years at risk for the given registry, sex, and age group. |

This is the data that serves as the input to the model **run.R**.
