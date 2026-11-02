# Spanish PBAT — validation and cross-country link community analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23892059.svg)](https://doi.org/10.5281/zenodo.23892059)

Analysis code for:

Goicoechea, C., Wallman-Jones, A., Ciarrochi, J., Hayes, S. C., Hofmann, S. G.,
& Perakakis, P. (2026). Spanish validation of the PBAT and cross-country mapping
of process-outcome networks using link community analysis. *Journal of Contextual
Behavioral Science, 41*, 101021. https://doi.org/10.1016/j.jcbs.2026.101021

Preprint: https://doi.org/10.31234/osf.io/sxhwt_v4 · OSF project (data and
materials): https://doi.org/10.17605/OSF.IO/UJ65Z

The study validates the Spanish version of the Process-Based Assessment Tool
(PBAT) in 602 Spanish adults, and compares predictive process-outcome networks
across Spain, the U.S., Sweden and Poland using Boruta feature selection and link
community detection.

## Running the analysis

Scripts are numbered and run in order from the repository root:

```r
source("code/auxFunctions.R")                    # helpers, loaded by the others
source("code/01_Demographics.R")                 # demographic summaries
source("code/02_Descriptives.R")                 # descriptives, t-tests, correlations -> tables/
source("code/03_Boruta_es.R")                    # Boruta on the Spanish data       -> results/
source("code/03b_Boruta_sw.R")                   # Boruta on the Swedish data       -> results/
source("code/04_Boruta_plot.R")                  # Boruta heatmap                   -> figures/
source("code/05_LinkCommunities_from_table.R")   # link community networks          -> figures/
source("code/06_Invariance.R")                   # measurement invariance (ES vs SW) -> results/
source("code/07_Reliability.R")                  # reliability indices (ES)         -> results/
source("code/08_Boruta_table_docx.R")            # Boruta table, plain              -> tables/
source("code/09_Boruta_table_styled_docx.R")     # Boruta table, styled             -> tables/
```

Most file-write calls in the scripts are commented out; uncomment them to
regenerate the saved outputs.

## R packages

`Boruta` (feature importance), `linkcomm` (link community detection), `lavaan`
(measurement invariance), `psych` (reliability), `officer`/`flextable` (Word
tables), plus base R for everything else. Analyses were run under R 4.4.2.

## Layout

- `code/` — numbered analysis scripts plus `auxFunctions.R`.
- `data/` — questionnaire data (Spain; Sweden item-level from Larsson & Sundström, 2024).
- `results/` — Boruta rankings, invariance fits, reliability, normality checks.
- `tables/` — the manuscript tables as CSV, and the Boruta tables as Word files.
- `figures/` — link community networks per country, and the Boruta heatmap.

Note that the figures here are regenerated output: their community colour
assignment and node layout differ from the figures in the published article,
which were produced by an earlier run.

## License

The code is released under the MIT licence (`LICENSE`). The Spanish data (`data/data_es.csv`) are released under CC BY 4.0; the Swedish data keep the terms of their original publication (`data/LICENSE.md`).
