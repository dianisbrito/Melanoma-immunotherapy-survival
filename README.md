# Melanoma Immunotherapy Survival Analysis

**Author:** Diana Brito Hoyos — Biologist & Biostatistician | Data Analyst

🔗 **Live dashboard:** [https://lxcbd2-diana-brito.shinyapps.io/melanoma-immunotherapy-survival/](https://lxcbd2-diana-brito.shinyapps.io/melanoma-immunotherapy-survival/)
📄 **Formal report (live):** [https://dianisbrito.github.io/Melanoma-immunotherapy-survival/report.html](https://dianisbrito.github.io/Melanoma-immunotherapy-survival/report.html)

A biostatistics portfolio piece demonstrating survival analysis methodology for immuno-oncology biomarker discovery: Kaplan-Meier estimation, Cox proportional hazards regression under four different covariate-engineering strategies, and cross-validated logistic regression for binary treatment-response prediction.

> ⚠️ **Note on data and origin:** This project is inspired by a real applied-biostatistics **team consulting project** (Statistical Consulting course, Master's in Applied Statistics, UNC) on melanoma immunotherapy response, conducted for two external consultants by a 4-person consulting team including myself, using real (de-identified) patient cohort data. **All data in this repository is entirely synthetic** — see [`data/generate_data.R`](./data/generate_data.R). The synthetic cohort is generated from a known statistical process (Dirichlet-distributed compositional immune-cell data, a Weibull proportional-hazards survival mechanism, and a logistic response model) that reproduces the *methodology and statistical structure* of the original study, without using any real patient, clinical, or genomic data. Immune cell-type names follow standard, publicly documented immunology nomenclature (CIBERSORT-style categories) and are not tied to any real patient or institution.

---

## Why this demonstrates biostatistics rigor

- **Correct handling of compositional data**: immune cell-type proportions sum to 1 by construction, which violates the independence assumptions of standard regression. This project implements the **centered log-ratio (CLR) transform** (Aitchison's approach) from scratch, alongside PCA and cluster-based dimensionality reduction, and compares all three against a naive "use raw proportions" baseline.
- **Distinguishing model stability from predictive validation**: the dashboard and report explicitly separate **bootstrap-resampled concordance** (a stability/variance check) from **k-fold cross-validated logistic performance** (genuine out-of-sample assessment) — a distinction that is often blurred in applied work but matters for correctly interpreting reported metrics.
- **Honest reporting of a null/modest result**: cross-validated predictive performance is modest across all four strategies given n≈98 and 22–26 covariates — reported and interpreted honestly (a demonstration of biomarker-discovery challenges in small cohorts) rather than cherry-picked.
- **No train/test leakage**: covariate-engineering transformers (PCA rotation, cluster assignment) are fit exclusively on training folds and applied — not refit — on held-out data during cross-validation.

## What's included

1. **[`report.Rmd`](./report.Rmd) / [`report.html`](./report.html)** — a formal, narrative biostatistics report: methodology, Kaplan-Meier curves with log-rank tests, all four Cox model specifications, bootstrap stability analysis, cross-validated logistic regression with ROC curves, and an honest discussion of results.
2. **[`app.R`](./app.R)** — an interactive Shiny dashboard to explore each covariate strategy live: adjustable stratification variable (Kaplan-Meier), strategy selector with live concordance comparison (Cox), and adjustable cross-validation fold count with ROC curves (logistic regression).
3. **[`src/analysis_functions.R`](./src/analysis_functions.R)** — the statistical engine: CLR transform, PCA/cluster covariate builders (with proper train-fit/test-apply separation), Cox and logistic model fitting, bootstrap concordance, manual ROC/AUC computation, and cross-validated evaluation.
4. **[`data/generate_data.R`](./data/generate_data.R)** — the synthetic cohort generator.

## Running locally

```r
# From the repository root, in R or RStudio:
source("data/generate_data.R")   # optional — app.R and report.Rmd call this automatically

# Render the formal report:
rmarkdown::render("report.Rmd")

# Launch the interactive dashboard:
shiny::runApp(".")
```

### Required R packages

```r
install.packages(c("shiny", "survival", "ggplot2", "DT", "rmarkdown", "knitr"))
```

## Tech stack

`R` · `survival` (Cox, Kaplan-Meier) · `Shiny` · `ggplot2` · `rmarkdown`/`knitr` · statistical methods: CLR/Aitchison transform, PCA, k-means clustering, bootstrap resampling, k-fold cross-validation, ROC/AUC

