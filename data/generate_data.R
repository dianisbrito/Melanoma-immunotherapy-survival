# ============================================================
# Synthetic melanoma immunotherapy survival cohort generator
#
# Generates a fully synthetic patient cohort (n = 98, matching the scale
# of the original consulting study) with:
#   - demographics (age, sex)
#   - tumor microenvironment immune-cell composition (22 cell subtypes,
#     CIBERSORT-style proportions, generated via a Dirichlet distribution
#     so they sum to 1 by construction — compositional data)
#   - two continuous biomarker signature scores + a checkpoint-gene
#     expression score (CTLA4), all generic/fictional
#   - treatment response (binary) and survival time/event, generated from
#     a known data-generating process so a handful of variables carry a
#     genuine (simulated) signal, and the rest are pure noise
#
# No real patient, clinical, or genomic data is used anywhere in this
# repository. Cell-type names follow standard, publicly documented
# immunology nomenclature (e.g. CIBERSORT categories) — they are generic
# scientific terms, not tied to any real patient or study record.
#
# Run with: Rscript data/generate_data.R
# ============================================================

set.seed(2024)

n_patients <- 98

CELL_TYPES <- c(
  "B_cells_naive", "B_cells_memory", "Plasma_cells",
  "T_cells_CD8", "T_cells_CD4_naive", "T_cells_CD4_memory_resting",
  "T_cells_CD4_memory_activated", "T_cells_follicular_helper",
  "T_cells_regulatory_Tregs", "T_cells_gamma_delta",
  "NK_cells_resting", "NK_cells_activated",
  "Monocytes", "Macrophages_M0", "Macrophages_M1", "Macrophages_M2",
  "Dendritic_cells_resting", "Dendritic_cells_activated",
  "Mast_cells_resting", "Mast_cells_activated",
  "Eosinophils", "Neutrophils"
)

# ------------------------------------------------------------------
# Dirichlet sampler (base R, no extra package needed)
# ------------------------------------------------------------------
rdirichlet <- function(n, alpha) {
  k <- length(alpha)
  x <- matrix(rgamma(n * k, shape = alpha), nrow = n, ncol = k, byrow = TRUE)
  x / rowSums(x)
}

# Cell-type "concentration" parameters — a few cell types get higher
# alpha (i.e. tend to be more abundant), mirroring a realistic
# CIBERSORT-style composition where T cells / macrophages dominate.
alpha_base <- c(
  2, 1.5, 1,      # B/plasma
  6, 3, 4, 3, 2,  # T cell subsets
  2, 0.8,         # Tregs, gamma-delta
  2, 1.5,         # NK
  3, 3, 2.5, 2,   # Monocytes/Macrophages
  1.5, 1,         # Dendritic
  1, 0.8,         # Mast
  1, 1            # Eosinophils, Neutrophils
)

generate_cohort <- function(n = n_patients) {
  age <- round(rnorm(n, mean = 60, sd = 14))
  age <- pmin(pmax(age, 20), 90)
  sex <- sample(c("M", "F"), n, replace = TRUE, prob = c(0.65, 0.35))

  # Compositional tumor microenvironment (sums to 1 per patient)
  composition <- rdirichlet(n, alpha_base)
  colnames(composition) <- CELL_TYPES
  composition <- as.data.frame(composition)

  # Generic continuous biomarker scores (independent of composition)
  signature_score_1 <- rnorm(n, mean = 5, sd = 3)          # analogous to "Fus"
  signature_score_2 <- rnorm(n, mean = 0, sd = 1)          # analogous to "hc"
  protein_score <- rnorm(n, mean = 0, sd = 1)              # analogous to "hc_prot"
  ctla4_expression <- rnorm(n, mean = 8, sd = 4)

  # ------------------------------------------------------------------
  # True (simulated) signal: a handful of variables genuinely drive
  # survival/response, mirroring the qualitative conclusion of the
  # original study (T follicular helper cells and dendritic cells
  # activated were significant across multiple model specifications).
  # ------------------------------------------------------------------
  linear_predictor <- (
    0.03 * signature_score_1
    - 10 * composition$T_cells_follicular_helper
    + 9 * composition$Dendritic_cells_activated
    + 4 * composition$T_cells_regulatory_Tregs
    - 0.06 * ctla4_expression
    + rnorm(n, 0, 0.35)
  )

  # Survival time via a simple proportional-hazards data-generating
  # mechanism (Weibull baseline hazard)
  shape <- 1.1
  scale_baseline <- 2300
  u <- runif(n)
  time_event <- (-log(u) / (exp(linear_predictor) / scale_baseline)) ^ (1 / shape)
  time_event <- pmin(time_event, 3000)

  # Administrative censoring at 1600 days (matches original study's
  # follow-up horizon), plus some random early censoring
  censor_time <- pmin(1600, rexp(n, rate = 1 / 1150))
  time <- pmin(time_event, censor_time)
  event <- as.integer(time_event <= censor_time)

  # Binary treatment response — correlated with the same signal but not
  # identical to the survival outcome (mirrors the original study's
  # separate logistic-regression response variable)
  response_prob <- plogis(1.6 * linear_predictor)
  response <- rbinom(n, 1, response_prob)
  response_label <- ifelse(response == 1, "Responde", "No responde")

  rna_seq_state <- "PRE"  # all patients pre-treatment biopsy, as in the original study

  cohort <- data.frame(
    patient_id = sprintf("PT-%03d", seq_len(n)),
    age = age,
    sex = sex,
    rna_seq_state = rna_seq_state,
    signature_score_1 = round(signature_score_1, 4),
    signature_score_2 = round(signature_score_2, 4),
    protein_score = round(protein_score, 4),
    ctla4_expression = round(ctla4_expression, 4),
    time = round(time, 1),
    event = event,
    response = response,
    response_label = response_label
  )

  cohort <- cbind(cohort, round(composition, 6))
  cohort
}

save_cohort <- function(output_dir = "data") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  cohort <- generate_cohort()
  write.csv(cohort, file.path(output_dir, "melanoma_cohort.csv"), row.names = FALSE)
  cohort
}

if (sys.nframe() == 0) {
  cohort <- save_cohort()
  cat(sprintf(
    "Generated synthetic cohort: %d patients, %d events (deaths), %d responders\n",
    nrow(cohort), sum(cohort$event), sum(cohort$response)
  ))
}
