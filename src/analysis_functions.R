# ============================================================
# Reusable statistical analysis functions for the melanoma
# immunotherapy survival analysis.
#
# Implements the four covariate-engineering strategies compared in the
# original study, applied to both Cox proportional-hazards and logistic
# regression models:
#   1. All raw covariates
#   2. Principal Component Analysis (PCA) scores
#   3. K-means cluster centroids
#   4. Centered log-ratio (CLR) transform — the standard implementation
#      of Aitchison's approach for compositional data
#
# Requires: survival (base R distribution / apt r-cran-survival)
# ============================================================

suppressMessages(library(survival))

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

BIOMARKER_VARS <- c("signature_score_1", "signature_score_2", "protein_score", "ctla4_expression")

# ------------------------------------------------------------------
# Covariate engineering: the four strategies
# ------------------------------------------------------------------

#' Centered log-ratio (CLR) transform — the standard Aitchison approach
#' for compositional data. A small pseudo-count avoids log(0).
clr_transform <- function(composition_matrix, pseudo_count = 1e-6) {
  x <- as.matrix(composition_matrix) + pseudo_count
  log_x <- log(x)
  geometric_mean_log <- rowMeans(log_x)
  sweep(log_x, 1, geometric_mean_log, "-")
}

#' PCA covariates: retain components explaining a target cumulative
#' variance (default 90%; the original study used 99%, but with far
#' fewer synthetic patients here we use a more conservative default
#' to keep the demo numerically stable).
build_pca_covariates <- function(data, variance_threshold = 0.90) {
  numeric_vars <- c(BIOMARKER_VARS, CELL_TYPES)
  X <- scale(data[, numeric_vars])
  pca <- prcomp(X, center = FALSE, scale. = FALSE)
  var_explained <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
  n_components <- min(which(var_explained >= variance_threshold))
  n_components <- max(n_components, 2)
  scores <- as.data.frame(pca$x[, seq_len(n_components), drop = FALSE])
  colnames(scores) <- paste0("pca", seq_len(n_components))
  list(covariates = scores, model = pca, n_components = n_components,
       variance_explained = var_explained[n_components])
}

#' K-means cluster covariates: cluster the cell-type variables into k
#' clusters, then use each cluster's mean proportion (per patient) as a
#' covariate — mirroring the original study's cluster-centroid approach.
build_cluster_covariates <- function(data, k = 5, seed = 42) {
  set.seed(seed)
  X <- scale(t(as.matrix(data[, CELL_TYPES])))  # cluster the VARIABLES, not patients
  km <- kmeans(X, centers = k, nstart = 20)

  cluster_assignment <- km$cluster
  cluster_scores <- sapply(seq_len(k), function(cl) {
    vars_in_cluster <- names(cluster_assignment)[cluster_assignment == cl]
    if (length(vars_in_cluster) == 1) {
      data[[vars_in_cluster]]
    } else {
      rowMeans(data[, vars_in_cluster, drop = FALSE])
    }
  })
  colnames(cluster_scores) <- paste0("cluster", seq_len(k))
  list(covariates = as.data.frame(cluster_scores), assignment = cluster_assignment)
}

#' CLR covariates as a data frame, ready to bind into a model formula.
build_clr_covariates <- function(data) {
  clr <- clr_transform(data[, CELL_TYPES])
  colnames(clr) <- paste0("clr_", CELL_TYPES)
  as.data.frame(clr)
}

# ------------------------------------------------------------------
# Model fitting per strategy
# ------------------------------------------------------------------

STRATEGIES <- c("Full covariates", "PCA", "Clusters", "Aitchison (CLR)")

#' Fit a covariate-engineering strategy on training data, returning both
#' the training covariates and a transformer function that applies the
#' SAME fitted transformation (PCA rotation, cluster assignment) to new
#' data. This matters methodologically: refitting PCA/clusters separately
#' on held-out data would leak information and produce incomparable
#' covariate spaces between train and test folds.
fit_strategy <- function(data, strategy) {
  if (strategy == "Full covariates") {
    vars <- c(BIOMARKER_VARS, CELL_TYPES)
    return(list(
      covariates = data[, vars],
      apply = function(newdata) newdata[, vars]
    ))
  }

  if (strategy == "PCA") {
    numeric_vars <- c(BIOMARKER_VARS, CELL_TYPES)
    train_means <- colMeans(data[, numeric_vars])
    train_sds <- apply(data[, numeric_vars], 2, sd)
    X <- scale(data[, numeric_vars], center = train_means, scale = train_sds)
    pca <- prcomp(X, center = FALSE, scale. = FALSE)
    var_explained <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
    n_components <- max(min(which(var_explained >= 0.90)), 2)

    project <- function(newdata) {
      Xn <- scale(newdata[, numeric_vars], center = train_means, scale = train_sds)
      scores <- as.data.frame(Xn %*% pca$rotation[, seq_len(n_components), drop = FALSE])
      colnames(scores) <- paste0("pca", seq_len(n_components))
      scores
    }
    return(list(covariates = project(data), apply = project, n_components = n_components))
  }

  if (strategy == "Clusters") {
    k <- 5
    set.seed(42)
    Xt <- scale(t(as.matrix(data[, CELL_TYPES])))
    km <- kmeans(Xt, centers = k, nstart = 20)
    assignment <- km$cluster

    project <- function(newdata) {
      scores <- sapply(seq_len(k), function(cl) {
        vars_in_cluster <- names(assignment)[assignment == cl]
        if (length(vars_in_cluster) == 1) newdata[[vars_in_cluster]]
        else rowMeans(newdata[, vars_in_cluster, drop = FALSE])
      })
      colnames(scores) <- paste0("cluster", seq_len(k))
      as.data.frame(scores)
    }
    return(list(covariates = project(data), apply = project, assignment = assignment))
  }

  if (strategy == "Aitchison (CLR)") {
    # CLR is applied row-wise (per patient) — no train-specific fitting
    # is needed, so train and test transform identically and independently.
    project <- function(newdata) {
      clr <- clr_transform(newdata[, CELL_TYPES])
      colnames(clr) <- paste0("clr_", CELL_TYPES)
      as.data.frame(clr)
    }
    return(list(covariates = project(data), apply = project))
  }

  stop("Unknown strategy: ", strategy)
}

#' Convenience wrapper matching the previous API (in-sample covariates only).
get_covariates <- function(data, strategy) fit_strategy(data, strategy)$covariates

fit_cox_model <- function(data, strategy) {
  covars <- get_covariates(data, strategy)
  model_data <- cbind(data[, c("time", "event")], covars)
  formula_str <- paste("Surv(time, event) ~", paste(colnames(covars), collapse = " + "))
  coxph(as.formula(formula_str), data = model_data, x = TRUE)
}

fit_logistic_model <- function(data, strategy) {
  covars <- get_covariates(data, strategy)
  model_data <- cbind(data.frame(response = data$response, time = data$time), covars)
  formula_str <- paste("response ~ time +", paste(colnames(covars), collapse = " + "))
  glm(as.formula(formula_str), data = model_data, family = binomial())
}

# ------------------------------------------------------------------
# Model evaluation
# ------------------------------------------------------------------

#' Leave-one-out cross-validated concordance index for a Cox model
#' specification (refits the model n times, leaving one patient out,
#' and evaluates concordance on the omitted patient's risk ranking
#' relative to the training set — approximated here via the standard
#' quantifies model stability/variability across resamples, not
#' out-of-sample predictive performance — see the README for the distinction.
bootstrap_concordance <- function(data, strategy, n_boot = 100, seed = 42) {
  set.seed(seed)
  n <- nrow(data)
  cindices <- numeric(n_boot)

  for (b in seq_len(n_boot)) {
    boot_idx <- sample(seq_len(n), n, replace = TRUE)
    boot_data <- data[boot_idx, ]
    fit <- tryCatch(fit_cox_model(boot_data, strategy), error = function(e) NULL)
    if (is.null(fit)) {
      cindices[b] <- NA
      next
    }
    conc <- tryCatch(summary(fit)$concordance[1], error = function(e) NA)
    cindices[b] <- conc
  }
  cindices[!is.na(cindices)]
}

#' Confusion-matrix-based metrics (accuracy, sensitivity, specificity,
#' precision) for a Cox model's survival-probability-based classification
#' at a given follow-up time and probability threshold.
cox_prediction_metrics <- function(data, strategy, eval_time = 800, threshold = 0.5) {
  fit <- fit_cox_model(data, strategy)
  surv_probs <- summary(survfit(fit, newdata = get_covariates(data, strategy)), times = eval_time)$surv
  predicted <- ifelse(surv_probs >= threshold, 1, 0)  # 1 = "Responde" (predicted survivor)
  observed <- as.integer(data$time >= eval_time | data$event == 0)

  confusion_metrics(observed, predicted)
}

#' Standard confusion-matrix metrics from binary observed/predicted vectors.
confusion_metrics <- function(observed, predicted) {
  tp <- sum(observed == 1 & predicted == 1)
  tn <- sum(observed == 0 & predicted == 0)
  fp <- sum(observed == 0 & predicted == 1)
  fn <- sum(observed == 1 & predicted == 0)

  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  sensitivity <- ifelse((tp + fn) > 0, tp / (tp + fn), NA)
  specificity <- ifelse((tn + fp) > 0, tn / (tn + fp), NA)
  precision <- ifelse((tp + fp) > 0, tp / (tp + fp), NA)

  list(tp = tp, tn = tn, fp = fp, fn = fn,
       accuracy = accuracy, sensitivity = sensitivity,
       specificity = specificity, precision = precision)
}

#' Manual ROC curve + AUC (trapezoidal rule) — avoids depending on pROC.
compute_roc <- function(observed, predicted_prob) {
  thresholds <- sort(unique(c(0, predicted_prob, 1)), decreasing = TRUE)
  tpr <- numeric(length(thresholds))
  fpr <- numeric(length(thresholds))

  n_pos <- sum(observed == 1)
  n_neg <- sum(observed == 0)

  for (i in seq_along(thresholds)) {
    pred_class <- as.integer(predicted_prob >= thresholds[i])
    tpr[i] <- ifelse(n_pos > 0, sum(pred_class == 1 & observed == 1) / n_pos, 0)
    fpr[i] <- ifelse(n_neg > 0, sum(pred_class == 1 & observed == 0) / n_neg, 0)
  }

  ord <- order(fpr, tpr)
  fpr <- fpr[ord]; tpr <- tpr[ord]
  auc <- sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2)

  list(fpr = fpr, tpr = tpr, auc = auc)
}

#' Cross-validated (k-fold) logistic regression deviance and predictive
#' metrics for a given covariate strategy.
evaluate_logistic_model <- function(data, strategy, k_folds = 5, seed = 42) {
  set.seed(seed)
  n <- nrow(data)
  folds <- sample(rep(seq_len(k_folds), length.out = n))

  deviances <- numeric(k_folds)
  all_observed <- c(); all_predicted_prob <- c()

  for (fold in seq_len(k_folds)) {
    train <- data[folds != fold, ]
    test <- data[folds == fold, ]

    strat_fit <- tryCatch(fit_strategy(train, strategy), error = function(e) NULL)
    if (is.null(strat_fit)) { deviances[fold] <- NA; next }

    train_model_data <- cbind(data.frame(response = train$response, time = train$time), strat_fit$covariates)
    formula_str <- paste("response ~ time +", paste(colnames(strat_fit$covariates), collapse = " + "))
    fit <- tryCatch(glm(as.formula(formula_str), data = train_model_data, family = binomial()),
                     error = function(e) NULL)
    if (is.null(fit)) { deviances[fold] <- NA; next }

    deviances[fold] <- fit$deviance

    test_covars <- strat_fit$apply(test)
    test_model_data <- cbind(data.frame(time = test$time), test_covars)
    pred_prob <- tryCatch(
      predict(fit, newdata = test_model_data, type = "response"),
      error = function(e) rep(NA, nrow(test))
    )

    all_observed <- c(all_observed, test$response)
    all_predicted_prob <- c(all_predicted_prob, pred_prob)
  }

  valid <- !is.na(all_predicted_prob)
  observed <- all_observed[valid]
  predicted_prob <- all_predicted_prob[valid]
  predicted_class <- as.integer(predicted_prob >= 0.5)

  metrics <- confusion_metrics(observed, predicted_class)
  roc <- compute_roc(observed, predicted_prob)

  list(deviances = deviances[!is.na(deviances)], metrics = metrics, roc = roc)
}
