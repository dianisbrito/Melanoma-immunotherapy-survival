# ============================================================
# Melanoma Immunotherapy Survival Analysis — Interactive Dashboard
#
# Companion Shiny app to report.Rmd: lets the user interactively explore
# the four covariate-engineering strategies (Full / PCA / Clusters /
# Aitchison-CLR) across Kaplan-Meier, Cox, and logistic regression
# analyses, on a fully synthetic melanoma immunotherapy cohort.
#
# Note: all data is synthetic — see data/generate_data.R.
# Run with: shiny::runApp()
# ============================================================

library(shiny)
library(survival)
library(ggplot2)
library(DT)

source("data/generate_data.R")
source("src/analysis_functions.R")

set.seed(2024)
cohort <- generate_cohort()
cohort$age_group <- cut(cohort$age, breaks = c(0, 50, 70, 100), labels = c("20-50", "50-70", "70+"))

ui <- fluidPage(
  titlePanel("🧬 Melanoma Immunotherapy Survival Analysis"),
  tags$p(
    "Demo dashboard — synthetic data, inspired by a real applied biostatistics team consulting project. ",
    "All patients, biomarker values, and outcomes are simulated — see the README for details.",
    style = "color: #666; margin-bottom: 20px;"
  ),

  tabsetPanel(
    # ------------------------------------------------------------
    tabPanel("Kaplan-Meier",
      sidebarLayout(
        sidebarPanel(
          selectInput("km_var", "Stratify by:", choices = c("Sex" = "sex", "Age group" = "age_group")),
          width = 3
        ),
        mainPanel(
          plotOutput("km_plot", height = "450px"),
          verbatimTextOutput("km_logrank"),
          width = 9
        )
      )
    ),

    # ------------------------------------------------------------
    tabPanel("Cox Models",
      sidebarLayout(
        sidebarPanel(
          selectInput("cox_strategy", "Covariate strategy:", choices = STRATEGIES),
          sliderInput("n_boot", "Bootstrap resamples:", min = 10, max = 100, value = 40, step = 10),
          width = 3
        ),
        mainPanel(
          h4("Model summary"),
          DTOutput("cox_table"),
          h4("Concordance comparison across all strategies"),
          plotOutput("cox_concordance_plot", height = "350px"),
          h4("Bootstrap concordance (stability) for selected strategy"),
          verbatimTextOutput("cox_boot_summary"),
          width = 9
        )
      )
    ),

    # ------------------------------------------------------------
    tabPanel("Logistic Regression",
      sidebarLayout(
        sidebarPanel(
          selectInput("logistic_strategy", "Covariate strategy:", choices = STRATEGIES),
          sliderInput("k_folds", "Cross-validation folds:", min = 3, max = 10, value = 5),
          width = 3
        ),
        mainPanel(
          h4("Cross-validated performance"),
          DTOutput("logistic_metrics_table"),
          h4("ROC curve (selected strategy)"),
          plotOutput("roc_plot", height = "450px"),
          width = 9
        )
      )
    ),

    # ------------------------------------------------------------
    tabPanel("Cohort data",
      DTOutput("cohort_table")
    )
  )
)

server <- function(input, output, session) {

  # -- Kaplan-Meier tab --
  output$km_plot <- renderPlot({
    formula_str <- paste("Surv(time, event) ~", input$km_var)
    fit <- survfit(as.formula(formula_str), data = cohort)
    plot(fit, col = c("#E64B35", "#00A087", "#3C5488"), lwd = 2, conf.int = TRUE,
         xlab = "Time (days)", ylab = "Survival probability",
         main = paste("Kaplan-Meier survival by", input$km_var))
    legend("bottomleft", legend = names(fit$strata), col = c("#E64B35", "#00A087", "#3C5488"), lwd = 2)
  })

  output$km_logrank <- renderPrint({
    formula_str <- paste("Surv(time, event) ~", input$km_var)
    lr <- survdiff(as.formula(formula_str), data = cohort)
    p_val <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
    cat(sprintf("Log-rank test: chi-squared = %.2f, df = %d, p-value = %.4f\n",
                lr$chisq, length(lr$n) - 1, p_val))
    if (p_val < 0.05) {
      cat("=> Statistically significant difference in survival between groups (alpha = 0.05)")
    } else {
      cat("=> No statistically significant difference in survival between groups (alpha = 0.05)")
    }
  })

  # -- Cox Models tab --
  cox_fit_selected <- reactive({
    fit_cox_model(cohort, input$cox_strategy)
  })

  output$cox_table <- renderDT({
    s <- summary(cox_fit_selected())$coefficients
    df <- as.data.frame(round(s, 4))
    df$Variable <- rownames(df)
    df <- df[, c("Variable", setdiff(colnames(df), "Variable"))]
    datatable(df, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$cox_concordance_plot <- renderPlot({
    conc_df <- do.call(rbind, lapply(STRATEGIES, function(strat) {
      fit <- fit_cox_model(cohort, strat)
      data.frame(strategy = strat, concordance = summary(fit)$concordance[1])
    }))
    conc_df$strategy <- factor(conc_df$strategy, levels = STRATEGIES)
    conc_df$selected <- conc_df$strategy == input$cox_strategy

    ggplot(conc_df, aes(x = strategy, y = concordance, fill = selected)) +
      geom_col(show.legend = FALSE) +
      scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#B0BEC5")) +
      labs(x = NULL, y = "Concordance index", title = "In-sample concordance by strategy (selected in red)") +
      ylim(0, 1) +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_text(angle = 15, hjust = 1))
  })

  output$cox_boot_summary <- renderPrint({
    cidx <- bootstrap_concordance(cohort, input$cox_strategy, n_boot = input$n_boot)
    cat(sprintf("Strategy: %s\n", input$cox_strategy))
    cat(sprintf("Bootstrap concordance: mean = %.3f, sd = %.3f, range = [%.3f, %.3f]\n",
                mean(cidx), sd(cidx), min(cidx), max(cidx)))
    cat(sprintf("(based on %d valid bootstrap resamples)\n", length(cidx)))
    cat("\nNote: this quantifies model STABILITY under resampling, not out-of-sample\n")
    cat("predictive validation (see the Logistic Regression tab for true held-out evaluation).")
  })

  # -- Logistic Regression tab --
  logistic_eval <- reactive({
    evaluate_logistic_model(cohort, input$logistic_strategy, k_folds = input$k_folds)
  })

  output$logistic_metrics_table <- renderDT({
    all_results <- do.call(rbind, lapply(STRATEGIES, function(strat) {
      res <- evaluate_logistic_model(cohort, strat, k_folds = input$k_folds)
      data.frame(
        Strategy = strat, Mean_Deviance = round(mean(res$deviances), 2),
        Accuracy = round(res$metrics$accuracy, 3), Sensitivity = round(res$metrics$sensitivity, 3),
        Specificity = round(res$metrics$specificity, 3), AUC = round(res$roc$auc, 3)
      )
    }))
    datatable(all_results, rownames = FALSE, options = list(dom = "t"))
  })

  output$roc_plot <- renderPlot({
    roc <- logistic_eval()$roc
    plot(roc$fpr, roc$tpr, type = "l", lwd = 2, col = "#E64B35",
         xlab = "False Positive Rate", ylab = "True Positive Rate",
         main = sprintf("ROC curve — %s (AUC = %.3f)", input$logistic_strategy, roc$auc),
         xlim = c(0, 1), ylim = c(0, 1))
    abline(0, 1, lty = 2, col = "gray50")
  })

  # -- Cohort data tab --
  output$cohort_table <- renderDT({
    datatable(cohort, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 15))
  })
}

shinyApp(ui, server)

