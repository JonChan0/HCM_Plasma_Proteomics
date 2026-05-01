#Script to output summary plots across different model types for performance in CV during training..
#Author: Jonathan Chan
#Date: 2025-02-24

library(tidyverse)
theme_set(theme_classic())

args <- commandArgs(trailingOnly = T)

input_folder <- args[1] #Defines the input folder of cv_results.csv e.g. '../OUTPUT/UKB/ML/2_models/1_hcm_cc_noprs/'
output_folder <- args[2] #Defines the output folder e.g. '../OUTPUT/UKB/ML/3_summary_plots/cv_performance/'
model_class_of_interest <- args[3]

#Import--------------------------------------------------------------------------
csv_importer <- function(input_folder) {
  csv_filenames <- list.files(input_folder, pattern = 'cv_results.csv')
  model_names <- str_match(csv_filenames, '(.+)_cv_results\\.csv')[, 2]

  cv <- map(str_c(input_folder, csv_filenames), ~ read_csv(.))
  names(cv) <- model_names

  return(cv)
}

#Per-Model evaluation----------------------------------------------------------------
## This evaluates model performance within a single model class, evaluating how different hyperparameter combinations perform

### Comparing the rank of each performance metric across each of the hyperparameter combinations for a model type
### Comparing via a grouped errorbar plot each of the hyperparameter combinations across each score metric

permodel_comp <- function(
  input_cv,
  model_class,
  output_path,
  scores_of_interest = c(
    "roc_auc",
    "f1",
    "average_precision",
    "balanced_accuracy"
  )
) {
  if (!dir.exists(output_path)) {
    dir.create(output_path)
  }

  plot_tb1 <- input_cv %>%
    select(params, contains('rank')) %>%
    pivot_longer(
      contains('rank'),
      names_to = 'score_metric',
      values_to = 'rank'
    ) %>%
    mutate(score_metric = str_match(score_metric, 'rank_test_(.+)')[, 2]) %>%
    filter(score_metric %in% scores_of_interest)

  output_rank_plot <- ggplot(
    plot_tb1,
    aes(y = rank, x = score_metric, group = params, color = params)
  ) +
    geom_point(size = 2) +
    geom_line() +
    xlab('Score Metric') +
    ylab('Model Rank') +
    scale_y_continuous(
      transform = 'reverse',
      n.breaks = length(unique(plot_tb1$params))
    ) +
    theme(axis.text.x = element_text(angle = 90)) +
    guides(col = guide_legend(position = 'right', ncol = 1)) +
    labs(col = 'Params')

  # print(output_rank_plot)

  plot_tb2 <- input_cv %>%
    select(-contains('time')) %>%
    select(params, contains('mean'), contains('std')) %>%
    pivot_longer(
      contains('mean'),
      names_to = 'score_metric',
      values_to = 'mean'
    ) %>%
    mutate(score_metric = str_match(score_metric, 'mean_test_(.+)')[, 2]) %>%
    pivot_longer(
      contains('std'),
      names_to = 'score_metric2',
      values_to = 'std'
    ) %>%
    mutate(score_metric2 = str_match(score_metric2, 'std_test_(.+)')[, 2]) %>%
    filter(score_metric == score_metric2) %>%
    select(-score_metric2) %>%
    filter(score_metric %in% scores_of_interest)

  output_grouped_errorbar_plot <- ggplot(
    plot_tb2,
    aes(y = score_metric, x = mean, group = params)
  ) +
    geom_point(aes(color = params), position = position_dodge(0.9)) +
    geom_errorbar(
      aes(xmin = mean - std, xmax = mean + std, color = params),
      position = position_dodge(0.9),
      width = 0.2
    ) +
    # theme(axis.text.x=element_text(angle=90))+
    scale_x_continuous(limits = c(0, 1), n.breaks = 10) +
    labs(col = 'Params') +
    guides(col = guide_legend(position = 'right', ncol = 1)) +
    labs(caption = 'Error bar = SD across 5 folds') +
    ylab('Score Metric') +
    xlab('Mean Score') +
    labs(
      title = str_wrap(
        'Mean score across metrics for each model hyperparameter combination'
      )
    )

  print(output_grouped_errorbar_plot)

  ggsave(
    str_c(output_path, model_class, '_intramodel_errorbar_plot.png'),
    output_grouped_errorbar_plot,
    dpi = 600,
    width = 12,
    height = 9
  )
  ggsave(
    str_c(output_path, model_class, '_intramodel_rank_plot.png'),
    output_rank_plot,
    dpi = 600,
    width = 12,
    height = 9
  )
}

# Cross-Model Comparison----------------------------------------------------------
## This compares model performance across different model classes.
## It allows you to specify which scoring metric to decide which is the best model hyperparameter setting to use.

# Define default output folder
DEFAULT_OUTPUT_FOLDER <- "crossmodel_comparison"

# Helper function to ensure output directory exists
ensure_output_dir <- function(output_folder) {
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
    message(sprintf("Created output directory: %s", output_folder))
  }
}

# Function to extract best models from each CV result
extract_best_models <- function(
  cv_results_list,
  names = NULL,
  metric = "mean_test_roc_auc"
) {
  if (is.null(names)) {
    names <- paste0("Model", 1:length(cv_results_list))
  }

  map_df(seq_along(cv_results_list), function(i) {
    df <- cv_results_list[[i]]

    # Find the index of the best model based on the selected metric
    best_idx <- which.max(df[[metric]])

    # Extract the row with best performance
    best_row <- df[best_idx, ]

    # Extract all mean test metrics and their standard deviations
    mean_metrics <- best_row %>%
      select(starts_with("mean_test_")) %>%
      pivot_longer(
        cols = everything(),
        names_to = "metric",
        values_to = "value"
      ) %>%
      mutate(metric = str_replace(metric, "mean_test_", ""))

    # Get corresponding std metrics
    std_metrics <- best_row %>%
      select(starts_with("std_test_")) %>%
      pivot_longer(
        cols = everything(),
        names_to = "metric",
        values_to = "std"
      ) %>%
      mutate(metric = str_replace(metric, "std_test_", ""))

    # Combine mean and std
    left_join(mean_metrics, std_metrics, by = "metric") %>%
      mutate(model_type = names[i], params = best_row$params)
  })
}

# Function to plot performance comparison across models with error bars
plot_model_comparison <- function(
  cv_results_list,
  model_names = NULL,
  metrics = c('roc_auc', 'balanced_accuracy', 'f1', 'average_precision'),
  output_folder = DEFAULT_OUTPUT_FOLDER,
  filename = "model_comparison.png",
  width = 10,
  height = 7,
  dpi = 600,
  save = TRUE,
  xlimits = c(0, 1),
  filter_models = ''
) {
  if (is.null(model_names)) {
    model_names <- paste0("Model", 1:length(cv_results_list))
  }

  # Add "mean_test_" prefix to metrics if not already present
  metrics_full <- ifelse(
    !str_detect(metrics, "^mean_test_"),
    paste0("mean_test_", metrics),
    metrics
  )

  # Extract best models for each CV result with standard deviations
  best_models <- extract_best_models(
    cv_results_list,
    model_names,
    metrics_full[1]
  )

  if (filter_models != '') {
    best_models <- best_models %>% filter(str_detect(model_type, filter_models))
  }

  # Filter for requested metrics and plot with error bars
  p <- best_models %>%
    filter(metric %in% metrics) %>%
    ggplot(aes(y = model_type, x = value, color = metric)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(
      aes(xmin = value - std, xmax = value + std),
      position = position_dodge(width = 0.5),
      width = 0.2
    ) +
    scale_x_continuous(limits = xlimits) +
    labs(
      title = "Validation Performance Comparison",
      caption = "Points show mean values, error bars show standard deviation",
      y = "Model Type",
      x = "Score",
      color = "Metric"
    ) +
    theme_classic() +
    theme(
      #axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )

  # Filter for requested metrics and plot with error bars
  p2 <- best_models %>%
    filter(metric %in% metrics) %>%
    ggplot(aes(y = model_type, x = value, color = metric)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(
      aes(xmin = value - std, xmax = value + std),
      position = position_dodge(width = 0.5),
      width = 0.2
    ) +
    geom_text(
      aes(label = sprintf("%.2f ± %.2f", value, std)),
      position = position_dodge(width = 0.5),
      vjust = -0.5,
      size = 3.5
    ) +
    scale_x_continuous(limits = xlimits) +
    labs(
      title = "Validation Performance Comparison (With Mean & SD)",
      caption = "Points show mean values, error bars show standard deviation",
      y = "Model Type",
      x = "Score",
      color = "Metric"
    ) +
    theme_classic() +
    theme(legend.position = "top")

  # Save the plot if requested
  if (save) {
    ensure_output_dir(output_folder)
    output_path <- file.path(output_folder, filename)
    ggsave(output_path, plot = p, width = width, height = height, dpi = dpi)
    message(sprintf("Plot saved to: %s", output_path))
    output_path <- str_c(output_folder, 'model_comparison_labelled.png')
    ggsave(output_path, plot = p2, width = width, height = height, dpi = dpi)
    message(sprintf("Plot saved to: %s", output_path))
  }

  return(p)
}

##--------------------------------------------------------------------------------
## Main

main <- function(input_folder, output_folder, model_class_of_interest) {
  cv <- csv_importer(input_folder)

  #Output the per-model class plots
  walk2(
    cv,
    names(cv),
    ~ permodel_comp(
      .x,
      .y,
      output_path = str_glue('{output_folder}permodel_comparison/')
    )
  )

  #Output the cross-model comparison plots
  plot_model_comparison(
    cv,
    names(cv),
    output_folder = str_glue('{output_folder}crossmodel_comparison/'),
    metrics = c('roc_auc'),
    xlimits = c(0.5, 1)
  )

  #Only perform cross-model comparison for a certain class of models
  plot_model_comparison(
    cv,
    names(cv),
    output_folder = str_c(
      str_glue('{output_folder}crossmodel_comparison/'),
      model_class_of_interest
    ),
    metrics = c('roc_auc'),
    filter_models = model_class_of_interest,
    xlimits = c(0.5, 1),
    width = 9,
    height = 3
  )
}

main(input_folder, output_folder, model_class_of_interest)
