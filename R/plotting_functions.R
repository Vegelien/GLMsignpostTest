# ============================================================================
# Plotting Functions for Power and Estimation Results
# ============================================================================

library(DBI)
library(RSQLite)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# ============================================================================
# POWER PLOTS
# ============================================================================

#' Create Power View in Database
#' 
#' Creates a view that calculates power (proportion of p < 0.05) for each
#' parameter setting and gamma value
create_power_view <- function(db_path = "power_simulations.db") {
  con <- dbConnect(SQLite(), db_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  dbExecute(con, "DROP VIEW IF EXISTS Power;")
  
  dbExecute(con, "
    CREATE VIEW Power AS
    SELECT
      p.*,
      COUNT(r.id) AS num_results,
      AVG(CASE WHEN r.p_value00 < 0.05 THEN 1.0 ELSE 0 END) AS power_gamma0,
      AVG(CASE WHEN r.p_value01 < 0.05 THEN 1.0 ELSE 0 END) AS power_gamma01,
      AVG(CASE WHEN r.p_value02 < 0.05 THEN 1.0 ELSE 0 END) AS power_gamma02,
      AVG(CASE WHEN r.p_value03 < 0.05 THEN 1.0 ELSE 0 END) AS power_gamma03,
      AVG(CASE WHEN r.p_value04 < 0.05 THEN 1.0 ELSE 0 END) AS power_gamma04,
      AVG(CASE WHEN r.p_value05 < 0.05 THEN 1.0 ELSE 0 END) AS power_gamma05
    FROM parameters p
    JOIN results r ON p.id = r.param_id
    GROUP BY p.id;
  ")
  
  message("✓ Power view created")
}


#' Plot Power Grid
#' 
#' Creates faceted power plot comparing test types and specifications
#' Rows = test types, Columns = specifications (3 columns only)
#' 
#' @param db_path Path to power database
#' @param model_filter Filter by GLM_model (default: "logistic")
#' @return ggplot object
plot_power_comparison <- function(db_path = "power_simulations.db",
                                  model_filter = "logistic") {
  
  con <- dbConnect(SQLite(), db_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  # Load power data
  power_data <- dbGetQuery(con, "SELECT * FROM Power")
  
  # Filter by model
  power_data <- power_data %>%
    filter(GLM_model == model_filter)

  # Create specification label - only 3 categories
  power_data <- power_data %>%
    mutate(
      spec_label = case_when(
        model_specification == "well_specified" & is.na(n_beta) ~ "Well-specified\n(Oracle)",
        model_specification == "well_specified" & !is.na(n_beta) ~ "Well-specified\n(Estimated)",
        model_specification == "misspecified" & !is.na(n_beta) ~ "Misspecified\n(Estimated)",
        TRUE ~ NA_character_
      ),
      spec_label = factor(
        spec_label,
        levels = c(
          "Well-specified\n(Oracle)",
          "Well-specified\n(Estimated)",
          "Misspecified\n(Estimated)"
        )
      )
    ) %>%
    filter(!is.na(spec_label))  # Remove misspecified oracle (not used)

  power_data <- power_data %>%
    mutate(
      test_type = recode(
        test_type,
        AS_SW_plugin = "Asymptotic",
        .default = test_type
      )
    )
  
  # Reshape to long format
  power_long <- power_data %>%
    pivot_longer(
      cols = starts_with("power_gamma"),
      names_to = "gamma",
      values_to = "power"
    ) %>%
    mutate(
      gamma = as.numeric(gsub("power_gamma0?", "0.", gamma))
    )
  
  # Compute breaks
  default_breaks <- pretty_breaks()(c(0, 1))
  all_breaks <- sort(unique(c(default_breaks, 0.05)))
  
  # Create plot - test_type as rows, spec_label as columns
  p <- ggplot(power_long, aes(x = n, y = power, color = factor(gamma))) +
    geom_line(size = 0.8) +
    geom_point(size = 1.5) +
    scale_y_continuous(limits = c(0, 1), breaks = all_breaks) +
    geom_hline(yintercept = 0.05, linetype = "dotted", color = "gray50") +
    facet_grid(test_type ~ spec_label) +
    labs(
      title = "Power Comparison Across Specifications",
      x = "Sample Size (n)",
      y = "Power",
      color = expression(gamma)
    ) +
    theme_minimal() +
    guides(color = guide_legend(nrow = 1, byrow = TRUE, title.position = "left")) +
    theme(
      text = element_text(size = 11),
      legend.position = "bottom",
      legend.box = "horizontal",
      strip.text = element_text(size = 10)
    )
  
  return(p)
}


# ============================================================================
# THETA_HAT DISTRIBUTION PLOTS
# ============================================================================

#' Load Raw Theta_hat Data
#' 
#' Loads all individual theta_hat results for boxplots
#' 
#' @param db_path Path to estimation database
#' @param model_filter Filter by GLM_model
#' @param lambda_filter Filter by lambda
#' @return Data frame in long format
load_theta_hat_data <- function(db_path = "estimation_simulations.db",
                                model_filter = "logistic",
                                lambda_filter = Inf) {
  
  con <- dbConnect(SQLite(), db_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  # Load raw data
  query <- "
    SELECT
      p.n, p.p, p.lambda, p.rho, p.GLM_model, p.model_specification, p.n_beta,
      r.replicate_id,
      r.theta_hat00, r.theta_hat01, r.theta_hat02, r.theta_hat03, r.theta_hat04,
      r.theta_hat05, r.theta_hat06, r.theta_hat07, r.theta_hat08, r.theta_hat09, r.theta_hat10
    FROM results r
    JOIN parameters p ON p.id = r.param_id
    WHERE p.GLM_model = ? AND p.lambda = ?
  "
  
  df <- dbGetQuery(con, query, params = list(model_filter, lambda_filter))
  
  # Extract theta_hat columns
  theta_cols <- grep("^theta_hat\\d\\d$", names(df), value = TRUE)
  
  # Pivot to long format
  theta_long <- df %>%
    pivot_longer(
      cols = all_of(theta_cols),
      names_to = "gamma",
      values_to = "theta_hat"
    ) %>%
    mutate(
      gamma = as.numeric(sub("theta_hat", "", gamma)) / 10
    )
  
  return(theta_long)
}


#' Plot Theta_hat Distribution
#' 
#' Creates boxplot of theta_hat distribution across specifications
#' Columns = 3 specifications only
#' 
#' @param db_path Path to estimation database
#' @param model_filter Filter by GLM_model
#' @param lambda_filter Filter by lambda
#' @param rho_boxplot Rho value to use for boxplots in the misspecified case
#'   (NULL defaults to the first available)
#' @param rho_median_lines Numeric vector of rho values to plot as median-only
#'   lines in the misspecified case
#' @return ggplot object
plot_theta_hat_distribution <- function(db_path = "estimation_simulations.db",
                                        model_filter = "logistic",
                                        lambda_filter = Inf,
                                        rho_boxplot = NULL,
                                        rho_median_lines = numeric(0)) {
  
  theta_data <- load_theta_hat_data(db_path, model_filter, lambda_filter)
  
  # Filter to only 3 specifications
  theta_data <- theta_data %>%
    mutate(
      spec_label = case_when(
        model_specification == "well_specified" & is.na(n_beta) ~ "Well-specified\n(Oracle)",
        model_specification == "well_specified" & !is.na(n_beta) ~ "Well-specified\n(Estimated)",
        model_specification == "misspecified" & !is.na(n_beta) ~ "Misspecified\n(Estimated)",
        TRUE ~ NA_character_
      ),
      spec_label = factor(
        spec_label,
        levels = c(
          "Well-specified\n(Oracle)",
          "Well-specified\n(Estimated)",
          "Misspecified\n(Estimated)"
        )
      )
    ) %>%
    filter(!is.na(spec_label))

  # Determine rho values to show for misspecified case
  misspecified_data <- theta_data %>%
    filter(spec_label == "Misspecified\n(Estimated)")

  available_rho <- sort(unique(misspecified_data$rho[!is.na(misspecified_data$rho)]))

  rho_boxplot <- if (is.null(rho_boxplot)) {
    if (length(available_rho) > 0) {
      available_rho[1]
    } else {
      NA_real_
    }
  } else {
    rho_boxplot
  }

  line_rhos <- rho_median_lines[!is.na(rho_median_lines)]
  line_rhos <- intersect(line_rhos, setdiff(available_rho, rho_boxplot))

  boxplot_data <- theta_data %>%
    filter(
      spec_label != "Misspecified\n(Estimated)" |
        (spec_label == "Misspecified\n(Estimated)" & (is.na(rho_boxplot) | rho == rho_boxplot))
    )

  line_data <- theta_data %>%
    filter(spec_label == "Misspecified\n(Estimated)", rho %in% line_rhos) %>%
    group_by(spec_label, rho, n, gamma) %>%
    summarise(theta_hat = median(theta_hat, na.rm = TRUE), .groups = "drop")

  dodge_width <- 0.75
  n_levels <- sort(unique(theta_data$n))
  gamma_levels <- sort(unique(theta_data$gamma))

  boxplot_data <- boxplot_data %>%
    mutate(
      n_pos = match(n, n_levels)
    )

  line_data <- line_data %>%
    arrange(n, rho, gamma) %>%
    mutate(
      n_pos = match(n, n_levels),
      x_pos = n_pos + dodge_width * ((match(gamma, gamma_levels) - 0.5) / length(gamma_levels) - 0.5)
    )

  p <- ggplot(boxplot_data, aes(x = n_pos, y = theta_hat)) +
    geom_boxplot(
      aes(fill = factor(gamma), group = interaction(n_pos, gamma)),
      outlier.shape = NA,
      position = position_dodge(width = dodge_width)
    ) +
    geom_line(
      data = line_data,
      aes(
        x = x_pos,
        group = interaction(n, rho),
        linetype = factor(rho)
      ),
      color = "black",
      size = 0.8,
      position = position_identity()
    ) +
    geom_point(
      data = line_data,
      aes(x = x_pos, color = factor(gamma), shape = factor(rho)),
      size = 1.5,
      position = position_identity()
    ) +
    facet_grid(. ~ spec_label) +
    scale_x_continuous(
      breaks = seq_along(n_levels),
      labels = n_levels
    ) +
    scale_y_continuous(breaks = pretty_breaks()) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = expression(paste("Distribution of ", hat(theta), " Estimates")),
      x = "Sample Size (n)",
      y = expression(hat(theta)),
      fill = expression(gamma),
      color = expression(gamma),
      linetype = expression(rho),
      shape = expression(rho)
    ) +
    theme_minimal() +
    guides(
      fill = guide_legend(nrow = 1, byrow = TRUE, title.position = "left"),
      color = guide_legend(nrow = 1, byrow = TRUE, title.position = "left"),
      linetype = guide_legend(nrow = 1, byrow = TRUE, title.position = "left"),
      shape = guide_legend(nrow = 1, byrow = TRUE, title.position = "left")
    ) +
    theme(
      text = element_text(size = 11),
      strip.text = element_text(size = 10),
      legend.position = "bottom",
      legend.box = "horizontal"
    )
  
  return(p)
}


# ============================================================================
# LOSS PLOTS
# ============================================================================

#' Load Loss Data
#' 
#' Loads target loss, null loss, and calculates relative improvement
#' 
#' @param db_path Path to estimation database
#' @param model_filter Filter by GLM_model
#' @param lambda_filter Filter by lambda
#' @return Data frame in long format
load_loss_data <- function(db_path = "estimation_simulations.db",
                           model_filter = "logistic",
                           lambda_filter = Inf) {
  
  con <- dbConnect(SQLite(), db_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  # Load raw data
  query <- "
    SELECT
      p.n, p.p, p.lambda, p.rho, p.GLM_model, p.model_specification, p.n_beta,
      e.lambda_est,
      r.replicate_id,
      r.target_loss00, r.target_loss01, r.target_loss02, r.target_loss03, r.target_loss04,
      r.target_loss05, r.target_loss06, r.target_loss07, r.target_loss08, r.target_loss09, r.target_loss10,
      r.null_loss00, r.null_loss01, r.null_loss02, r.null_loss03, r.null_loss04,
      r.null_loss05, r.null_loss06, r.null_loss07, r.null_loss08, r.null_loss09, r.null_loss10
    FROM results r
    JOIN parameters p ON p.id = r.param_id
    JOIN estimation_settings e ON e.id = r.estimation_id
    WHERE p.GLM_model = ? AND p.lambda = ?
  "
  
  df <- dbGetQuery(con, query, params = list(model_filter, lambda_filter))
  
  # Extract columns
  target_cols <- grep("^target_loss\\d\\d$", names(df), value = TRUE)
  null_cols <- grep("^null_loss\\d\\d$", names(df), value = TRUE)
  
  # Pivot to long format
  target_long <- df %>%
    select(n, rho, model_specification, n_beta, lambda_est, replicate_id, all_of(target_cols)) %>%
    pivot_longer(
      cols = all_of(target_cols),
      names_to = "gamma",
      values_to = "target_loss"
    ) %>%
    mutate(gamma = as.numeric(sub("target_loss", "", gamma)) / 10)

  null_long <- df %>%
    select(n, rho, model_specification, n_beta, lambda_est, replicate_id, all_of(null_cols)) %>%
    pivot_longer(
      cols = all_of(null_cols),
      names_to = "gamma",
      values_to = "null_loss"
    ) %>%
    mutate(gamma = as.numeric(sub("null_loss", "", gamma)) / 10)

  # Merge and calculate relative improvement
  merged <- target_long %>%
    left_join(
      null_long,
      by = c(
        "n",
        "rho",
        "model_specification",
        "n_beta",
        "lambda_est",
        "replicate_id",
        "gamma"
      ),
      na_matches = "na"
    ) %>%
    mutate(
      loss_diff = null_loss - target_loss,
      relative_improvement = loss_diff / null_loss
    )
  
  return(merged)
}


#' Plot Relative Loss Improvement
#' 
#' Creates boxplot of relative loss improvement across specifications
#' Columns = 3 specifications only
#' 
#' @param db_path Path to estimation database
#' @param model_filter Filter by GLM_model
#' @param lambda_filter Filter by lambda
#' @param y_limits Numeric vector of length 2 giving y-axis limits for the
#'   relative improvement plot
#' @param rho_boxplot Rho value to use for boxplots in the misspecified case
#'   (NULL defaults to the first available)
#' @param rho_median_lines Numeric vector of rho values to plot as median-only
#'   lines in the misspecified case
#' @return ggplot object
plot_relative_loss_improvement <- function(db_path = "estimation_simulations.db",
                                           model_filter = "logistic",
                                           lambda_filter = Inf,
                                           y_limits = c(-1, 1),
                                           rho_boxplot = NULL,
                                           rho_median_lines = numeric(0)) {
  
  loss_data <- load_loss_data(db_path, model_filter, lambda_filter)
  
  # Filter to only 3 specifications
  loss_data <- loss_data %>%
    mutate(
      spec_label = case_when(
        model_specification == "well_specified" & is.na(n_beta) ~ "Well-specified\n(Oracle)",
        model_specification == "well_specified" & !is.na(n_beta) ~ "Well-specified\n(Estimated)",
        model_specification == "misspecified" & !is.na(n_beta) ~ "Misspecified\n(Estimated)",
        TRUE ~ NA_character_
      ),
      spec_label = factor(
        spec_label,
        levels = c(
          "Well-specified\n(Oracle)",
          "Well-specified\n(Estimated)",
          "Misspecified\n(Estimated)"
        )
      )
    ) %>%
    filter(!is.na(spec_label))

  misspecified_data <- loss_data %>%
    filter(spec_label == "Misspecified\n(Estimated)")

  available_rho <- sort(unique(misspecified_data$rho[!is.na(misspecified_data$rho)]))

  rho_boxplot <- if (is.null(rho_boxplot)) {
    if (length(available_rho) > 0) {
      available_rho[1]
    } else {
      NA_real_
    }
  } else {
    rho_boxplot
  }

  line_rhos <- rho_median_lines[!is.na(rho_median_lines)]
  line_rhos <- intersect(line_rhos, setdiff(available_rho, rho_boxplot))

  boxplot_data <- loss_data %>%
    filter(
      spec_label != "Misspecified\n(Estimated)" |
        (spec_label == "Misspecified\n(Estimated)" & (is.na(rho_boxplot) | rho == rho_boxplot))
    )

  line_data <- loss_data %>%
    filter(spec_label == "Misspecified\n(Estimated)", rho %in% line_rhos) %>%
    group_by(spec_label, rho, n, gamma) %>%
    summarise(relative_improvement = median(relative_improvement, na.rm = TRUE), .groups = "drop")

  dodge_width <- 0.75
  n_levels <- sort(unique(loss_data$n))
  gamma_levels <- sort(unique(loss_data$gamma))

  boxplot_data <- boxplot_data %>%
    mutate(
      n_pos = match(n, n_levels)
    )

  line_data <- line_data %>%
    arrange(n, rho, gamma) %>%
    mutate(
      n_pos = match(n, n_levels),
      x_pos = n_pos + dodge_width * ((match(gamma, gamma_levels) - 0.5) / length(gamma_levels) - 0.5)
    )

  p <- ggplot(boxplot_data, aes(x = n_pos, y = relative_improvement)) +
    geom_boxplot(
      aes(fill = factor(gamma), group = interaction(n_pos, gamma)),
      outlier.shape = NA,
      position = position_dodge(width = dodge_width)
    ) +
    geom_line(
      data = line_data,
      aes(
        x = x_pos,
        group = interaction(n, rho),
        linetype = factor(rho)
      ),
      color = "black",
      size = 0.8,
      position = position_identity()
    ) +
    geom_point(
      data = line_data,
      aes(x = x_pos, color = factor(gamma), shape = factor(rho)),
      size = 1.5,
      position = position_identity()
    ) +
    facet_grid(. ~ spec_label) +
    scale_x_continuous(
      breaks = seq_along(n_levels),
      labels = n_levels
    ) +
    scale_y_continuous(breaks = pretty_breaks()) +
    coord_cartesian(ylim = y_limits) +
    labs(
      title = "Relative Loss Improvement",
      x = "Sample Size (n)",
      y = "Relative Improvement",
      fill = expression(gamma),
      color = expression(gamma),
      linetype = expression(rho),
      shape = expression(rho)
    ) +
    theme_minimal() +
    guides(
      fill = guide_legend(nrow = 1, byrow = TRUE, title.position = "left"),
      color = guide_legend(nrow = 1, byrow = TRUE, title.position = "left"),
      linetype = guide_legend(nrow = 1, byrow = TRUE, title.position = "left"),
      shape = guide_legend(nrow = 1, byrow = TRUE, title.position = "left")
    ) +
    theme(
      text = element_text(size = 11),
      strip.text = element_text(size = 10),
      legend.position = "bottom",
      legend.box = "horizontal"
    )
  
  return(p)
}


# ============================================================================
# CONVENIENCE FUNCTIONS
# ============================================================================

#' Generate All Three Main Plots
#' 
#' @param power_db Path to power database
#' @param estimation_db Path to estimation database
#' @param save_plots Logical, whether to save plots as PNG
#' @param theta_rho_boxplot Rho value to use for theta_hat boxplots in the
#'   misspecified case (NULL defaults to the first available)
#' @param theta_rho_median_lines Numeric vector of rho values to display as
#'   median-only lines for theta_hat in the misspecified case
#' @param loss_rho_boxplot Rho value to use for loss boxplots in the misspecified
#'   case (NULL defaults to the first available)
#' @param loss_rho_median_lines Numeric vector of rho values to display as
#'   median-only lines for loss in the misspecified case
#' @return List of ggplot objects
generate_all_plots <- function(power_db = "power_simulations.db",
                               estimation_db = "estimation_simulations.db",
                               save_plots = TRUE,
                               theta_rho_boxplot = NULL,
                               theta_rho_median_lines = numeric(0),
                               loss_rho_boxplot = NULL,
                               loss_rho_median_lines = numeric(0)) {
  
  message("Creating power plot...")
  p1 <- plot_power_comparison(power_db)
  
  message("Creating theta_hat distribution plot...")
  p2 <- plot_theta_hat_distribution(
    estimation_db,
    rho_boxplot = theta_rho_boxplot,
    rho_median_lines = theta_rho_median_lines
  )

  message("Creating relative loss improvement plot...")
  p3 <- plot_relative_loss_improvement(
    estimation_db,
    rho_boxplot = loss_rho_boxplot,
    rho_median_lines = loss_rho_median_lines
  )
  
  if (save_plots) {
    ggsave("plot_power_comparison.png", p1, width = 12, height = 8, dpi = 300)
    ggsave("plot_theta_hat_distribution.png", p2, width = 12, height = 6, dpi = 300)
    ggsave("plot_relative_loss_improvement.png", p3, width = 12, height = 8, dpi = 300)
    message("✓ Plots saved as PNG files")
  }
  
  return(list(
    power = p1,
    theta_hat = p2,
    loss = p3
  ))
}


# ============================================================================
# USAGE EXAMPLE
# ============================================================================

if (FALSE) {
  # Create power view first
  create_power_view("power_simulations.db")
  
  # Generate all plots
  plots <- generate_all_plots(
    power_db = "power_simulations.db",
    estimation_db = "estimation_simulations.db",
    save_plots = TRUE
  )
  
  # View individual plots
  print(plots$power)
  print(plots$theta_hat)
  print(plots$loss)
}