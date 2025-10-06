# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(tidyr)
  library(cowplot)
  library(R.matlab)
  library(fGarch)
})

set.seed(423) # Set seed for reproducibility

# ==============================================================================
# Configuration
# ==============================================================================

# File paths
OUTPUT_DIR <- "/data/p_02865/data/ds003838/output/hep_anxiety_sim/"
PLOT_OUTPUT_DIR <- file.path(OUTPUT_DIR, "plots")
CHAN_NAMES_DIR <- "/data/p_02865/code/supplement/orig_chan_labels.mat"
BASE_DIR <- "/data/p_02865/data/ds003838/derivatives/"
PCA_DATA_PATH <- file.path("/data/p_02865/code/pwr_simulations", "sim_res_reg.csv")

# Analysis parameters
DEFAULT_SAMPLE_SIZE <- 3000
DEFAULT_RHO <- 0.3
ROI_ELECTRODES <- c("Cz", "C3", "C4", "Fz", "F3", "F4", "FC1", "FC2")
TIME_WINDOW <- c(200, 500)

# Create output directory if it doesn't exist
if (!dir.exists(PLOT_OUTPUT_DIR)) {
  dir.create(PLOT_OUTPUT_DIR, recursive = TRUE)
}

# Generate HEP file paths
generate_hep_files <- function(base_dir) {
  expand.grid(
    hp = c("0.1", "0.5"),
    bl = c("[-150 -50]", "[]"),
    ic = c("0.0", "1.5", "3.5", "NaN"),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      key = paste0("hp_", hp, "_bl_", gsub("\\[|\\]|-", "", bl), "_ic_", ic),
      key = gsub("150 50", "m150_m50", key),
      key = gsub(" ", "_", key),
      path = file.path(base_dir, paste0("hp_", hp, "_bl_", bl, "_ic_", ic, "_all_HEP.mat"))
    ) %>%
    select(key, path) %>%
    deframe()
}

HEP_FILES <- generate_hep_files(BASE_DIR)

# ==============================================================================
# Data Loading Functions
# ==============================================================================

#' Load and prepare simulation results
load_simulation_data <- function() {
  cat("Loading simulation results...\n")
  
  # Load correlation simulation results
  rho_continuous_path <- file.path(OUTPUT_DIR, "sim_results_rho_continuous.RData")
  rho_grouped_path <- file.path(OUTPUT_DIR, "sim_results_rho_grouped.RData")
  
  if (!file.exists(rho_continuous_path) || !file.exists(rho_grouped_path)) {
    stop("Simulation result files not found. Please run simulations first.")
  }
  
  load(rho_continuous_path, envir = .GlobalEnv)
  load(rho_grouped_path, envir = .GlobalEnv)
  
  # Extract and combine power data
  power_data_continuous <- sim_results_rho_continuous$power_data %>% 
    filter(sample_size == DEFAULT_SAMPLE_SIZE) %>%
    mutate(analysis_type = "Continuous")
    
  power_data_grouped <- sim_results_rho_grouped$power_data %>% 
    filter(sample_size == DEFAULT_SAMPLE_SIZE) %>%
    mutate(analysis_type = "Grouped")
  
  return(bind_rows(power_data_continuous, power_data_grouped))
}

#' Load PCA simulation results
load_pca_data <- function() {
  cat("Loading PCA simulation results...\n")
  
  if (!file.exists(PCA_DATA_PATH)) {
    warning("PCA data file not found: ", PCA_DATA_PATH)
    return(NULL)
  }
  
  return(read.csv(PCA_DATA_PATH))
}

#' Extract channel names from MATLAB file
extract_chan_names <- function(chan_names) {
  labels <- sapply(chan_names$orig.chan.labels, function(x) x[[1]][1, 1])
  return(labels)
}

#' Load HEP data from MATLAB file
load_hep_data <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("HEP file not found: ", file_path)
  }
  return(readMat(file_path))
}

# ==============================================================================
# Analysis Functions
# ==============================================================================

#' Generate anxiety values correlated with HEP
#' @param hep_values Vector of HEP amplitude values
#' @param rho Correlation coefficient (-1 < rho < 1)
#' @return Vector of simulated GAD-7 scores
generate_anxiety_corr <- function(hep_values, rho) {
  if (!is.finite(rho) || abs(rho) >= 1) {
    stop("rho must be finite and strictly between -1 and 1.")
  }
  
  n <- length(hep_values)
  hep_z <- as.numeric(scale(hep_values))
  z_noise <- rnorm(n, 0, 1)
  latent <- rho * hep_z + sqrt(1 - rho^2) * z_noise
  
  # Apply skewed normal transformation
  u <- pnorm(latent)
  anxiety_cont <- qsnorm(u, mean = 5, sd = 4, xi = 3)
  
  # Clamp to valid GAD-7 range [0, 21]
  anxiety <- pmax(0, pmin(21, anxiety_cont))
  
  return(anxiety)
}

#' Extract ROI-averaged HEP values
#' @param hep_data Loaded HEP data structure
#' @param electrodes Vector of electrode names to include in ROI
#' @param time_window Time window in ms [start, end]
#' @param time_vector Time vector for the data
#' @param channel_names Vector of all channel names
#' @return Vector of ROI-averaged HEP amplitudes
extract_roi_hep <- function(hep_data, 
                            electrodes = ROI_ELECTRODES,
                            time_window = TIME_WINDOW, 
                            time_vector = NULL, 
                            channel_names = NULL) {
  
  # Validate inputs
  if (is.null(hep_data) || is.null(time_vector) || is.null(channel_names)) {
    warning("Missing required inputs in extract_roi_hep")
    return(NULL)
  }
  
  hep_values <- hep_data$HEP.all
  
  if (is.null(dim(hep_values)) || length(dim(hep_values)) < 3) {
    warning("HEP data is not 3D")
    return(NULL)
  }
  
  # Find electrode and time indices
  electrode_indices <- match(electrodes, channel_names)
  time_indices <- which(time_vector >= time_window[1] & time_vector <= time_window[2])
  
  if (anyNA(electrode_indices) || length(time_indices) == 0) {
    warning("Invalid electrode or time indices")
    return(NULL)
  }
  
  if (any(electrode_indices > dim(hep_values)[1]) || any(time_indices > dim(hep_values)[2])) {
    warning("Index out of bounds")
    return(NULL)
  }
  
  # Extract and average ROI data
  roi_hep <- hep_values[electrode_indices, time_indices, , drop = FALSE]
  avg_hep <- apply(roi_hep, 3, function(x) mean(x, na.rm = TRUE))
  
  if (all(is.na(avg_hep))) {
    warning("All values NA after averaging")
    return(NULL)
  }
  
  return(avg_hep)
}

# ==============================================================================
# Plotting Functions
# ==============================================================================

#' Custom theme for publication-ready plots
theme_power <- function(base_size = 8) { 
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", size = 0.25),
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0.5), 
      plot.subtitle = element_text(size = base_size - 1, hjust = 0.5, color = "grey40"),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 1),  
      legend.title = element_text(size = base_size, face = "bold"),  
      legend.text = element_text(size = base_size - 1),  
      legend.position = "bottom",
      strip.text = element_text(size = base_size, face = "bold"),
      strip.background = element_rect(fill = "grey95", color = "grey80")
    )
}

#' Create power curve plot by effect size
#' @param data Data frame with power analysis results
#' @param title_suffix Optional suffix for plot title
#' @return ggplot object
plot_power_by_effect <- function(data, title_suffix = "") {
  ggplot(data, aes(x = mean_effect, y = power, color = analysis_type, linetype = analysis_type)) +
    geom_line(size = 0.75) +
    geom_point(size = 1.5, alpha = 0.8) +
    geom_hline(yintercept = 0.95, linetype = "dashed", alpha = 0.6, color = "grey50", size = 0.5) +
    scale_y_continuous(
      limits = c(0, 1),
      labels = percent_format(accuracy = 1),
      breaks = seq(0, 1, 0.2)
    ) +
    scale_x_continuous(
      breaks = pretty,
      labels = function(x) sprintf("%.3f", x)
    ) +
    scale_color_manual(
      name = "Analysis Type",
      values = c("Continuous" = "#666666", "Grouped" = "#666666"),
      labels = c("Continuous" = "Continuous Anxiety", "Grouped" = "Low vs High Groups")
    ) +
    scale_linetype_manual(
      name = "Analysis Type", 
      values = c("Continuous" = "solid", "Grouped" = "solid"),
      labels = c("Continuous" = "Continuous Anxiety", "Grouped" = "Low vs High Groups")
    ) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = "Mean Effect Size (regression coefficient)",
      y = "Statistical Power",
      caption = NULL
    ) +
    annotate("text", x = max(data$mean_effect) * 0.02, y = 0.95, 
             label = "95%", 
             hjust = 0, vjust = -0.3, size = 2.5, color = "grey50") +
    theme_power()
}

#' Create histogram of simulated anxiety scores
#' @param n_subjects Number of subjects to simulate
#' @param rho Correlation between HEP and anxiety
#' @return ggplot object
plot_anxiety_histogram <- function(n_subjects = DEFAULT_SAMPLE_SIZE, rho = DEFAULT_RHO) {
  # Load channel names
  chan_names <- readMat(CHAN_NAMES_DIR)
  channel_names <- extract_chan_names(chan_names)
  
  # Load real HEP data for realistic amplitude distributions
  first_file <- HEP_FILES[[1]]
  cat("Loading real HEP data from:", basename(first_file), "\n")
    
  hep_data <- load_hep_data(first_file)
  time_vector <- seq(-200, 800, length.out = dim(hep_data$HEP.all)[2])
  real_hep_amplitudes <- extract_roi_hep(hep_data, 
                                         time_vector = time_vector, 
                                         channel_names = channel_names)
  
  # Sample with replacement to get desired number of subjects
  hep_values <- sample(real_hep_amplitudes, n_subjects, replace = TRUE)
  
  # Generate correlated anxiety values
  anxiety_values <- generate_anxiety_corr(hep_values, rho)
  
  # Create histogram
  ggplot(data.frame(anxiety = anxiety_values), aes(x = anxiety)) +
    geom_histogram(bins = 22, fill = "#666666", alpha = 0.7, color = "white", size = 0.2) +
    scale_x_continuous(
      breaks = seq(0, 21, 3),
      limits = c(0, 21)
    ) +
    scale_y_continuous(labels = comma_format()) +
    labs(
      title = "C. Simulated GAD-7 Scores",
      x = "GAD-7 Score",
      y = "Frequency"
    ) +
    theme_power() +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
}

#' Create individual PCA control plots
#' @param pca_data PCA simulation results data frame
#' @param control_type Type of control analysis ("no", "all", "pca")
#' @return ggplot object
plot_pca_individual <- function(pca_data, control_type) {
  if (is.null(pca_data)) {
    warning("PCA data is NULL")
    return(ggplot() + theme_void())
  }
  
  # Process PCA data for specific control type
  pca_power_data <- pca_data %>%
    filter(controls == control_type) %>%
    group_by(n_sub, sim_value) %>%
    summarise(
      power = mean(sig, na.rm = TRUE),
      mean_coef = mean(result_coef, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(n_sub_factor = factor(n_sub)) %>%
    filter(!is.na(power), !is.infinite(power), power >= 0, power <= 1)
    
  # Set title based on control type
  title_text <- case_when(
    control_type == "no" ~ "D. No Controls",
    control_type == "all" ~ "E. All Controls", 
    control_type == "pca" ~ "F. PCA Controls",
    TRUE ~ paste("PCA:", control_type)
  )
  
  if (nrow(pca_power_data) == 0) {
    warning("No data available for control type: ", control_type)
    return(ggplot() + theme_void())
  }
  
  # Create the plot
  ggplot(pca_power_data, aes(x = sim_value, y = power, color = n_sub_factor, group = n_sub_factor)) +
    geom_line(size = 0.75) +
    geom_point(size = 1.2, alpha = 0.7) +
    geom_hline(yintercept = 0.95, linetype = "dashed", alpha = 0.6, color = "grey50", size = 0.5) +
    scale_y_continuous(
      limits = c(0, 1),
      labels = percent_format(accuracy = 1),
      breaks = seq(0, 1, 0.2)
    ) +
    scale_x_continuous(
      limits = c(0.1, NA),
      breaks = pretty,
      labels = function(x) sprintf("%.2f", x)
    ) +
    scale_color_grey(
      name = "Subjects",
      start = 0.2, end = 0.8
    ) +
    labs(
      title = title_text,
      x = "Simulated association value",
      y = "Statistical Power"
    ) +
    annotate("text", x = max(pca_power_data$sim_value) * 0.02, y = 0.95, 
             label = "95%", hjust = 0, vjust = -0.3, size = 2.5, color = "grey50") +
    theme_power() +
    theme(
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      legend.position = "bottom",
      legend.key.size = unit(0.3, "cm"),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 7)
    )
}

# ==============================================================================
# Main Execution
# ==============================================================================

main <- function() {
  cat("Starting power curve plotting...\n")
  
  # Load data
  all_power_data <- load_simulation_data()
  pca_results <- load_pca_data()
  
  # Separate data by analysis type
  continuous_data <- filter(all_power_data, analysis_type == "Continuous")
  grouped_data <- filter(all_power_data, analysis_type == "Grouped")
  
  # Generate individual plots
  cat("Generating power curve plots...\n")
  
  # 1. Power by effect size plots
  p1a <- plot_power_by_effect(continuous_data)
  p1b <- plot_power_by_effect(grouped_data)
  
  # 2. Anxiety histogram
  p_hist <- plot_anxiety_histogram()
  
  # 3. PCA control plots
  p_pca_no <- plot_pca_individual(pca_results, "no")
  p_pca_all <- plot_pca_individual(pca_results, "all")
  p_pca_pca <- plot_pca_individual(pca_results, "pca")
  
  # Save individual plots
  cat("Saving individual plots...\n")
  plot_specs <- list(width = 10, height = 6, units = "in", dpi = 300)
  
  ggsave(file.path(PLOT_OUTPUT_DIR, "power_curves_continuous_by_effect.svg"), 
         p1a, width = plot_specs$width, height = plot_specs$height, 
         units = plot_specs$units, dpi = plot_specs$dpi)
         
  ggsave(file.path(PLOT_OUTPUT_DIR, "power_curves_grouped_by_effect.svg"), 
         p1b, width = plot_specs$width, height = plot_specs$height, 
         units = plot_specs$units, dpi = plot_specs$dpi)
         
  ggsave(file.path(PLOT_OUTPUT_DIR, "power_curves_pca_no_controls.svg"), 
         p_pca_no, width = plot_specs$width, height = plot_specs$height, 
         units = plot_specs$units, dpi = plot_specs$dpi)
         
  ggsave(file.path(PLOT_OUTPUT_DIR, "power_curves_pca_all_controls.svg"), 
         p_pca_all, width = plot_specs$width, height = plot_specs$height, 
         units = plot_specs$units, dpi = plot_specs$dpi)
         
  ggsave(file.path(PLOT_OUTPUT_DIR, "power_curves_pca_controls.svg"), 
         p_pca_pca, width = plot_specs$width, height = plot_specs$height, 
         units = plot_specs$units, dpi = plot_specs$dpi)
  
  # Create combined plot
  cat("Creating combined plot...\n")
  combined_power <- plot_grid(
    p_hist + theme(legend.position = "none") + 
      labs(title = "A. Simulated GAD-7 Scores") + 
      theme(plot.title = element_text(hjust = 0)),
    p1a + theme(legend.position = "none") + 
      labs(title = "B. Continuous Analysis") + 
      theme(plot.title = element_text(hjust = 0)),
    p1b + theme(legend.position = "none") + 
      labs(title = "C. Grouped Analysis") + 
      theme(plot.title = element_text(hjust = 0)),
    p_pca_no + theme(legend.position = "none") + 
      theme(plot.title = element_text(hjust = 0)),
    p_pca_all + theme(legend.position = "none") + 
      theme(plot.title = element_text(hjust = 0)),
    p_pca_pca + theme(legend.position = "none") + 
      theme(plot.title = element_text(hjust = 0)),
    ncol = 3, nrow = 2,
    align = "hv",
    rel_widths = c(1, 1, 1),
    rel_heights = c(1, 1)
  )
  
  # Save combined plots
  combined_specs <- list(width = 270, height = 120, units = "mm", dpi = 300)
  
  ggsave(file.path(PLOT_OUTPUT_DIR, "power_curves_combined.svg"), 
         combined_power, width = combined_specs$width, height = combined_specs$height, 
         units = combined_specs$units, dpi = combined_specs$dpi)
         
  ggsave(file.path(PLOT_OUTPUT_DIR, "power_curves_combined.png"), 
         combined_power, width = combined_specs$width, height = combined_specs$height, 
         units = combined_specs$units, dpi = 600, device = "png", bg = "white", type = "cairo")
  
  cat("Power curve plotting completed successfully!\n")
  cat("Plots saved to:", PLOT_OUTPUT_DIR, "\n")
}

# Run the main function
if (!interactive()) {
  main()
}