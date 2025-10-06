# Equivalence Margin Analysis for HEP-Anxiety Power Simulations
# Finds the smallest detectable effect (SESOI) based on power threshold
# Uses empirical coefficients from simulations as equivalence bounds

library(dplyr)
library(ggplot2)

# ---------- Core Functions ----------

find_equivalence_bound <- function(power_data, target_sample_size, power_threshold) {
  # Filter to target sample size and sort by effect strength
  sample_data <- power_data %>%
    filter(sample_size == target_sample_size) %>%
    arrange(strength)
  
  if (nrow(sample_data) == 0) {
    stop("No data for sample size ", target_sample_size)
  }
  
  # Find smallest effect that reaches power threshold
  adequate <- sample_data %>% filter(power >= power_threshold)
  
  if (nrow(adequate) == 0) {
    stop("No effects reach power threshold ", power_threshold)
  }
  
  boundary_strength <- min(adequate$strength)
  effect_row <- sample_data %>% filter(strength == boundary_strength) %>% slice(1)
  
  list(
    boundary_strength = boundary_strength,
    equivalence_bound = abs(effect_row$mean_effect),
    boundary_power = effect_row$power,
    mean_effect = effect_row$mean_effect,
    sample_size = target_sample_size
  )
}

create_power_plot <- function(plot_data, result, label, variant, power_threshold) {
  ggplot(plot_data, aes(strength, power)) +
    geom_line(color = "steelblue", linewidth = 1) +
    geom_point(color = "steelblue", size = 2) +
    geom_hline(yintercept = power_threshold, linetype = "dashed", color = "red") +
    geom_vline(xintercept = result$boundary_strength, linetype = "dotted", color = "darkgreen") +
    annotate("text", 
             x = result$boundary_strength, 
             y = min(1, power_threshold) + 0.05,
             label = sprintf("SESOI ≈ %.4f", result$equivalence_bound),
             color = "darkgreen", angle = 90, vjust = -0.4, size = 3.5) +
    labs(
      title = sprintf("Power Curve: %s (%s), N=%d", label, variant, result$sample_size),
      subtitle = sprintf("Smallest effect with ≥%.0f%% power", power_threshold * 100),
      x = "Effect Strength (ρ)",
      y = "Power"
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"))
}

analyze_variant <- function(power_data, label, variant, target_sample_size, 
                           power_threshold, output_dir) {
  cat(sprintf("  Analyzing variant: %s\n", variant))
  
  # Find equivalence bound
  result <- find_equivalence_bound(power_data, target_sample_size, power_threshold)
  
  # Create and save plot
  plot_data <- power_data %>% filter(sample_size == target_sample_size)
  p <- create_power_plot(plot_data, result, label, variant, power_threshold)
  
  plot_file <- file.path(
    output_dir, 
    sprintf("power_curve_%s_%s_n%d.png", label, variant, target_sample_size)
  )
  ggsave(plot_file, p, width = 8, height = 5, dpi = 300)
  
  # Extract grid type if available
  grid_type <- if ("grid_type" %in% names(power_data)) {
    unique(power_data$grid_type)[1]
  } else {
    NA_character_
  }
  
  # Return summary
  list(
    label = label,
    variant = variant,
    grid_type = grid_type,
    sample_size = target_sample_size,
    power_threshold = power_threshold,
    boundary_strength = result$boundary_strength,
    equivalence_bound = result$equivalence_bound,
    boundary_power = result$boundary_power,
    mean_effect = result$mean_effect,
    plot_file = plot_file,
    timestamp = Sys.time()
  )
}

analyze_simulation_object <- function(sim_object, label, target_sample_size, 
                                     power_threshold, output_dir) {
  cat(sprintf("\n--- Analyzing: %s ---\n", label))
  
  power_data <- sim_object$power_data
  
  # Split by anxiety grouping if present
  if ("anx_group" %in% names(power_data)) {
    variants <- list(
      continuous = power_data %>% filter(anx_group == FALSE),
      grouped = power_data %>% filter(anx_group == TRUE)
    )
  } else {
    variants <- list(combined = power_data)
  }
  
  # Analyze each variant
  results <- list()
  for (vname in names(variants)) {
    vdata <- variants[[vname]]
    if (nrow(vdata) == 0) next
    
    result <- analyze_variant(vdata, label, vname, target_sample_size, 
                             power_threshold, output_dir)
    
    # Save individual result
    save_file <- file.path(
      output_dir,
      sprintf("equivalence_bound_%s_%s_n%d.RData", label, vname, target_sample_size)
    )
    save(result, file = save_file)
    
    result$save_file <- save_file
    results[[length(results) + 1]] <- result
  }
  
  results
}

# ---------- Main Analysis Function ----------

analyze_equivalence_margins <- function(result_files, target_sample_size = 100,
                                       power_threshold = 0.80, output_dir) {
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("\n=== Equivalence Margin Analysis ===\n")
  cat(sprintf("Target sample size: %d\n", target_sample_size))
  cat(sprintf("Power threshold: %.2f\n", power_threshold))
  cat(sprintf("Output directory: %s\n\n", output_dir))
  
  # Analyze each file
  all_results <- list()
  for (f in result_files) {
    if (!file.exists(f)) {
      warning("File not found: ", f)
      next
    }
    
    cat(sprintf("Loading: %s\n", basename(f)))
    
    # Load simulation results
    loaded_names <- load(f)
    obj_name <- loaded_names[1]
    obj <- get(obj_name)
    
    # Extract label from object name
    label <- gsub("^sim_results_(rho_)?", "", obj_name)
    
    # Analyze object
    results <- analyze_simulation_object(obj, label, target_sample_size, 
                                        power_threshold, output_dir)
    
    # Add source file info
    for (r in results) {
      r$source_file <- f
      all_results[[length(all_results) + 1]] <- r
    }
  }
  
  if (length(all_results) == 0) {
    stop("No valid results found in provided files")
  }
  
  # Create summary table
  summary_df <- do.call(rbind, lapply(all_results, as.data.frame))
  
  # Save summary
  csv_file <- file.path(
    output_dir,
    sprintf("equivalence_summary_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S"))
  )
  write.csv(summary_df, csv_file, row.names = FALSE)
  
  cat("\n=== Analysis Complete ===\n")
  cat(sprintf("Summary saved to: %s\n\n", csv_file))
  
  print(summary_df[, c("label", "variant", "sample_size", "boundary_strength", 
                       "equivalence_bound", "boundary_power")])
  
  invisible(list(
    summary = summary_df,
    results = all_results,
    csv_file = csv_file
  ))
}

# ---------- Run Analysis ----------

cat("\n=== HEP-Anxiety Equivalence Margin Analysis ===\n")

# Configuration
target_sample_size <- 300
power_threshold <- 0.80
output_dir <- "/data/p_02865/data/ds003838/output/hep_anxiety_sim_test"

# Input files (from sim_roi_cor_and_group.R with legacy naming)
result_files <- c(
  file.path(output_dir, "sim_results_rho_continuous.RData"),
  file.path(output_dir, "sim_results_rho_grouped.RData")
)

# Run analysis
analysis_results <- analyze_equivalence_margins(
  result_files = result_files,
  target_sample_size = target_sample_size,
  power_threshold = power_threshold,
  output_dir = output_dir
)

cat("\nEquivalence bounds determined!\n")
