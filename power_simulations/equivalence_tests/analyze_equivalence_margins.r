# SIMPLIFIED Analysis: Find largest effect with <80% power and use EMPIRICAL coefficient as equivalence bound
# This approach is more precise because it uses the actual observed coefficients from simulation

# Load required libraries
library(dplyr)
library(ggplot2)


find_empirical_equivalence_bound <- function(power_data,
                                             target_sample_size = 100,
                                             power_threshold = 0.95) {
  sample_data <- power_data %>%
    filter(sample_size == target_sample_size) %>%
    arrange(strength)
  if (nrow(sample_data) == 0) stop("No rows for target sample size ", target_sample_size)
  adequate <- sample_data %>% filter(power >= power_threshold)
  if (nrow(adequate) == 0) stop("No effects reach the power threshold ", power_threshold)
  boundary_strength <- min(adequate$strength)
  effect_row <- sample_data %>% filter(strength == boundary_strength) %>% slice(1)
  empirical_bound <- abs(effect_row$mean_effect)
  mode <- "smallest_detectable"
  list(
    theoretical_boundary_strength = boundary_strength,
    empirical_equivalence_bound = empirical_bound,
    boundary_power = effect_row$power,
    sample_size = target_sample_size,
    empirical_mean_effect = effect_row$mean_effect,
    mode = mode
  )
}

analyze_equivalence_from_object <- function(sim_object,
                                            object_label,
                                            target_sample_size = 100,
                                            power_threshold = 0.95,
                                            output_dir) {
  power_data <- sim_object$power_data

  # Determine variants to analyze separately (continuous vs grouped)
  variants <- list()
  if ("anx_group" %in% names(power_data)) {
    vals <- sort(unique(power_data$anx_group))
    for (v in vals) {
      vlab <- if (isTRUE(v)) "grouped" else if (identical(v, FALSE)) "continuous" else paste0("anx_group=", v)
      variants[[vlab]] <- power_data %>% dplyr::filter(anx_group == v)
    }
  } else {
    variants[["unknown"]] <- power_data
  }
  
  out <- list()
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  cat("\n--- Analyzing:", object_label, " (", length(variants), " variant(s)) ---\n", sep="")
  for (vname in names(variants)) {
    pdv <- variants[[vname]]
    if (nrow(pdv) == 0) next
    cat("Variant:", vname, "\n")
  res <- find_empirical_equivalence_bound(pdv, target_sample_size, power_threshold)
    if (is.null(res)) next
  plot_df <- pdv %>% dplyr::filter(sample_size == res$sample_size)
    p <- ggplot(plot_df, aes(strength, power)) +
      geom_line(color="steelblue", linewidth=1) +
      geom_point(color="steelblue", size=2) +
      geom_hline(yintercept = power_threshold, linetype="dashed", color="red") +
      geom_vline(xintercept = res$theoretical_boundary_strength, linetype="dotted", color="darkgreen") +
  annotate("text", x = res$theoretical_boundary_strength, y = 0.05 + min(1, power_threshold),
     label = paste0("SESOI ≈ ", sprintf("%.4f", res$empirical_equivalence_bound), " (coef)"),
         color="darkgreen", angle=90, vjust=-0.4, size=3) +
      labs(title = paste0("Power Curve (", object_label, " – ", vname, ") N=", res$sample_size),
       subtitle = paste0("Smallest effect with ≥ ", power_threshold*100, "% power"),
           x = "Effect Strength (rho / strength grid)", y = "Power") +
      scale_y_continuous(limits=c(0,1)) + theme_minimal()
  plot_file <- file.path(output_dir, sprintf("power_curve_empirical_BOUND_%s_%s_n%d.png", object_label, vname, res$sample_size))
    ggsave(plot_file, p, width=8, height=5, dpi=300)
    # Capture grid_type if available to clarify that strength corresponds to rho in correlation sims
    grid_type_here <- if ("grid_type" %in% names(pdv)) unique(as.character(na.omit(pdv$grid_type))) else NA_character_
    if (length(grid_type_here) != 1) grid_type_here <- paste(na.omit(grid_type_here), collapse=",")

    out[[length(out)+1]] <- list(
      label = object_label,
      variant = vname,
      grid_type = grid_type_here,
      target_sample_size = target_sample_size,
      power_threshold = power_threshold,
      theoretical_boundary_strength = res$theoretical_boundary_strength,
      empirical_equivalence_bound = res$empirical_equivalence_bound,
      bound_unit = "coefficient",
      boundary_power = res$boundary_power,
      empirical_mean_effect = res$empirical_mean_effect,
      mode = res$mode,
      plot_file = plot_file,
      timestamp = Sys.time()
    )
  }
  if (length(out) == 0) return(NULL)
  out
}

analyze_equivalence_files <- function(result_files,
                                      target_sample_size = 100,
                                      power_threshold = 0.95,
                                      output_dir = file.path(getwd(), "pwr_simulations_output")) {
  summaries <- list()
  for (f in result_files) {
    cat("\nLoading:", f, "\n")
  loaded_names <- load(f)
    o <- loaded_names[1]
    obj <- get(o, inherits = TRUE)
    label <- gsub("^sim_results_", "", o)
    s_list <- analyze_equivalence_from_object(obj, label, target_sample_size, power_threshold, output_dir)
    if (!is.null(s_list)) {
      for (s in s_list) {
        save_file <- file.path(output_dir, sprintf("equivalence_bound_%s_%s_n%d.RData", label, s$variant, target_sample_size))
        save(s, file = save_file)
        s$save_file <- save_file
        s$source_file <- f
        summaries[[length(summaries)+1]] <- s
      }
    }
  }
  if (length(summaries) == 0) {
    stop("No suitable simulation result objects found in supplied files.")
  }
  summary_df <- do.call(rbind, lapply(summaries, function(x) {
    df <- as.data.frame(x)
    # keep columns order readable
  df[, c("label","variant","grid_type","target_sample_size","power_threshold",
           "theoretical_boundary_strength","empirical_equivalence_bound","bound_unit",
           "boundary_power","empirical_mean_effect","mode","plot_file","source_file","timestamp")]
  }))
  csv_file <- file.path(output_dir, sprintf("equivalence_empirical_summary_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S")))
  write.csv(summary_df, csv_file, row.names = FALSE)
  cat("\nSummary written to:", csv_file, "\n")
  invisible(list(summary_table = summary_df, files = result_files, target_sample_size = target_sample_size, power_threshold = power_threshold, csv = csv_file))
}

## ================= MAIN EXECUTION (edit parameters as needed) =================
cat("\n=== Empirical Equivalence Bound Extraction (sim_roi_cor outputs) ===\n")
target_sample_size <- 100       
power_threshold <- 0.80         
output_dir <- "/data/p_02865/data/ds003838/output/hep_anxiety_sim_tmp/"

# Correlation-based simulation result files produced by sim_roi_cor_and_group.R
result_files <- c(
  file.path(output_dir, "sim_results_continuous.RData"),
  file.path(output_dir, "sim_results_grouped.RData")
)

analysis_summary <- analyze_equivalence_files(result_files,
                                              target_sample_size = target_sample_size,
                                              power_threshold = power_threshold,
                                              output_dir = output_dir)

cat("\nCompleted empirical equivalence bound extraction.\n")
print(analysis_summary$summary_table)