# HEP-Anxiety Power Simulations (Correlation-Based Only)
# Streamlined version - no fallbacks, correlation approach only

library(R.matlab)
library(lme4)
library(lmerTest)
library(dplyr)
library(parallel)
library(fGarch)

set.seed(423)

# ---------- Configuration ----------

chan_names_dir <- "/data/p_02865/code/supplement/orig_chan_labels.mat"
base_dir <- "/data/p_02865/data/ds003838/derivatives/"

hep_files <- list(
  "hp_0.1_bl_m150_m50_ic_0.0" = paste0(base_dir, "hp_0.1_bl_[-150 -50]_ic_0.0_all_HEP.mat"),
  "hp_0.1_bl_m150_m50_ic_1.5" = paste0(base_dir, "hp_0.1_bl_[-150 -50]_ic_1.5_all_HEP.mat"),
  "hp_0.1_bl_m150_m50_ic_3.5" = paste0(base_dir, "hp_0.1_bl_[-150 -50]_ic_3.5_all_HEP.mat"),
  "hp_0.1_bl_m150_m50_ic_NaN" = paste0(base_dir, "hp_0.1_bl_[-150 -50]_ic_NaN_all_HEP.mat"),
  "hp_0.1_bl_none_ic_0.0" = paste0(base_dir, "hp_0.1_bl_[]_ic_0.0_all_HEP.mat"),
  "hp_0.1_bl_none_ic_1.5" = paste0(base_dir, "hp_0.1_bl_[]_ic_1.5_all_HEP.mat"),
  "hp_0.1_bl_none_ic_3.5" = paste0(base_dir, "hp_0.1_bl_[]_ic_3.5_all_HEP.mat"),
  "hp_0.1_bl_none_ic_NaN" = paste0(base_dir, "hp_0.1_bl_[]_ic_NaN_all_HEP.mat"),
  "hp_0.5_bl_m150_m50_ic_0.0" = paste0(base_dir, "hp_0.5_bl_[-150 -50]_ic_0.0_all_HEP.mat"),
  "hp_0.5_bl_m150_m50_ic_1.5" = paste0(base_dir, "hp_0.5_bl_[-150 -50]_ic_1.5_all_HEP.mat"),
  "hp_0.5_bl_m150_m50_ic_3.5" = paste0(base_dir, "hp_0.5_bl_[-150 -50]_ic_3.5_all_HEP.mat"),
  "hp_0.5_bl_m150_m50_ic_NaN" = paste0(base_dir, "hp_0.5_bl_[-150 -50]_ic_NaN_all_HEP.mat"),
  "hp_0.5_bl_none_ic_0.0" = paste0(base_dir, "hp_0.5_bl_[]_ic_0.0_all_HEP.mat"),
  "hp_0.5_bl_none_ic_1.5" = paste0(base_dir, "hp_0.5_bl_[]_ic_1.5_all_HEP.mat"),
  "hp_0.5_bl_none_ic_3.5" = paste0(base_dir, "hp_0.5_bl_[]_ic_3.5_all_HEP.mat"),
  "hp_0.5_bl_none_ic_NaN" = paste0(base_dir, "hp_0.5_bl_[]_ic_NaN_all_HEP.mat")
)

# ---------- Helper Functions ----------

extract_chan_names <- function(chan_names) {
  labels <- sapply(chan_names$orig.chan.labels, function(x) x[[1]][1, 1])
  return(labels)
}

extract_roi_hep <- function(hep_data, electrodes = c("Cz","C3","C4","Fz","F3","F4","FC1","FC2"),
                            time_window = c(200, 500), time_vector, channel_names) {
  hep_values <- hep_data$HEP.all
  electrode_indices <- match(electrodes, channel_names)
  time_indices <- which(time_vector >= time_window[1] & time_vector <= time_window[2])
  
  roi_hep <- hep_values[electrode_indices, time_indices, , drop = FALSE]
  avg_hep <- apply(roi_hep, 3, function(x) mean(x, na.rm = TRUE))
  
  return(avg_hep)
}

# ---------- Anxiety Generation (Correlation-Based) ----------

generate_anxiety <- function(hep_values, rho) {
  if (!is.finite(rho) || abs(rho) >= 1) {
    stop("rho must be finite and strictly between -1 and 1")
  }
  
  n <- length(hep_values)
  hep_z <- as.numeric(scale(hep_values))
  z_noise <- rnorm(n, 0, 1)
  latent <- rho * hep_z + sqrt(1 - rho^2) * z_noise
  
  # Skew transform to GAD-7 distribution
  u <- pnorm(latent)
  anxiety_cont <- qsnorm(u, mean = 5, sd = 4, xi = 3)
  
  # Clamp to valid GAD-7 range [0, 21]
  anxiety <- pmax(0, pmin(21, anxiety_cont))
  
  return(anxiety)
}

# ---------- Model Fitting & Testing ----------

run_mixed_model <- function(data) {
  data$multiverse <- factor(data$multiverse)
  contrasts(data$multiverse) <- contr.sum(nlevels(data$multiverse))
  
  if (is.factor(data$anxiety)) {
    data$anxiety <- factor(data$anxiety, levels = c("low", "high"))
    contrasts(data$anxiety) <- contr.treatment(2, base = 1)
  }
  
  model <- lmer(hep_amp ~ anxiety + multiverse + anxiety:multiverse + (1 | subject), data = data)
  list(model = model, summary = summary(model))
}

find_anxiety_row <- function(ms) {
  rows <- rownames(ms$coefficients)
  candidates <- rows[grepl("^anxiety", rows) & !grepl(":", rows)]
  return(candidates[1])
}

run_equivalence_test <- function(model, equivalence_margin) {
  low_bound <- -equivalence_margin
  high_bound <- equivalence_margin
  
  ms <- summary(model)
  rowname <- find_anxiety_row(ms)
  coef_anx <- ms$coefficients[rowname, "Estimate"]
  se_anx <- ms$coefficients[rowname, "Std. Error"]
  df_anx <- ms$coefficients[rowname, "df"]
  
  t_upper <- (coef_anx - high_bound) / se_anx
  p_upper <- pt(t_upper, df = df_anx, lower.tail = TRUE)
  t_lower <- (coef_anx - low_bound) / se_anx
  p_lower <- pt(t_lower, df = df_anx, lower.tail = FALSE)
  equivalent <- p_upper < 0.05 && p_lower < 0.05
  
  list(
    equivalent = equivalent,
    p1 = p_lower,
    p2 = p_upper,
    equivalence_p = max(p_upper, p_lower),
    margin = equivalence_margin,
    coefficient = coef_anx,
    coef_row = rowname
  )
}

# ---------- Data Pre-loading ----------

preload_hep_data <- function(hep_files, chan_names_data) {
  cat("Pre-loading HEP data...\n")
  channel_names <- extract_chan_names(chan_names_data)
  
  hep_data_cache <- list()
  for (i in seq_along(hep_files)) {
    mv_name <- names(hep_files)[i]
    file_path <- hep_files[[i]]
    
    cat(sprintf("Loading %d/%d: %s\n", i, length(hep_files), mv_name))
    
    hep_data <- readMat(file_path)
    time_vector <- seq(-200, 800, length.out = dim(hep_data$HEP.all)[2])
    hep_data_cache[[mv_name]] <- list(hep_data = hep_data, time_vector = time_vector)
  }
  
  first_key <- names(hep_data_cache)[1]
  orig_n_subjects <- dim(hep_data_cache[[first_key]]$hep_data$HEP.all)[3]
  
  list(cache = hep_data_cache, channel_names = channel_names, orig_n_subjects = orig_n_subjects)
}

# ---------- Single Iteration ----------

run_single_iteration <- function(iter, sample_size, rho_values, anx_group, preloaded_data,
                                equivalence_margin, seed_base = 12345) {
  set.seed(seed_base + iter + sample_size * 1000)
  
  hep_cache <- preloaded_data$cache
  channel_names <- preloaded_data$channel_names
  orig_n_subjects <- preloaded_data$orig_n_subjects
  
  # Bootstrap subjects
  boot_indices <- sample(1:orig_n_subjects, size = sample_size, replace = (sample_size > orig_n_subjects))
  
  # Build long data across multiverses
  all_data <- data.frame()
  for (mv in names(hep_cache)) {
    cached_item <- hep_cache[[mv]]
    boot_hep <- cached_item$hep_data$HEP.all[, , boot_indices, drop = FALSE]
    
    avg_hep <- extract_roi_hep(
      list(HEP.all = boot_hep),
      time_vector = cached_item$time_vector,
      channel_names = channel_names
    )
    
    tmp <- data.frame(subject = 1:length(avg_hep), multiverse = mv, hep_amp = avg_hep)
    tmp <- tmp[!is.na(tmp$hep_amp), ]
    all_data <- rbind(all_data, tmp)
  }
  
  iter_results <- list()
  
  for (rho in rho_values) {
    # Per-subject HEP average for anxiety generation
    subject_avg_hep <- tapply(all_data$hep_amp, all_data$subject, mean, na.rm = TRUE)
    
    # Generate anxiety scores
    anxiety_scores <- generate_anxiety(subject_avg_hep, rho = rho)
    
    n_subjects_initial <- length(subject_avg_hep)
    
    if (anx_group) {
      # Grouping: low = 0-4, high = 10-21, exclude 5-9
      anxiety_numeric <- as.numeric(anxiety_scores)
      anxiety_binary <- ifelse(anxiety_numeric >= 10, "high",
                              ifelse(anxiety_numeric <= 4, "low", NA))
      
      valid_subjects <- !is.na(anxiety_binary)
      anxiety_var <- factor(anxiety_binary[valid_subjects], levels = c("low","high"))
      subject_names_valid <- names(subject_avg_hep)[valid_subjects]
      
      n_subjects_retained <- sum(valid_subjects)
      n_low_anxiety <- sum(anxiety_binary == "low", na.rm = TRUE)
      n_high_anxiety <- sum(anxiety_binary == "high", na.rm = TRUE)
      
      # Skip if only one group present
      if (nlevels(anxiety_var) < 2 || any(table(anxiety_var) == 0)) {
        next
      }
    } else {
      anxiety_var <- as.numeric(anxiety_scores)
      n_subjects_retained <- n_subjects_initial
      n_low_anxiety <- NA
      n_high_anxiety <- NA
    }
    
    # Merge anxiety with HEP data
    if (anx_group) {
      sim_data <- all_data %>%
        left_join(
          data.frame(subject = as.integer(subject_names_valid), anxiety = anxiety_var),
          by = "subject"
        ) %>%
        filter(!is.na(anxiety))
    } else {
      sim_data <- all_data %>%
        left_join(
          data.frame(subject = as.integer(names(subject_avg_hep)), anxiety = anxiety_var),
          by = "subject"
        )
    }
    
    # Fit model and test
    mm <- run_mixed_model(sim_data)
    model_obj <- mm$model
    model_sum <- mm$summary
    anx_row <- find_anxiety_row(model_sum)
    
    equiv_test <- run_equivalence_test(model_obj, equivalence_margin = equivalence_margin)
    
    # Sample anxiety for diagnostics
    samp_n <- min(200L, length(anxiety_scores))
    samp_idx <- sample(seq_along(anxiety_scores), size = samp_n)
    anxiety_sample <- data.frame(
      value = as.numeric(anxiety_scores[samp_idx]),
      rho = rho,
      anx_group = anx_group,
      sample_size = sample_size,
      iteration = iter
    )
    
    iter_results[[as.character(rho)]] <- list(
      rho = rho,
      p_value = model_sum$coefficients[anx_row, "Pr(>|t|)"],
      significant = model_sum$coefficients[anx_row, "Pr(>|t|)"] < 0.05,
      coefficient = model_sum$coefficients[anx_row, "Estimate"],
      equivalent = equiv_test$equivalent,
      equiv_p1 = equiv_test$p1,
      equiv_p2 = equiv_test$p2,
      equiv_p_max = equiv_test$equivalence_p,
      equiv_margin = equiv_test$margin,
      coef_row = anx_row,
      anx_group = anx_group,
      anxiety_sample = anxiety_sample,
      n_subjects_initial = n_subjects_initial,
      n_subjects_retained = n_subjects_retained,
      n_low_anxiety = n_low_anxiety,
      n_high_anxiety = n_high_anxiety,
      retention_rate = n_subjects_retained / n_subjects_initial
    )
  }
  
  list(sample_size = sample_size, iteration = iter, associations = iter_results)
}

# ---------- Main Simulation Function ----------

simulate_hep_anxiety <- function(rhos,
                                 anx_group,
                                 n_iterations,
                                 sample_sizes,
                                 equivalence_margin,
                                 output_dir,
                                 n_cores = NULL,
                                 parallel_mode = "fork",
                                 max_fork_workers = 32) {
  
  # Validate inputs
  if (any(!is.finite(rhos) | abs(rhos) >= 1)) {
    stop("All rhos must be finite and in (-1, 1)")
  }
  
  # Setup cores
  total_cores <- parallel::detectCores(logical = TRUE)
  if (is.null(n_cores)) {
    n_cores <- max(1, floor(total_cores * 0.6))
  }
  
  n_cores_env <- suppressWarnings(as.integer(Sys.getenv("HEP_MAX_WORKERS", "")))
  if (!is.na(n_cores_env) && n_cores_env > 0) {
    n_cores <- min(n_cores, n_cores_env)
  }
  
  if (parallel_mode == "fork") {
    n_cores <- min(n_cores, max_fork_workers)
  }
  if (parallel_mode == "sequential") {
    n_cores <- 1
  }
  
  cat(sprintf("Using %d cores in %s mode\n", n_cores, parallel_mode))
  
  # Pre-load HEP data
  preloaded_data <- preload_hep_data(hep_files, readMat(chan_names_dir))
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  run_iterations_for_size <- function(sample_size) {
    cat(sprintf("\n---------------- SAMPLE SIZE %d ----------------\n", sample_size))
    start_time <- Sys.time()
    
    worker_fun <- function(iter) {
      run_single_iteration(iter, sample_size, rhos, anx_group,
                          preloaded_data, equivalence_margin, seed_base = 12345)
    }
    
    if (parallel_mode == "sequential" || n_iterations == 1) {
      res_list <- lapply(1:n_iterations, worker_fun)
    } else if (parallel_mode == "fork") {
      res_list <- parallel::mclapply(1:n_iterations, worker_fun, mc.cores = n_cores)
    }
    
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    cat(sprintf("Completed: %d iterations in %.2f min (%.2f sec/iter)\n",
                n_iterations, elapsed, (elapsed*60)/n_iterations))
    res_list
  }
  
  all_results <- list()
  all_plot_data_list <- list()
  all_anxiety_samples_list <- list()
  
  for (ss in sample_sizes) {
    sample_results <- run_iterations_for_size(ss)
    valid_results <- sample_results[!sapply(sample_results, is.null)]
    cat(sprintf("Valid iterations: %d/%d\n", length(valid_results), length(sample_results)))
    
    cat("Processing results...\n")
    for (result in valid_results) {
      lst <- lapply(result$associations, function(a) {
        data_plot <- data.frame(
          sample_size = ss,
          iteration = result$iteration,
          rho = a$rho,
          anx_group = a$anx_group,
          coefficient = a$coefficient,
          significant = a$significant,
          equivalent = a$equivalent,
          equiv_p_max = a$equiv_p_max,
          n_subjects_initial = a$n_subjects_initial,
          n_subjects_retained = a$n_subjects_retained,
          n_low_anxiety = a$n_low_anxiety,
          n_high_anxiety = a$n_high_anxiety,
          retention_rate = a$retention_rate
        )
        list(plot = data_plot, anx = a$anxiety_sample)
      })
      all_plot_data_list <- c(all_plot_data_list, lapply(lst, `[[`, "plot"))
      all_anxiety_samples_list <- c(all_anxiety_samples_list, lapply(lst, `[[`, "anx"))
    }
    all_results[[as.character(ss)]] <- valid_results
    gc()
  }
  
  cat("Combining results...\n")
  all_plot_data <- do.call(rbind, all_plot_data_list)
  all_anxiety_samples <- do.call(rbind, all_anxiety_samples_list)
  
  cat("Calculating summaries...\n")
  power_data <- all_plot_data %>%
    group_by(sample_size, anx_group, rho) %>%
    summarise(
      power = mean(significant, na.rm = TRUE),
      equivalence_rate = mean(equivalent, na.rm = TRUE),
      mean_effect = mean(coefficient, na.rm = TRUE),
      se_effect = sd(coefficient, na.rm = TRUE) / sqrt(n()),
      n_iterations = n(),
      mean_n_subjects_initial = mean(n_subjects_initial, na.rm = TRUE),
      mean_n_subjects_retained = mean(n_subjects_retained, na.rm = TRUE),
      mean_retention_rate = mean(retention_rate, na.rm = TRUE),
      mean_n_low_anxiety = mean(n_low_anxiety, na.rm = TRUE),
      mean_n_high_anxiety = mean(n_high_anxiety, na.rm = TRUE),
      .groups = 'drop'
    )
  
  list(
    results = all_results,
    plot_data = all_plot_data,
    power_data = power_data,
    anxiety_samples = all_anxiety_samples
  )
}

# ---------- Run Simulations ----------

cat("Running HEP-Anxiety power simulations...\n\n")

output_dir <- "/data/p_02865/data/ds003838/output/hep_anxiety_sim_test"
n_iter <- 1000

cat("Continuous anxiety analysis...\n")
sim_results_continuous <- simulate_hep_anxiety(
  rhos = c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3),
  anx_group = FALSE,
  n_iterations = n_iter,
  sample_sizes = c(3000),
  equivalence_margin = 0.012,
  output_dir = output_dir
)

save(sim_results_continuous, file = file.path(output_dir, "sim_results_continuous.RData"))
write.csv(sim_results_continuous$power_data,
          file = file.path(output_dir, "power_data_continuous.csv"),
          row.names = FALSE)
cat("Continuous analysis saved!\n\n")

cat("Grouped anxiety analysis...\n")
sim_results_grouped <- simulate_hep_anxiety(
  rhos = c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3),
  anx_group = TRUE,
  n_iterations = n_iter,
  sample_sizes = c(3000),
  equivalence_margin = 0.16,
  output_dir = output_dir
)

save(sim_results_grouped, file = file.path(output_dir, "sim_results_grouped.RData"))
write.csv(sim_results_grouped$power_data,
          file = file.path(output_dir, "power_data_grouped.csv"),
          row.names = FALSE)
cat("Grouped analysis saved!\n\n")

# Create legacy-compatible aliases for plot_power_curves.R
cat("Creating legacy-compatible output...\n")

# Rename 'rho' column to 'strength' for backward compatibility
sim_results_continuous$power_data <- sim_results_continuous$power_data %>%
  rename(strength = rho) %>%
  mutate(grid_type = "rho")

sim_results_grouped$power_data <- sim_results_grouped$power_data %>%
  rename(strength = rho) %>%
  mutate(grid_type = "rho")

# Save with legacy names
sim_results_rho_continuous <- sim_results_continuous
sim_results_rho_grouped <- sim_results_grouped

save(sim_results_rho_continuous, file = file.path(output_dir, "sim_results_rho_continuous.RData"))
save(sim_results_rho_grouped, file = file.path(output_dir, "sim_results_rho_grouped.RData"))

cat("Legacy-compatible results saved!\n")

# Summary
cat("All simulations complete!\n")
cat("Results saved to:", output_dir, "\n")
