library(R.matlab)   # readMat
library(lme4)       # mixed models
library(lmerTest)   # p-values for lmer
library(dplyr)      # data wrangling
library(ggplot2)    # plotting
library(tidyr)      # (kept for compatibility)
library(gridExtra)  # (kept for compatibility)
library(parallel)   # parallel
library(fGarch)     # qsnorm for skew-normal transform
library(pbapply)
library(scales)     # for scales::percent in plots

set.seed(423) # fix rng

# ---------- IO & helpers ----------

load_hep_data <- function(file_path) readMat(file_path)

default_output_dir <- function() {
  "/data/p_02865/data/ds003838/output/hep_anxiety_sim"
}

chan_names_dir <- "/data/p_02865/code/matlab/HEP/supplement/orig_chan_labels.mat"
chan_names <- readMat(chan_names_dir)
base_dir <- "/data/p_02865/data/ds003838/derivatives/"

extract_chan_names <- function(chan_names) {
  labels <- sapply(chan_names$orig.chan.labels, function(x) x[[1]][1, 1])
  return(labels)
}

hep_files <- list(
  "hp_0.1_bl_m150_m50_ic_0.0"    = paste0(base_dir, "hp_0.1_bl_[-150 -50]_ic_0.0_all_HEP.mat"),
  "hp_0.1_bl_m150_m50_ic_1.5"    = paste0(base_dir, "hp_0.1_bl_[-150 -50]_ic_1.5_all_HEP.mat"),
  "hp_0.1_bl_m150_m50_ic_3.5"    = paste0(base_dir, "hp_0.1_bl_[-150 -50]_ic_3.5_all_HEP.mat"),
  "hp_0.1_bl_m150_m50_ic_NaN"    = paste0(base_dir, "hp_0.1_bl_[-150 -50]_ic_NaN_all_HEP.mat"),
  "hp_0.1_bl_none_ic_0.0"        = paste0(base_dir, "hp_0.1_bl_[]_ic_0.0_all_HEP.mat"),
  "hp_0.1_bl_none_ic_1.5"        = paste0(base_dir, "hp_0.1_bl_[]_ic_1.5_all_HEP.mat"),
  "hp_0.1_bl_none_ic_3.5"        = paste0(base_dir, "hp_0.1_bl_[]_ic_3.5_all_HEP.mat"),
  "hp_0.1_bl_none_ic_NaN"        = paste0(base_dir, "hp_0.1_bl_[]_ic_NaN_all_HEP.mat"),
  "hp_0.5_bl_m150_m50_ic_0.0"    = paste0(base_dir, "hp_0.5_bl_[-150 -50]_ic_0.0_all_HEP.mat"),
  "hp_0.5_bl_m150_m50_ic_1.5"    = paste0(base_dir, "hp_0.5_bl_[-150 -50]_ic_1.5_all_HEP.mat"),
  "hp_0.5_bl_m150_m50_ic_3.5"    = paste0(base_dir, "hp_0.5_bl_[-150 -50]_ic_3.5_all_HEP.mat"),
  "hp_0.5_bl_m150_m50_ic_NaN"    = paste0(base_dir, "hp_0.5_bl_[-150 -50]_ic_NaN_all_HEP.mat"),
  "hp_0.5_bl_none_ic_0.0"        = paste0(base_dir, "hp_0.5_bl_[]_ic_0.0_all_HEP.mat"),
  "hp_0.5_bl_none_ic_1.5"        = paste0(base_dir, "hp_0.5_bl_[]_ic_1.5_all_HEP.mat"),
  "hp_0.5_bl_none_ic_3.5"        = paste0(base_dir, "hp_0.5_bl_[]_ic_3.5_all_HEP.mat"),
  "hp_0.5_bl_none_ic_NaN"        = paste0(base_dir, "hp_0.5_bl_[]_ic_NaN_all_HEP.mat")
)

extract_roi_hep <- function(hep_data, electrodes = c("Cz","C3","C4","Fz","F3","F4","FC1","FC2"),
                            time_window = c(200, 500), time_vector = NULL, channel_names = NULL) {
  hep_values <- hep_data$HEP.all
  electrode_indices <- match(electrodes, channel_names)
  time_indices <- which(time_vector >= time_window[1] & time_vector <= time_window[2])
  roi_hep <- hep_values[electrode_indices, time_indices, , drop = FALSE]
  avg_hep <- apply(roi_hep, 3, function(x) mean(x, na.rm = TRUE))
  return(avg_hep)
}

# ---------- Anxiety simulation ----------

# (A) Coefficient-based 
generate_anxiety_coef <- function(hep_values, target_coef, intercept = 0, noise_sd = 1,
                                  max_iterations = 10, tolerance = 0.001) {
  n_subjects <- length(hep_values)
  if (abs(target_coef) < 1e-6) {
    anxiety_gad7 <- round(runif(n_subjects, min = 0, max = 21))
    if (length(hep_values) > 10) {
      achieved_coef <- coef(lm(hep_values ~ anxiety_gad7))[2]
      cat(sprintf("Target=%.4f: Achieved ordinal coef=%.4f (null case)\n", target_coef, achieved_coef))
    }
    return(anxiety_gad7)
  }
  latent_coef_guess <- target_coef * 3
  best_anxiety <- NULL; best_error <- Inf; best_coef <- NA
  for (iteration in 1:max_iterations) {
    latent_L <- intercept + latent_coef_guess * hep_values + rnorm(n_subjects, 0, noise_sd)
    cut_points <- quantile(latent_L, probs = seq(0, 1, length.out = 23))
    anxiety_gad7 <- as.numeric(cut(latent_L, breaks = cut_points, labels = 0:21, include.lowest = TRUE))
    anxiety_gad7[is.na(anxiety_gad7)] <- round(mean(anxiety_gad7, na.rm = TRUE))
    anxiety_gad7 <- pmax(0, pmin(21, anxiety_gad7))
    achieved_coef <- coef(lm(hep_values ~ anxiety_gad7))[2]
    error <- abs(achieved_coef - target_coef)
    if (error < best_error) { best_error <- error; best_anxiety <- anxiety_gad7; best_coef <- achieved_coef }
    if (error < tolerance) {
      if (length(hep_values) > 10) {
        cat(sprintf("Target=%.4f: Achieved ordinal coef=%.4f (converged in %d iterations)\n",
                    target_coef, achieved_coef, iteration))
      }
      return(anxiety_gad7)
    }
    den <- if (abs(achieved_coef) < 1e-8) sign(achieved_coef) * 1e-8 else achieved_coef
    ratio <- target_coef / den
    latent_coef_guess <- latent_coef_guess * ratio * 0.8
  }
  if (length(hep_values) > 10) {
    cat(sprintf("Target=%.4f: Achieved ordinal coef=%.4f (best after %d iterations)\n",
                target_coef, best_coef, max_iterations))
  }
  return(best_anxiety)
}

# (B) Correlation-based (clamp to [0,21] and round)
generate_anxiety_corr <- function(hep_values, rho) {
  if (!is.finite(rho) || abs(rho) >= 1) stop("rho must be finite and strictly between -1 and 1.")
  n <- length(hep_values)
  hep_z <- as.numeric(scale(hep_values))
  z_noise <- rnorm(n, 0, 1)
  latent <- rho * hep_z + sqrt(1 - rho^2) * z_noise
  # Skew transform
  u <- pnorm(latent)
  anxiety_cont <- qsnorm(u, mean = 5, sd = 4, xi = 3)
  # Clamp and discretize to real GAD-7 range
  anxiety <- pmax(0, pmin(21, anxiety_cont))
  # anxiety <- as.integer(round(anxiety)) # optional, if needed 
  return(anxiety)
}

# ---------- Model fitting & tests ----------

# Fit once; sum contrasts for multiverse so effect is grand mean across multiverses
# If anxiety is a 2-level factor (low/high), use treatment coding with 'low' as baseline so coef = high - low
run_mixed_model <- function(data) {
  data$multiverse <- factor(data$multiverse)
  contrasts(data$multiverse) <- contr.sum(nlevels(data$multiverse))
  
  if (is.factor(data$anxiety)) {
    data$anxiety <- factor(data$anxiety)  # drop empty levels if any
    if (nlevels(data$anxiety) < 2) stop("anxiety factor has fewer than 2 levels")
    data$anxiety <- factor(data$anxiety, levels = c("low", "high"))
    contrasts(data$anxiety) <- contr.treatment(2, base = 1)
  }
  
  model <- lmer(hep_amp ~ anxiety + multiverse + anxiety:multiverse + (1 | subject), data = data)
  list(model = model, summary = summary(model))
}

# Helper: find the right anxiety coefficient row name in summary (continuous or factor)
find_anxiety_row <- function(ms) {
  rows <- rownames(ms$coefficients)
  # exclude interactions
  candidates <- rows[grepl("^anxiety", rows) & !grepl(":", rows)]
  if (length(candidates) < 1) stop("'anxiety' coefficient not found in model summary")
  # For factor (low/high), this will be 'anxietyhigh'; for continuous it's 'anxiety'
  candidates[1]
}

# Manual TOST on the anxiety effect row (slope or high-low difference)
run_equivalence_test <- function(model, equivalence_margin = 0.005) {
  low_bound <- -equivalence_margin; high_bound <- equivalence_margin
  ms <- summary(model)
  rowname <- find_anxiety_row(ms)
  coef_anx <- ms$coefficients[rowname, "Estimate"]
  se_anx   <- ms$coefficients[rowname, "Std. Error"]
  df_anx   <- ifelse(inherits(model, "lmerMod"), ms$coefficients[rowname, "df"], df.residual(model))
  if (any(is.na(c(coef_anx, se_anx, df_anx)))) stop("NA in coef/SE/df for TOST")
  t_upper <- (coef_anx - high_bound) / se_anx; p_upper <- pt(t_upper, df = df_anx, lower.tail = TRUE)
  t_lower <- (coef_anx - low_bound) / se_anx;  p_lower <- pt(t_lower, df = df_anx, lower.tail = FALSE)
  equivalent <- p_upper < 0.05 && p_lower < 0.05
  list(equivalent = equivalent,
       p1 = p_lower, p2 = p_upper,
       equivalence_p = max(p_upper, p_lower),
       margin = equivalence_margin,
       method = ifelse(inherits(model,"lmerMod"), "lmer_manual_tost","lm_manual_tost"),
       coefficient = coef_anx,
       coef_row = rowname)
}

# Pre-load all HEP data once to avoid repeated I/O in parallel workers
preload_hep_data <- function(hep_files, chan_names_data) {
  cat("Pre-loading HEP data for parallel processing...\n")
  channel_names <- extract_chan_names(chan_names_data)
  
  hep_data_cache <- list()
  total_files <- length(hep_files)
  loaded_any <- FALSE
  
  for (i in seq_along(hep_files)) {
    mv_name <- names(hep_files)[i]
    file_path <- hep_files[[i]]
    
    cat(sprintf("Loading %d/%d: %s\n", i, total_files, mv_name))
    
    hep_data <- load_hep_data(file_path)
    time_vector <- seq(-200, 800, length.out = dim(hep_data$HEP.all)[2])
    hep_data_cache[[mv_name]] <- list(hep_data = hep_data, time_vector = time_vector)
    loaded_any <- TRUE
  }
  # Determine original subject count from first cached entry
  first_key <- names(hep_data_cache)[1]
  orig_n_subjects <- dim(hep_data_cache[[first_key]]$hep_data$HEP.all)[3]
  list(cache = hep_data_cache, channel_names = channel_names, orig_n_subjects = orig_n_subjects)
}

# Optimized single iteration using pre-loaded data
run_single_iteration_optimized <- function(iter, sample_size, grid_values, simulate_by = c("coefficient","correlation"),
                                          anx_group = FALSE, preloaded_data,
                                          equivalence_margin = 0.005, seed_base = 12345) {
  simulate_by <- match.arg(simulate_by)
  set.seed(seed_base + iter + sample_size * 1000)
  
  # Use pre-loaded data
  hep_cache <- preloaded_data$cache
  channel_names <- preloaded_data$channel_names
  orig_n_subjects <- preloaded_data$orig_n_subjects
  
  # Bootstrap subjects
  boot_indices <- if (sample_size > orig_n_subjects) {
    sample(1:orig_n_subjects, size = sample_size, replace = TRUE)
  } else {
    sample(1:orig_n_subjects, size = sample_size, replace = FALSE)
  }
    
  # Build long data across multiverses using cached data
  all_data <- data.frame()
  for (mv in names(hep_cache)) {
    cached_item <- hep_cache[[mv]]
    boot_hep <- cached_item$hep_data$HEP.all[, , boot_indices, drop = FALSE]
    
    avg_hep <- extract_roi_hep(list(HEP.all = boot_hep), 
                               time_vector = cached_item$time_vector, 
                               channel_names = channel_names)
    if (is.null(avg_hep) || all(is.na(avg_hep))) next
    
    tmp <- data.frame(subject = 1:length(avg_hep), multiverse = mv, hep_amp = avg_hep)
    tmp <- tmp[!is.na(tmp$hep_amp), ]
    if (nrow(tmp) > 0) all_data <- rbind(all_data, tmp)
  }
  if (nrow(all_data) == 0) return(NULL)
  
  iter_results <- list()
  
  for (g in grid_values) {
    # Per-subject HEP average across multiverses for anxiety generation
    subject_avg_hep <- tapply(all_data$hep_amp, all_data$subject, mean, na.rm = TRUE)
    
    if (simulate_by == "coefficient") {
      anxiety_scores <- generate_anxiety_coef(subject_avg_hep, target_coef = g)
      passed_value <- g
      passed_type  <- "beta_target"
      passed_rho   <- NA_real_
    } else {
      anxiety_scores <- generate_anxiety_corr(subject_avg_hep, rho = g)
      passed_value <- g
      passed_type  <- "rho"
      passed_rho   <- g
    }
    
    # Optional: convert to groups (low/high)
    n_subjects_initial <- length(subject_avg_hep)
    n_subjects_retained <- n_subjects_initial
    n_low_anxiety <- NA
    n_high_anxiety <- NA
    
    if (anx_group) {

      # Grouping rule: low = 0–4 (inclusive), high = 10–21 (inclusive)
      anxiety_numeric <- as.numeric(anxiety_scores)
      anxiety_binary <- ifelse(anxiety_numeric >= 10, "high", 
                              ifelse(anxiety_numeric <= 4, "low", NA))
      # exclude 5-9
      valid_subjects <- !is.na(anxiety_binary)
      anxiety_var <- factor(anxiety_binary[valid_subjects], levels = c("low","high"))
      subject_names_valid <- names(subject_avg_hep)[valid_subjects]
      
      # Track retention numbers
      n_subjects_retained <- sum(valid_subjects)
      n_low_anxiety <- sum(anxiety_binary == "low", na.rm = TRUE)
      n_high_anxiety <- sum(anxiety_binary == "high", na.rm = TRUE)
      
      # Guard — if only one group, skip this grid value for this iteration
      if (nlevels(anxiety_var) < 2 || any(table(anxiety_var) == 0)) {
        next
      }
    } else {
      anxiety_var <- as.numeric(anxiety_scores)
    }
    
    # Merge back
    if (anx_group && exists("subject_names_valid")) {
      # For grouped analysis with excluded subjects
      sim_data <- all_data %>%
        left_join(
          data.frame(subject = as.integer(subject_names_valid),
                     anxiety = anxiety_var),
          by = "subject"
        ) %>%
        filter(!is.na(anxiety))  # Keep only subjects with valid anxiety groups
    } else {
      # For continuous analysis
      sim_data <- all_data %>%
        left_join(
          data.frame(subject = as.integer(names(subject_avg_hep)),
                     anxiety = anxiety_var),
          by = "subject"
        )
    }
    
    # Fit model and run tests
    tryCatch({
      mm <- run_mixed_model(sim_data)
      model_obj <- mm$model
      model_sum <- mm$summary
      anx_row <- find_anxiety_row(model_sum)
      
    # Manual TOST
      equiv_test <- run_equivalence_test(model_obj, equivalence_margin = equivalence_margin)
      
      # Small sample of anxiety values for histogram sanity check
      samp_n <- min(200L, length(anxiety_scores))
      samp_idx <- sample(seq_along(anxiety_scores), size = samp_n)
      anxiety_sample <- data.frame(
        value = as.numeric(if (anx_group) as.numeric(anxiety_scores[samp_idx]) else anxiety_scores[samp_idx]),
        grid_value = passed_value,
        grid_type = passed_type,
        simulate_by = simulate_by,
        anx_group = anx_group,
        sample_size = sample_size,
        iteration = iter
      )
      
      iter_results[[as.character(passed_value)]] <- list(
        grid_value = passed_value,
        grid_type = passed_type,
        rho_passed = passed_rho,
        p_value = model_sum$coefficients[anx_row, "Pr(>|t|)"],
        significant = model_sum$coefficients[anx_row, "Pr(>|t|)"] < 0.05,
      coefficient = model_sum$coefficients[anx_row, "Estimate"],  # slope or high-low diff
        equivalent = equiv_test$equivalent,
        equiv_p1 = equiv_test$p1,
        equiv_p2 = equiv_test$p2,
        equiv_p_max = equiv_test$equivalence_p,
        equiv_margin = equiv_test$margin,
        equiv_method = equiv_test$method,
        coef_row = anx_row,
        anx_group = anx_group,
        anxiety_sample = anxiety_sample,
        # Subject retention information
        n_subjects_initial = n_subjects_initial,
        n_subjects_retained = n_subjects_retained,
        n_low_anxiety = n_low_anxiety,
        n_high_anxiety = n_high_anxiety,
        retention_rate = n_subjects_retained / n_subjects_initial
      )
    }, error = function(e) {
      warning(sprintf("Error in iteration %d, grid value %s: %s", iter, as.character(g), e$message))
    })
  }
  
  list(sample_size = sample_size, iteration = iter, associations = iter_results)
}

# ---------- Quick test ----------

test_contest_function <- function() {
  cat("Testing equivalence test function...\n")
  set.seed(123)
  test_data <- data.frame(
    subject = rep(1:20, each = 5),
    anxiety = rep(rnorm(20), each = 5),
    hep_amp = rnorm(100),
    multiverse = rep(c("A","B","C","D","E"), 20)
  )
  mm <- run_mixed_model(test_data)
  test_model <- mm$model
  cat("Fixed effects:", names(fixef(test_model)), "\n")
  cat("Anxiety coefficient:", fixef(test_model)["anxiety"], "\n")
  eq <- run_equivalence_test(test_model, equivalence_margin = 0.005)
  print(eq)
  cat("P1 =", eq$p1, ", P2 =", eq$p2, " | Equivalent =", eq$equivalent, "\n")
  TRUE
}

#test_contest_function()

# ---------- Main driver ----------

# simulate_by: "coefficient" (use association_strengths) or "correlation" (use rhos)
simulate_hep_anxiety <- function(association_strengths = c(0, 0.001, 0.002, 0.005, 0.008, 0.010, 0.011, 0.012, 0.013, 0.014, 0.015, 0.02),
                                 rhos = NULL,
                                 simulate_by = c("coefficient","correlation"),
                                 anx_group = FALSE,
                                 n_iterations = 1000,
                                 sample_sizes = c(10, 50, 100, 1000, 3000),
                                 equivalence_margin = 0.05,
                                 n_cores = NULL,
                                 use_progress_bar = TRUE,
                                 output_dir = NULL,
                                 parallel_mode = c("fork","sequential"),
                                 max_fork_workers = 32,
                                 plot_prefix = NULL) {

  
  simulate_by   <- match.arg(simulate_by)
  parallel_mode <- match.arg(parallel_mode)

  grid_values <- if (simulate_by == "coefficient") {
    association_strengths
  } else {
    if (is.null(rhos)) stop("Provide rhos for simulate_by='correlation'.")
    if (any(!is.finite(rhos) | abs(rhos) >= 1)) stop("All rhos must be finite & in (-1,1).")
    rhos
  }

  total_cores <- parallel::detectCores(logical = TRUE)
  if (is.null(n_cores)) n_cores <- max(1, floor(total_cores * 0.6))
  n_cores_env <- suppressWarnings(as.integer(Sys.getenv("HEP_MAX_WORKERS","")))
  if (!is.na(n_cores_env) && n_cores_env > 0) n_cores <- min(n_cores, n_cores_env)
  if (parallel_mode == "fork")  n_cores <- min(n_cores, max_fork_workers)
  if (parallel_mode == "sequential") n_cores <- 1

  preloaded_data <- preload_hep_data(hep_files, readMat(chan_names_dir))

  if (is.null(output_dir) || !nzchar(output_dir)) output_dir <- default_output_dir()
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  run_iterations_for_size <- function(sample_size) {
    cat(sprintf("\n---------------- SAMPLE SIZE %d ----------------\n", sample_size))
    start_time <- Sys.time()

    worker_fun <- function(iter) {
      run_single_iteration_optimized(iter, sample_size, grid_values, simulate_by, anx_group,
                                     preloaded_data, equivalence_margin, seed_base = 12345)
    }

  if (parallel_mode == "sequential" || n_iterations == 1) {
      res_list <- lapply(1:n_iterations, worker_fun)
  } else if (parallel_mode == "fork") {
      res_list <- parallel::mclapply(1:n_iterations, worker_fun, mc.cores = n_cores)
    }

    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    cat(sprintf("Sample size %d completed: %d iterations in %.2f min (%.2f sec/iter)\n",
                sample_size, n_iterations, elapsed, (elapsed*60)/max(1,n_iterations)))
    res_list
  }

  all_results <- list()
  all_plot_data_list <- list()
  all_anxiety_samples_list <- list()

  for (ss in sample_sizes) {
    sample_results <- run_iterations_for_size(ss)
    valid_results <- sample_results[!sapply(sample_results, is.null)]
    cat(sprintf("Valid iteration results: %d/%d\n", length(valid_results), length(sample_results)))

    cat("Processing results for plotting data...\n")
    for (result in valid_results) {
      lst <- lapply(result$associations, function(a) {
        if (is.null(a)) return(NULL)
        data_plot <- data.frame(
          sample_size  = ss,
          iteration    = result$iteration,
          strength     = as.numeric(a$grid_value),
          grid_type    = a$grid_type,
          anx_group    = a$anx_group,
          coefficient  = a$coefficient,
          significant  = a$significant,
          equivalent   = a$equivalent,
          equiv_p_max  = a$equiv_p_max,
          equiv_method = ifelse(is.null(a$equiv_method),"unknown",a$equiv_method),
          # Add retention information
          n_subjects_initial = ifelse(is.null(a$n_subjects_initial), ss, a$n_subjects_initial),
          n_subjects_retained = ifelse(is.null(a$n_subjects_retained), ss, a$n_subjects_retained),
          n_low_anxiety = ifelse(is.null(a$n_low_anxiety), NA, a$n_low_anxiety),
          n_high_anxiety = ifelse(is.null(a$n_high_anxiety), NA, a$n_high_anxiety),
          retention_rate = ifelse(is.null(a$retention_rate), 1.0, a$retention_rate)
        )
        list(plot = data_plot, anx = a$anxiety_sample)
      })
      lst <- Filter(Negate(is.null), lst)
      if (length(lst)) {
        all_plot_data_list <- c(all_plot_data_list, lapply(lst, `[[`, "plot"))
        all_anxiety_samples_list <- c(all_anxiety_samples_list, lapply(lst, `[[`, "anx"))
      }
    }
    all_results[[as.character(ss)]] <- valid_results
    gc()
  }

  cat("Combining all data (single rbind operation)...\n")
  all_plot_data <- if(length(all_plot_data_list) > 0) do.call(rbind, all_plot_data_list) else data.frame()
  all_anxiety_samples <- if(length(all_anxiety_samples_list) > 0) do.call(rbind, all_anxiety_samples_list) else data.frame()

  cat("Calculating summaries...\n")

  power_data <- all_plot_data %>%
    group_by(sample_size, grid_type, anx_group, strength) %>%
    summarise(
      power = mean(significant, na.rm = TRUE),
      equivalence_rate = mean(equivalent, na.rm = TRUE),
      mean_effect = mean(coefficient, na.rm = TRUE),
      se_effect = sd(coefficient, na.rm = TRUE) / sqrt(n()),
      n_iterations = n(),
      # Add retention statistics
      mean_n_subjects_initial = mean(n_subjects_initial, na.rm = TRUE),
      mean_n_subjects_retained = mean(n_subjects_retained, na.rm = TRUE),
      mean_retention_rate = mean(retention_rate, na.rm = TRUE),
      mean_n_low_anxiety = mean(n_low_anxiety, na.rm = TRUE),
      mean_n_high_anxiety = mean(n_high_anxiety, na.rm = TRUE),
      .groups = 'drop'
    )
  
  cat("Skipping plot generation for faster completion...\n")
  combined_plot <- NULL
  power_plot <- NULL  
  equivalence_plot <- NULL
  anxiety_hist <- NULL
  
  list(
    results = all_results,
    plot_data = all_plot_data,
    power_data = power_data,
    anxiety_samples = all_anxiety_samples,
    combined_plot = combined_plot,
    power_plot = power_plot,
    equivalence_plot = equivalence_plot,
    anxiety_hist = anxiety_hist
  )
}

# ---------- Power Analysis for rsHEP Registered Report ----------

cat("Running power simulations for rsHEP registered report...\n")

cat("Running production simulations for continuous anxiety (Hypotheses 1a, 1c)...\n")
n_iter = 1000
# sim_results_continuous <- simulate_hep_anxiety(
#   association_strengths = c(0, 0.002, 0.005, 0.008, 0.010, 0.012, 0.015, 0.018, 0.020),
#   simulate_by = "coefficient",
#   anx_group = FALSE,
#   n_iterations = n_iter,
#   sample_sizes = c(100, 300, 500, 1000, 1500, 2000, 3000),
#   equivalence_margin = 0.01,
#   n_cores = NULL,
#   use_progress_bar = TRUE
# )

cat("Running production simulations for group comparison (Hypotheses 1b, 1d)...\n")
# sim_results_grouped <- simulate_hep_anxiety(
#   association_strengths = c(0, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030),
#   simulate_by = "coefficient",
#   anx_group = TRUE,
#   n_iterations = n_iter,
#   sample_sizes = c(50, 100, 200, 300, 500),
#   equivalence_margin = 0.01,
#   n_cores = NULL,
#   use_progress_bar = TRUE
# )

cat("Running production correlation-based simulations (ONLY mode)...\n")

# Setup output directory
output_dir <- default_output_dir()
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("Running continuous anxiety simulation...\n")
sim_results_rho_continuous <- simulate_hep_anxiety(
  rhos = c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11, 0.12, 0.13, 0.14, 0.15),
  simulate_by = "correlation",
  anx_group = FALSE,
  n_iterations = n_iter,
  sample_sizes = c(300, 3000),
  equivalence_margin = 0.012,
  n_cores = NULL,
  use_progress_bar = TRUE
)

# Save continuous results immediately
cat("Saving continuous simulation results...\n")
save(sim_results_rho_continuous, file = file.path(output_dir, "sim_results_rho_continuous.RData"))
write.csv(sim_results_rho_continuous$power_data, 
          file = file.path(output_dir, "power_data_rho_continuous.csv"), 
          row.names = FALSE)
cat("Continuous simulation completed and saved!\n\n")

cat("Running grouped anxiety simulation...\n")
sim_results_rho_grouped <- simulate_hep_anxiety(
  # rhos = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.32, 0.34, 0.36, 0.38, 0.40, 0.42, 0.44, 0.46, 0.48, 0.50, 0.55, 0.60, 0.65),
  rhos = c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11, 0.12, 0.13, 0.14, 0.15),
  simulate_by = "correlation",
  anx_group = TRUE,
  n_iterations = n_iter,
  sample_sizes = c(300, 3000),
  equivalence_margin = 0.19,
  n_cores = NULL,
  use_progress_bar = TRUE
)

# Save grouped results immediately
cat("Saving grouped simulation results...\n")
save(sim_results_rho_grouped, file = file.path(output_dir, "sim_results_rho_grouped.RData"))
write.csv(sim_results_rho_grouped$power_data, 
          file = file.path(output_dir, "power_data_rho_grouped.csv"), 
          row.names = FALSE)
cat("Grouped simulation completed and saved!\n\n")

# Create legacy-compatible aliases
sim_results_continuous <- sim_results_rho_continuous
sim_results_grouped    <- sim_results_rho_grouped

# Provide a unified object similar to earlier sim_roi_1 output (optional helper)
sim_results_combined <- list(
  continuous = sim_results_continuous,
  grouped = sim_results_grouped,
  created = Sys.time(),
  note = "Combined container holding continuous and grouped correlation simulation results"
)

# Save optional combined result (but main results already saved above)
save(sim_results_combined, file = file.path(output_dir, "sim_results_combined.RData"))
cat("Combined results container saved!\n\n")

# Create a detailed effect size table for 95% power
create_detailed_effect_size_table <- function() {
  # Combine both continuous and grouped results
  all_power_data <- rbind(
    sim_results_continuous$power_data %>% mutate(source = "continuous"),
    sim_results_grouped$power_data %>% mutate(source = "grouped")
  )
  
  effect_size_summary <- data.frame()
  
  for (sample_size in unique(all_power_data$sample_size)) {
    for (anx_group in c(TRUE, FALSE)) {
      # Get data for this configuration
      subset_data <- all_power_data %>%
        filter(sample_size == !!sample_size, anx_group == !!anx_group) %>%
        arrange(strength)
      
      # Find minimum effect size for 95% power
      min_effect <- subset_data %>% filter(power >= 0.95) %>% slice(1)
      
      if (nrow(min_effect) > 0) {
        effect_size_summary <- rbind(effect_size_summary, data.frame(
          Sample_Size = sample_size,
          Analysis_Type = ifelse(anx_group, "Group_Comparison", "Continuous"),
          Min_Effect_Size = min_effect$strength,
          Achieved_Power = round(min_effect$power * 100, 1),
          Mean_Effect = round(min_effect$mean_effect, 4),
          Equivalence_Rate = round(min_effect$equivalence_rate * 100, 1)
        ))
      }
    }
  }
  
  return(effect_size_summary)
}

# Generate detailed effect size table for 95% power
effect_sizes_95_power <- create_detailed_effect_size_table()

write.csv(effect_sizes_95_power, file = file.path(output_dir, "effect_sizes_95_percent_power.csv"), row.names = FALSE)

build_legacy_sim_results <- function() {
  # Ensure aliases exist
  mk_alias <- function(pd) {
    if ("mean_coefficient" %in% names(pd) && !("mean_effect" %in% names(pd))) pd$mean_effect <- pd$mean_coefficient
    if ("se_coefficient" %in% names(pd) && !("se_effect" %in% names(pd)))   pd$se_effect   <- pd$se_coefficient
    pd
  }

  pd_cont <- mk_alias(sim_results_continuous$power_data)
  pd_group <- mk_alias(sim_results_grouped$power_data)

  # Bind rows but keep markers for analysis type to aid filtering later
  power_data_legacy <- bind_rows(
    mutate(pd_cont, Analysis_Type = "Continuous"),
    mutate(pd_group, Analysis_Type = "Group_Comparison")
  )

  # Compose a compact legacy-style object
  list(
    power_data   = power_data_legacy,
    plot_data    = bind_rows(
      mutate(sim_results_continuous$plot_data, Analysis_Type = "Continuous"),
      mutate(sim_results_grouped$plot_data,    Analysis_Type = "Group_Comparison")
    ),
    results      = list(
      continuous = sim_results_continuous$results,
      grouped    = sim_results_grouped$results
    ),
    created      = Sys.time(),
    note         = "Legacy-compatible container built from correlation-based simulations"
  )
}