#!/usr/bin/env Rscript
# rsHEP / LIFE: Join the 16 pipeline HEP summaries AFTER common exclusions.
# Run after LIFE_matching.R -> LIFE_generate_data.R -> exclusion.R:
#   Rscript LIFE_pipeline_matching.R [life_directory]
# Default: /data/p_02865/data/LIFE
# Requires: install.packages("readODS")
#
# INPUTS:
#   derivatives/data_after_exclusion.csv
#   derivatives/group/*_roi_HEP.csv   (16 files; the source of 5 fields each)
#   tabular/preprocessing_logs.ods    (hp_0.1 and hp_0.5 final status)
# OUTPUTS:
#   derivatives/data_with_pipeline_flags.csv  (ALL retained rows; none dropped)
#     Per pipeline: five original ROI/epoch fields + HEP z-score,
#     HEP 3-SD outlier flag, and combined pipeline exclusion flag.
#     Preprocessing history and other intermediate flags are QC-only.
#   derivatives/matching_and_summaries/pipeline_exclusion/*  (diagnostics)
#
# Exclusion coding:
#   Preprocessing status: "excluded" => excluded; "ok" or missing => not.
#   <pipeline>_n_epochs_kept < 300 => excluded; missing => not.
#   |<pipeline>_roi_HEP_SD| > 3 => excluded; missing => not.
#   <pipeline>_exclude = OR of the above three reasons, 1/0 only.
#   <pipeline>_roi_HEP_SD is a z-SCORE, not the group SD itself.
#   Missing source values are separately marked for QC, NEVER treated as an
#   automatic exclusion. In particular, no overall exclusion is computed.
#
# The 3-SD reference uses unique life_id participants who survived common
# exclusions and, for this pipeline, are not preprocessing failures and have
# >=300 epochs (when an epoch count is available). Unknown preprocessing or
# epoch information does not remove a participant from this reference.
# This avoids weighting unresolved duplicate LIFE rows or clearly failed
# recordings. Its reference sample size and mean/SD are saved for each pipeline.
#
# For repeated source IDs, a field is copied only when its nonmissing values
# agree. Conflicts leave that field blank and create a diagnostic event.
# Joined values are copied identically to all surviving rows of each LIFE ID.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: Rscript LIFE_pipeline_matching.R [life_directory]")
life <- if (length(args)) args[[1L]] else "/data/p_02865/data/LIFE"
deriv <- file.path(life, "derivatives")
out <- file.path(deriv, "matching_and_summaries", "pipeline_exclusion")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing required input: ", path)
  x <- read.csv(path, colClasses = "character", check.names = FALSE,
                stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM",
                na.strings = "")
  names(x) <- trimws(names(x))
  if (anyDuplicated(names(x))) stop("Duplicate headers in ", path)
  x
}
write_csv <- function(x, name) write.csv(x, file.path(out, name), na = "", row.names = FALSE)
norm <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | !nzchar(x)] <- NA_character_
  x
}
ids_norm <- function(x) {
  x <- toupper(norm(x))
  x[!is.na(x) & x == "NA"] <- NA_character_
  x
}
joined <- function(x, sep = " | ") {
  x <- unique(norm(x))
  paste(x[!is.na(x)], collapse = sep)
}
num <- function(x) {
  x <- suppressWarnings(as.numeric(gsub(",", ".", norm(x), fixed = TRUE)))
  x[!is.finite(x)] <- NA_real_
  x
}
empty <- function(cols) as.data.frame(
  setNames(rep(list(character()), length(cols)), cols),
  stringsAsFactors = FALSE, check.names = FALSE)
bind <- function(rows, cols) {
  if (!length(rows)) return(empty(cols))
  x <- do.call(rbind, rows)
  rownames(x) <- NULL
  x[, cols, drop = FALSE]
}

input_path <- file.path(deriv, "data_after_exclusion.csv")
data <- read_csv(input_path)
if (!nrow(data) || !"life_id" %in% names(data))
  stop("data_after_exclusion.csv needs nonempty rows and life_id")
data$life_id <- ids_norm(data$life_id)
if (anyNA(data$life_id)) stop("Missing life_id in data_after_exclusion.csv")
if ("overall_exclude" %in% names(data) && any(num(data$overall_exclude) == 1, na.rm = TRUE))
  stop("data_after_exclusion.csv still has overall_exclude=1 rows. Run exclusion.R first.")
# Avoid accidentally mixing previously generated HEP columns with fresh results.
old_pipe <- grep("(^preproc_hp_|_roi_(HEP|ECG)(_|$)|_n_epochs_(total|rejected|kept)(_|$))",
                 names(data), value = TRUE, ignore.case = TRUE, perl = TRUE)
if (length(old_pipe)) stop("Input already contains pipeline columns (",
  joined(head(old_pipe, 8L), ", "), "). Re-run the NEW LIFE_generate_data.R, then exclusion.R.")

uid <- unique(data$life_id)
nids <- length(uid)
row_ix <- seq_len(nrow(data)) + 1L
uid_row <- match(uid, data$life_id)
canonical <- function(x) {
  x <- ids_norm(x)
  x <- sub("(?i)_(Ruhe|Pilot)(?:_[0-9]+)?$", "", x, perl = TRUE)
  x <- sub("[_-]+2$", "", x, perl = TRUE)
  x <- sub("[_-]+$", "", x, perl = TRUE)
  sub("^LI", "", x)
}
canonical_uid <- canonical(uid)
# Optional context for QC only: distinguish IDs excluded at the common stage
# from genuinely unmatched identifiers in external HEP/preprocessing sources.
prior_path <- file.path(deriv, "data.csv")
prior_ids <- if (file.exists(prior_path)) {
  prior <- read_csv(prior_path)
  if ("life_id" %in% names(prior)) unique(ids_norm(prior$life_id)) else character()
} else character()
canonical_prior <- canonical(prior_ids)
skipped <- list()
issues <- list()
issue <- function(source, id, type, detail = "") {
  issues[[length(issues) + 1L]] <<- data.frame(
    source = source, life_id = ifelse(is.na(id), "", id),
    issue = type, detail = detail, stringsAsFactors = FALSE)
}
resolve_ids <- function(x, source) {
  ids <- ids_norm(x)
  ans <- match(ids, uid)
  for (i in which(!is.na(ids) & is.na(ans))) {
    candidate <- which(!is.na(canonical_uid) & canonical_uid == canonical(ids[i]))
    if (length(candidate) == 1L) ans[i] <- candidate
    else if (length(candidate)) {
      issue(source, ids[i], "ambiguous_alias", joined(uid[candidate]))
    } else {
      expected_excluded <- canonical(ids[i]) %in% canonical_prior
      skipped[[length(skipped) + 1L]] <<- data.frame(
        source = source, source_id = ids[i],
        reason = if (expected_excluded) "excluded_before_pipeline_join"
                 else "not_in_data_or_unrecognized_ID",
        stringsAsFactors = FALSE)
      if (!expected_excluded) issue(source, ids[i], "unmatched_ID",
                                   "Source participant not found among generated IDs")
    }
  }
  ans
}
add_column <- function(label, values_by_id) {
  if (label %in% names(data)) stop("Duplicate output column: ", label)
  if (length(values_by_id) != nids) stop("Internal error: unexpected participant vector length")
  data[[label]] <<- values_by_id[match(data$life_id, uid)]
}
coverage <- list()
add_coverage <- function(source, col, matched, available, repeated) {
  coverage[[length(coverage) + 1L]] <<- data.frame(
    source = source, column = col,
    n_matched_IDs = matched, n_IDs_with_value = available,
    n_repeated_source_IDs = repeated, stringsAsFactors = FALSE)
}

# ---- 1. Import the 16 HEP CSVs; keep the original five fields ------------
group_dir <- file.path(deriv, "group")
if (!dir.exists(group_dir)) stop("Missing directory: ", group_dir)
hep_files <- sort(list.files(group_dir, pattern = "roi_HEP[.]csv$",
                             full.names = TRUE, recursive = FALSE,
                             ignore.case = TRUE))
if (length(hep_files) != 16L)
  stop("Expected exactly 16 *_roi_HEP.csv files in ", group_dir,
       "; found ", length(hep_files))
stems <- sub("(?i)_roi_HEP[.]csv$", "", basename(hep_files), perl = TRUE)
if (anyDuplicated(stems) || any(!nzchar(stems))) stop("Repeated or blank pipeline stems")

# The HP filter must come from the filename; never infer it from source values.
detect_hp <- function(stem) {
  s <- tolower(stem)
  hp01 <- grepl("(^|[_-])hp[_-]?(01|0[.]1|0p1)(?=$|[_-])", s, perl = TRUE)
  hp05 <- grepl("(^|[_-])hp[_-]?(05|0[.]5|0p5)(?=$|[_-])", s, perl = TRUE)
  if (xor(hp01, hp05)) {
    if (hp01) "hp_0.1" else "hp_0.5"
  } else {
    NA_character_
  }
}
hp_for <- vapply(stems, detect_hp, character(1L))
if (anyNA(hp_for)) stop("Cannot infer HP (01 vs 05) from: ",
                      joined(stems[is.na(hp_for)], "; "))
if (!all(table(factor(hp_for, levels = c("hp_0.1", "hp_0.5"))) == 8L))
  stop("Expected 8 HEP files per HP level. Found: ",
       paste(names(table(hp_for)), table(hp_for), collapse = ", "))
roi_fields <- c("roi_HEP", "roi_ECG", "n_epochs_total",
                "n_epochs_rejected", "n_epochs_kept")

for (i in seq_along(hep_files)) {
  source <- basename(hep_files[i])
  tab <- read_csv(hep_files[i])
  missing_fields <- setdiff(c("subject_id", roi_fields), names(tab))
  if (length(missing_fields)) stop("Missing source fields in ", source, ": ",
                                   joined(missing_fields, ", "))
  ids <- resolve_ids(tab$subject_id, source)
  matched <- unique(ids[!is.na(ids)])
  duplicate_count <- sum(table(ids[!is.na(ids)]) > 1L)
  for (field in roi_fields) {
    vals <- rep(NA_character_, nids)
    raw <- norm(tab[[field]])
    for (k in matched) {
      candidates <- unique(raw[ids == k])
      candidates <- candidates[!is.na(candidates)]
      if (length(candidates) == 1L) vals[k] <- candidates
      if (length(candidates) > 1L)
        issue(source, uid[k], paste0("conflicting_", field),
              paste("Disagreeing values:", joined(candidates, "; "),
                    "(field left missing)"))
    }
    name <- paste0(stems[i], "_", field)
    add_column(name, vals)
    add_coverage(source, name, length(matched), sum(!is.na(vals)), duplicate_count)
  }
}

# ---- 2. Import HP preprocessing status and retain stage-level evidence ----
if (!requireNamespace("readODS", quietly = TRUE))
  stop("Install readODS: install.packages('readODS')")
pre_path <- file.path(life, "tabular", "preprocessing_logs.ods")
if (!file.exists(pre_path)) stop("Missing preprocessing log: ", pre_path)
pre <- as.data.frame(readODS::read_ods(pre_path, sheet = 1),
                     check.names = FALSE, stringsAsFactors = FALSE)
names(pre) <- trimws(names(pre))
needed_pre <- c("subject_id", "step", "hp_0.1", "hp_0.5")
if (!all(needed_pre %in% names(pre)))
  stop("preprocessing_logs.ods missing columns: ", joined(setdiff(needed_pre, names(pre))))
if (!"detail" %in% names(pre)) pre$detail <- rep(NA_character_, nrow(pre))
pre[] <- lapply(pre, norm)
pre_ids <- resolve_ids(pre$subject_id, basename(pre_path))
pre_steps <- tolower(pre$step)
pre_events <- list()
# Store status/step/history privately for each HP filter. These drive flags
# and QC reports but should not appear in the analysis-wide CSV.
pre_info <- list()
for (hp in c("hp_0.1", "hp_0.5")) {
  status <- stage <- detail <- rep(NA_character_, nids)
  history_stages <- history_details <- rep(NA_character_, nids)
  for (k in unique(pre_ids[!is.na(pre_ids)])) {
    valid <- which(pre_ids == k & !is.na(pre[[hp]]))
    if (!length(valid)) next
    # Mirror the original LIFE_generate_data.R precedence exactly:
    # last nonmissing a5, otherwise last nonmissing row.
    a5 <- valid[!is.na(pre_steps[valid]) & pre_steps[valid] == "a5"]
    chosen <- if (length(a5)) tail(a5, 1L) else tail(valid, 1L)
    status[k] <- pre[[hp]][chosen]
    stage[k] <- pre$step[chosen]
    detail[k] <- pre$detail[chosen]
    # Capture all historical "excluded" stages as audit information.
    failed_rows <- valid[tolower(pre[[hp]][valid]) == "excluded"]
    if (length(failed_rows)) {
      history_stages[k] <- joined(pre$step[failed_rows], "; ")
      history_details[k] <- joined(paste0(
        ifelse(is.na(pre$step[failed_rows]), "<step unknown>", pre$step[failed_rows]),
        ": ", ifelse(is.na(pre$detail[failed_rows]), "<no detail>",
                       pre$detail[failed_rows])), "; ")
      pre_events[[length(pre_events) + 1L]] <- data.frame(
        life_id = rep(uid[k], length(failed_rows)), high_pass = hp,
        ods_row = failed_rows + 1L, step = pre$step[failed_rows],
        status = pre[[hp]][failed_rows], detail = pre$detail[failed_rows],
        selected_as_final = failed_rows == chosen, stringsAsFactors = FALSE)
    }
    if (length(failed_rows) && tolower(status[k]) == "ok")
      issue("preprocessing_logs.ods", uid[k], "earlier_excluded_superseded_by_ok",
            paste(hp, "chosen stage", stage[k]))
  }
  pre_info[[hp]] <- list(status = status, step = stage, detail = detail,
                         excluded_stages = history_stages,
                         excluded_details = history_details)
  out_name <- paste0("preproc_", hp)
  add_coverage("preprocessing_logs.ods", out_name,
               length(unique(pre_ids[!is.na(pre_ids)])), sum(!is.na(status)),
               sum(table(pre_ids[!is.na(pre_ids)]) > 1L))
  observed <- tolower(status[!is.na(status)])
  unexpected <- unique(observed[!observed %in% c("excluded", "ok")])
  if (length(unexpected)) issue("preprocessing_logs.ods", "",
     "unrecognized_preprocessing_status", paste(hp, joined(unexpected, "; ")))
}
write_csv(bind(pre_events,
  c("life_id", "high_pass", "ods_row", "step", "status", "detail",
    "selected_as_final")), "pipeline_preprocessing_stages.csv")

# ---- 3. Per-pipeline rejection flags and reference distributions ----------
counts <- list()
references <- list()
events <- list()
missing_events <- list()
summary_by_anxiety <- list()
meta_group <- if ("gad7_group" %in% names(data)) {
  tolower(norm(data$gad7_group[uid_row]))
} else {
  rep(NA_character_, nids)
}
meta_group[!meta_group %in% c("high", "low")] <- "other_or_unknown"

for (i in seq_along(stems)) {
  stem <- stems[i]
  hp <- hp_for[i]
  hp_col <- paste0("preproc_", hp)
  h_col <- paste0(stem, "_roi_HEP")
  ep_col <- paste0(stem, "_n_epochs_kept")
  h <- num(data[[h_col]][uid_row])
  ep <- num(data[[ep_col]][uid_row])
  ep_total <- num(data[[paste0(stem, "_n_epochs_total")]][uid_row])
  ep_rejected <- num(data[[paste0(stem, "_n_epochs_rejected")]][uid_row])
  inconsistent <- !is.na(ep) & !is.na(ep_total) & !is.na(ep_rejected) &
    abs(ep_total - ep_rejected - ep) > 1e-8
  if (any(inconsistent)) for (k in which(inconsistent))
    issue(basename(hep_files[i]), uid[k], "epoch_accounting_mismatch",
          paste("total=", ep_total[k], "rejected=", ep_rejected[k],
                "kept=", ep[k]))
  invalid_ep <- !is.na(ep) & (ep < 0 | ep != floor(ep))
  if (any(invalid_ep)) {
    for (k in which(invalid_ep)) issue(basename(hep_files[i]), uid[k],
      "invalid_epoch_count", paste(ep_col, ep[k], "(left unknown)"))
    ep[invalid_ep] <- NA_real_
  }
  status <- tolower(norm(pre_info[[hp]]$status))
  pre_fail <- !is.na(status) & status == "excluded"
  ep_fail <- !is.na(ep) & ep < 300
  h_missing <- is.na(h)
  ep_missing <- is.na(ep)
  status_unknown <- is.na(status) | !status %in% c("excluded", "ok")
  # Exclude only with positive, observed evidence. Missing values remain QC.
  reference <- !pre_fail & !ep_fail & !h_missing
  nh <- sum(reference)
  mu <- if (nh) mean(h[reference]) else NA_real_
  sd_h <- if (nh > 1L) sd(h[reference]) else NA_real_
  z <- if (is.finite(sd_h) && sd_h > 0) (h - mu) / sd_h else
    rep(NA_real_, nids)
  h_outlier <- !is.na(z) & abs(z) > 3  # prereg: strictly outside +/- 3 SD
  excluded <- pre_fail | ep_fail | h_outlier
  # Three analysis fields only. Other individual reasons and missing-data
  # flags remain available in pipeline_exclusion_events.csv and
  # pipeline_missing_data.csv, plus the aggregate QC summary files.
  add_column(paste0(stem, "_roi_HEP_SD"), z)
  add_column(paste0(stem, "_roi_HEP_outliers"),
             ifelse(is.na(z), NA_integer_, as.integer(h_outlier)))
  add_column(paste0(stem, "_exclude"), as.integer(excluded))

  # Each event appears once per participant rather than once per unresolved row.
  add_event <- function(criterion, source_col, flag, values, explanation = NULL) {
    ks <- which(flag)
    if (!length(ks)) return(invisible(NULL))
    if (is.null(explanation)) explanation <- rep("", nids)
    events[[length(events) + 1L]] <<- data.frame(
      pipeline = stem, life_id = uid[ks], criterion = criterion,
      source_column = source_col, value = as.character(values[ks]),
      detail = as.character(explanation[ks]), stringsAsFactors = FALSE)
  }
  add_event("preprocessing_excluded", hp_col, pre_fail, status,
            paste0("selected_step=", pre_info[[hp]]$step,
                   "; history=", pre_info[[hp]]$excluded_details))
  add_event("fewer_than_300_epochs", ep_col, ep_fail, ep)
  add_event("roi_HEP_outside_3SD", h_col, h_outlier, h,
            paste0("z=", signif(z, 6)))
  missing_flags <- list(roi_HEP_missing = h_missing, epochs_missing = ep_missing,
                        preprocessing_unknown = status_unknown)
  for (kind in names(missing_flags)) {
    ks <- which(missing_flags[[kind]])
    if (length(ks)) missing_events[[length(missing_events) + 1L]] <- data.frame(
      pipeline = stem, life_id = uid[ks], missing_type = kind,
      source_column = switch(kind, roi_HEP_missing = h_col,
        epochs_missing = ep_col, preprocessing_unknown = hp_col),
      stringsAsFactors = FALSE)
  }
  if (any(ep_fail & !pre_fail))
    issue(basename(pre_path), "", "epochs_below_300_but_preprocessing_not_excluded",
          paste(stem, sum(ep_fail & !pre_fail), "IDs"))
  if (!is.finite(sd_h) || sd_h <= 0)
    issue(basename(hep_files[i]), "", "cannot_estimate_HEP_SD",
          paste(stem, "n_reference=", nh, "sd=", sd_h))
  references[[i]] <- data.frame(pipeline = stem, high_pass = hp,
    n_reference_IDs = nh, mean_roi_HEP = mu, sd_roi_HEP = sd_h,
    cutoff_lower = mu - 3 * sd_h, cutoff_upper = mu + 3 * sd_h)
  counts[[i]] <- data.frame(pipeline = stem, high_pass = hp,
    n_total_IDs = nids, n_preproc_excluded = sum(pre_fail),
    n_below_300_epochs = sum(ep_fail), n_roi_HEP_outliers = sum(h_outlier),
    n_any_pipeline_exclusion = sum(excluded),
    n_available_for_analysis = sum(!excluded & !h_missing & !ep_missing),
    n_HEP_missing = sum(h_missing), n_epoch_count_missing = sum(ep_missing),
    n_epoch_accounting_mismatch = sum(inconsistent),
    n_preproc_status_unknown = sum(status_unknown),
    n_preproc_and_epochs = sum(pre_fail & ep_fail),
    n_preproc_and_outlier = sum(pre_fail & h_outlier),
    n_epochs_and_outlier = sum(ep_fail & h_outlier))
  for (group in c("high", "low", "other_or_unknown")) {
    g <- meta_group == group
    summary_by_anxiety[[length(summary_by_anxiety) + 1L]] <- data.frame(
      pipeline = stem, gad7_group = group, n_group = sum(g),
      n_excluded = sum(g & excluded),
      n_analysis_available = sum(g & !excluded & !h_missing & !ep_missing))
  }
}

# ---- 4. QC reports and unresolved duplicate-ID navigation ------------------
write_csv(do.call(rbind, counts), "pipeline_exclusion_counts.csv")
write_csv(do.call(rbind, references), "pipeline_HEP_reference_stats.csv")
write_csv(do.call(rbind, summary_by_anxiety), "pipeline_anxiety_group_counts.csv")
write_csv(bind(events, c("pipeline", "life_id", "criterion", "source_column",
                         "value", "detail")), "pipeline_exclusion_events.csv")
write_csv(bind(missing_events, c("pipeline", "life_id", "missing_type",
                                 "source_column")), "pipeline_missing_data.csv")
write_csv(bind(coverage, c("source", "column", "n_matched_IDs",
                            "n_IDs_with_value", "n_repeated_source_IDs")),
          "pipeline_matching_coverage.csv")
write_csv(bind(issues, c("source", "life_id", "issue", "detail")),
          "pipeline_matching_issues.csv")
write_csv(bind(skipped, c("source", "source_id", "reason")),
          "pipeline_skipped_source_IDs.csv")

repeated <- uid[table(data$life_id)[uid] > 1L]
nav <- lapply(repeated, function(id) {
  pos <- which(data$life_id == id)
  one <- pos[1L]
  rejected <- stems[vapply(stems, function(stem)
    data[[paste0(stem, "_exclude")]][one] == 1L, logical(1L))]
  data.frame(life_id = id, n_rows = length(pos),
    data_with_pipeline_flags_rows = joined(pos + 1L),
    data_after_exclusion_csv_rows = joined(row_ix[pos]),
    pipelines_with_exclusion = joined(rejected, "; "),
    n_pipelines_with_exclusion = length(rejected),
    stringsAsFactors = FALSE)
})
write_csv(bind(nav, c("life_id", "n_rows", "data_with_pipeline_flags_rows",
                      "data_after_exclusion_csv_rows", "pipelines_with_exclusion",
                      "n_pipelines_with_exclusion")), "pipeline_nonunique_IDs.csv")
if (nrow(data) != length(row_ix) || length(unique(data$life_id)) != nids)
  stop("Internal consistency error: participant rows changed during pipeline join")

output_path <- file.path(deriv, "data_with_pipeline_flags.csv")
write.csv(data, output_path, na = "", row.names = FALSE)
report <- capture.output({
  cat("========== LIFE pipeline HEP matching and exclusion flags =========\n")
  cat("Input:  ", input_path, "\n", sep = "")
  cat("Output: ", output_path, "\n", sep = "")
  cat("Rows retained (unchanged): ", nrow(data), "\n", sep = "")
  cat("Unique LIFE IDs:           ", nids, "\n", sep = "")
  cat("Non-unique LIFE IDs:       ", length(repeated), "\n", sep = "")
  cat("HEP pipelines matched:     ", length(stems), " (8 for each HP)\n", sep = "")
  cat("No overall_exclude recalculated; no data rows removed.\n")
  cat("Per pipeline output: 5 input fields + HEP z, HEP outlier, combined exclude.\n")
  cat("Detailed exclusion/missing/preprocessing information is QC-only.\n")
  cat("Exclusion conditions: preprocessing='excluded', epochs<300, |HEP z|>3.\n")
  cat("Preprocessing 'ok'/missing: not excluded. Missing HEP/epochs: QC only.\n")
  cat("HEP z-score references computed on post-common-exclusion participants.\n")
  cat("\nPer-pipeline participant counts:\n")
  print(do.call(rbind, counts), row.names = FALSE)
  cat("\nMatching issues: ", length(issues), "\n", sep = "")
  cat("Diagnostics: ", out, "\n", sep = "")
})
writeLines(report, file.path(out, "pipeline_summary.txt"))
cat(paste(report, collapse = "\n"), "\n")
