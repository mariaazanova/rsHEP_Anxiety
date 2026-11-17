#!/usr/bin/env Rscript
# rsHEP / LIFE: preregistered exclusion audit, AFTER 04_LIFE_generate_data.R
# Usage: Rscript 05_exclusion.R [life_directory]
# Default life_directory: /data/p_02865/data/LIFE
# Requires: readODS (install.packages("readODS"))
# Input:  derivatives/data.csv, tabular/preprocessing_logs.ods,
#         derivatives/matching_and_summaries/medications/{LIFE_medications_by_ID.csv,
#         medications_cns_yes_qc.csv}
# Output: derivatives/data_after_exclusion.csv
#         derivatives/matching_and_summaries/exclusion/*
#
# IMPORTANT DESIGN CHOICES:
# - Each participant is dropped from the whole multiverse only if they fail a
#   common exclusion criterion OR all 16 pipelines. One bad pipeline is not a
#   whole-participant exclusion. Downstream analyses MUST filter on the selected
#   <pipeline>_exclude flag, not only on overall_exclude.
# - Multiple unresolved rows for the same life_id are not silently adjudicated:
#   a positive exclusion on ANY row is propagated to ALL rows of that ID.
# - We NEVER interpret missing clinical disease history as a positive response.
#   Unexpected nonempty codes are reported for review but do not block export.
# - Clinical yes/no coding: 1=yes, 0=no, 98=unknown. Blank/NA is also
#   unknown, never a positive exclusion. Any other value is flagged for review,
#   not silently treated as "no". All observed values are audited.
# - No additional exclusion for missing control covariates (IDS, BP, BMI, etc.);
#   these are handled within the respective control analyses

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: Rscript 05_exclusion.R [life_directory]")
life <- if (length(args)) args[1L] else "/data/p_02865/data/LIFE"
deriv <- file.path(life, "derivatives")
source_dir <- file.path(deriv, "matching_and_summaries")
out_dir <- file.path(source_dir, "exclusion")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Remaining clinical variables needing separate assessment:
# ecg_rhythm, ecg_other_rhythm, ecg_sa_block, ecg_av_block,
# herzkatheter, herzkranzgefaess_oder_stent.

read_csv <- function(path) {
  if (!file.exists(path)) stop("Required input missing: ", path)
  x <- read.csv(path, colClasses = "character", check.names = FALSE,
                stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM", na.strings = "")
  if (anyDuplicated(names(x))) stop("Duplicate headers in ", path)
  x
}
save_csv <- function(x, filename) {
  write.csv(x, file.path(out_dir, filename), row.names = FALSE, na = "")
}
norm <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | !nzchar(x) | toupper(x) %in% c("NA", "N/A", "NULL")] <- NA_character_
  x
}
ids_norm <- function(x) toupper(norm(x))
join <- function(x, sep = " | ") {
  x <- unique(norm(x))
  paste(x[!is.na(x)], collapse = sep)
}
number <- function(x) suppressWarnings(as.numeric(gsub(",", ".", norm(x), fixed = TRUE)))
empty <- function(cols) as.data.frame(setNames(rep(list(character()), length(cols)), cols),
                                        check.names = FALSE, stringsAsFactors = FALSE)
bind <- function(parts, cols) {
  if (!length(parts)) return(empty(cols))
  x <- do.call(rbind, parts)
  rownames(x) <- NULL
  x[, cols, drop = FALSE]
}

input <- read_csv(file.path(deriv, "data.csv"))
if (!nrow(input)) stop("data.csv contains zero rows")
if (!"life_id" %in% names(input)) stop("data.csv requires life_id")
input$life_id <- ids_norm(input$life_id)
if (anyNA(input$life_id)) stop("data.csv has empty life_id")
n <- nrow(input)
row_number <- seq_len(n) + 1L  # physical CSV line, including header
uid <- unique(input$life_id)
groups <- split(seq_len(n), factor(input$life_id, levels = uid))
first <- !duplicated(input$life_id)
if ("non-unique" %in% names(input) &&
    !identical(as.integer(number(input[["non-unique"]])),
               as.integer(table(input$life_id)[input$life_id] > 1L)))
  stop("Input non-unique column is not consistent with repeated life_id")

# A row-level yes in any repeated source row excludes the whole LIFE ID.
by_id_any <- function(flag) {
  stopifnot(length(flag) == n)
  ans <- rep(FALSE, n)
  for (g in groups) ans[g] <- any(flag[g], na.rm = TRUE)
  ans
}
by_id_missing_any <- function(flag) {
  ans <- rep(FALSE, n)
  for (g in groups) ans[g] <- any(flag[g])
  ans
}

# Record questionable source coding, contradictions, and missing mappings.
issues <- list()
add_issue <- function(type, source, detail, life_id = "") {
  issues[[length(issues) + 1L]] <<- data.frame(
    issue_type = type, life_id = life_id, source = source, detail = detail,
    stringsAsFactors = FALSE)
}
code_audit <- list()
unknown_cells <- list()
record_codes <- function(source, values, parsed, unknown = c("98")) {
  observed <- toupper(norm(values))
  tab <- table(ifelse(is.na(observed), "<missing>", observed))
  for (v in names(tab)) {
    label <- if (v == "<missing>") "missing_unknown" else if (v %in% toupper(unknown))
      "coded_unknown" else {
        p <- unique(parsed[!is.na(observed) & observed == v])
        if (length(p) == 1L && !is.na(p)) if (p) "yes" else "no" else "UNRECOGNIZED"
      }
    code_audit[[length(code_audit) + 1L]] <<- data.frame(
      source_column = source, value = v, n_rows = as.integer(tab[[v]]),
      interpretation = label, stringsAsFactors = FALSE)
  }
}

# The default for exclusion-source binary variables has been confirmed:
# 1=yes, 0=no, 98=unknown, blank/NA=unknown. A different nonempty code needs
# explicit per-variable adjudication; notably, 2 is NOT implicitly "no".
YES <- c("1")
NO <- c("0")
UNKNOWN <- c("98")
parse_binary <- function(x, col, yes = YES, no = NO, unknown = UNKNOWN) {
  y <- toupper(norm(x))
  out <- rep(NA, length(y))
  out[!is.na(y) & y %in% toupper(yes)] <- TRUE
  out[!is.na(y) & y %in% toupper(no)] <- FALSE
  is_unknown <- is.na(y) | y %in% toupper(unknown)
  bad <- !is_unknown & is.na(out)
  if (any(bad)) {
    add_issue("unrecognized_binary_code", col,
              paste(paste(unique(y[bad]), collapse = "; "),
                    "(rows:", sum(bad), ")"))
  }
  missing_or_unknown <- which(is_unknown)
  if (length(missing_or_unknown)) {
    unknown_cells[[length(unknown_cells) + 1L]] <<- data.frame(
      life_id = input$life_id[missing_or_unknown],
      data_csv_row = row_number[missing_or_unknown],
      source_column = col,
      source_value = ifelse(is.na(y[missing_or_unknown]), "", y[missing_or_unknown]),
      unknown_type = ifelse(is.na(y[missing_or_unknown]), "missing", "coded_unknown"),
      stringsAsFactors = FALSE)
  }
  record_codes(col, x, out, unknown)
  out
}
require_col <- function(col) {
  if (!col %in% names(input)) stop("Missing required data.csv column: ", col)
}

# Criterion events capture EXACT source columns, original values and optional
# clinical/preprocessing explanations. One event per positive input row/column.
events <- list()
add_events <- function(criterion, source, positive, category = "common",
                       pipeline = "", detail = NULL) {
  if (!any(positive)) return(invisible(NULL))
  ix <- which(positive)
  if (is.null(detail)) detail <- rep("", n)
  if (!source %in% names(input)) values <- rep("", n) else values <- input[[source]]
  events[[length(events) + 1L]] <<- data.frame(
    life_id = input$life_id[ix], data_csv_row = row_number[ix],
    type = category, pipeline = pipeline, criterion = criterion,
    source_column = source, source_value = ifelse(is.na(values[ix]), "", values[ix]),
    detail = ifelse(is.na(detail[ix]), "", detail[ix]),
    stringsAsFactors = FALSE)
}

# ---- 1. Pipeline-INDEPENDENT criteria (preregistration section 2.2.4) ------
# Keep individual criterion decisions to diagnose overlap and sequential loss.
common <- list()
register_common <- function(name, source, positive, detail = NULL) {
  # Store the triggering ORIGINAL cells in the audit; propagate the DECISION
  # to all unresolved rows of the same participant separately.
  triggered_cells <- positive %in% TRUE
  participant_positive <- by_id_any(triggered_cells)
  common[[name]] <<- if (name %in% names(common))
    common[[name]] | participant_positive else participant_positive
  add_events(name, source, triggered_cells, "common", detail = detail)
}

# GAD-7 must have all seven valid items. This is *not* GAD7 <= 4/ >=10:
# scores 5-9 remain in the continuous-anxiety analyses.
require_col("all_GAD7")
gad_complete <- parse_binary(input$all_GAD7, "all_GAD7", yes = c("1"), no = c("0"))
# all_GAD7=0 specifically establishes invalidity; an absent flag remains
# unknown rather than becoming an automatic exclusion (original GAD items
# were assessed upstream by LIFE_generate_data.R).
invalid_gad <- gad_complete %in% FALSE
register_common("invalid_GAD7", "all_GAD7", invalid_gad,
                detail = rep("less than seven valid GAD-7 items", n))

require_col("eeg_alcohol_today")
# EEG alcohol: 0=no, 1=yes, 98=unknown; unknown never excludes.
alcohol <- parse_binary(input$eeg_alcohol_today, "eeg_alcohol_today")
register_common("alcohol_on_EEG_day", "eeg_alcohol_today", alcohol)

# Known brain conditions likely to affect EEG/ECG. The data owner confirmed
# that meningitis, head trauma, and skull fracture are NOT exclusion criteria.
neurological <- c("Intrakranielle_Blutung", "head_brain_op", "parkinson",
                  "epilepsy", "multiple_sclerosis", "schlaganfall")
for (col in neurological) {
  require_col(col)
  p <- parse_binary(input[[col]], col)
  register_common("neurological_history", col, p)
}

# Cardiac history: named history, echo_arrhythmia, and both implanted
# pacemaker fields (explicitly confirmed by the data owner).
# Other rhythm categories with ambiguous coding are left for separate QC.
cardiac <- c("infarct", "bypass", "rhythm_disturbances",
             "coronary_disease_angina_pectoris", "echo_arrhythmia",
             "ecg_pacemaker_atrial", "ecg_pacemaker_chamber")
for (col in cardiac) {
  require_col(col)
  p <- parse_binary(input[[col]], col)
  register_common("cardiac_history", col, p)
}

require_col("medications_cns_exclude_definite")
med_positive <- parse_binary(input$medications_cns_exclude_definite,
                             "medications_cns_exclude_definite")

# ---- 2. Medication category details (not a second exclusion decision) -----
med_path <- file.path(source_dir, "medications", "LIFE_medications_by_ID.csv")
qc_path <- file.path(source_dir, "medications", "medications_cns_yes_qc.csv")
meds <- read_csv(med_path)
qc <- read_csv(qc_path)
if (!all(c("life_id", "cns_yes_atc_codes") %in% names(meds)))
  stop("LIFE_medications_by_ID.csv needs life_id and cns_yes_atc_codes")
if (!all(c("atc_code", "cns_final_reason") %in% names(qc)))
  stop("medications_cns_yes_qc.csv needs atc_code and cns_final_reason")
meds$life_id <- ids_norm(meds$life_id)
qc$atc_code <- toupper(norm(qc$atc_code))
qc$cns_final_reason <- norm(qc$cns_final_reason)
# Generation used exact match, then unique canonical match; preserve the same
# contract here, leaving ambiguous canonical IDs unresolved.
canon <- function(x) {
  x <- ids_norm(x)
  x <- sub("(?i)_(Ruhe|Pilot)(?:_[0-9]+)?$", "", x, perl = TRUE)
  x <- sub("[_-]+2$", "", x)
  x <- sub("[_-]+$", "", x)
  sub("^LI", "", x)
}
resolve_id <- function(x, source) {
  raw <- ids_norm(x)
  result <- raw
  for (i in which(!is.na(raw) & !raw %in% uid)) {
    candidates <- uid[canon(uid) == canon(raw[i]) & !is.na(canon(uid))]
    if (length(candidates) == 1L) result[i] <- candidates
    else {
      result[i] <- NA_character_
      add_issue(if (length(candidates)) "ambiguous_identifier" else "unmatched_identifier",
                source, raw[i])
    }
  }
  result[is.na(raw) | !(result %in% uid)] <- NA_character_
  result
}
meds$matched_life_id <- resolve_id(meds$life_id, basename(med_path))
extract_atc <- function(text) {
  x <- toupper(join(text, sep = " "))
  hits <- regmatches(x, gregexpr("[A-Z][0-9]{2}[A-Z]{2}[0-9]{2}", x,
                                    perl = TRUE))[[1L]]
  unique(hits)
}
med_detail_by_id <- setNames(rep("", length(uid)), uid)
unmatched_codes <- list()
for (id in uid) {
  raw_codes <- meds$cns_yes_atc_codes[!is.na(meds$matched_life_id) &
                                     meds$matched_life_id == id]
  codes <- extract_atc(raw_codes)
  if (!length(codes)) next
  items <- character()
  for (atc in codes) {
    reason <- unique(qc$cns_final_reason[!is.na(qc$atc_code) & qc$atc_code == atc])
    reason <- reason[!is.na(reason)]
    if (!length(reason)) {
      unmatched_codes[[length(unmatched_codes) + 1L]] <- data.frame(
        life_id = id, atc_code = atc, issue = "ATC has no cns_final_reason",
        stringsAsFactors = FALSE)
      reason <- "category UNMAPPED"
    }
    if (length(reason) > 1L)
      add_issue("multiple_ATC_categories", "medications_cns_yes_qc.csv",
                paste(atc, join(reason)), id)
    items <- c(items, paste0(atc, " [", join(reason, "; "), "]"))
  }
  med_detail_by_id[id] <- join(items, "; ")
}
med_details <- unname(med_detail_by_id[input$life_id])
for (id in uid) {
  idx <- which(input$life_id == id)
  pos <- any(med_positive[idx], na.rm = TRUE)
  codes_present <- nzchar(med_detail_by_id[id])
  if (pos && !codes_present)
    add_issue("positive_CNS_medication_missing_ATC", basename(med_path),
              "exclusion is retained, but individual ATC categories unavailable", id)
  if (!pos && codes_present)
    add_issue("CNS_ATC_without_positive_flag", basename(med_path),
              "confirm cns_exclude_definite versus cns_yes_atc_codes", id)
}
register_common("CNS_medications", "medications_cns_exclude_definite",
                med_positive, detail = med_details)

# Precomputed, preregistered extrasystole threshold: count > mean + 4 SD.
require_col("extrasystolic_exclusion")
extra <- parse_binary(input$extrasystolic_exclusion, "extrasystolic_exclusion",
                      yes = "1", no = "0")
register_common("extrasystoles_above_4SD", "extrasystolic_exclusion", extra,
                detail = if ("extrasystolic_N_SD" %in% names(input))
                  paste0("z=", input$extrasystolic_N_SD) else NULL)

# VIGALL B+C >900 one-second epochs, not >=900.
require_col("drowsy_15min")
drowsy <- parse_binary(input$drowsy_15min, "drowsy_15min",
                       yes = "1", no = "0")
register_common("drowsy_over_15min", "drowsy_15min", drowsy,
                detail = if (all(c("vigall_b", "vigall_c") %in% names(input)))
                  paste0("vigall_b=", input$vigall_b, ", vigall_c=", input$vigall_c)
                else NULL)

# ---- 3. Read preprocessing log and preserve stage + detail ----------------
if (!requireNamespace("readODS", quietly = TRUE))
  stop("Install readODS before running: install.packages('readODS')")
pre_path <- file.path(life, "tabular", "preprocessing_logs.ods")
if (!file.exists(pre_path)) stop("Missing: ", pre_path)
pre <- as.data.frame(readODS::read_ods(pre_path, sheet = 1),
                     check.names = FALSE, stringsAsFactors = FALSE)
names(pre) <- trimws(names(pre))
if (!all(c("subject_id", "step", "hp_0.1", "hp_0.5") %in% names(pre)))
  stop("preprocessing_logs.ods needs subject_id, step, hp_0.1, hp_0.5")
pre[] <- lapply(pre, norm)
pre$life_id <- resolve_id(pre$subject_id, basename(pre_path))
pre$step <- tolower(pre$step)
pre$detail_for_audit <- if ("detail" %in% names(pre)) pre$detail else
  if ("details" %in% names(pre)) pre$details else rep(NA_character_, nrow(pre))
if (!"detail" %in% names(pre) && !"details" %in% names(pre))
  add_issue("missing_optional_detail_column", basename(pre_path),
            "Neither detail nor details exists; stage is still reported")

# 'excluded' may appear at the stage where a subject failed; otherwise an
# empty hp cell alone is not positive evidence of exclusion. A subject WITH
# log rows and no excluded event is provisionally treated as non-excluded.
pre_status <- function(x, col) {
  y <- tolower(norm(x))
  out <- rep(NA, length(y))
  out[!is.na(y) & grepl("^(excluded|rejected|failed|fail)(\\b|$)", y, perl = TRUE)] <- TRUE
  out[!is.na(y) & y %in% c("included", "not excluded", "passed", "pass",
                                    "ok", "valid", "false", "0")] <- FALSE
  bad <- !is.na(y) & is.na(out)
  if (any(bad))
    add_issue("unrecognized_preprocessing_status", col, join(unique(y[bad]), "; "))
  out
}
pre_events <- list()
pre_flag_by_hp <- list()
pre_detail_by_hp <- list()
for (hp in c("hp_0.1", "hp_0.5")) {
  h <- pre_status(pre[[hp]], paste0("preprocessing_logs.ods:", hp))
  positives <- which(!is.na(h) & h & !is.na(pre$life_id))
  for (j in positives) {
    pre_events[[length(pre_events) + 1L]] <- data.frame(
      life_id = pre$life_id[j], hp = hp, step = ifelse(is.na(pre$step[j]), "", pre$step[j]),
      detail = ifelse(is.na(pre$detail_for_audit[j]), "", pre$detail_for_audit[j]),
      log_ods_row = j + 1L, status = pre[[hp]][j], stringsAsFactors = FALSE)
  }
  f <- setNames(rep(NA, length(uid)), uid)
  descriptions <- setNames(rep("", length(uid)), uid)
  for (id in uid) {
    at <- which(!is.na(pre$life_id) & pre$life_id == id)
    if (length(at)) {
      # Any stage exclusion takes precedence over any subsequent summary cell.
      f[id] <- any(h[at], na.rm = TRUE)
      ev <- at[!is.na(h[at]) & h[at]]
      if (length(ev)) {
        descriptions[id] <- join(vapply(ev, function(j) {
          stage <- if (is.na(pre$step[j])) "unknown-step" else pre$step[j]
          detail <- pre$detail_for_audit[j]
          if (is.na(detail)) stage else paste0(stage, " (", detail, ")")
        }, character(1L)), sep = "; ")
      }
    }
  }
  pre_flag_by_hp[[hp]] <- unname(f[input$life_id])
  pre_detail_by_hp[[hp]] <- unname(descriptions[input$life_id])
  summary_name <- paste0("preproc_", hp)
  require_col(summary_name)
  final_summary <- pre_status(input[[summary_name]], summary_name)
  for (id in uid) {
    at <- which(input$life_id == id)
    if (any(final_summary[at], na.rm = TRUE) && !isTRUE(f[id])) {
      # Summary can preserve an exclusion not visible in a corrupted/incomplete log.
      pre_flag_by_hp[[hp]][at] <- TRUE
      pre_detail_by_hp[[hp]][at] <- "summary=excluded (stage not available)"
      add_issue("preprocessing_summary_without_log_event", summary_name,
                "summary indicates excluded; detailed log lacks excluded event", id)
    }
    if (isTRUE(f[id]) && any(!is.na(final_summary[at]) & !final_summary[at]))
      add_issue("preprocessing_summary_disagrees_with_log", summary_name,
                "earlier stage says excluded; later summary says included", id)
  }
}
pre_log_events <- bind(pre_events,
                       c("life_id", "hp", "step", "detail", "log_ods_row", "status"))
save_csv(pre_log_events, "preprocessing_exclusion_events.csv")

# ---- 4. Per-pipeline flags: common + matching high-pass + ROI HEP outlier --
# Derive pipeline names from the real columns created in LIFE_generate_data.R,
# not from manually guessed ICA/baseline naming.
roi_flags <- grep("_roi_HEP_outliers$", names(input), value = TRUE)
if (length(roi_flags) != 16L)
  stop("Expected 16 '*_roi_HEP_outliers' columns, got ", length(roi_flags))
stems <- sub("_roi_HEP_outliers$", "", roi_flags)
if (anyDuplicated(stems)) stop("Repeated pipeline names")

# Detect hp_01 and hp_05, allowing hp_0.1/hp_0.5, hp01/hp05 and 0p1/0p5.
detect_hp <- function(stem) {
  s <- tolower(stem)
  h01 <- grepl("(^|[_-])hp[_-]?(01|0[.]1|0p1)(?=$|[_-])", s, perl = TRUE)
  h05 <- grepl("(^|[_-])hp[_-]?(05|0[.]5|0p5)(?=$|[_-])", s, perl = TRUE)
  if (xor(h01, h05)) return(if (h01) "hp_0.1" else "hp_0.5")
  NA_character_
}
assigned_hp <- vapply(stems, detect_hp, character(1L))
if (anyNA(assigned_hp))
  stop("Cannot infer unique high-pass filter for pipeline stems: ",
       join(stems[is.na(assigned_hp)], "; "),
       ". Update detect_hp() using the real filenames.")
if (!all(table(factor(assigned_hp, levels = c("hp_0.1", "hp_0.5"))) == 8L))
  stop("Expected 8 preprocessing pipelines for each high-pass filter")

common_exclude <- Reduce(`|`, common)
# As set above, each common flag already propagates over all rows of an ID.
input$common_exclude <- as.integer(common_exclude)
# Build participant-level source reasons, not just one row's partial evidence.
source_reason <- function(e) {
  if (!nzchar(e$detail)) return(e$source_column)
  paste0(e$source_column, " [", e$detail, "]")
}
all_event_cols <- c("life_id", "data_csv_row", "type", "pipeline", "criterion",
                    "source_column", "source_value", "detail")
common_event_table <- bind(events, all_event_cols)
common_reason_by_id <- setNames(rep("", length(uid)), uid)
for (id in uid) {
  ee <- common_event_table[common_event_table$life_id == id, , drop = FALSE]
  if (nrow(ee)) common_reason_by_id[id] <- join(vapply(seq_len(nrow(ee)),
                                                     function(i) source_reason(ee[i, ]),
                                                     character(1L)))
}
input$exclusion_reasons <- unname(common_reason_by_id[input$life_id])

pipe_matrix <- matrix(FALSE, nrow = n, ncol = length(stems),
                      dimnames = list(NULL, stems))
pipe_event_list <- list()
pre_log_missing_reported <- character()
for (i in seq_along(stems)) {
  stem <- stems[i]
  hp <- assigned_hp[i]
  pc <- paste0("preproc_", hp)
  pre_excluded <- pre_flag_by_hp[[hp]]
  # No matching preprocessing record is a *review* item, not a silently passed
  # subject; subsequent missing HEP is independently an unusable pipeline.
  if (anyNA(pre_excluded) && !(pc %in% pre_log_missing_reported)) {
    add_issue("no_preprocessing_log", pc,
              paste(sum(is.na(pre_excluded[first])), "unique IDs have no matching log"))
    pre_log_missing_reported <- c(pre_log_missing_reported, pc)
  }
  hep_col <- paste0(stem, "_roi_HEP")
  epochs_col <- paste0(stem, "_n_epochs_kept")
  require_col(hep_col)
  require_col(epochs_col)
  outlier <- parse_binary(input[[roi_flags[i]]], roi_flags[i], yes = "1", no = "0")
  missing_hep <- !is.finite(number(input[[hep_col]]))
  epochs <- number(input[[epochs_col]])
  too_few_epochs <- !is.na(epochs) & epochs < 300
  # Missing HEP makes this pipeline unanalyzable irrespective of missing log.
  # Check the preregistered 300-epoch cutoff independently to catch any
  # discrepancy with the imported preprocessing summary.
  # Do not recalculate the preregistered 3-SD outlier flags.
  pipe_matrix[, i] <- by_id_any(common_exclude | (pre_excluded %in% TRUE) |
                               (outlier %in% TRUE) | missing_hep | too_few_epochs)

  # Record each distinct pipeline-specific source as its own event.
  ev_add <- function(crit, source, pos, desc = NULL) {
    # Log the source cell that actually failed, before participant propagation.
    pos <- pos %in% TRUE
    if (any(pos)) {
      ix <- which(pos)
      val <- if (source %in% names(input)) input[[source]][ix] else rep("", length(ix))
      if (is.null(desc)) desc <- rep("", n)
      pipe_event_list[[length(pipe_event_list) + 1L]] <<- data.frame(
        life_id = input$life_id[ix], data_csv_row = row_number[ix], type = "pipeline",
        pipeline = stem, criterion = crit, source_column = source,
        source_value = ifelse(is.na(val), "", val),
        detail = ifelse(is.na(desc[ix]), "", desc[ix]),
        stringsAsFactors = FALSE)
    }
  }
  ev_add("preprocessing_failed", pc, pre_excluded %in% TRUE,
         pre_detail_by_hp[[hp]])
  ev_add("ROI_HEP_outlier_over_3SD", roi_flags[i], outlier %in% TRUE,
         if (paste0(stem, "_roi_HEP_SD") %in% names(input))
           paste0("z=", input[[paste0(stem, "_roi_HEP_SD")]]) else NULL)
  ev_add("ROI_HEP_unavailable", hep_col, missing_hep)
  ev_add("fewer_than_300_epochs", epochs_col, too_few_epochs)
  possible_mismatch <- too_few_epochs & !(pre_excluded %in% TRUE)
  if (any(possible_mismatch))
    add_issue("epoch_count_disagrees_with_preproc_flag", epochs_col,
              paste(sum(by_id_any(possible_mismatch)[first]),
                    "IDs have <300 epochs but are not marked excluded by hp log"))
  if (any(is.na(epochs) & !missing_hep))
    add_issue("HEP_available_but_epochs_missing", epochs_col,
              paste(sum(by_id_any(is.na(epochs) & !missing_hep)[first]),
                    "unique IDs have a valid HEP but unknown epoch count"))
  input[[paste0(stem, "_exclude")]] <- as.integer(pipe_matrix[, i])
}
pipe_event_table <- bind(pipe_event_list, all_event_cols)
event_table <- rbind(common_event_table, pipe_event_table)
save_csv(event_table, "exclusion_events.csv")

# Compact per-pipeline reasons, repeated identically across unresolved rows.
pipe_reason_by_id <- setNames(rep("", length(uid)), uid)
for (id in uid) {
  e <- pipe_event_table[pipe_event_table$life_id == id, , drop = FALSE]
  if (!nrow(e)) next
  phrases <- unique(paste0(e$pipeline, ": ", vapply(seq_len(nrow(e)),
                                                   function(j) source_reason(e[j, ]),
                                                   character(1L))))
  pipe_reason_by_id[id] <- join(phrases)
}
input$pipeline_exclusion_reasons <- unname(pipe_reason_by_id[input$life_id])

# Overall exclusion = failed common criteria OR failed EVERY multiverse branch.
# A participant who remains eligible in even one branch stays in the final
# table; downstream models apply <pipeline>_exclude for their own branch.
all_pipelines_failed <- rowSums(pipe_matrix) == ncol(pipe_matrix)
overall <- common_exclude | all_pipelines_failed
input$overall_exclude <- as.integer(overall)
# Include common source reasons and a terse collective branch explanation.
input$overall_exclusion_reasons <- ifelse(
  overall,
  ifelse(common_exclude, input$exclusion_reasons,
         "all_16_pipelines_excluded: see pipeline_exclusion_reasons"), "")

# Preserve a stable and human-readable column order, independent of input names.
new_cols <- c("common_exclude", "exclusion_reasons", paste0(stems, "_exclude"),
              "pipeline_exclusion_reasons", "overall_exclude",
              "overall_exclusion_reasons")
input <- input[, c(setdiff(names(input), new_cols), new_cols), drop = FALSE]

# ---- 5. Diagnostic counts: individual, sequential, overlap, pipeline ------
# Always count *both* physical data rows and independent participants.
meta <- data.frame(life_id = uid, group = NA_character_, stringsAsFactors = FALSE)
if ("gad7_group" %in% names(input)) {
  g <- tolower(norm(input$gad7_group[first]))
  for (i in seq_along(uid)) {
    if (g[i] %in% c("high", "low")) meta$group[i] <- g[i]
  }
}
if ("GAD7" %in% names(input)) {
  s <- number(input$GAD7[first])
  ok <- gad_complete[first] %in% TRUE
  middle <- which(is.na(meta$group) & ok & !is.na(s) & s >= 5 & s <= 9)
  if (length(middle)) meta$group[middle] <- "middle_5_9"
}
meta$group[is.na(meta$group)] <- "invalid_or_unknown"
group_row <- meta$group[match(input$life_id, meta$life_id)]

criteria <- names(common)
criterion_rows <- list()
previous <- rep(FALSE, n)
for (crit in criteria) {
  f <- common[[crit]]
  only_new <- f & !previous
  crit_events <- common_event_table[common_event_table$criterion == crit, , drop = FALSE]
  criterion_rows[[length(criterion_rows) + 1L]] <- data.frame(
    criterion = crit, level = "common", pipeline = "",
    rows_flagged = sum(f), participants_flagged = sum(f[first]),
    newly_excluded_rows = sum(only_new),
    newly_excluded_participants = sum(only_new[first]),
    high_anxiety_participants = sum(f[first] & group_row[first] == "high"),
    low_anxiety_participants = sum(f[first] & group_row[first] == "low"),
    source_columns = join(crit_events$source_column, "; "),
    stringsAsFactors = FALSE)
  previous <- previous | f
}
for (i in seq_along(stems)) {
  stem <- stems[i]
  f <- pipe_matrix[, i]
  pe <- pipe_event_table[pipe_event_table$pipeline == stem, , drop = FALSE]
  criterion_rows[[length(criterion_rows) + 1L]] <- data.frame(
    criterion = "pipeline_total_including_common", level = "pipeline", pipeline = stem,
    rows_flagged = sum(f), participants_flagged = sum(f[first]),
    newly_excluded_rows = sum(f & !common_exclude),
    newly_excluded_participants = sum(f[first] & !common_exclude[first]),
    high_anxiety_participants = sum(f[first] & group_row[first] == "high"),
    low_anxiety_participants = sum(f[first] & group_row[first] == "low"),
    source_columns = join(pe$source_column, "; "), stringsAsFactors = FALSE)
}
criterion_summary <- bind(criterion_rows, c("criterion", "level", "pipeline",
  "rows_flagged", "participants_flagged", "newly_excluded_rows",
  "newly_excluded_participants", "high_anxiety_participants",
  "low_anxiety_participants", "source_columns"))
save_csv(criterion_summary, "exclusion_criterion_counts.csv")

# Individual source-column counts: overlapping criteria are deliberately NOT
# represented as additive participant losses.
if (nrow(event_table)) {
  key <- unique(event_table[, c("type", "pipeline", "criterion", "source_column")])
  counts <- lapply(seq_len(nrow(key)), function(j) {
    q <- key[j, ]
    x <- event_table[event_table$type == q$type & event_table$pipeline == q$pipeline &
                     event_table$criterion == q$criterion &
                     event_table$source_column == q$source_column, , drop = FALSE]
    data.frame(q, rows_flagged = length(unique(x$data_csv_row)),
               participants_flagged = length(unique(x$life_id)))
  })
  save_csv(do.call(rbind, counts), "exclusion_source_counts.csv")
} else save_csv(data.frame(type=character(), pipeline=character(), criterion=character(),
                           source_column=character(), rows_flagged=integer(),
                           participants_flagged=integer()), "exclusion_source_counts.csv")

anxiety <- lapply(c("high", "low", "middle_5_9", "invalid_or_unknown"), function(level) {
  i <- meta$group == level
  data.frame(group = level,
    participants_before = sum(i), participants_common_excluded = sum(i & common_exclude[first]),
    participants_overall_excluded = sum(i & overall[first]),
    participants_retained = sum(i & !overall[first]),
    rows_before = sum(group_row == level),
    rows_overall_excluded = sum(group_row == level & overall),
    rows_retained = sum(group_row == level & !overall))
})
save_csv(do.call(rbind, anxiety), "exclusion_anxiety_groups.csv")

# Full all-row decision inventory (includes physically excluded rows). The
# original upstream row number remains attached for cross-referencing.
decision_cols <- c("life_id", "common_exclude", "exclusion_reasons",
                   paste0(stems, "_exclude"), "pipeline_exclusion_reasons",
                   "overall_exclude", "overall_exclusion_reasons")
decision_file <- cbind(data.frame(data_csv_row = row_number, gad7_group = group_row),
                       input[, decision_cols, drop = FALSE])
save_csv(decision_file, "exclusion_decisions_all_rows.csv")
save_csv(decision_file[overall, , drop = FALSE], "exclusion_removed_rows.csv")

# The same numbers as the MLM would see after applying each branch's flag.
pipeline_flow <- lapply(seq_along(stems), function(i) {
  failed <- pipe_matrix[, i]
  eligible <- !failed
  data.frame(pipeline = stems[i], high_pass = assigned_hp[i],
    rows_eligible = sum(eligible), participants_eligible = sum(eligible[first]),
    rows_excluded = sum(failed), participants_excluded = sum(failed[first]),
    eligible_high_anxiety_IDs = sum(eligible[first] & group_row[first] == "high"),
    eligible_low_anxiety_IDs = sum(eligible[first] & group_row[first] == "low"))
})
save_csv(do.call(rbind, pipeline_flow), "exclusion_pipeline_counts.csv")

# Count overlapping common causes without falsely adding overlapping losses.
common_burden <- rowSums(do.call(cbind, common))
overlap <- table(common_burden[first])
save_csv(data.frame(n_common_criteria = as.integer(names(overlap)),
                    n_participants = as.integer(overlap)),
         "exclusion_common_overlap.csv")

# ---- 6. Filter rows, refresh duplicated-ID flags, repeat matching audit ----
remaining <- input[!overall, , drop = FALSE]
original_rows_remaining <- row_number[!overall]
if ("non-unique" %in% names(remaining)) {
  remaining[["non-unique"]] <- as.integer(
    table(remaining$life_id)[remaining$life_id] > 1L)
}
if (any(remaining$overall_exclude != 0L)) stop("Internal error: excluded row retained")
if (nrow(remaining) + sum(overall) != n)
  stop("Internal error: filtering did not partition the original rows")
if (any(!overall & rowSums(!pipe_matrix) == 0L))
  stop("Internal error: participant retained without an eligible pipeline")
if (nrow(remaining) && !identical(remaining[["non-unique"]],
   as.integer(table(remaining$life_id)[remaining$life_id] > 1L)))
  stop("Internal error: duplicated-ID flag not updated")

# Retain existing manual-navigation context from LIFE_generate_data.R.
old_prov_path <- file.path(source_dir, "data_row_provenance.csv")
old_prov <- if (file.exists(old_prov_path)) read_csv(old_prov_path) else NULL
old_unresolved_path <- file.path(source_dir, "data_unresolved_rows.csv")
old_unresolved <- if (file.exists(old_unresolved_path)) read_csv(old_unresolved_path) else NULL
nav <- list()
conflict_rows <- list()
for (id in unique(remaining$life_id[duplicated(remaining$life_id) |
                                       duplicated(remaining$life_id, fromLast = TRUE)])) {
  i <- which(remaining$life_id == id)
  original <- original_rows_remaining[i]
  prov <- if (!is.null(old_prov) && "data_csv_row" %in% names(old_prov))
    old_prov[match(as.character(original), old_prov$data_csv_row), , drop = FALSE] else NULL
  old <- if (!is.null(old_unresolved))
    old_unresolved[old_unresolved$life_id == id, , drop = FALSE] else NULL
  pick <- function(col, frame) if (!is.null(frame) && nrow(frame) && col %in% names(frame))
    join(frame[[col]]) else ""
  values_at <- function(col) if (col %in% names(remaining))
    join(paste0(i + 1L, ": ", ifelse(is.na(norm(remaining[[col]][i])),
                                     "<missing>", remaining[[col]][i]))) else ""
  cols <- setdiff(names(remaining), c("life_id", "non-unique"))
  state <- vapply(cols, function(col) {
    x <- norm(remaining[[col]][i])
    unique_nonmissing <- unique(x[!is.na(x)])
    if (length(unique_nonmissing) > 1L) "conflict" else
      if (length(unique_nonmissing) == 1L && anyNA(x)) "complementary" else "same"
  }, character(1L))
  nav[[length(nav) + 1L]] <- data.frame(
    life_id = id, n_rows = length(i),
    data_after_exclusion_csv_rows = join(i + 1L),
    original_data_csv_rows = join(original),
    primary_input_csv_rows = pick("primary_input_csv_row", prov),
    all_input_csv_rows = pick("all_input_csv_rows", prov),
    selection_rules = pick("selection_rule", prov),
    eeg_dates = values_at("date_eeg"), eeg_groups = values_at("group"),
    recording_paths = pick("recording_paths", old),
    selected_recording_path = pick("selected_recording_path", old),
    n_conflicting_columns = sum(state == "conflict"),
    n_complementary_columns = sum(state == "complementary"),
    differing_columns = join(cols[state != "same"], "; "),
    stringsAsFactors = FALSE)
  for (col in cols[state != "same"]) {
    x <- norm(remaining[[col]][i])
    conflict_rows[[length(conflict_rows) + 1L]] <- data.frame(
      life_id = id, column = col, comparison = state[col],
      values_by_after_exclusion_csv_row =
        join(paste0(i + 1L, ": ", ifelse(is.na(x), "<missing>", x))),
      stringsAsFactors = FALSE)
  }
}
nav_cols <- c("life_id", "n_rows", "data_after_exclusion_csv_rows",
  "original_data_csv_rows", "primary_input_csv_rows", "all_input_csv_rows",
  "selection_rules", "eeg_dates", "eeg_groups", "recording_paths",
  "selected_recording_path", "n_conflicting_columns", "n_complementary_columns",
  "differing_columns")
conflict_cols <- c("life_id", "column", "comparison", "values_by_after_exclusion_csv_row")
save_csv(bind(nav, nav_cols), "exclusion_unresolved_rows.csv")
save_csv(bind(conflict_rows, conflict_cols), "exclusion_unresolved_columns.csv")

# Report cells where multiple unmerged rows for one life_id disagree about
# exclusion-relevant source columns. The positive answer wins conservatively.
disagreement <- list()
relevant <- unique(c("all_GAD7", "eeg_alcohol_today", neurological, cardiac,
  "medications_cns_exclude_definite", "extrasystolic_exclusion", "drowsy_15min",
  "preproc_hp_0.1", "preproc_hp_0.5", roi_flags))
for (id in uid[table(input$life_id)[uid] > 1L]) {
  ix <- which(input$life_id == id)
  for (col in relevant) {
    x <- norm(input[[col]][ix])
    if (length(unique(x[!is.na(x)])) > 1L)
      disagreement[[length(disagreement) + 1L]] <- data.frame(
        life_id = id, column = col,
        values_by_original_data_csv_row = join(paste0(row_number[ix], ": ",
                                            ifelse(is.na(x), "<missing>", x))),
        stringsAsFactors = FALSE)
  }
}
save_csv(bind(disagreement, c("life_id", "column", "values_by_original_data_csv_row")),
         "exclusion_duplicate_row_disagreements.csv")

# ---- 7. Non-blocking review items and final export -------------------------
review <- data.frame(
  topic = c("source coding", "neurological variable coverage", "cardiac rhythm coding",
            "implanted pacemaker definition", "brain tumour / Alzheimer's",
            "non-excluding neurologic histories", "duplicate rows",
            "missing preprocessing records", "overall logic"),
  action = c(
    "Confirmed default: yes=1, no=0, unknown=98/blank. Investigate any unexpected codes in exclusion_code_audit.csv; no automatic interpretation of 2.",
    "Confirm epilepsy, multiple_sclerosis, schlaganfall and intracranial haemorrhage are the intended selected neurological exclusions.",
    "Review ECG/echo rhythm categories and whether ecg_sa_block/ecg_av_block should count; their coding is not yet specified. Secondary MEDANAM cardiac fields were already cross-checked with KARDANAM and intentionally removed.",
    "CONFIRMED: ecg_pacemaker_atrial and ecg_pacemaker_chamber identify implanted pacemakers; either positive indicator excludes.",
    "Current selected-variable mapping has no explicit brain tumour / Alzheimer indicators; locate reliable additional fields if available.",
    "CONFIRMED: meningitis, schaedel_hirn_trauma and schaedelfraktur do not trigger exclusion.",
    "Review exclusion_duplicate_row_disagreements.csv; exclusion from any repeated row propagates to all rows.",
    "Check participant coverage in preprocessing_logs.ods; missing records remain unknown, not positive exclusions.",
    "Overall excludes common failures or participants excluded from every pipeline, preserving the preregistered multiverse."),
  stringsAsFactors = FALSE)
save_csv(review, "exclusion_review_items.csv")
save_csv(bind(unmatched_codes, c("life_id", "atc_code", "issue")),
         "exclusion_medication_unmatched_ATC.csv")
save_csv(bind(code_audit, c("source_column", "value", "n_rows", "interpretation")),
         "exclusion_code_audit.csv")
save_csv(bind(unknown_cells,
              c("life_id", "data_csv_row", "source_column", "source_value", "unknown_type")),
         "exclusion_unknown_source_cells.csv")
issue_table <- bind(issues, c("issue_type", "life_id", "source", "detail"))
save_csv(issue_table, "exclusion_issues.csv")

summary_lines <- capture.output({
  cat("========== Preregistered exclusion audit =========\n")
  cat("Input: ", file.path(deriv, "data.csv"), "\n", sep="")
  cat("Input rows / unique IDs: ", n, " / ", length(uid), "\n", sep="")
  cat("Failed any common criterion: ", sum(common_exclude), " rows / ",
      sum(common_exclude[first]), " IDs\n", sep="")
  cat("Failed every pipeline: ", sum(all_pipelines_failed), " rows / ",
      sum(all_pipelines_failed[first]), " IDs\n", sep="")
  cat("Overall excluded: ", sum(overall), " rows / ",
      sum(overall[first]), " IDs\n", sep="")
  cat("Overall retained: ", nrow(remaining), " rows / ",
      length(unique(remaining$life_id)), " IDs\n", sep="")
  cat("Retained IDs still having >1 row: ", sum(table(remaining$life_id) > 1L), "\n")
  cat("Retained rows with non-unique ID: ", sum(remaining[["non-unique"]] == 1), "\n")
  cat("Anxiety groups (unique IDs):\n"); print(anxiety)
  cat("\nCommon criteria and pipeline summary:\n")
  print(criterion_summary[, c("criterion", "pipeline", "participants_flagged",
    "newly_excluded_participants")], row.names = FALSE)
  cat("\nPreprocessing events by stage and hp:\n")
  if (nrow(pre_log_events)) print(with(pre_log_events, table(hp, step)))
  cat("\nPotential coding/source issues: ", nrow(issue_table), "\n", sep="")
  if (nrow(issue_table)) print(sort(table(issue_table$issue_type), decreasing = TRUE))
  cat("\nSource clinical history: positive indicates exclusion regardless of conflicting duplicates.\n")
  cat("GAD 5-9 remains eligible for continuous models, excluded only from the high-vs-low model.\n")
  cat("Unmeasured control covariates do not exclude from main analysis.\n")
  cat("All pipeline-specific flags MUST be applied when reshaping to long format.\n")
  cat("Source/coding issues are retained in diagnostics; export is not gated.\n")
})
writeLines(summary_lines, file.path(out_dir, "exclusion_summary.txt"))
cat(paste(summary_lines, collapse="\n"), "\n")

# QC findings are informational. Only missing required inputs or failed
# internal consistency checks above may interrupt the processing run.
final_path <- file.path(deriv, "data_after_exclusion.csv")
write.csv(remaining, final_path, row.names = FALSE, na = "")
cat("Saved final table: ", final_path, "\n", sep="")
cat("Saved exclusion diagnostics to: ", out_dir, "\n", sep="")
