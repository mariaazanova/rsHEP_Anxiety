#!/usr/bin/env Rscript
# rsHEP / LIFE: PIPELINE-INDEPENDENT exclusions, after LIFE_generate_data.R
# Usage: Rscript exclusion.R [life_directory]
# Default life_directory: /data/p_02865/data/LIFE
# Required inputs:
#   derivatives/data.csv
#   derivatives/matching_and_summaries/medications/LIFE_medications_by_ID.csv
#   derivatives/matching_and_summaries/medications/medications_cns_yes_qc.csv
# Outputs:
#   derivatives/data_after_exclusion.csv
#   derivatives/matching_and_summaries/exclusion/*
#
# Scope: ONLY questionnaire validity, alcohol on EEG day, relevant medical
# history, CNS-active medications, preregistered extrasystole exclusion and
# excessive drowsiness. These are participant-level, NOT pipeline-specific.
#
# This script DOES NOT use preprocessing HP status, per-pipeline epoch counts,
# ROI HEP values/outliers, or preprocessing_logs.ods for decisions or audits.
# Any pipeline columns already present in data.csv are passed through unchanged.
# Participant-level exclusion flags/reasons and the non-unique indicator are
# retained in audit outputs, but omitted from the filtered analysis dataset.
#
# Clinical flag coding: 1=yes/exclude, 0=no, 98=unknown, missing=unknown.
# all_GAD7 is a validity flag: 0=invalid/exclude, 1=valid/retain;
# missing/98=unknown/retain. Other unexpected codes are audited but not
# interpreted. A positive exclusion on ANY row for a repeated life_id is
# propagated to ALL rows of that participant; conflicts are documented.
#
# Preregistration: section 2.2.4 and corresponding upstream generated flags.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: Rscript exclusion.R [life_directory]")
life <- if (length(args)) args[1L] else "/data/p_02865/data/LIFE"
deriv <- file.path(life, "derivatives")
source_dir <- file.path(deriv, "matching_and_summaries")
out_dir <- file.path(source_dir, "exclusion")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# The final dataset is always written after internal consistency checks.
# Clinical-coverage questions and unrecognized source codes appear in the
# diagnostics and do not impose an additional release gate.
# Confirmed by the data owner: BOTH ECG pacemaker indicators reflect
# implanted pacemakers. A positive response in either is an exclusion.
# No additional exclusion for meningitis, head trauma or skull fracture.
# Secondary MEDANAM cardiac fields were already checked against KARDANAM,
# found equivalent, and deliberately dropped in LIFE_generate_data.R.
# Other possible clinical variables can be reviewed through non-blocking QC;
# no unregistered criteria are added automatically.

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

# Criterion events capture source columns and original values, with optional
# explanations. One event per positive input row/source column.
events <- list()
add_events <- function(criterion, source, positive, detail = NULL) {
  if (!any(positive)) return(invisible(NULL))
  ix <- which(positive)
  if (is.null(detail)) detail <- rep("", n)
  if (!source %in% names(input)) values <- rep("", n) else values <- input[[source]]
  events[[length(events) + 1L]] <<- data.frame(
    life_id = input$life_id[ix], data_csv_row = row_number[ix],
    criterion = criterion, source_column = source,
    source_value = ifelse(is.na(values[ix]), "", values[ix]),
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
  add_events(name, source, triggered_cells, detail = detail)
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
# Only echo_arrhythmia is used for arrhythmia exclusion; rhythm_disturbances
# is deliberately not evaluated. Other rhythm categories are left for QC.
cardiac <- c("infarct", "bypass",
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

# ---- 3. Overall decisions: ONLY pipeline-independent criteria ------------
common_exclude <- Reduce(`|`, common)
if (anyNA(common_exclude)) stop("Internal error: NA in common_exclude")
# Maintain both public flag names for compatibility. Their values are identical
# because pipeline-dependent exclusions will be applied by a separate script.
input$common_exclude <- as.integer(common_exclude)
input$overall_exclude <- as.integer(common_exclude)

# Preserve *triggering* raw source rows separately from propagated decisions.
event_cols <- c("life_id", "data_csv_row", "criterion", "source_column",
                "source_value", "detail")
event_table <- bind(events, event_cols)
save_csv(event_table, "exclusion_events.csv")
source_reason <- function(e) {
  detail <- norm(e$detail)
  if (is.na(detail)) e$source_column else
    paste0(e$source_column, " [", detail, "]")
}
reason_by_id <- setNames(rep("", length(uid)), uid)
for (id in uid) {
  ev <- event_table[event_table$life_id == id, , drop = FALSE]
  if (nrow(ev)) {
    reason_by_id[id] <- join(vapply(seq_len(nrow(ev)), function(j)
      source_reason(ev[j, , drop = FALSE]), character(1L)))
  }
}
input$exclusion_reasons <- unname(reason_by_id[input$life_id])
if (any(common_exclude & !nzchar(input$exclusion_reasons)))
  stop("Internal error: an excluded participant has no source reason")

new_cols <- c("common_exclude", "overall_exclude", "exclusion_reasons")
# Any original pipeline columns are only copied, with no inspection of values.
input <- input[, c(setdiff(names(input), new_cols), new_cols), drop = FALSE]

# ---- 4. Common exclusion diagnostics -------------------------------------
meta <- data.frame(life_id = uid, group = NA_character_, stringsAsFactors = FALSE)
if ("gad7_group" %in% names(input)) {
  g <- tolower(norm(input$gad7_group[first]))
  meta$group[g %in% c("high", "low")] <- g[g %in% c("high", "low")]
}
if ("GAD7" %in% names(input)) {
  s <- number(input$GAD7[first])
  valid_gad <- gad_complete[first] %in% TRUE
  middle <- is.na(meta$group) & valid_gad & !is.na(s) & s >= 5 & s <= 9
  meta$group[middle] <- "middle_5_9"
}
meta$group[is.na(meta$group)] <- "invalid_or_unknown"
group_row <- meta$group[match(input$life_id, meta$life_id)]

criteria <- names(common)
criterion_rows <- list()
previous <- rep(FALSE, n)
for (crit in criteria) {
  failed <- common[[crit]]
  only_new <- failed & !previous
  crit_events <- event_table[event_table$criterion == crit, , drop = FALSE]
  criterion_rows[[length(criterion_rows) + 1L]] <- data.frame(
    criterion = crit,
    rows_flagged = sum(failed),
    participants_flagged = sum(failed[first]),
    newly_excluded_rows = sum(only_new),
    newly_excluded_participants = sum(only_new[first]),
    high_anxiety_participants = sum(failed[first] & meta$group == "high"),
    low_anxiety_participants = sum(failed[first] & meta$group == "low"),
    source_columns = join(crit_events$source_column, "; "),
    stringsAsFactors = FALSE)
  previous <- previous | failed
}
criterion_cols <- c("criterion", "rows_flagged", "participants_flagged",
                    "newly_excluded_rows", "newly_excluded_participants",
                    "high_anxiety_participants", "low_anxiety_participants",
                    "source_columns")
criterion_summary <- bind(criterion_rows, criterion_cols)
save_csv(criterion_summary, "exclusion_criterion_counts.csv")

if (nrow(event_table)) {
  key <- unique(event_table[, c("criterion", "source_column"), drop = FALSE])
  source_counts <- lapply(seq_len(nrow(key)), function(j) {
    q <- key[j, , drop = FALSE]
    ev <- event_table[event_table$criterion == q$criterion &
                      event_table$source_column == q$source_column, , drop = FALSE]
    data.frame(q, rows_flagged = length(unique(ev$data_csv_row)),
               participants_flagged = length(unique(ev$life_id)))
  })
  save_csv(do.call(rbind, source_counts), "exclusion_source_counts.csv")
} else {
  save_csv(data.frame(criterion = character(), source_column = character(),
                      rows_flagged = integer(), participants_flagged = integer()),
           "exclusion_source_counts.csv")
}

anxiety <- lapply(c("high", "low", "middle_5_9", "invalid_or_unknown"), function(level) {
  i <- meta$group == level
  data.frame(group = level,
    participants_before = sum(i),
    participants_excluded = sum(i & common_exclude[first]),
    participants_retained = sum(i & !common_exclude[first]),
    rows_before = sum(group_row == level),
    rows_excluded = sum(group_row == level & common_exclude),
    rows_retained = sum(group_row == level & !common_exclude))
})
save_csv(do.call(rbind, anxiety), "exclusion_anxiety_groups.csv")

burden <- rowSums(do.call(cbind, common))
overlap <- table(burden[first])
save_csv(data.frame(n_common_criteria = as.integer(names(overlap)),
                    n_participants = as.integer(overlap)),
         "exclusion_common_overlap.csv")

# Full decision inventory, including rows subsequently removed.
decision_cols <- c("life_id", "common_exclude", "overall_exclude",
                   "exclusion_reasons")
decision_file <- cbind(data.frame(data_csv_row = row_number,
                                  gad7_group = group_row),
                       input[, decision_cols, drop = FALSE])
save_csv(decision_file, "exclusion_decisions_all_rows.csv")
save_csv(decision_file[common_exclude, , drop = FALSE], "exclusion_removed_rows.csv")

# ---- 4a. ADDITIONAL evidence audit for common exclusions ONLY ------------
# Entirely diagnostic. No additional criteria, changes to the decisions,
# changes to GAD-7 validity, or changes to the main data exports.
#
# Three companion files capture the evidence at source, criterion and ID level:
#   exclusion_evidence_by_source.csv      one LIFE ID x exclusion source column
#   exclusion_evidence_by_criterion.csv   one LIFE ID x common criterion
#   exclusion_evidence_by_ID.csv          one LIFE ID (overall assessment)
#   exclusion_evidence_summary.txt        interpretation and overall counts
#
# "Information complete" is conservative: EVERY original data.csv row must
# have a recognized definite code for EVERY applicable source column. Thus
# duplicates with 0 and blank are incomplete, while 1 anywhere still triggers
# the existing exclusion (and is recorded as confirmed exclusion). This audit
# assesses the IMPLEMENTED criteria only, not conditions without mapped columns
# (such as the clinical coverage gaps listed in exclusion_review_items.csv).

criterion_sources <- list(
  invalid_GAD7 = "all_GAD7",
  alcohol_on_EEG_day = "eeg_alcohol_today",
  neurological_history = neurological,
  cardiac_history = cardiac,
  CNS_medications = "medications_cns_exclude_definite",
  extrasystoles_above_4SD = "extrasystolic_exclusion",
  drowsy_over_15min = "drowsy_15min"
)
if (!identical(names(criterion_sources), names(common)))
  stop("Evidence audit mapping does not match the implemented common criteria")
if (!all(unlist(criterion_sources, use.names = FALSE) %in% names(input)))
  stop("Evidence audit cannot find an implemented exclusion source")

# This reproduces the source-code interpretation of parse_binary(), with the
# single pre-existing GAD-7 reversal: all_GAD7=0 means invalid (exclude).
# In particular, neither 98/missing nor UNRECOGNIZED constitutes exclusion.
source_evidence_parts <- list()
for (crit in names(criterion_sources)) {
  for (col in criterion_sources[[crit]]) {
    raw <- norm(input[[col]])
    observed <- toupper(raw)
    yes_code <- if (crit == "invalid_GAD7") "0" else "1"
    no_code <- if (crit == "invalid_GAD7") "1" else "0"
    yes <- !is.na(observed) & observed == yes_code
    no <- !is.na(observed) & observed == no_code
    missing <- is.na(observed)
    coded_unknown <- !is.na(observed) & observed == "98"
    unrecognized <- !is.na(observed) &
      !(observed %in% c(yes_code, no_code, "98"))
    count_per_id <- function(flag) as.integer(vapply(groups, function(ix)
      sum(flag[ix]), integer(1L)))
    n_yes <- count_per_id(yes)
    n_no <- count_per_id(no)
    n_missing <- count_per_id(missing)
    n_unknown <- count_per_id(coded_unknown)
    n_unrecognized <- count_per_id(unrecognized)
    n_rows <- as.integer(lengths(groups))
    complete <- n_missing == 0L & n_unknown == 0L & n_unrecognized == 0L
    status <- ifelse(n_yes > 0L, "confirmed_exclusion",
                     ifelse(n_no == n_rows, "confirmed_no_exclusion",
                            "incomplete_assessment"))
    values <- vapply(groups, function(ix) {
      x <- raw[ix]
      join(ifelse(is.na(x), "<missing>", x), "; ")
    }, character(1L))
    values_by_row <- vapply(groups, function(ix) {
      x <- raw[ix]
      paste(paste0(row_number[ix], ":",
                   ifelse(is.na(x), "<missing>", x)), collapse = " | ")
    }, character(1L))
    source_evidence_parts[[length(source_evidence_parts) + 1L]] <- data.frame(
      life_id = uid, criterion = crit, source_column = col,
      n_input_rows = n_rows,
      source_values = unname(values),
      values_by_data_csv_row = unname(values_by_row),
      n_positive_rows = n_yes, n_negative_rows = n_no,
      n_missing_rows = n_missing, n_coded_98_rows = n_unknown,
      n_unrecognized_rows = n_unrecognized,
      confirmed_exclusion = as.integer(n_yes > 0L),
      information_complete = as.integer(complete),
      conflicting_definite_codes = as.integer(n_yes > 0L & n_no > 0L),
      evidence_status = status, stringsAsFactors = FALSE)
  }
}
source_evidence <- do.call(rbind, source_evidence_parts)
rownames(source_evidence) <- NULL
save_csv(source_evidence, "exclusion_evidence_by_source.csv")

# Group clinical source columns by the exact criteria used in the script;
# a confirmed positive in ANY source establishes that criterion's exclusion.
criterion_evidence_parts <- list()
for (crit in names(criterion_sources)) {
  sub <- source_evidence[source_evidence$criterion == crit, , drop = FALSE]
  sub_groups <- split(seq_len(nrow(sub)), factor(sub$life_id, levels = uid))
  confirmed <- vapply(sub_groups, function(ix)
    any(sub$confirmed_exclusion[ix] == 1L), logical(1L))
  complete <- vapply(sub_groups, function(ix)
    all(sub$information_complete[ix] == 1L), logical(1L))
  # Audit only: ensure our independent reconstruction exactly matches the
  # established exclusion flags; do not override any decisions.
  if (!identical(unname(confirmed), unname(common[[crit]][first])))
    stop("Evidence audit does not reproduce common exclusion: ", crit)
  source_names <- function(predicate) vapply(sub_groups, function(ix)
    join(sub$source_column[ix][predicate[ix]], "; "), character(1L))
  has_missing <- sub$n_missing_rows > 0L
  has_unknown <- sub$n_coded_98_rows > 0L
  has_unrecognized <- sub$n_unrecognized_rows > 0L
  incomplete <- sub$information_complete == 0L
  conflict <- sub$conflicting_definite_codes == 1L
  criterion_evidence_parts[[length(criterion_evidence_parts) + 1L]] <- data.frame(
    life_id = uid, criterion = crit,
    confirmed_exclusion = as.integer(confirmed),
    information_complete = as.integer(complete),
    evidence_status = ifelse(confirmed, "confirmed_exclusion",
                             ifelse(complete, "confirmed_no_exclusion",
                                    "incomplete_assessment")),
    n_source_columns = as.integer(length(criterion_sources[[crit]])),
    n_incomplete_source_columns = as.integer(vapply(sub_groups, function(ix)
      sum(incomplete[ix]), integer(1L))),
    positive_source_columns = unname(source_names(sub$confirmed_exclusion == 1L)),
    source_columns_with_missing = unname(source_names(has_missing)),
    source_columns_with_98 = unname(source_names(has_unknown)),
    source_columns_with_unrecognized = unname(source_names(has_unrecognized)),
    source_columns_with_conflicting_codes = unname(source_names(conflict)),
    stringsAsFactors = FALSE)
}
criterion_evidence <- do.call(rbind, criterion_evidence_parts)
rownames(criterion_evidence) <- NULL
save_csv(criterion_evidence, "exclusion_evidence_by_criterion.csv")

# Overall classification keeps confirmed exclusion separate from missingness:
# an excluded person can ALSO have incomplete information for other sources.
# Only a RETAINED person with no missing/98/unrecognized codes for all mapped
# exclusion sources can be described as fully assessed under these criteria.
crit_groups <- split(seq_len(nrow(criterion_evidence)),
                     factor(criterion_evidence$life_id, levels = uid))
source_groups <- split(seq_len(nrow(source_evidence)),
                       factor(source_evidence$life_id, levels = uid))
criterion_list <- function(predicate) vapply(crit_groups, function(ix)
  join(criterion_evidence$criterion[ix][predicate[ix]], "; "), character(1L))
source_list <- function(predicate) vapply(source_groups, function(ix)
  join(source_evidence$source_column[ix][predicate[ix]], "; "), character(1L))
crit_missing <- criterion_evidence$information_complete == 0L
source_missing <- source_evidence$n_missing_rows > 0L
source_unknown <- source_evidence$n_coded_98_rows > 0L
source_unrecognized <- source_evidence$n_unrecognized_rows > 0L
source_conflict <- source_evidence$conflicting_definite_codes == 1L
n_incomplete_criteria <- as.integer(vapply(crit_groups, function(ix)
  sum(crit_missing[ix]), integer(1L)))
confirmed_by_id <- as.logical(common_exclude[first])
if (!identical(confirmed_by_id, as.logical(vapply(crit_groups, function(ix)
  any(criterion_evidence$confirmed_exclusion[ix] == 1L), logical(1L)))))
  stop("Evidence audit overall classification disagrees with common_exclude")

id_evidence <- data.frame(
  life_id = uid,
  n_input_rows = as.integer(lengths(groups)),
  original_data_csv_rows = vapply(groups, function(ix)
    paste(row_number[ix], collapse = "; "), character(1L)),
  common_exclude = as.integer(confirmed_by_id),
  evidence_status = ifelse(confirmed_by_id, "confirmed_excluded",
                           ifelse(n_incomplete_criteria == 0L,
                                  "assessed_no_exclusion",
                                  "retained_incomplete_assessment")),
  information_complete = as.integer(n_incomplete_criteria == 0L),
  n_incomplete_criteria = n_incomplete_criteria,
  confirmed_exclusion_criteria = unname(criterion_list(
    criterion_evidence$confirmed_exclusion == 1L)),
  criteria_with_incomplete_information = unname(criterion_list(crit_missing)),
  source_columns_with_missing = unname(source_list(source_missing)),
  source_columns_with_98 = unname(source_list(source_unknown)),
  source_columns_with_unrecognized = unname(source_list(source_unrecognized)),
  source_columns_with_conflicting_codes = unname(source_list(source_conflict)),
  stringsAsFactors = FALSE)
save_csv(id_evidence, "exclusion_evidence_by_ID.csv")

# Add a self-contained key and counts without changing exclusion_summary.txt.
evidence_summary <- c(
  "Common-exclusion EVIDENCE AUDIT (diagnostic only; no changed exclusions)",
  "Scope: implemented pipeline-INDEPENDENT exclusion criteria only.",
  "Confirmed exclusion: at least one definitive positive source code, even if other information is missing.",
  "Assessed no exclusion: every mapped source has definite negative information on every data.csv row.",
  "Retained incomplete assessment: no positive exclusion, but at least one source has blank/NA, 98, or unrecognized codes.",
  "Information complete: no blanks/NA, code 98, or unrecognized codes in any mapped source, including duplicate rows.",
  "GAD-7 is unchanged: all_GAD7=0 means excluded; all_GAD7=1 means valid; unknown codes do not exclude.",
  "Unmapped clinical conditions are not covered by this completeness indicator (see exclusion_review_items.csv).",
  "", "Participants by overall evidence status:",
  capture.output(print(table(id_evidence$evidence_status))),
  "", paste0("Confirmed excluded with incomplete source information: ",
             sum(confirmed_by_id & n_incomplete_criteria > 0L)),
  paste0("Retained with incomplete source information: ",
         sum(!confirmed_by_id & n_incomplete_criteria > 0L)),
  "", "Incomplete assessments by criterion among retained participants:",
  capture.output(print(table(criterion_evidence$criterion[
    criterion_evidence$information_complete == 0L &
      !(criterion_evidence$life_id %in% uid[confirmed_by_id])]))))
writeLines(evidence_summary,
           file.path(out_dir, "exclusion_evidence_summary.txt"))

# ---- 5. Filter participants, then refresh non-unique ID diagnostics -------
remaining <- input[!common_exclude, , drop = FALSE]
original_rows_remaining <- row_number[!common_exclude]
if ("non-unique" %in% names(remaining)) {
  remaining[["non-unique"]] <- as.integer(
    table(remaining$life_id)[remaining$life_id] > 1L)
}
if (nrow(remaining) + sum(common_exclude) != n)
  stop("Internal error: filtering did not partition the original data")
if (any(remaining$overall_exclude != 0L))
  stop("Internal error: excluded participant retained")
if ("non-unique" %in% names(remaining) &&
    !identical(remaining[["non-unique"]],
               as.integer(table(remaining$life_id)[remaining$life_id] > 1L)))
  stop("Internal error: non-unique flag was not refreshed correctly")

# Use the original provenance and matching diagnostics for manual navigation.
old_prov_path <- file.path(source_dir, "data_row_provenance.csv")
old_unresolved_path <- file.path(source_dir, "data_unresolved_rows.csv")
old_prov <- if (file.exists(old_prov_path)) read_csv(old_prov_path) else NULL
old_unresolved <- if (file.exists(old_unresolved_path)) read_csv(old_unresolved_path) else NULL

# Only check non-pipeline columns for differences. Names from upstream HEP/ECG
# and HP pipelines are identified without reading their values. Those values
# play no role in participant exclusion or matching QC at this stage.
is_pipeline_column <- function(col) {
  grepl("(^preproc_hp_|_roi_(HEP|ECG)(_|$)|_n_epochs_(total|rejected|kept)$|(^|_)hp_(01|05)(_|$))",
        col, ignore.case = TRUE, perl = TRUE)
}
comparison_cols <- setdiff(names(remaining)[!vapply(names(remaining),
    is_pipeline_column, logical(1L))],
    c("life_id", "non-unique", "common_exclude", "overall_exclude",
      "exclusion_reasons"))

nav <- list()
conflicts <- list()
repeated_ids <- unique(remaining$life_id[duplicated(remaining$life_id) |
                                           duplicated(remaining$life_id, fromLast = TRUE)])
for (id in repeated_ids) {
  ix <- which(remaining$life_id == id)
  original <- original_rows_remaining[ix]
  prov <- if (!is.null(old_prov) && "data_csv_row" %in% names(old_prov))
    old_prov[match(as.character(original), old_prov$data_csv_row), , drop = FALSE] else NULL
  old <- if (!is.null(old_unresolved) && "life_id" %in% names(old_unresolved))
    old_unresolved[ids_norm(old_unresolved$life_id) == id, , drop = FALSE] else NULL
  pick <- function(col, frame) {
    if (!is.null(frame) && nrow(frame) && col %in% names(frame))
      join(frame[[col]]) else ""
  }
  values_at <- function(col) {
    if (!(col %in% names(remaining))) return("")
    x <- norm(remaining[[col]][ix])
    join(paste0(ix + 1L, ": ", ifelse(is.na(x), "<missing>", x)))
  }
  state <- vapply(comparison_cols, function(col) {
    x <- norm(remaining[[col]][ix])
    known <- unique(x[!is.na(x)])
    if (length(known) > 1L) "conflict" else
      if (length(known) == 1L && anyNA(x)) "complementary" else "same"
  }, character(1L))
  nav[[length(nav) + 1L]] <- data.frame(
    life_id = id, n_rows = length(ix),
    data_after_exclusion_csv_rows = join(ix + 1L),
    original_data_csv_rows = join(original),
    primary_input_csv_rows = pick("primary_input_csv_row", prov),
    all_input_csv_rows = pick("all_input_csv_rows", prov),
    selection_rules = pick("selection_rule", prov),
    eeg_dates = values_at("date_eeg"), eeg_groups = values_at("group"),
    recording_paths = pick("recording_paths", old),
    selected_recording_path = pick("selected_recording_path", old),
    n_conflicting_columns = sum(state == "conflict"),
    n_complementary_columns = sum(state == "complementary"),
    differing_columns = join(comparison_cols[state != "same"], "; "),
    stringsAsFactors = FALSE)
  for (col in comparison_cols[state != "same"]) {
    x <- norm(remaining[[col]][ix])
    conflicts[[length(conflicts) + 1L]] <- data.frame(
      life_id = id, column = col, comparison = state[col],
      values_by_after_exclusion_csv_row =
        join(paste0(ix + 1L, ": ", ifelse(is.na(x), "<missing>", x))),
      stringsAsFactors = FALSE)
  }
}
nav_cols <- c("life_id", "n_rows", "data_after_exclusion_csv_rows",
  "original_data_csv_rows", "primary_input_csv_rows", "all_input_csv_rows",
  "selection_rules", "eeg_dates", "eeg_groups", "recording_paths",
  "selected_recording_path", "n_conflicting_columns", "n_complementary_columns",
  "differing_columns")
conflict_cols <- c("life_id", "column", "comparison",
                   "values_by_after_exclusion_csv_row")
save_csv(bind(nav, nav_cols), "exclusion_unresolved_rows.csv")
save_csv(bind(conflicts, conflict_cols), "exclusion_unresolved_columns.csv")

# Contradictions between unresolved rows of the same participant are reported;
# the positive response determines participant-level eligibility.
disagreement <- list()
relevant <- unique(c("all_GAD7", "eeg_alcohol_today", neurological, cardiac,
                     "medications_cns_exclude_definite",
                     "extrasystolic_exclusion", "drowsy_15min"))
for (id in uid[table(input$life_id)[uid] > 1L]) {
  ix <- which(input$life_id == id)
  for (col in relevant) {
    x <- norm(input[[col]][ix])
    if (length(unique(x[!is.na(x)])) > 1L) {
      disagreement[[length(disagreement) + 1L]] <- data.frame(
        life_id = id, column = col,
        values_by_original_data_csv_row = join(paste0(row_number[ix], ": ",
                                              ifelse(is.na(x), "<missing>", x))),
        stringsAsFactors = FALSE)
    }
  }
}
save_csv(bind(disagreement,
  c("life_id", "column", "values_by_original_data_csv_row")),
  "exclusion_duplicate_row_disagreements.csv")

# ---- 6. Informational review, audits and final data export ---------------
review <- data.frame(
  topic = c("source coding", "neurological history", "cardiac rhythm coverage",
            "implanted pacemakers", "other neurological histories",
            "brain tumour / Alzheimer's", "duplicate participant rows",
            "pipeline-specific exclusions"),
  action = c(
    "Confirmed: yes=1, no=0, unknown=98/blank; all_GAD7=0 means invalid. Other codes are only audited.",
    "The selected exclusions include epilepsy, multiple sclerosis, stroke, hemorrhage, brain surgery and Parkinson's; coverage of additional diseases can be reviewed separately.",
    "Additional ECG rhythm indicators have not been adjudicated; no new cardiac criteria were silently added.",
    "Both ECG pacemaker indicators identify implanted pacemakers; either positive one triggers exclusion.",
    "Meningitis, head trauma and skull fracture are not exclusions, per data owner.",
    "No explicit brain-tumour or Alzheimer's columns were identified in the selected data; does not block export.",
    "Positive evidence in any unresolved row applies to every row of that participant; see duplicate-row disagreements.",
    "Handled by a later script; no preprocessing or HEP/epoch decisions are performed here."),
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
  cat("========== LIFE pipeline-independent exclusion audit =========\n")
  cat("Input: ", file.path(deriv, "data.csv"), "\n", sep = "")
  cat("Input: ", n, " rows / ", length(uid), " unique IDs\n", sep = "")
  cat("Common exclusions: ", sum(common_exclude), " rows / ",
      sum(common_exclude[first]), " unique IDs\n", sep = "")
  cat("Retained: ", nrow(remaining), " rows / ",
      length(unique(remaining$life_id)), " unique IDs\n", sep = "")
  cat("Retained IDs with >1 row: ", length(repeated_ids), "\n", sep = "")
  cat("Retained rows with duplicated ID: ",
      sum(duplicated(remaining$life_id) |
          duplicated(remaining$life_id, fromLast = TRUE)), "\n", sep = "")
  cat("\nExclusions per preregistered common criterion:\n")
  print(criterion_summary[, c("criterion", "participants_flagged",
                              "newly_excluded_participants")], row.names = FALSE)
  cat("\nAnxiety-group retention:\n"); print(do.call(rbind, anxiety), row.names = FALSE)
  cat("\nPotential source/coding issues: ", nrow(issue_table), "\n", sep = "")
  if (nrow(issue_table)) print(sort(table(issue_table$issue_type), decreasing = TRUE))
  cat("\nNo HP preprocessing, ROI HEP/outlier, or per-pipeline epoch exclusions performed.\n")
  cat("GAD-7 5–9 is retained for continuous-anxiety analyses.\n")
  cat("Missing control covariates do not exclude from the main analysis.\n")
  cat("QC warnings are diagnostic only, not export gates.\n")
})
writeLines(summary_lines, file.path(out_dir, "exclusion_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

# The retained participants have no common exclusions. Do not carry their
# constant-zero flags, empty exclusion reasons, or matching-only duplicate
# marker into analysis data. All decisions and duplicate-ID details were
# written to the dedicated QC files above. This does NOT change row filtering.
analysis_data <- remaining[, !names(remaining) %in%
  c("common_exclude", "overall_exclude", "exclusion_reasons", "non-unique"),
  drop = FALSE]
if (nrow(analysis_data) != nrow(remaining) ||
    !identical(analysis_data$life_id, remaining$life_id))
  stop("Internal error: cleaning the export changed the retained participants")

final_path <- file.path(deriv, "data_after_exclusion.csv")
write.csv(analysis_data, final_path, row.names = FALSE, na = "")
cat("Saved final dataset: ", final_path, "\n", sep = "")
cat("Excluded QC-only columns: common_exclude, overall_exclude, ",
    "exclusion_reasons, non-unique\n", sep = "")
cat("Saved exclusion diagnostics: ", out_dir, "\n", sep = "")
