# LIFE: prepare analysis data with row reconciliation, HEP/ECG pipelines, HR/HRV, scores.
# Run after LIFE_matching_main.R
# Requires readODS and readxl:
#   install.packages(c("readODS", "readxl"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: Rscript LIFE_generate_data_extended_v3.R [matching_output_dir]")
life <- "/data/p_02865/data/LIFE"
out <- if (length(args)) args[[1L]] else file.path(life, "derivatives/matching_and_summaries")
if (!dir.exists(out)) stop("Missing matching output directory: ", out)

read_csv <- function(path) read.csv(path, colClasses = "character",
  check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM",
  na.strings = "")
optional_csv <- function(name) {
  path <- file.path(out, name)
  if (file.exists(path)) read_csv(path) else NULL
}
write_csv <- function(x, name) write.csv(x, file.path(out, name),
                                          na = "", row.names = FALSE)
normalize <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | !nzchar(x)] <- NA_character_
  x
}
normalize_id <- function(x) {
  x <- toupper(normalize(x))
  x[!is.na(x) & x == "NA"] <- NA_character_
  x
}
show <- function(x) if (length(x) == 1L && !is.na(x)) as.character(x) else ""
joined <- function(x, sep = " | ") {
  x <- unique(as.character(x[!is.na(x)]))
  paste(x[nzchar(x)], collapse = sep)
}
empty_table <- function(cols) {
  z <- setNames(rep(list(character()), length(cols)), cols)
  as.data.frame(z, check.names = FALSE, stringsAsFactors = FALSE)
}
bind_events <- function(lst, cols) {
  if (!length(lst)) return(empty_table(cols))
  result <- do.call(rbind, lst)
  rownames(result) <- NULL
  result[, cols, drop = FALSE]
}
conflicting <- function(a, b, cols) {
  cols <- intersect(cols, intersect(names(a), names(b)))
  cols[!is.na(a[cols]) & !is.na(b[cols]) & a[cols] != b[cols]]
}
changes <- function(a, b, cols) {
  cols <- intersect(cols, intersect(names(a), names(b)))
  cols[(is.na(a[cols]) != is.na(b[cols])) |
       (!is.na(a[cols]) & !is.na(b[cols]) & a[cols] != b[cols])]
}
contains <- function(x, pattern) grepl(pattern, x, perl = TRUE, ignore.case = TRUE)
yes <- function(x) !is.na(x) & toupper(trimws(as.character(x))) %in%
  c("TRUE", "T", "YES", "Y", "1")

# ---- 1. Mapping: keep only requested variables, plus required core fields --
raw <- optional_csv("IDmatch_EEG.csv")
manual <- optional_csv("IDmatch_variables_manual.csv")
if (is.null(raw) || is.null(manual))
  stop("Missing IDmatch_EEG.csv or idmatch_variables_manual.csv in ", out)
# Some spreadsheet programs save underscores as UTF-7 +AF8- strings.
decode <- function(x) gsub("+AF8-", "_", x, fixed = TRUE)
names(manual) <- decode(names(manual))
manual[] <- lapply(manual, decode)
required_map <- c("original_column", "keep", "new_name")
if (!all(required_map %in% names(manual)))
  stop("Mapping lacks: ", joined(setdiff(required_map, names(manual))))
manual[] <- lapply(manual, normalize)
manual$keep[is.na(manual$keep)] <- "0"
if (any(!manual$keep %in% c("0", "1"))) stop("keep must be 0 or 1")
if (anyDuplicated(names(raw))) stop("Duplicate columns in IDmatch_EEG.csv")
identifier_sources <- c("Column1", "Column2", "Table3.TEILNEHMER_SIC")
if (!all(identifier_sources %in% names(raw)))
  stop("Missing one or more required identifier columns: ",
       joined(setdiff(identifier_sources, names(raw))))
forced <- c("Table3.TEILNEHMER_SIC", "Table3.EEG_GRUPPE", "Table3.EEG_DATUM")
selection <- manual[manual$keep == "1" |
                      (!is.na(manual$original_column) &
                         manual$original_column %in% forced), , drop = FALSE]
if (anyNA(selection$original_column) || anyNA(selection$new_name))
  stop("Selected variables must have both original_column and new_name")
if (anyDuplicated(selection$original_column) || anyDuplicated(selection$new_name))
  stop("Selected original_column and new_name must each be unique")
if (any(selection$original_column %in% c("Column1", "Column2")))
  stop("Remove Column1/Column2 from keep=1 mapping; they are always retained")
if (any(selection$new_name %in% c("life_id", "Column2")))
  stop("life_id and Column2 are reserved output names")
if (!all(forced %in% selection$original_column))
  stop("SIC, EEG_GRUPPE and EEG_DATUM must exist in manual mapping")
missing_sources <- setdiff(selection$original_column, names(raw))
if (length(missing_sources))
  stop("Variables not in IDmatch_EEG.csv: ", joined(missing_sources))
name_for <- function(original) {
  idx <- match(original, selection$original_column)
  if (is.na(idx)) NA_character_ else selection$new_name[idx]
}
sic <- name_for("Table3.TEILNEHMER_SIC")
grp <- name_for("Table3.EEG_GRUPPE")
dt <- name_for("Table3.EEG_DATUM")
age <- name_for("Table3.age_EEG")
if (any(grepl(paste0("^(?:", grp, "|", dt, ")[2-9][0-9]*$"),
              selection$new_name)))
  stop("Names such as group2/date_eeg2 are reserved for alternative sessions")

values <- data.frame(life_id = normalize_id(raw$Column1),
  Column2 = normalize(raw$Column2), raw[selection$original_column],
  check.names = FALSE, stringsAsFactors = FALSE)
names(values)[-(1:2)] <- selection$new_name
values[] <- lapply(values, normalize)
if (anyNA(values$life_id)) stop("Blank participant life_id: cannot reconcile")
if (!nrow(values)) stop("IDmatch_EEG.csv has no rows")
values <- values[, c("life_id", "Column2", sic,
                     setdiff(names(values), c("life_id", "Column2", sic))), drop = FALSE]
orig <- setNames(selection$original_column, selection$new_name)

# Core EEG evidence: actual session date/age and EEG_RUHE/Befund fields.
# EEG_SUBSTANZ variables are acquisition-related but do not by themselves
# prove which EEG recording a tabular row belongs to.
core_eeg <- unique(c(dt, age[!is.na(age)], names(orig)[
  contains(orig, "^Table3[.]EEG_RUHE[_]|^Table3[.]BEFUND_EEG_")]))
core_eeg <- intersect(core_eeg, names(values))
core_eeg <- setdiff(core_eeg, grp)
# All selected EEG-associated measures must follow the primary EEG row.
eeg_session <- unique(c(grp, dt, age[!is.na(age)], names(orig)[
  contains(orig, "^Table3[.]EEG_|^Table3[.]BEFUND_EEG_")]))
eeg_session <- intersect(eeg_session, names(values))
# VIGALL must come ONLY from the primary tabular row, even for close matches.
vigall <- names(orig)[contains(orig, "VIGALL") |
                        contains(names(orig), "^vigall")]
vigall <- intersect(vigall, names(values))
eeg_session <- setdiff(eeg_session, vigall)
identifiers <- c("life_id", "Column2", sic)

# Columns empty in the entire input may be dropped from the final CSV.
all_empty <- setdiff(names(values)[vapply(values, function(x)
  all(is.na(x)), logical(1))], identifiers)


# ---- 2. Prior matching diagnostics and EEG recording paths (NO dates) ------
inventory <- optional_csv("EEG_recordings_inventory.csv")
participants <- optional_csv("EEG_participants_inventory.csv")
prior_duplicates <- optional_csv("IDmatch_duplicates.csv")
prior_issues <- optional_csv("EEG_file_issues.csv")
if (is.null(inventory)) {
  warning("EEG_recordings_inventory.csv missing; recording paths unavailable")
  inventory <- data.frame(life_id = character(), vhdr_file = character(),
                          vhdr_path = character(), selected = character(),
                          ready = character(), suffix = character())
}
if (!all(c("life_id", "vhdr_file", "vhdr_path") %in% names(inventory)))
  stop("EEG_recordings_inventory.csv must have life_id, vhdr_file, vhdr_path")
for (field in c("selected", "ready", "suffix"))
  if (!field %in% names(inventory))
    inventory[[field]] <- rep(NA_character_, nrow(inventory))
inventory$life_id <- normalize_id(inventory$life_id)
inventory$is_2 <- contains(inventory$vhdr_file, "_2[.]vhdr$")
inventory$is_selected <- yes(inventory$selected)
inventory$is_ready <- yes(inventory$ready)

# Previously rejected _2 filenames are included in the *review table* only.
# They are not silently promoted to a validated, selected recording.
extra_files <- data.frame(life_id = character(), vhdr_file = character(),
  vhdr_path = character(), selected = character(), ready = character(),
  suffix = character(), is_2 = logical(), is_selected = logical(),
  is_ready = logical(), stringsAsFactors = FALSE)
if (!is.null(prior_issues) &&
    all(c("issue", "filename", "path") %in% names(prior_issues))) {
  unresolved_files <- prior_issues[prior_issues$issue == "unrecognized_vhdr" &
    contains(prior_issues$filename, "_2[.]vhdr$"), , drop = FALSE]
  for (i in seq_len(nrow(unresolved_files))) {
    f <- unresolved_files$filename[i]
    # e.g. LI12345678_Pilot_2.vhdr -> LI12345678
    proposed <- normalize_id(sub("(?i)_([^_]+)_2[.]vhdr$", "", f, perl = TRUE))
    if (!(proposed %in% values$life_id)) next
    if (any(inventory$life_id == proposed & inventory$vhdr_file == f)) next
    extra_files[nrow(extra_files) + 1L, ] <- list(
      proposed, f, unresolved_files$path[i], "FALSE", "FALSE", "2",
      TRUE, FALSE, FALSE)
  }
}
recordings <- rbind(inventory[, names(extra_files), drop = FALSE], extra_files)


# Manual resolution is only needed for otherwise ambiguous rows.
# No association between _2 and a particular tabular row is inferred from the
# filename alone; the user can explicitly specify that row here.
# Older date-derived overrides are deliberately ignored to preserve unresolved rows.
# Set TRUE only when data_primary_row_overrides.csv has been MANUALLY checked.
use_manual_overrides <- FALSE
overrides <- if (use_manual_overrides) optional_csv("data_primary_row_overrides.csv") else NULL
if (!is.null(overrides)) {
  needed <- c("life_id", "input_csv_row")
  if (!all(needed %in% names(overrides)))
    stop("data_primary_row_overrides.csv must have life_id,input_csv_row")
  overrides$life_id <- normalize_id(overrides$life_id)
  overrides$input_csv_row <- suppressWarnings(as.integer(overrides$input_csv_row))
  if (anyNA(overrides$input_csv_row) || anyNA(overrides$life_id) ||
      anyDuplicated(overrides$life_id))
    stop("Override IDs and row numbers must be present and unique")
  if (any(overrides$input_csv_row < 2L | overrides$input_csv_row > nrow(values) + 1L))
    stop("Invalid input_csv_row in override; first data row is CSV row 2")
  if (any(values$life_id[overrides$input_csv_row - 1L] != overrides$life_id))
    stop("An override input_csv_row does not belong to its life_id")
}
override_row <- function(id) {
  if (is.null(overrides)) return(NA_integer_)
  i <- match(id, overrides$life_id)
  if (is.na(i)) NA_integer_ else overrides$input_csv_row[i]
}

# ---- 3. Row reconciliation -------------------------------------------------
cell_actions <- integer()
record <- function(i) list(id = values$life_id[i],
  vals = unlist(values[i, -1L, drop = FALSE], use.names = TRUE),
  primary = as.integer(i + 1L), sources = as.integer(i + 1L),
  alts = list(), reason = "single_input_row")
# Only count decisions; do not save columns/values in verbose audit tables.
log_cell <- function(main, other, col, before, supplied, action, rule) {
  if (!action %in% names(cell_actions)) cell_actions[action] <<- 1L
  else cell_actions[action] <<- cell_actions[action] + 1L
}
source_ids <- function(r) joined(sort(r$sources))

has_core <- function(r) any(!is.na(r$vals[core_eeg]))
core_count <- function(r) sum(!is.na(r$vals[core_eeg]))
choose_best <- function(rows, id) {
  ov <- override_row(id)
  if (!is.na(ov)) {
    w <- which(vapply(rows, function(r) ov %in% r$sources, logical(1)))
    if (length(w) == 1L) return(w)
    stop("Reviewed override not found in candidate records for ", id)
  }
  # Only used when every row is already mutually compatible outside
  # group/date; tie-breaking here cannot discard a conflicting EEG/VIGALL value.
  available_vig <- vigall[vapply(vigall, function(col)
    any(vapply(rows, function(r) !is.na(r$vals[[col]]), logical(1L))), logical(1L))]
  joint <- function(r) {
    has_core(r) && !is.na(r$vals[[dt]]) &&
      (is.na(age) || !is.na(r$vals[[age]])) &&
      length(available_vig) > 0L && all(!is.na(r$vals[available_vig]))
  }
  scores <- t(vapply(rows, function(r) c(
    joint = as.integer(joint(r)),
    has_eeg = as.integer(has_core(r)),
    date_age = as.integer(!is.na(r$vals[[dt]]) &&
                            (is.na(age) || !is.na(r$vals[[age]]))),
    vig = sum(!is.na(r$vals[vigall])),
    core = core_count(r),
    all = sum(!is.na(r$vals)),
    input = -r$primary), numeric(7L)))
  order(-scores[, "joint"], -scores[, "has_eeg"],
        -scores[, "date_age"], -scores[, "vig"],
        -scores[, "core"], -scores[, "all"], -scores[, "input"])[1L]
}

# Choose a primary EEG/VIGALL row only when it has a *clear* advantage.
# No _2 filename or recording metadata is used to select a tabular row.
select_primary <- function(rows, id) {
  ov <- override_row(id)
  if (!is.na(ov)) return(list(index = choose_best(rows, id),
                              reason = "reviewed_primary_override"))
  has_eeg <- vapply(rows, has_core, logical(1L))
  if (sum(has_eeg) == 1L)
    return(list(index = which(has_eeg), reason = "only_row_with_core_eeg_information"))

  available_vig <- vigall[vapply(vigall, function(col)
    any(vapply(rows, function(r) !is.na(r$vals[[col]]), logical(1L))), logical(1L))]
  full_joint <- vapply(rows, function(r)
    has_core(r) && !is.na(r$vals[[dt]]) &&
      (is.na(age) || !is.na(r$vals[[age]])) &&
      length(available_vig) > 0L && all(!is.na(r$vals[available_vig])), logical(1L))
  if (sum(full_joint) == 1L)
    return(list(index = which(full_joint),
                reason = "only_row_with_complete_eeg_and_vigall"))

  # If one record covers *every* nonmissing EEG/VIGALL field found in the
  # other candidates (and at least one more), it wins. Conflicting VALUES in
  # fields covered by both are logged, never copied from the secondary row.
  signal <- unique(c(core_eeg, vigall))
  cover <- lapply(rows, function(r) !is.na(r$vals[signal]))
  dominated <- vapply(seq_along(rows), function(i) {
    if (sum(has_eeg) && !has_eeg[i]) return(FALSE)
    all(vapply(seq_along(rows), function(j)
      all(!cover[[j]] | cover[[i]]), logical(1L))) &&
      any(vapply(seq_along(rows), function(j)
        any(cover[[i]] & !cover[[j]]), logical(1L)))
  }, logical(1L))
  if (sum(dominated) == 1L)
    return(list(index = which(dominated),
                reason = "dominant_eeg_vigall_completeness"))

  if (!any(has_eeg)) {
    with_vig <- vapply(rows, function(r)
      any(!is.na(r$vals[vigall])), logical(1L))
    if (sum(with_vig) == 1L)
      return(list(index = which(with_vig),
                  reason = "only_row_with_vigall_no_eeg"))
    return(list(index = NA_integer_, reason = "no_clear_eeg_vigall_primary"))
  }
  list(index = NA_integer_, reason = "competing_eeg_vigall_rows_no_clear_primary")
}

# Promote all alternate date/group pairs to separate columns: group2/date_eeg2,
# group3/date_eeg3, ... . Keep GROUP and DATE from a single same source row.
add_alt <- function(main, other) {
  g <- other$vals[[grp]]
  d <- other$vals[[dt]]
  # Only a genuine nonmissing disagreement needs an alternative visit pair.
  # A second row with missing group/date is NOT another recorded session.
  if (!length(conflicting(main$vals, other$vals, c(grp, dt)))) return(main)
  dup <- vapply(main$alts, function(a)
    identical(g, a$group) && identical(d, a$date), logical(1))
  if (any(dup)) return(main)
  main$alts[[length(main$alts) + 1L]] <-
    list(group = g, date = d, source_rows = source_ids(other))
  log_cell(main, other, paste0(grp, "/", dt),
           paste(show(main$vals[[grp]]), show(main$vals[[dt]]), sep = " / "),
           paste(show(g), show(d), sep = " / "),
           "saved_alternative_group_date", "preserve_secondary_visit")
  main
}

# Absorb selected rows without creating synthetic EEG/VIGALL combinations.
# Different nonmissing values are retained from primary and logged.
absorb <- function(main, other, rule, keep_alt = TRUE) {
  if (length(conflicting(main$vals, other$vals, c("Column2", sic))))
    stop("Internal error: merging different identifiers for ", main$id)
  if (keep_alt) main <- add_alt(main, other)
  different_visits <- length(conflicting(main$vals, other$vals, c(grp, dt, age))) > 0L
  for (col in names(main$vals)) {
    a <- main$vals[[col]]
    b <- other$vals[[col]]
    if (is.na(a) && is.na(b)) next
    if (!is.na(a) && !is.na(b)) {
      if (a != b) log_cell(main, other, col, a, b,
                           "kept_primary_conflict", rule)
      next
    }
    if (!is.na(a)) next
    if (col %in% vigall) {
      log_cell(main, other, col, a, b, "skipped_secondary_vigall", rule)
    } else if (col %in% eeg_session && different_visits) {
      log_cell(main, other, col, a, b, "skipped_secondary_eeg_session", rule)
    } else if (col %in% c(grp, dt) && length(main$alts)) {
      # Do not mix alternate visit date/group into the primary paired fields.
      log_cell(main, other, col, a, b, "kept_primary_visit_pair", rule)
    } else {
      main$vals[[col]] <- b
      log_cell(main, other, col, a, b, "filled_from_secondary", rule)
    }
  }
  main$sources <- c(main$sources, other$sources)
  main
}

# Distinct SICs (or conflicting Column2) cannot automatically share a row.
# If there are several known values, missing-identifier rows remain separate.
partition_on <- function(rows, col) {
  found <- unique(vapply(rows, function(r) show(r$vals[[col]]), character(1)))
  found <- found[nzchar(found)]
  if (length(found) < 2L) return(list(rows))
  keys <- vapply(rows, function(r) {
    a <- show(r$vals[[col]])
    if (nzchar(a)) paste0("known:", a) else "MISSING:REVIEW"
  }, character(1))
  unname(split(rows, factor(keys, levels = unique(keys))))
}

# All rows of a partition can be merged via (a) ONLY when no pair has
# conflicting nonmissing cells outside EEG group/date. Missing vs observed
# is complementary, not a conflict. Also covers identical duplicates.
can_merge_group_date <- function(rows) {
  if (length(rows) < 2L) return(TRUE)
  cols <- setdiff(names(rows[[1L]]$vals), c(grp, dt))
  for (i in seq_len(length(rows) - 1L))
    for (j in seq.int(i + 1L, length(rows)))
      if (length(conflicting(rows[[i]]$vals, rows[[j]]$vals, cols)))
        return(FALSE)
  TRUE
}

ids <- unique(values$life_id)
final <- list()
for (id in ids) {
  source_records <- lapply(which(values$life_id == id), record)
  sic_parts <- partition_on(source_records, sic)
  for (sp in sic_parts) for (part in partition_on(sp, "Column2")) {
    if (length(part) == 1L) {
      final[[length(final) + 1L]] <- part[[1L]]
      next
    }
    # Rule A: no conflicts other than group/date, even when some cells are NA.
    # Pick the EEG/VIGALL-complete row; NEVER backfill VIGALL from another row.
    if (can_merge_group_date(part)) {
      best <- choose_best(part, id)
      main <- part[[best]]
      main$reason <- if (is.na(override_row(id)))
        "only_group_date_disagreement_or_complementary" else "reviewed_primary_override"
      for (j in setdiff(seq_along(part), best)) {
        other <- part[[j]]
        rule <- if (!length(changes(main$vals, other$vals, names(main$vals))))
          "exact_duplicate" else "group_date_or_complementary"
        main <- absorb(main, other, rule)
      }
      final[[length(final) + 1L]] <- main
      next
    }
    # Rule B: prefer a row with EEG and VIGALL acquired together.
    # If no uniquely superior row exists, retain the competing records.
    decision <- select_primary(part, id)
    preferred <- decision$index
    reason <- decision$reason
    if (is.na(preferred)) {
      for (r in part) {
        r$reason <- paste0("UNRESOLVED:", reason)
        final[[length(final) + 1L]] <- r
      }
      next
    }
    main <- part[[preferred]]
    main$reason <- reason
    for (j in setdiff(seq_along(part), preferred))
      main <- absorb(main, part[[j]], reason)
    final[[length(final) + 1L]] <- main
  }
}

# Exactly once: each input CSV data row must appear in the final provenance.
final <- final[order(vapply(final, function(r) min(r$sources), integer(1)))]
coverage <- unlist(lapply(final, `[[`, "sources"), use.names = FALSE)
if (!identical(sort(as.integer(coverage)), seq_len(nrow(values)) + 1L))
  stop("AUDIT ERROR: an original input row was lost or used more than once")

# ---- 4. Assemble selected rows (no output yet) -----------------------------
output_cols <- setdiff(names(values), all_empty)
if (!all(identifiers %in% output_cols)) stop("Identifiers must never be dropped")
max_alt <- max(vapply(final, function(r) length(r$alts), integer(1L)))
alt_cols <- if (max_alt) unlist(lapply(seq_len(max_alt), function(k)
  c(paste0(grp, k + 1L), paste0(dt, k + 1L)))) else character()
if (any(alt_cols %in% output_cols)) stop("Alternative names collide with existing fields")
result_data <- as.data.frame(setNames(rep(list(rep(NA_character_, length(final))),
                                length(c(output_cols, alt_cols))),
                                c(output_cols, alt_cols)), check.names = FALSE)
for (i in seq_along(final)) {
  r <- final[[i]]
  result_data$life_id[i] <- r$id
  for (col in setdiff(output_cols, "life_id"))
    result_data[[col]][i] <- r$vals[[col]]
  if (length(r$alts)) for (j in seq_along(r$alts)) {
    result_data[[paste0(grp, j + 1L)]][i] <- r$alts[[j]]$group
    result_data[[paste0(dt, j + 1L)]][i] <- r$alts[[j]]$date
  }
}
# Remove fully empty alternates, too.
if (length(alt_cols)) {
  alt_empty <- alt_cols[vapply(result_data[alt_cols], function(x)
    all(is.na(x)), logical(1L))]
  result_data <- result_data[, setdiff(names(result_data), alt_empty), drop = FALSE]
}

# Remove the entire requested EP-EW interval in the manually defined order.
# Include any variables newly selected between these two endpoints in future maps.
start <- which(manual$new_name == "coronary_disease_angina_pectoris_2")
end <- which(manual$new_name == "diabetes_question_2")
if (length(start) != 1L || length(end) != 1L || end < start)
  stop("Cannot identify EP-EW boundaries in idmatch_variables_manual.csv")
block_sources <- manual$original_column[seq.int(start, end)]
block <- selection$new_name[selection$original_column %in% block_sources]
result_data <- result_data[, !names(result_data) %in% block, drop = FALSE]

# non-unique describes the FINAL table, not repetitions in the initial IDmatch.
# This is a deliberately non-syntactic CSV header; use [[...]] in R.
id_counts <- table(result_data$life_id)
result_data[["non-unique"]] <- as.integer(id_counts[result_data$life_id] > 1L)
result_data <- result_data[, c("life_id", "Column2", sic, "non-unique",
  setdiff(names(result_data), c("life_id", "Column2", sic, "non-unique"))), drop = FALSE]

# ---- 5. Safe matching of participant IDs across external data sources -----
# Exact matching first. For unmatched IDs, allow LIFE's documented LI-prefix /
# trailing suffix corrections, but NEVER guess when a fallback is ambiguous.
main_ids <- unique(result_data$life_id)
canon_id <- function(x) {
  x <- normalize_id(x)
  x <- sub("(?i)_(Ruhe|Pilot)(?:_[0-9]+)?$", "", x, perl = TRUE)
  x <- sub("[_-]+2$", "", x, perl = TRUE)
  x <- sub("[_-]+$", "", x, perl = TRUE)
  sub("^LI", "", x)
}
canonical_main <- canon_id(main_ids)
issues <- list()
add_issue <- function(source, participant, type, detail) {
  issues[[length(issues) + 1L]] <<- data.frame(
    source = source, life_id = participant, issue = type,
    details = detail, stringsAsFactors = FALSE)
}
resolve_ids <- function(x, source) {
  raw_id <- normalize_id(x)
  ans <- match(raw_id, main_ids)
  missing <- which(!is.na(raw_id) & is.na(ans))
  for (i in missing) {
    candidates <- which(!is.na(canonical_main) &
                        canonical_main == canon_id(raw_id[i]))
    if (length(candidates) == 1L) ans[i] <- candidates
    if (length(candidates) > 1L)
      add_issue(source, raw_id[i], "ambiguous_id_alias",
                paste(main_ids[candidates], collapse = " | "))
  }
  ans
}
add_column <- function(name, indexed_values) {
  if (name %in% names(result_data)) stop("New column already exists: ", name)
  result_data[[name]] <<- indexed_values[match(result_data$life_id, main_ids)]
}
coverage <- list()
add_coverage <- function(source, column, matched, nonmissing, duplicates, note = "") {
  coverage[[length(coverage) + 1L]] <<- data.frame(
    source = source, column = column, participants_found = matched,
    participants_with_value = nonmissing, source_duplicate_ids = duplicates,
    note = note, stringsAsFactors = FALSE)
}

# ---- 6. Group CSVs: HEP required, ECG optional ---------------------------
# Filenames, rather than total CSV count, distinguish the 16 HEP pipelines
# from their optional ECG counterparts. Other CSVs in the folder are ignored.
group_dir <- file.path(life, "derivatives/group")
if (!dir.exists(group_dir)) stop("Group directory missing: ", group_dir)
group_files <- sort(list.files(group_dir, pattern = "[.]csv$", full.names = TRUE,
                               ignore.case = TRUE, recursive = FALSE))
hep_files <- group_files[grepl("(?i)roi_HEP[.]csv$", basename(group_files), perl = TRUE)]
ecg_files <- group_files[grepl("(?i)roi_ECG[.]csv$", basename(group_files), perl = TRUE)]
if (!length(hep_files)) stop("No *roi_HEP.csv files in ", group_dir)
if (length(hep_files) != 16L)
  warning("Expected 16 roi_HEP.csv files, found ", length(hep_files))
if (!length(ecg_files))
  message("No roi_ECG.csv files found; importing HEP only.")
if (length(ecg_files) && length(ecg_files) != length(hep_files))
  warning("Found ", length(ecg_files), " ECG files and ", length(hep_files),
          " HEP files; importing the available files independently.")

hep_names <- sub("(?i)[.]csv$", "", basename(hep_files), perl = TRUE)
ecg_names <- sub("(?i)[.]csv$", "", basename(ecg_files), perl = TRUE)
group_names <- c(hep_names, ecg_names)
if (anyDuplicated(group_names) || any(group_names %in% names(result_data)))
  stop("Group filename-derived column names collide with existing columns")

read_group_roi <- function(filename, column_name, label) {
  tab <- read_csv(filename)
  if (!all(c("subject_id", column_name) %in% names(tab)))
    stop("Group CSV ", basename(filename), " needs subject_id and ", column_name)
  keys <- resolve_ids(tab$subject_id, label)
  roi <- normalize(tab[[column_name]])
  vals <- rep(NA_character_, length(main_ids))
  matched <- unique(keys[!is.na(keys)])
  repeated <- 0L
  for (k in matched) {
    inds <- which(keys == k)
    if (length(inds) > 1L) repeated <- repeated + 1L
    nonmissing <- unique(roi[inds][!is.na(roi[inds])])
    if (length(nonmissing) == 1L) vals[k] <- nonmissing
    if (length(nonmissing) > 1L)
      add_issue(label, main_ids[k], paste0("conflicting_", column_name, "_duplicates"),
                paste(length(inds), "source rows; left missing"))
  }
  add_column(label, vals)
  add_coverage(label, column_name, length(matched), sum(!is.na(vals)), repeated)
}
for (i in seq_along(hep_files))
  read_group_roi(hep_files[i], "roi_HEP", hep_names[i])
for (i in seq_along(ecg_files))
  read_group_roi(ecg_files[i], "roi_ECG", ecg_names[i])

# ---- 7. Preprocessing log ODS: later a5 supersedes a1/a2 -------------------
if (!requireNamespace("readODS", quietly = TRUE))
  stop("Please install readODS before running: install.packages('readODS')")
pre_path <- file.path(life, "tabular/preprocessing_logs.ods")
if (!file.exists(pre_path)) stop("Missing preprocessing log: ", pre_path)
pre <- as.data.frame(readODS::read_ods(pre_path, sheet = 1),
                     check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("subject_id", "step", "hp_0.1", "hp_0.5") %in% names(pre)))
  stop("preprocessing_logs.ods must have subject_id, step, hp_0.1, hp_0.5; found: ",
       paste(names(pre), collapse = ", "))
pre[] <- lapply(pre, normalize)
pre_keys <- resolve_ids(pre$subject_id, "preprocessing_logs.ods")
pre_step <- tolower(pre$step)
pre_matched <- unique(pre_keys[!is.na(pre_keys)])
pre_repeated <- sum(table(pre_keys[!is.na(pre_keys)]) > 1L)
for (src in c("hp_0.1", "hp_0.5")) {
  v <- rep(NA_character_, length(main_ids))
  for (k in pre_matched) {
    valid <- which(pre_keys == k & !is.na(pre[[src]]))
    if (!length(valid)) next
    # Prefer the last nonmissing a5 result; otherwise the last nonmissing row.
    a5 <- valid[!is.na(pre_step[valid]) & pre_step[valid] == "a5"]
    chosen <- if (length(a5)) tail(a5, 1L) else tail(valid, 1L)
    v[k] <- pre[[src]][chosen]
  }
  new <- paste0("preproc_", src)
  add_column(new, v)
  add_coverage("preprocessing_logs.ods", new, length(pre_matched),
               sum(!is.na(v)), pre_repeated, "last available; a5 takes precedence")
}

# ---- 8. CER-S: read the correctly separated XLSX, never repeated headers ---
if (!requireNamespace("readxl", quietly = TRUE))
  stop("Please install readxl: install.packages('readxl')")
ecga_path <- file.path(life, "tabular/AMPS_CER-S/ECGAna_allSubj.xlsx")
if (!file.exists(ecga_path)) stop("Missing CER-S Excel table: ", ecga_path)
ecga <- as.data.frame(readxl::read_excel(ecga_path, sheet = 1,
  col_types = "text", .name_repair = "minimal"),
  check.names = FALSE, stringsAsFactors = FALSE)
names(ecga) <- trimws(names(ecga))
ecga[] <- lapply(ecga, normalize)
beats <- c("N Beats", "V Beats", "S Beats")
if (!all(beats %in% names(ecga)))
  stop("CER-S Excel missing beat columns: ", paste(setdiff(beats, names(ecga)), collapse = ", "))

# An occasional embedded header ("N Beats", "V Beats", "S Beats") or
# 'Excluded' status is NOT a beat count. Accept only nonnegative integers.
clean_beat <- function(x) {
  x <- normalize(x)
  numeric_cell <- !is.na(x) & grepl("^[0-9]+(?:\\.0+)?$", x, perl = TRUE)
  out <- rep(NA_character_, length(x))
  out[numeric_cell] <- sub("\\.0+$", "", x[numeric_cell])
  out
}
invalid_beats <- setNames(vapply(beats, function(col)
  sum(!is.na(ecga[[col]]) & is.na(clean_beat(ecga[[col]]))), integer(1L)), beats)
for (col in beats) ecga[[col]] <- clean_beat(ecga[[col]])
# Automatically locate the filename column; preserve its original header.
file_col <- names(ecga)[vapply(ecga, function(x) {
  x <- basename(normalize(x))
  any(!is.na(x) & grepl("(?i)^ECG_.+[.]hea[.]csv$", x, perl = TRUE))
}, logical(1L))]
if (length(file_col) != 1L)
  stop("Expected exactly one CER-S filename column with ECG_LI..._Ruhe.hea.csv; ",
       "candidates: ", paste(file_col, collapse = ", "))
ecga_files <- basename(normalize(ecga[[file_col]]))
ecga_ids <- sub("(?i)[.]hea[.]csv$", "", ecga_files, perl = TRUE)
ecga_ids <- sub("(?i)^ECG_", "", ecga_ids, perl = TRUE)
ecga_ids <- sub("(?i)_(Ruhe|Pilot)(?:_[0-9]+)?$", "", ecga_ids, perl = TRUE)
valid_ecga <- !is.na(ecga_files) &
  grepl("(?i)^ECG_.+[.]hea[.]csv$", ecga_files, perl = TRUE)
ecga_keys <- rep(NA_integer_, nrow(ecga))
ecga_keys[valid_ecga] <- resolve_ids(ecga_ids[valid_ecga], "ECGAna_allSubj.xlsx")
matched_ecga <- unique(ecga_keys[!is.na(ecga_keys)])
ecga_repeated <- sum(table(ecga_keys[!is.na(ecga_keys)]) > 1L)
beat_values <- as.data.frame(setNames(lapply(beats, function(x)
  rep(NA_character_, length(main_ids))), beats), check.names = FALSE)
for (k in matched_ecga) {
  inds <- which(ecga_keys == k)
  # Skip embedded header / excluded-only records when a numeric record exists.
  usable <- inds[rowSums(!is.na(ecga[inds, beats, drop = FALSE])) > 0L]
  if (!length(usable)) {
    add_issue("ECGAna_allSubj.xlsx", main_ids[k], "no_numeric_beat_counts",
              joined(ecga_files[inds]))
    next
  }
  inds <- usable
  # When several files exist, prefer _Ruhe_2, then _Ruhe; do not use arbitrary
  # source row ordering if equally preferred files report different results.
  is2 <- grepl("(?i)_Ruhe_2[.]hea[.]csv$", ecga_files[inds], perl = TRUE)
  isruhe <- grepl("(?i)_Ruhe(?:_[0-9]+)?[.]hea[.]csv$",
                  ecga_files[inds], perl = TRUE)
  priority <- as.integer(is2) * 10L + as.integer(isruhe)
  chosen <- inds[priority == max(priority)]
  if (length(chosen) > 1L) {
    v <- ecga[chosen, beats, drop = FALSE]
    v[] <- lapply(v, normalize)
    if (any(vapply(v, function(x)
      length(unique(x[!is.na(x)])) > 1L, logical(1L)))) {
      add_issue("ECGAna_allSubj.xlsx", main_ids[k], "conflicting_equal_priority_ecg_files",
                paste(basename(ecga_files[chosen]), collapse = " | "))
      next
    }
    # Equal-priority duplicates with compatible values: choose the most complete.
    completeness <- rowSums(!is.na(v))
    chosen <- chosen[which.max(completeness)]
  }
  for (col in beats) beat_values[[col]][k] <- normalize(ecga[[col]][chosen])
}
for (j in seq_along(beats)) {
  new <- c("ecg_n_beats", "ecg_v_beats", "ecg_s_beats")[j]
  vals <- beat_values[[beats[j]]]
  add_column(new, vals)
  add_coverage("ECGAna_allSubj.xlsx", new, length(matched_ecga),
               sum(!is.na(vals)), ecga_repeated,
               "_Ruhe_2 preferred for duplicated ECG filenames")
}

# ---- 9. GAD-7: match questionnaire SIC to Column2/sic, then IDmatch alias --
if (!requireNamespace("readxl", quietly = TRUE))
  stop("Please install readxl: install.packages('readxl')")
gad_path <- file.path(life, "tabular/Behavior/data/PV0491_T00019_NODUP.xlsx")
if (!file.exists(gad_path)) stop("Missing GAD-7 workbook: ", gad_path)
# Read identifiers as TEXT to preserve leading zeroes and exact ID spelling.
gad <- as.data.frame(readxl::read_excel(gad_path, sheet = 1,
  col_types = "text", .name_repair = "minimal"),
  check.names = FALSE, stringsAsFactors = FALSE)
names(gad) <- trimws(names(gad))
gad_items <- paste0("GAD7_GAD7_", seq_len(7L))
gad_columns <- c(gad_items, "GAD7_DATUM")
if (!all(c("GAD7_SIC", gad_columns) %in% names(gad)))
  stop("GAD-7 workbook missing columns: ",
       paste(setdiff(c("GAD7_SIC", gad_columns), names(gad)), collapse = ", "))
gad[] <- lapply(gad, normalize)

# Some Excel date cells imported as text are actually Excel serial numbers,
# e.g. 40604.5060416667 (the fractional part is time of day). Export as
# YYYY-MM-DD using the conventional Excel 1900 date system.
format_gad_date <- function(x) {
  x <- normalize(x)
  output <- x
  serial <- !is.na(x) & grepl("^[0-9]+(?:\\.[0-9]+)?$", x, perl = TRUE)
  if (any(serial)) {
    ii <- which(serial)
    n <- suppressWarnings(as.numeric(x[ii]))
    valid <- is.finite(n) & n >= 20000 & n <= 80000
    output[ii[valid]] <- format(as.Date(floor(n[valid]),
                                   origin = "1899-12-30"), "%Y-%m-%d")
    output[ii[!valid]] <- NA_character_
  }
  other <- which(!is.na(x) & !serial)
  for (fmt in c("%Y-%m-%d", "%d.%m.%Y", "%d/%m/%Y", "%Y/%m/%d")) {
    if (!length(other)) break
    # The first ten characters also cover ISO datetime strings.
    trial <- suppressWarnings(as.Date(substr(x[other], 1L, 10L), format = fmt))
    good <- which(!is.na(trial))
    if (length(good)) output[other[good]] <- format(trial[good], "%Y-%m-%d")
    other <- other[-good]
  }
  if (length(other)) warning("Unrecognized GAD7_DATUM value(s): ",
                              paste(unique(x[other]), collapse = " | "))
  output[other] <- NA_character_
  output
}
gad$GAD7_DATUM <- format_gad_date(gad$GAD7_DATUM)
gad_key <- normalize_id(gad$GAD7_SIC)

# Include identifiers from ALL ORIGINAL tabular rows, even if they were later
# absorbed into one selected row. Never use life_id itself as a GAD7_SIC alias.
direct_aliases <- lapply(main_ids, function(id) {
  x <- values[values$life_id == id, c("Column2", sic), drop = FALSE]
  unique(na.omit(normalize_id(unlist(x, use.names = FALSE))))
})

# Direct matches to either Column2 or SIC have priority over any fallback.
gad_hits <- lapply(direct_aliases, function(keys)
  which(!is.na(gad_key) & gad_key %in% keys))
gad_method <- ifelse(lengths(gad_hits) > 0L, "direct_Column2_or_sic", "unmatched")

# The IDmatch sheet may contain an additional alias in its second column.
# Use it only for LIFE IDs without a direct questionnaire match.
fallback_aliases <- replicate(length(main_ids), character(), simplify = FALSE)
need_fallback <- which(lengths(gad_hits) == 0L)
if (length(need_fallback)) {
  idmatch_path <- file.path(life, "tabular/IDmatch.xlsx")
  if (file.exists(idmatch_path) && "IDmatch" %in% readxl::excel_sheets(idmatch_path)) {
    alt <- as.data.frame(readxl::read_excel(idmatch_path,
      sheet = "IDmatch", col_types = "text", .name_repair = "minimal"),
      check.names = FALSE, stringsAsFactors = FALSE)
    if (ncol(alt) < 2L) stop("IDmatch.xlsx / IDmatch must have at least two columns")
    alias_life <- normalize_id(alt[[1L]])
    alias_value <- normalize_id(alt[[2L]])
    # Earlier LIFE_matching used exact IDs, then unique LI-prefix/suffix cleanup.
    # Apply the same restricted fallback only when its result is unambiguous.
    pos <- match(alias_life, main_ids)
    for (j in which(!is.na(alias_life) & is.na(pos))) {
      hits <- which(!is.na(canonical_main) &
                    canonical_main == canon_id(alias_life[j]))
      if (length(hits) == 1L) pos[j] <- hits
    }
    for (j in which(!is.na(pos) & !is.na(alias_value))) {
      k <- pos[j]
      fallback_aliases[[k]] <- unique(c(fallback_aliases[[k]], alias_value[j]))
    }
    for (k in need_fallback) {
      hits <- which(!is.na(gad_key) & gad_key %in% fallback_aliases[[k]])
      if (length(hits)) {
        gad_hits[[k]] <- hits
        gad_method[k] <- "IDmatch_sheet_alias"
      }
    }
  } else {
    warning("GAD-7 fallback unavailable: missing IDmatch.xlsx or 'IDmatch' sheet")
  }
}

# A questionnaire row must not be claimed by multiple LIFE participants.
claimed <- unlist(gad_hits, use.names = FALSE)
claim_counts <- tabulate(claimed, nbins = nrow(gad))
gad_selected <- rep(NA_integer_, length(main_ids))
for (k in seq_along(main_ids)) {
  hit <- unique(gad_hits[[k]])
  if (!length(hit)) next
  if (any(claim_counts[hit] > 1L)) {
    add_issue("GAD7", main_ids[k], "shared_questionnaire_identifier",
              "At least one GAD7_SIC matches multiple LIFE participants")
    next
  }
  if (length(hit) > 1L) {
    add_issue("GAD7", main_ids[k], "multiple_questionnaire_rows",
              paste("Matching Excel rows:", paste(hit + 1L, collapse = " | ")))
    next
  }
  gad_selected[k] <- hit
}

for (col in gad_columns) {
  matched_values <- rep(NA_character_, length(main_ids))
  has_match <- which(!is.na(gad_selected))
  if (length(has_match)) matched_values[has_match] <- gad[[col]][gad_selected[has_match]]
  add_column(col, matched_values)  # Same value on EVERY final row of a LIFE ID.
  add_coverage("PV0491_T00019_NODUP.xlsx", col, length(has_match),
               sum(!is.na(matched_values)), sum(table(gad_key[!is.na(gad_key)]) > 1L),
               "GAD7_SIC via Column2/sic; IDmatch alias only if unmatched")
}
gad_direct_n <- sum(!is.na(gad_selected) & gad_method == "direct_Column2_or_sic")
gad_fallback_n <- sum(!is.na(gad_selected) & gad_method == "IDmatch_sheet_alias")
gad_ambiguous_n <- sum(lengths(gad_hits) > 0L & is.na(gad_selected))

# ---- 10. HR/HRV: import all measurements by participant ID ---------------
# Keep the names from the source file; only prefix when a name already exists.
hrv_path <- file.path(life, "derivatives/HR_HRV_all_participants.csv")
if (!file.exists(hrv_path)) stop("Missing HR/HRV file: ", hrv_path)
hrv <- read_csv(hrv_path)
names(hrv) <- trimws(names(hrv))
if (anyDuplicated(names(hrv))) stop("Duplicate column names in HR/HRV file")
hrv_id_candidates <- c("subject_id", "life_id", "participant_id", "SubjectID", "subject", "Subject", "ID", "id")
hrv_id_col <- intersect(hrv_id_candidates, names(hrv))
if (!length(hrv_id_col)) {
  warning("No recognized HR/HRV ID header; assuming the FIRST column is the ID: ",
          names(hrv)[1L])
  hrv_id_col <- names(hrv)[1L]
}
hrv_id_col <- hrv_id_col[1L]
hrv_keys <- resolve_ids(hrv[[hrv_id_col]], "HR_HRV_all_participants.csv")
hrv_matched <- unique(hrv_keys[!is.na(hrv_keys)])
if (!length(hrv_matched))
  warning("No HR/HRV identifiers matched LIFE participant IDs. Check ID format and source.")
hrv_repeated <- sum(table(hrv_keys[!is.na(hrv_keys)]) > 1L)
hrv_fields <- setdiff(names(hrv), hrv_id_col)
if (!length(hrv_fields)) warning("HR/HRV table has no measurement columns")
hrv_names <- setNames(character(length(hrv_fields)), hrv_fields)
for (src in hrv_fields) {
  label <- if (src %in% names(result_data)) paste0("HR_HRV_", src) else src
  if (label %in% names(result_data) || label %in% hrv_names)
    stop("Ambiguous HR/HRV output column: ", label)
  hrv_names[src] <- label
  vals <- rep(NA_character_, length(main_ids))
  original <- normalize(hrv[[src]])
  for (k in hrv_matched) {
    matches <- which(hrv_keys == k)
    nonmissing <- unique(original[matches][!is.na(original[matches])])
    if (length(nonmissing) == 1L) vals[k] <- nonmissing
    if (length(nonmissing) > 1L)
      add_issue("HR_HRV_all_participants.csv", main_ids[k],
                "conflicting_duplicate_measurements",
                paste(src, ": ", length(nonmissing), " different values; left missing"))
  }
  add_column(label, vals)
  add_coverage("HR_HRV_all_participants.csv", label, length(hrv_matched),
               sum(!is.na(vals)), hrv_repeated)
}

# ---- 11. Derived participant-level values ---------------------------------
# Scores are calculated for each output row; external measurements already
# repeat identically across unresolved duplicate rows of the same LIFE ID.
# All questionnaire items must be coded 0,1,2,3. Unexpected values become NA.
item_numeric <- function(x) {
  x <- normalize(x)
  n <- suppressWarnings(as.numeric(gsub(",", ".", x, fixed = TRUE)))
  n[!is.finite(n) | is.na(n) | !(n %in% 0:3)] <- NA_real_
  n
}
positive_number <- function(x) {
  n <- suppressWarnings(as.numeric(gsub(",", ".", normalize(x), fixed = TRUE)))
  n[!is.finite(n)] <- NA_real_
  n
}
nonmissing_row_sum <- function(x) {
  out <- rowSums(x, na.rm = TRUE)
  out[rowSums(!is.na(x)) == 0L] <- NA_real_
  out
}

# GAD-7: requested available-item sum, with a separate preregistration validity
# flag. A partially observed score is NOT valid for the registered models.
if (!all(gad_items %in% names(result_data)))
  stop("GAD-7 items missing after questionnaire matching")
gad_mat <- as.data.frame(lapply(result_data[gad_items], item_numeric))
gad_valid_count <- rowSums(!is.na(gad_mat))
result_data$all_GAD7 <- as.integer(gad_valid_count == 7L)
result_data$GAD7 <- nonmissing_row_sum(gad_mat)
result_data$gad7_group <- ifelse(result_data$all_GAD7 == 1L &
                                  result_data$GAD7 <= 4, "low",
                           ifelse(result_data$all_GAD7 == 1L &
                                    result_data$GAD7 >= 10, "high", NA_character_))

# IDS-SR30: score ONE value from each alternative appetite and weight pair,
# using max(item 11, item 12), max(item 13, item 14).
# LIFE-specific mapping verified by the data owner: IDS_I09A is scored item 9
# (severity of mood variation, 0-3). IDS_I09B is NOT included in IDS total.
# Valid IDS-SR total requires all 28 scored items (range 0-84).
ids_regular <- sprintf("IDS_I%02d", c(1:10, 15:30))  # 26 ordinary items
ids_slots <- c(ids_regular, "IDS_appetite", "IDS_weight") # 28 scored slots
ids_mat <- matrix(NA_real_, nrow = nrow(result_data), ncol = length(ids_slots),
                  dimnames = list(NULL, ids_slots))
for (field in intersect(setdiff(ids_regular, "IDS_I09"), names(result_data)))
  ids_mat[, field] <- item_numeric(result_data[[field]])
# The input field IDS_I09A represents the standard scored item 9 in LIFE.
if ("IDS_I09A" %in% names(result_data)) ids_mat[, "IDS_I09"] <- item_numeric(result_data[["IDS_I09A"]]) else warning("IDS_I09A is missing")
ids_pair <- function(a, b) {
  a <- if (a %in% names(result_data)) item_numeric(result_data[[a]]) else
    rep(NA_real_, nrow(result_data))
  b <- if (b %in% names(result_data)) item_numeric(result_data[[b]]) else
    rep(NA_real_, nrow(result_data))
  out <- pmax(a, b, na.rm = TRUE)
  out[is.infinite(out)] <- NA_real_
  out
}
ids_mat[, "IDS_appetite"] <- ids_pair("IDS_I11", "IDS_I12")
ids_mat[, "IDS_weight"] <- ids_pair("IDS_I13", "IDS_I14")
result_data$IDS_n_items <- rowSums(!is.na(ids_mat))
result_data$all_IDS <- as.integer(result_data$IDS_n_items == 28L)
result_data$IDS <- nonmissing_row_sum(ids_mat)

# Extrasystoles: exactly (S + V) and (S + V) / N, as requested.
# Do not treat missing beat counts as zero or divide by zero.
required_beats <- c("ecg_s_beats", "ecg_v_beats", "ecg_n_beats")
if (!all(required_beats %in% names(result_data)))
  stop("Cannot calculate extrasystoles; beat columns are missing")
s <- positive_number(result_data$ecg_s_beats)
v <- positive_number(result_data$ecg_v_beats)
n <- positive_number(result_data$ecg_n_beats)
result_data$extrasystolic_N <- ifelse(
  !is.na(s) | !is.na(v),
  rowSums(cbind(s, v), na.rm = TRUE),
  NA_real_
)
result_data$extrasystolic_proportion <- ifelse(!is.na(result_data$extrasystolic_N) &
                                               !is.na(n) & n > 0,
                                               result_data$extrasystolic_N / n, NA_real_)

# Extrasystole exclusion: > 4 SD ABOVE the sample mean (one-sided), as preregistered.
# Count each participant once when estimating the sample mean and SD;
# assign the same score and exclusion flag to all unresolved rows of that ID.
extra_unique <- !duplicated(result_data$life_id)
extra_reference <- result_data$extrasystolic_N[extra_unique]
extra_mean <- if (any(!is.na(extra_reference)))
  mean(extra_reference, na.rm = TRUE) else NA_real_
extra_sd <- if (sum(!is.na(extra_reference)) >= 2L)
  sd(extra_reference, na.rm = TRUE) else NA_real_
result_data$extrasystolic_N_SD <- if (is.finite(extra_sd) && extra_sd > 0)
  (result_data$extrasystolic_N - extra_mean) / extra_sd else
  rep(NA_real_, nrow(result_data))
result_data$extrasystolic_exclusion <- ifelse(
  is.na(result_data$extrasystolic_N_SD), NA_integer_,
  as.integer(result_data$extrasystolic_N_SD > 4)
)

# HEP z-score in SD units and |z| > 3 outlier flag for each pipeline.
# Unresolved duplicated LIFE IDs do NOT get extra statistical weight.
# Derived values are replicated onto every row of the same participant.
hep_outlier_counts <- setNames(integer(length(hep_names)), hep_names)
for (label in hep_names) {
  src <- positive_number(result_data[[label]])
  unique_pos <- !duplicated(result_data$life_id)
  reference <- src[unique_pos]
  mu <- if (any(!is.na(reference))) mean(reference, na.rm = TRUE) else NA_real_
  sigma <- if (sum(!is.na(reference)) >= 2L) sd(reference, na.rm = TRUE) else NA_real_
  z <- if (is.finite(sigma) && sigma > 0) (src - mu) / sigma else
    rep(NA_real_, length(src))
  z_name <- paste0(label, "_SD")
  flag_name <- paste0(label, "_outliers")
  if (z_name %in% names(result_data) || flag_name %in% names(result_data))
    stop("Derived HEP column name collision: ", label)
  result_data[[z_name]] <- z
  result_data[[flag_name]] <- ifelse(is.na(z), NA_integer_, as.integer(abs(z) > 3))
  hep_outlier_counts[label] <- sum(result_data[[flag_name]][unique_pos] == 1L,
                                  na.rm = TRUE)
}

# ---- 12. Save ONLY useful compact outputs ------------------------------
if (anyDuplicated(names(result_data))) stop("Duplicate output column names")
if (nrow(result_data) != length(final)) stop("Internal error: row count changed")
if (!identical(result_data[["non-unique"]],
               as.integer(table(result_data$life_id)[result_data$life_id] > 1L)))
  stop("Internal error: non-unique flag does not match final IDs")
# Analysis data lives at derivatives/data.csv; diagnostics stay in out/.
data_path <- file.path(life, "derivatives", "data.csv")
write.csv(result_data, data_path, na = "", row.names = FALSE)

# Row provenance: original source row numbers + choices, without all cell values.
provenance <- do.call(rbind, lapply(seq_along(final), function(i) {
  r <- final[[i]]
  data.frame(data_csv_row = i + 1L, life_id = r$id,
    non_unique = result_data[["non-unique"]][i],
    primary_input_csv_row = r$primary, all_input_csv_rows = source_ids(r),
    n_input_rows_merged = length(r$sources), selection_rule = r$reason,
    n_alternative_eeg_dates = length(r$alts), stringsAsFactors = FALSE)
}))
write_csv(provenance, "data_row_provenance.csv")

file_info <- function(id) {
  f <- recordings[!is.na(recordings$life_id) & recordings$life_id == id, , drop = FALSE]
  list(paths = joined(f$vhdr_path),
       selected = joined(f$vhdr_path[f$is_selected]))
}
unresolved <- list()
for (id in unique(result_data$life_id[result_data[["non-unique"]] == 1L])) {
  positions <- which(result_data$life_id == id)
  f <- file_info(id)
  unresolved[[length(unresolved) + 1L]] <- data.frame(
    life_id = id, n_rows = length(positions),
    data_csv_rows = joined(positions + 1L),
    primary_input_csv_rows = joined(provenance$primary_input_csv_row[positions]),
    eeg_dates = joined(paste0(positions + 1L, ": ",
                       ifelse(is.na(result_data[[dt]][positions]), "<missing>",
                              result_data[[dt]][positions]))),
    eeg_groups = joined(paste0(positions + 1L, ": ",
                        ifelse(is.na(result_data[[grp]][positions]), "<missing>",
                               result_data[[grp]][positions]))),
    recording_paths = f$paths, selected_recording_path = f$selected,
    selection_rules = joined(provenance$selection_rule[positions]),
    stringsAsFactors = FALSE)
}
unresolved_cols <- c("life_id", "n_rows", "data_csv_rows", "primary_input_csv_rows",
  "eeg_dates", "eeg_groups", "recording_paths", "selected_recording_path",
  "selection_rules")
write_csv(bind_events(unresolved, unresolved_cols), "data_unresolved_rows.csv")
# Dates and candidate SICs/Excel row numbers help manual matching, without
# dumping full source rows into the audit. All values reflect original matches.
issue_cols <- c("source", "life_id", "issue", "details")
issue_table <- bind_events(issues, issue_cols)
issue_table$eeg_date <- rep(NA_character_, nrow(issue_table))
issue_table$gad7_date <- rep(NA_character_, nrow(issue_table))
issue_table$gad7_sic <- rep(NA_character_, nrow(issue_table))
issue_table$gad7_excel_rows <- rep(NA_character_, nrow(issue_table))
for (i in seq_len(nrow(issue_table))) {
  id <- issue_table$life_id[i]
  k <- match(id, main_ids)
  if (!is.na(k)) {
    # Include dates across unresolved LIFE rows, plus saved alternative visits.
    positions <- which(result_data$life_id == id)
    date_fields <- intersect(unique(c(dt, grep(paste0("^", dt, "[2-9][0-9]*$"),
                                       names(result_data), value = TRUE))),
                             names(result_data))
    issue_table$eeg_date[i] <- joined(unlist(result_data[positions, date_fields,
                                                  drop = FALSE], use.names = FALSE))
    candidates <- unique(gad_hits[[k]])
  } else if (issue_table$source[i] == "GAD7") {
    # For a questionnaire-side identifier that could not resolve to LIFE.
    candidates <- which(!is.na(gad_key) & gad_key == normalize_id(id))
  } else {
    candidates <- integer()
  }
  if (length(candidates)) {
    issue_table$gad7_date[i] <- joined(gad$GAD7_DATUM[candidates])
    issue_table$gad7_sic[i] <- joined(gad_key[candidates])
    issue_table$gad7_excel_rows[i] <- joined(candidates + 1L)
  }
}
write_csv(issue_table, "data_match_issues.csv")

summary_coverage <- bind_events(coverage, c("source", "column", "participants_found",
                                            "participants_with_value", "source_duplicate_ids", "note"))
count_cells <- function(action) if (action %in% names(cell_actions))
  unname(cell_actions[[action]]) else 0L
n_original_repeated <- sum(table(values$life_id) > 1L)
n_final_repeated <- sum(table(result_data$life_id) > 1L)
report <- capture.output({
  cat("========== LIFE analysis data preparation ==========\n")
  cat("Matched IDmatch input rows:         ", nrow(values), "\n")
  cat("Unique participants:                ", length(main_ids), "\n")
  cat("Participants with repeated inputs: ", n_original_repeated, "\n")
  cat("Final rows in data.csv:             ", nrow(result_data), "\n")
  cat("Participants with one final row:    ", sum(table(result_data$life_id) == 1L), "\n")
  cat("Participants still non-unique:      ", n_final_repeated, "\n")
  cat("Original rows absorbed:            ", nrow(values) - nrow(result_data), "\n")
  cat("Complementary cells filled:        ", count_cells("filled_from_secondary"), "\n")
  cat("Primary conflicting cells kept:    ", count_cells("kept_primary_conflict"), "\n")
  cat("Secondary VIGALL cells NOT used:   ", count_cells("skipped_secondary_vigall"), "\n")
  cat("Other session EEG cells NOT used:  ", count_cells("skipped_secondary_eeg_session"), "\n")
  cat("Alternative group/date pairs saved:", count_cells("saved_alternative_group_date"), "\n")
  cat("Fully empty selected columns dropped:", length(all_empty), "\n")
  cat("MEDANAM EP-EW columns removed:       ", paste(block, collapse = ", "), "\n")
  cat("Final output columns:               ", ncol(result_data), "\n")
  cat("Group HEP CSVs imported:             ", length(hep_files), "\n")
  cat("Group ECG CSVs imported:             ", length(ecg_files), "\n")
  cat("HR/HRV columns imported:             ", length(hrv_fields), "\n")
  cat("GAD-7 complete (unique participants):", sum(result_data$all_GAD7[!duplicated(result_data$life_id)] == 1L), "\n")
  cat("IDS complete (unique participants):  ", sum(result_data$all_IDS[!duplicated(result_data$life_id)] == 1L), "\n")
  cat("GAD-7 groups high/low (unique participants):\n")
  print(table(result_data$gad7_group[!duplicated(result_data$life_id)], useNA = "ifany"))
  cat("Extrasystole count mean (unique participants):", extra_mean, "\n")
  cat("Extrasystole count SD (unique participants):  ", extra_sd, "\n")
  cat("Extrasystole >4-SD exclusions (unique participants):",
      sum(result_data$extrasystolic_exclusion[extra_unique] == 1L, na.rm = TRUE), "\n")
  cat("HEP 3-SD outliers, each pipeline (unique participants):\n")
  print(hep_outlier_counts)
  cat("Prior EEG recording inventory loaded:", !is.null(inventory), "\n")
  if (!is.null(participants)) cat("Earlier LIFE matching participants: ", nrow(participants), "\n")
  if (!is.null(prior_duplicates)) cat("Earlier LIFE repeated-row cases:    ", nrow(prior_duplicates), "\n")
  cat("\n---------- External data coverage ----------\n")
  print(summary_coverage[, c("column", "participants_found", "participants_with_value",
                             "source_duplicate_ids"), drop = FALSE], row.names = FALSE)
  cat("\nCER-S filenames interpreted from column:", file_col, "\n")
  cat("CER-S invalid/non-numeric count cells ignored:",
      paste(paste(names(invalid_beats), invalid_beats, sep = "="), collapse = "; "), "\n")
  cat("Final rows flagged non-unique:", sum(result_data[["non-unique"]] == 1L), "\n")
  cat("GAD-7 matched directly:              ", gad_direct_n, "\n")
  cat("GAD-7 matched via IDmatch fallback:  ", gad_fallback_n, "\n")
  cat("GAD-7 ambiguous cases left missing:  ", gad_ambiguous_n, "\n")
  cat("GAD-7 participants without a match:  ", sum(lengths(gad_hits) == 0L), "\n")
  cat("External matching issues needing review:", length(issues), "\n")
  if (length(issues)) {
    print(table(vapply(issues, function(x) x$source[1L], character(1L))))
  }
  cat("\nNon-unique participants and EEG dates: see data_unresolved_rows.csv\n")
  cat("EEG/VIGALL primary row is preserved; VIGALL is NEVER backfilled.\n")
  cat("External source values are replicated across all surviving rows of an ID.\n")
})
writeLines(report, file.path(out, "data_summary.txt"))
cat(paste(report, collapse = "\n"), "\n")
cat("\nSaved analysis data to: ", data_path, "\n", sep = "")
cat("Saved data_row_provenance.csv, data_unresolved_rows.csv, ",
    "data_match_issues.csv, data_summary.txt to: ", out, "\n", sep = "")
